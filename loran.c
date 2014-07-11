#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <alloca.h>

#include "loran.h"
#include "loran_reference.h"

/* Phase codes for GRI A and B */

int lc_pc[LC_PC_CNT][LC_GRCNT * LC_STA_GRIWCNT] = {
	{ 1, 1,-1,-1, 1,-1, 1,-1,  1,-1,-1, 1, 1, 1, 1, 1},	/* MasA */
	{ 1,-1,-1, 1, 1, 1, 1, 1,  1, 1,-1,-1, 1,-1, 1,-1},	/* MasB	*/
	{ 1, 1, 1, 1, 1,-1,-1, 1,  1,-1, 1,-1, 1, 1,-1,-1},	/* SlaveA */
	{ 1,-1, 1,-1, 1, 1,-1,-1,  1, 1, 1, 1, 1,-1,-1, 1},	/* SlaveB */
	{ 1, 0, 0, 0,-1, 0, 0, 0,  1, 0, 0, 0,-1, 0, 0, 0}	/* Synthetic */
};

int lc_pc_acorr[LC_PC_CNT] = {16, 16, 16, 16, 4};

/*
 * Comb filter delays(minimal FRI count) for output level (gamma) = 0.7
 */
size_t lc_comb_delay[LC_COMB_KP_MAX + 1] = {
	1, 2, 5, 9, 19, 38, 77, 154, 308, 616, 1233
     /* 0  1  2  3  4   5   6   7    8    9    10 */
};

static void	lc_calc_free_intervals(struct lc_chain *chain);
static _Bool	lc_is_offset_free(struct lc_chain *chain, size_t offset);
static void	lc_calc_seek_ints(struct lc_chain *chain);
static void	lc_start_seek(struct lc_chain *chain);

void lc_init(struct lc_chain *chain)
{
	assert(chain->gri >= 4000 && chain->gri < 10000);
	assert(chain->station_cnt > 0 && chain->station_cnt <= LC_STA_MAXCNT);

	chain->grin = (LC_FS / 100000) * chain->gri;
	chain->frin = LC_GRCNT * chain->grin;

	chain->next_sample = 0;

	chain->seek_cand_cnt = 0;
	chain->free_int_cnt = 0;

	for(size_t i = 0; i < chain->station_cnt; ++i) {
		struct lc_station *sta = &chain->sta[i];

		sta->idx = i;
		sta->offset = 0;
		sta->fris_passed = 0;
		sta->last_gri_idx = 0;
		sta->state = LC_STAST_IDLE;

		sta->comb_kp = 5;

		for(size_t j = 0; j < LC_PC_CNT; ++j) {
			memset(sta->xcorr_tmp[j], 0,
			       LC_STA_WINSZ * sizeof(lc_type_data_s));
			memset(sta->xcorr[j], 0,
			       LC_STA_WINSZ * sizeof(lc_type_comb_s));
			sta->xcorr_var[j] = 0;
		}
		memset(sta->win0, 0,
		       LC_STA_WINSZ * sizeof(lc_type_comb_s));
		sta->win0_var = 0;
	}

	printf("------------------------\n"
	       "LORAN-C CONSTANTS\n"
	       "- FS: %u\n"
	       "- WININT: %u\n"
	       "- WINSZ: %u\n"
	       "- WINSHIFT: %u\n"
	       "- STA_SPAN: %u\n"
	       "- sizeof(lc_chain): %u\n"
	       "- sizeof(lc_station): %u\n"
	       "------------------------\n",
	       LC_FS, LC_STA_WININT, LC_STA_WINSZ,
	       LC_STA_WINSHIFT, LC_STA_SPAN,
	       sizeof(struct lc_chain), sizeof(struct lc_station));

	lc_calc_free_intervals(chain);
	lc_start_seek(chain);
}

/*
 * Calculate spanzones given station offset
 */

static void lc_calc_spans(struct lc_chain *chain, struct lc_station *sta)
{
	assert(sta->offset < chain->grin);

	for(size_t i = 0; i < LC_GRCNT; ++i) {
		sta->spans[i].begin = sta->offset + chain->grin * i;
		sta->spans[i].end = sta->spans[i].begin + LC_STA_SPAN;
		sta->spans[i].win_idx = LC_STA_GRIWCNT * i;
		sta->spans[i].win_offset = 0;
	}

	sta->spans_cnt = LC_GRCNT;

	if(sta->spans[LC_GRCNT-1].end > chain->frin) {
		struct lc_spanzone *nsz = &sta->spans[LC_GRCNT];
		nsz->begin = 0;
		nsz->end = sta->spans[LC_GRCNT-1].end - chain->frin;
		/* Calculate window index and offset into it for new spanzone */
		size_t sploff = chain->frin - sta->spans[LC_GRCNT-1].begin;
		size_t wcnt = sploff / LC_STA_WININT;
		size_t woff = sploff - wcnt * LC_STA_WININT;
		nsz->win_idx = sta->spans[LC_GRCNT-1].win_idx + wcnt;
		nsz->win_offset = woff;

		sta->spans[LC_GRCNT-1].end = chain->frin;
		sta->spans_cnt = LC_GRCNT + 1;
	}

	sta->span_next = 0;
	sta->span0_open = 0;
}

/*
 * Loran-C timing constraints (USCG Loran-C signal spec, 1994)
 *
 * Minimum TD: 10.9 ms
 * Minimum difference between TDs: 9.9 ms
 * Maximum TD: GRI - 9.9 ms
 *
 * These constraints are rather simple so can be directly included in calcula-
 * tion of "allowed" ("free") intervals in GRI
 */

/*
 * Calculates allowed ("free") intervals boundaries on GRI
 *
 * Takes into account timing constraints described above
 */
static void lc_calc_free_intervals(struct lc_chain *chain)
{
	/* Sort station offsets in ascending order */
	int		*offs;
	offs = alloca(sizeof(int) * chain->station_cnt);

	size_t busy_cnt = 0;
	for(size_t i = 0; i < chain->station_cnt; ++i) {
		if(!LC_STA_IS_BUSY(&chain->sta[i]))
			continue;

		size_t j;
		int tmp = chain->sta[i].offset;

		for(j = busy_cnt; j > 0 && offs[j-1] > tmp; j--)
			offs[j] = offs[j-1];
		offs[j] = tmp;
		busy_cnt++;
	}

	printf("Offs: ");
	for(size_t i = 0; i < busy_cnt; ++i)
		printf("%d ", offs[i]);
	printf("\n");

	/*
	 * Each station applies following constraints:
	 *  - ending of preceding "free" interval: sta->offset - 9.9ms
	 *  - start of following "free" interval: sta->offset + 9.9ms
	 */
	chain->free_int_cnt = 0;
	int intbegin, intend;
	if(busy_cnt)
		intbegin = offs[busy_cnt - 1] + LC_INT_99MS;
	else
		intbegin = 0;
	intbegin = intbegin > (int) chain->grin ? intbegin - chain->grin : 0;
	for(size_t i = 0; i < busy_cnt; ++i) {
		intend = offs[i] - LC_INT_99MS + 1;
		if(intend > intbegin) {
			chain->free_int[chain->free_int_cnt].begin = intbegin;
			chain->free_int[chain->free_int_cnt].end = intend;
			chain->free_int_cnt++;
		}
		intbegin = offs[i] + LC_INT_99MS;
	}

	if(busy_cnt)
		intend = offs[0] - LC_INT_99MS + 1;
	else
		intend = chain->grin;
	intend = intend < 0 ? intend + chain->grin : chain->grin;

	if(intend > intbegin) {
		chain->free_int[chain->free_int_cnt].begin = intbegin;
		chain->free_int[chain->free_int_cnt].end = intend;
		chain->free_int_cnt++;
	}
}

/*
 * Checks, whether offset is in one of "free" intervals, or not
 * Can be used to test, whether station position meets Loran-C constraints
 */
static _Bool lc_is_offset_free(struct lc_chain *chain, size_t offset)
{
	for(size_t i = 0; i < chain->free_int_cnt; ++i) {
		struct lc_int	*intv = &chain->free_int[i];
		if(offset >= intv->begin && offset < intv->end)
			return 1;
	}
	return 0;
}

/*
 * Distribute (more or less equally) intervals to seek through among IDLE
 *  stations
 */
static void lc_calc_seek_ints(struct lc_chain *chain)
{
	size_t		free_size = 0;
	for(size_t i = 0; i < chain->free_int_cnt; ++i)
		free_size+= chain->free_int[i].end - chain->free_int[i].begin;

	size_t		idle_stations = 0;
	for(size_t i = 0; i < chain->station_cnt; ++i)
		if(LC_STA_IS_IDLE(&chain->sta[i]))
			idle_stations++;

	printf("Free size: %u, idle stations: %u\n", free_size, idle_stations);
	if(!idle_stations)
		return;

	size_t		sta_intsz = free_size / idle_stations;

	struct lc_station *sta;
	struct lc_int	*intv;
	struct lc_int	*sta_intv;

	size_t		assigned_total = 0;
	size_t		assigned_sta;
	size_t		intv_current;

	/* Find first non-busy station */
	sta = &chain->sta[0];
	while(!LC_STA_IS_IDLE(sta) && sta < chain->sta + chain->station_cnt)
		sta++;
	assert(sta < chain->sta + chain->station_cnt);
	idle_stations--;
	sta->seek_ints_cnt = 0;
	sta_intv = &sta->seek_ints[0];
	assigned_sta = 0;

	/* Get first interval */
	intv = &chain->free_int[0];
	intv_current = intv->begin;

	while(assigned_total < free_size) {
		if(assigned_sta >= sta_intsz) {
			/* Get next station */
			do {
				sta++;
			} while(sta < chain->sta + chain->station_cnt &&
				!LC_STA_IS_IDLE(sta));

			assert(sta < chain->sta + chain->station_cnt);

			idle_stations--;
			sta->seek_ints_cnt = 0;
			sta_intv = &sta->seek_ints[0];

			assigned_sta = 0;
			printf("STA CHANGE: %d\n", sta - chain->sta);
		}
		if(intv_current >= intv->end) {
			/* Get next interval */
			intv++;
			intv_current = intv->begin;
			printf("INTV CHANGE\n");
		}

		/*
		 * demand is amount of space that station can get
		 * offer is amount of space that current interval offers to
		 *  station(s)
		 */
		size_t demand = sta_intsz - assigned_sta;
		size_t offer = intv->end - intv_current;

		/* Last free station handles rest of intervals */
		if(!idle_stations)
			demand = free_size - assigned_total;

		size_t assigned;
		assigned = offer > demand ? demand : offer;

		sta_intv->begin = intv_current;
		sta_intv->end = sta_intv->begin + assigned;

		printf("ASSIGNED: %u .. %u\n", sta_intv->begin, sta_intv->end);

		intv_current+= assigned;
		assigned_sta+= assigned;
		assigned_total+= assigned;

		sta->seek_ints_cnt++;
		sta_intv++;
	}

	/*
	 * If there are "free" stations and no intervals assigned to them,
	 * don't forget to write down that fact
	 */
	while(++sta < chain->sta + chain->station_cnt) {
		if(LC_STA_IS_IDLE(sta))
			sta->seek_ints_cnt = 0;
	}
}

/*
 * Perform candidate selection based on set of candidates, Loran timing
 *  constraints and points that are already in LOCKING/TRACKING state
 */
static void lc_select_candidates(struct lc_chain *chain)
{
	if(!chain->seek_cand_cnt)
		return;
	assert(chain->seek_cand_cnt <= LC_CAND_MAXCNT);

	uint8_t *idxs = alloca(sizeof(uint8_t) * chain->seek_cand_cnt);

	for(size_t i = 0; i < chain->seek_cand_cnt; ++i) {
		uint32_t sta_off = chain->seek_cand[i].offset;

		size_t j;
		for(j = i; j > 0 &&
		    chain->seek_cand[idxs[j - 1]].offset > sta_off; --j)
			idxs[j] = idxs[j - 1];
		idxs[j] = i;
	}

	printf("\tSORTED: ");
	for(size_t i = 0; i < chain->seek_cand_cnt; ++i)
		printf("%u ", idxs[i]);
	printf("\n");
	/*
	 * Task: leave candidate set that satisfies Loran-C timing constaints
	 *  i.e. there can be two (or more) candidates distributed in a way that
	 *  beginning of one falls into protected 9.9ms interval after beginning
	 *  of other.
	 * These "overlappings" can form chain-like sets on a GRI "wheel" (GRI
	 *  resembles a wheel, because the end of GRI interval is followed
	 *  immediately by its beginning thanks to periodic nature of Loran
	 *  signal).
	 * To accomplish the task simple algorithm is presented: even if there
	 *  are lot of candidates, and they overlap with each other forming
	 *  large "chain" of overlappings on GRI, there is also interval free
	 *  of any candidates, so that "chain" of overlappings is not a cycle.
	 * First we try to find such position, that is also a true "beginning"
	 *  of chained overlaps, and then go through chain, removing certain
	 *  candidates, thus satisfying Loran constraints.
	 *
	 * If algorithm is unable to find such position, it means that
	 *  overlappings "chain" is cyclic and there is no single solution that
	 *  is based only on candidates' offsets...
	 */

	uint32_t restr_end = 0;
	size_t start_idx = 0;
	for(; start_idx < chain->seek_cand_cnt; ++start_idx) {
		uint32_t off = chain->seek_cand[idxs[start_idx]].offset;

		if(off >= restr_end && off >= LC_INT_99MS) {
			/*
			 * We've found candidate not overlapped by any previous
			 */
			break;
		}
		restr_end = off + LC_INT_99MS;
	}

	/*
	 * Iterate from start_idx over all candidates (in cycle), and set
	 * passed/not passed flag for each of them
	 */
	restr_end = 0;
	size_t i = start_idx;
	while(1) {
		struct lc_cand *cand = &chain->seek_cand[idxs[i]];

		if(cand->offset >= restr_end) {
			/* Current candidate passes selection.. */
			cand->passed = 1;
			restr_end = off + LC_INT_99MS;
		} else {
			/* Doesn't pass - discard */
			cand->passed = 0;
		}

		++i;
		if(i >= chain->seek_cand_cnt) {
			i = 0;
			if(restr_end > chain->grin)
				restr_end -= chain->grin;
			else
				restr_end = 0;
		}
		if(i == start_idx)
			break;
	}
}


/*
 * Initiate seek process for all IDLE stations
 *
 * Assumes that "free" intervals are already calculated
 */
static void lc_start_seek(struct lc_chain *chain)
{
	lc_calc_seek_ints(chain);

	chain->seek_sta_cnt = 0;

	for(size_t i = 0; i < chain->station_cnt; ++i) {
		struct lc_station *sta = &chain->sta[i];

		assert(sta->state != LC_STAST_SEEK);

		if(!LC_STA_IS_IDLE(sta))
			continue;
		if(!sta->seek_ints_cnt)
			continue;

		printf("#SEEK# STA %u begins seeking through %u intervals: ",
		       i, sta->seek_ints_cnt);
		for(size_t j = 0; j < sta->seek_ints_cnt; ++j)
			printf("%u .. %u, ", sta->seek_ints[j].begin,
			       sta->seek_ints[j].end);
		printf("\n");

		sta->seek_cur_int = &sta->seek_ints[0];
		sta->seek_complete = 0;
		sta->comb_kp = 6;	/* Default SEEK value */
		sta->offset = sta->seek_cur_int->begin;

		sta->fris_passed = 0;
		sta->last_gri_idx = 0;
		sta->ssr_call_period = lc_comb_delay[sta->comb_kp];

		sta->lock_enabled = 0;

		sta->state = LC_STAST_SEEK;

		lc_calc_spans(chain, sta);

		chain->seek_sta_cnt++;
	}
}

/*
 * Called when all SEEK stations have processed their seek intervals and went
 *  into IDLE state
 *
 * Task: process results of seek - filter set of candidates and assign tasks
 *  to IDLE stations, put them into LOCKING phase, and restart SEEK process
 *  for stations that were left IDLE, possibly adjusting KP and other params.
 */

static void lc_seek_completed(struct lc_chain *chain)
{
	printf("#SEEK# COMPLETED! We have %u candidates\n",
	       chain->seek_cand_cnt);
	for(size_t i = 0; i < chain->seek_cand_cnt; ++i) {
		printf("\t@ %5u: PC %u, ACC %6u\n",
		       chain->seek_cand[i].offset, chain->seek_cand[i].pc,
		       chain->seek_cand[i].accuracy);
	}

	lc_select_candidates(chain);

	/*
	 * Here we should assign tasks to IDLE stations and put them into
	 * LOCKING state
	 */

	lc_calc_free_intervals(chain);
	lc_calc_seek_ints(chain);

	assert(0);
}

/*******************************************************************************
 * Station service routines - called every ssr_call_period FRIs
 */

static void lc_sta_service_locking(struct lc_chain *chain,
				   struct lc_station *sta)
{
	assert(0);	/* Die */
}

/*
 * Station service routine: seek process
 *
 * Tries to proceed with seeking at new offset, maybe on new interval
 *  If no more seek intervals are available, sets "seek_complete" flag and
 *  changes state to "IDLE"
 */
static void lc_sta_service_seek(struct lc_chain *chain, struct lc_station *sta)
{
	assert(sta->state == LC_STAST_SEEK);
	assert(chain->seek_sta_cnt != 0);

	/* Process gathered data */
	lc_process_seek(chain, sta);

	/* Try to advance station offset */
	size_t		new_offset;
	new_offset = sta->offset + LC_STA_WINSHIFT;

	if(new_offset >= sta->seek_cur_int->end) {
		/* Try to proceed within next interval */
		if(++(sta->seek_cur_int) < sta->seek_ints + sta->seek_ints_cnt){
			new_offset = sta->seek_cur_int->begin;
		} else {
			/* That was last interval available: seek is over.. */
			sta->seek_complete = 1;
			sta->state = LC_STAST_IDLE;
			sta->ssr_call_period = 0;

			printf("#SEEK# STA %d completed seek\n",
			       sta - chain->sta);
			chain->seek_sta_cnt--;
			if(!chain->seek_sta_cnt) {
				/*
				 * All SEEK stations have completed seek at the
				 * moment
				 */
				lc_seek_completed(chain);
			}

//			if(sta->lock_enabled) {
//				if(sta->lock_point < LC_STA_WINSHIFT)
//					sta->offset = sta->lock_point + chain->grin - LC_STA_WINSHIFT;
//				else
//					sta->offset = sta->lock_point - LC_STA_WINSHIFT;
//				sta->ssr_call_period = 1;
//				sta->state = LC_STAST_LOCKING;
//				printf("#LOCK# STA %d locking at %u\n",
//				       sta - chain->sta, sta->offset);
//			}

			return;
		}
	}

	/* Clear windows :3 */
	memset(sta->win0, 0, sizeof(lc_type_comb_s) * LC_STA_WINSZ);
	for(size_t i = 0; i < LC_PC_CNT; ++i)
		memset(sta->xcorr[i], 0, sizeof(lc_type_comb_s) * LC_STA_WINSZ);

	sta->fris_passed = 0;
	sta->offset = new_offset;
	lc_calc_spans(chain, sta);
	sta->span_next = 0;
}

/*
 * Station service routine: called each ssr_call_period FRIs
 */
static void lc_sta_service(struct lc_chain *chain, struct lc_station *sta)
{
	switch(sta->state) {
	case LC_STAST_IDLE:
		/* Disable service routine calls for IDLE stations */
		printf("Hush! SSR called for IDLE STA\n");
		sta->ssr_call_period = 0;
		break;
	case LC_STAST_SEEK:
		lc_sta_service_seek(chain, sta);
		break;
	case LC_STAST_BUSY:
		break;
	case LC_STAST_LOCKING:
		lc_sta_service_locking(chain, sta);
		break;
	}
}

static void lc_sta_xcorr_tmp_filter(struct lc_station *sta)
{
	for(size_t pc = 0; pc < LC_PC_CNT; ++pc) {
		lc_type_data_s	*src, *src_end;
		lc_type_comb_s	*dst;

		src = sta->xcorr_tmp[pc];
		src_end = sta->xcorr_tmp[pc] + LC_STA_WINSZ;
		dst = sta->xcorr[pc];

		while(src < src_end) {
			lc_type_comb_s cs = (lc_type_comb_s) *src << LC_BITS_COMBF;
			cs = cs / lc_pc_acorr[pc];
			cs = (*dst * ((1 << sta->comb_kp) - 1) + cs) / (1 << sta->comb_kp);
			*dst++ = cs;
			*src++ = 0;	/* Filling xcorr_tmp with zeros */
		}
	}
}

/*
 * Process new samples targeted at particular integration window
 *
 * Almost all immediate signal filtering happens here
 */
static void lc_win_new_samples(struct lc_station *sta, size_t win_idx,
			       size_t data_offset,
			       lc_type_sample *data, size_t count)
{
	lc_type_comb_s		*wptr;
	lc_type_sample		*data_end;

	assert(win_idx < LC_STA_WINCNT);
	assert(data_offset < LC_STA_WINSZ);

//	printf("STA %u WINDOW %2u OFF %5u SIZE %4u\n",
//	       sta->idx, win_idx, data_offset, count);
	wptr = sta->win0 + data_offset;
	data_end = data + count;

	lc_type_sample	sd = sta->dcrem_buf;

	while(data < data_end) {
		lc_type_sample	s = *data++;		/* Take sample */
		lc_type_data_s	ss;			/* DC notch output */
		/* DC removing notch */
		s+= sd - (sd >> LC_DCREM_KP);		/* Add delayed sample
							   multiplied by (1-k)*/
		ss = s - sd;				/* Subtract delayed
							   sample */
		sd = s;					/* Write to delay buf */
		ss-= ss >> (LC_DCREM_KP + 1);		/* Mult by 'b' */

		/*
		 * Now, we should do the following:
		 * -accumulate sample to sta->xcorr_tmp windows, adding or sub-
		 *  tracting it based on value of corresponding PC bit
		 * -if it is window 0 (first pulse of station, first bit of PC)
		 *  then update comb-filtered version of 1st pulse in sta->win0
		 */

		for(size_t pc = 0; pc < LC_PC_CNT; ++pc) {
			int pcbit = lc_pc[pc][win_idx];
			if(pcbit > 0)
				sta->xcorr_tmp[pc][data_offset] += ss;
			else if(pcbit < 0)
				sta->xcorr_tmp[pc][data_offset] -= ss;
		}

		if(win_idx == 0) {
			/* Comb filter - integrator over FRI */
			lc_type_comb_s cs;
			cs = (lc_type_comb_s) ss << LC_BITS_COMBF; /* Convert to fixed */

			//ss = ((*wptr << sta->comb_kp) - *wptr + ss) >> sta->comb_kp;
			cs = (*wptr * ((1 << sta->comb_kp) - 1) + cs) / (1 << sta->comb_kp);
			*wptr++ = cs;
		}
		data_offset++;
	}
	sta->dcrem_buf = sd;
}

/*
 * Handle samples targeted at particular station, at particular GRI window group
 *
 * Redirects these samples to individual windows...
 */
static void lc_sta_new_samples(struct lc_chain *chain, struct lc_station *sta,
			       size_t data_offset, lc_type_sample *data,
			       size_t count)
{
	struct lc_spanzone *sz = &sta->spans[sta->span_next];

	size_t		startw, widx;
	size_t		win_offset, data_end;

	if(!count)
		return;

	assert(data_offset + count <= sz->end - sz->begin);

//	if(data_offset == 0 && sta->span_next == 0)
//		sta->span0_open = 1;
//	if(!sta->span0_open)
//		return;

//	printf("\nSTA %u SPANZONE %u OFF %5u SIZE %4u\n",
//	       sta->idx,  sta->span_next, data_offset, count);
	/* Here data_offset is offset from beginning of _spanzone_ */
	win_offset = data_offset + sz->win_offset;
	startw = win_offset / LC_STA_WININT;
	win_offset-= startw * LC_STA_WININT;
	assert(win_offset < LC_STA_WININT);
	startw+= sz->win_idx;
//	printf("\tSTARTW %u WOFF %u\n", startw, win_offset);
	/*
	 * Now win_offset is offset from beginning of _window_ startw,
	 * and startw is first window that new data belongs to
	 */

	data_end = win_offset + count;

	for(widx = startw;; ++widx) {
		if(win_offset < LC_STA_WINSZ) {
			size_t ov_end = LC_MIN(data_end, LC_STA_WINSZ);
			lc_win_new_samples(sta, widx, win_offset,
					   data, ov_end - win_offset);
		}

		if(data_end <= LC_STA_WININT)
			break;

		size_t readcount = LC_STA_WININT - win_offset;

		data_end-= LC_STA_WININT;
		data+= readcount;
		win_offset = 0;
	}

	/* If data fills spanzone completely, switch to next spanzone */
	if(data_offset + count == sz->end - sz->begin) {
		if(sta->span_next == sta->spans_cnt - 1) {
			sta->span_next = 0;
			/*
			 * Increment passed FRIs counter, cause we've just
			 * filled last spanzone
			 */
			sta->fris_passed++;
			/* Perform comb-filtering of gathered xcorr_tmp wins */
			lc_sta_xcorr_tmp_filter(sta);
		} else {
			sta->span_next++;
		}
	}

	/* Call SSR if enough FRIs have passed through */
	if(sta->ssr_call_period && sta->fris_passed >= sta->ssr_call_period) {
		lc_sta_service(chain, sta);
		sta->fris_passed = 0;
	}
}

/*
 * Handle new raw samples for chain
 *
 * Tasks:
 *  -place supplied samples at the right place (redirect to correct stations)
 *  -handle station that is overlappping end of FRI
 *  -handle too large (several FRIs) sample buffers
 */
void lc_new_samples(struct lc_chain *chain, lc_type_sample *buf, size_t count)
{
	size_t buf_begin, buf_end;

	if(!count)
		return;

	while(count) {
		buf_begin = chain->next_sample;
		buf_end = buf_begin + count;
		if(buf_end > chain->frin)
			buf_end = chain->frin;

		for(size_t i = 0; i < chain->station_cnt; ++i) {
			struct lc_station *sta = &chain->sta[i];

			struct lc_spanzone *sz = &sta->spans[sta->span_next];

			size_t ov_begin = LC_MAX(sz->begin, buf_begin);
			size_t ov_end = LC_MIN(sz->end, buf_end);

			if(ov_begin >= ov_end)
				continue;

			/* We have buffer overlapping with spanzone.. */
			lc_sta_new_samples(chain, sta, ov_begin - sz->begin,
					   buf + ov_begin - buf_begin,
					   ov_end - ov_begin);
		}

		size_t read_cnt;
		read_cnt = buf_end - buf_begin;
		buf+= read_cnt;
		count-= read_cnt;
		chain->next_sample = 0;
	}
	chain->next_sample = buf_end;

	if(chain->next_sample >= chain->frin)
		chain->next_sample-= chain->frin;
}
