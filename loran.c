#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <alloca.h>

#include "loran.h"
#include "loran_reference.h"

/* Phase codes for GRI A and B */

_Bool lc_pc[LC_PC_CNT][LC_GRCNT * LC_STA_GRIWCNT] = {
	{1, 1, 0, 0, 1, 0, 1, 0,  1, 0, 0, 1, 1, 1, 1, 1},	/* MasA */
	{1, 0, 0, 1, 1, 1, 1, 1,  1, 1, 0, 0, 1, 0, 1, 0},	/* MasB	*/
	{1, 1, 1, 1, 1, 0, 0, 1,  1, 0, 1, 0, 1, 1, 0, 0},	/* SlaveA */
	{1, 0, 1, 0, 1, 1, 0, 0,  1, 1, 1, 1, 1, 0, 0, 1}	/* SlaveB */
};

/*
 * Comb filter delays(minimal FRI count) for output level (gamma) = 0.7
 */
size_t lc_comb_delay[LC_COMB_KP_MAX + 1] = {
	0, 2, 5, 9, 19, 38, 77, 154, 308, 616, 1233
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

	chain->free_int_cnt = 0;
	chain->next_sample = 0;

	for(size_t i = 0; i < chain->station_cnt; ++i) {
		struct lc_station *sta = &chain->sta[i];

		sta->offset = 0;
		sta->fris_passed = 0;
		sta->last_gri_idx = 0;
		sta->state = LC_STAST_IDLE;

		sta->comb_kp = 5;

		for(size_t j = 0; j < LC_GRCNT; ++j)
			for(size_t k = 0; k < LC_STA_GRIWCNT; ++k) {
				size_t w = k + j * LC_STA_GRIWCNT;
				sta->wnd[w].offset = k * LC_STA_WININT;
				sta->wnd[w].offset+= j * chain->grin;

				sta->wnd[w].dcrem_buf = 0;
				memset(sta->wnd[w].data, 0,
				       LC_STA_WINSZ * sizeof(lc_type_data_s));
			}
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
		sta->comb_kp = 7;	/* Default SEEK value */
		sta->offset = sta->seek_cur_int->begin;

		sta->fris_passed = 0;
		sta->last_gri_idx = 0;
		sta->ssr_call_period = lc_comb_delay[sta->comb_kp + 1];

		sta->state = LC_STAST_SEEK;

		chain->seek_sta_cnt++;
	}
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

	/* Calculate cross-correlation */
	lc_corr_sta(sta);

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

			printf("#SEEK# STA %u completed seek\n",
			       sta - chain->sta);
			chain->seek_sta_cnt--;
			if(!chain->seek_sta_cnt) {
				/*
				 * All SEEK stations have completed seek at the
				 * moment
				 */
				printf("#SEEK# COMPLETED!\n");
			}

			return;
		}
	}

	sta->fris_passed = 0;
	sta->last_gri_idx = 0;
	sta->offset = new_offset;
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
	struct lc_subwin *win;
	lc_type_data_s	*wptr;
	lc_type_sample	*data_end;

	assert(win_idx < LC_GRCNT * LC_STA_GRIWCNT);
	assert(data_offset < LC_STA_WINSZ);

	win = &sta->wnd[win_idx];
	wptr = win->data + data_offset;
	data_end = data + count;

	lc_type_sample	sd = win->dcrem_buf;
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
		/* Comb filter - integrator over FRI */
		ss<<= LC_BITS_COMBF;			/* Convert to fixed */

		//ss = ((*wptr << sta->comb_kp) - *wptr + ss) >> sta->comb_kp;
		ss = (*wptr * ((1 << sta->comb_kp) - 1) + ss) / (1 << sta->comb_kp);
		*wptr++ = ss;
	}
	win->dcrem_buf = sd;
}

/*
 * Handle samples targeted at particular station, at particular GRI window group
 *
 * Redirects these samples to individual windows...
 */
static void lc_sta_new_samples(struct lc_chain *chain, struct lc_station *sta,
			       size_t gri_idx, size_t data_offset,
			       lc_type_sample *data, size_t count)
{
	size_t		griw_base, startw_idx, widx;
	size_t		startw_offset;
	size_t		data_end;

	if(!count)
		return;
	assert((data_offset + count) <= LC_STA_SPAN);

	griw_base = gri_idx * LC_STA_GRIWCNT;
	startw_idx = data_offset / LC_STA_WININT;
	assert(startw_idx < LC_STA_GRIWCNT);

	if(sta->offset == 33356 && 0)
		printf("Loran station samples:\n"
		       "\tSTA offset: %u\n"
		       "\tData GRI index: %u\n"
		       "\tData offset: %u\n"
		       "\tData length: %u\n"
		       "\tWindow begin: %u\n",
		       sta->offset, gri_idx,
		       data_offset,
		       count, startw_idx + griw_base);


	startw_offset = startw_idx * LC_STA_WININT;
	assert(data_offset >= startw_offset);
	data_offset-= startw_offset;
	assert(data_offset < LC_STA_WININT);
	data_end = data_offset + count;

	for(widx = (startw_idx + griw_base);; ++widx) {
		/*
		printf("\tDATA for window %u\n"
		       "\t\tStart of data: %u\n"
		       "\t\tEnd of data: %u\n",
		       i, data_offset, data_end);
		       */
		if(data_offset < LC_STA_WINSZ) {
			size_t spansize = data_end > LC_STA_WINSZ ?
						(LC_STA_WINSZ - data_offset):
						(data_end - data_offset);
			lc_win_new_samples(sta, widx, data_offset,
					   data, spansize);
		}

		size_t readcount;
		if(data_end <= LC_STA_WININT)
			//readcount = data_end - data_offset;
			break;
		readcount = LC_STA_WININT - data_offset;

		//printf("\tRead: %u samples\n", readcount);

		//if(data_end <= LC_STA_WININT)
		//	break;

		data_end-= LC_STA_WININT;
		data+= readcount;
		data_offset = 0;
	}

	sta->last_gri_idx = gri_idx;
	/*
	 * If station has yet filled in last window in FRI, check for this
	 * and call service routine if needed
	 */
	if(widx == LC_STA_WINCNT-1 && data_end >= LC_STA_WINSZ) {
		sta->fris_passed++;
		if(sta->ssr_call_period &&
		   sta->fris_passed >= sta->ssr_call_period) {
			lc_sta_service(chain, sta);
			sta->fris_passed = 0;
		}
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
			for(size_t j = 0; j < LC_GRCNT; ++j) {
				/* Station region boundaries */
				size_t sta_begin = sta->offset +
						   j * chain->grin;
				size_t sta_end = sta_begin + LC_STA_SPAN;

				/**printf("##STA %u\tGRI %u\tBEGIN: %u\tEND: %u###\n",
				       i, j, sta_begin, sta_end);*/

				/* Overlap region boundaries */
				size_t ov_begin, ov_end;
				/*
				 * Length of station region that goes beyond
				 * FRI boundary
				 */
				size_t ovi_span = sta_end - chain->frin;
				/*
				 * Remap data from near beginning of FRI beyond
				 * the end of FRI and pass it to corresponding
				 * "overlapping" station
				 */
				if(sta_end > chain->frin &&
				   sta->last_gri_idx == LC_GRCNT-1 &&
				   buf_begin < ovi_span) {
					/**printf("##OVERLAPPING##\n");*/
					ov_begin = buf_begin;
					ov_end = buf_end > ovi_span ?
							 ovi_span : buf_end;
					lc_sta_new_samples(chain, sta, j,
							   chain->frin + ov_begin - sta_begin,
							   buf,
							   ov_end - ov_begin);
				}
				/*
				 * Process non-"overlapping" station data
				 */
				ov_begin = sta_begin > buf_begin ?
						   sta_begin : buf_begin;
				ov_end = buf_end > sta_end ?
						 sta_end : buf_end;
				if(ov_end > ov_begin) {
					/**printf("##NORMAL##\n");*/
					lc_sta_new_samples(chain, sta, j,
							   ov_begin - sta_begin,
							   buf + ov_begin - buf_begin,
							   ov_end - ov_begin);
				}
			}
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

