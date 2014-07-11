#include "loran.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <float.h>
#include <assert.h>

#define LC_HWAVE_SZT	5U
#define LC_HWAVE_SZ	(LC_HWAVE_SZT * LC_FS / 1000000)
//#define LC_WIN_HWAVESCNT (LC_STA_WINSZ * 2 * 100000 / LC_FS)
#define LC_WIN_HWAVESCNT LC_STA_WINSZ

#define LC_SGN(x) ((x) < 0 ? -1 : 1)

struct lc_halfwave
{
	/* Last sample index belonging to halfwave */
	uint16_t	last;
	/* Magnitude of halfwave */
	lc_type_data_s	mag;
};

struct lc_halfwave_buffer
{
	size_t first_sample;
	size_t hw_count;
	struct lc_halfwave hwaves[LC_WIN_HWAVESCNT];
};

static struct lc_halfwave_buffer pulwinhw;

/*
 * Support for partial window processing: generates halfwaves for samples
 *  [first, last], probably including those before start and after end belonging
 *  to first/last halfwave to correctly determine their magnitude
 */
static void lc_hwaves_win_gen(lc_type_comb_s *win, size_t first, size_t last,
			      struct lc_halfwave_buffer *hwbuf)
{
	size_t send, sstart;
	struct lc_halfwave *hwcur, *hwend;

	assert(last >= first);

	/* Determine first sample */
	for(sstart = first; sstart > 0; --sstart)
		if(LC_SGN(win[sstart]) == LC_SGN(win[sstart-1]))
			continue;

	hwbuf->first_sample = sstart;

	hwcur = hwbuf->hwaves;
	hwend = hwcur + LC_WIN_HWAVESCNT;

	for(size_t i = sstart; i < (LC_STA_WINSZ - 1) && hwcur < hwend; ++i) {
		if(LC_SGN(win[i]) == LC_SGN(win[i + 1]))
			continue;

		lc_type_comb_s hwamp = 0;
		send = i;
		/* Determine halfwave amplitude - interpolation is welcome! */
		for(size_t j = sstart; j <= send; ++j)
			if(abs(win[j]) > abs(hwamp))
				hwamp = win[j];

		hwcur->mag = hwamp >> LC_BITS_COMBF;
		hwcur->last = send;

		if(send >= last)
			break;
		hwcur++;
		sstart = i + 1;
	}
	/* Finish halfwaves list */
	hwbuf->hw_count = hwcur - hwbuf->hwaves + 1;
}

static size_t lc_hwave_find_peak(size_t hwidx, lc_type_comb_s *win,
				 struct lc_halfwave_buffer *hwbuf)
{
	assert(hwidx < LC_WIN_HWAVESCNT);

	size_t	start = hwidx ? (hwbuf->hwaves[hwidx-1].last + 1U) :
			hwbuf->first_sample;
	size_t	last = hwbuf->hwaves[hwidx].last;

	assert(last >= start);

	int	maxabs = 0;
	size_t	maxidx = 0;

	for(; start <= last; ++start)
		if(abs(win[start]) > maxabs) {
			maxabs = abs(win[start]);
			maxidx = start;
		}

	return maxidx;
}

/*
 * Tries to find pulse beginning by successive trial-and-error method,
 *  continuously approximating signal envelope by reference envelope at every
 *  point using least-squares method, computing normalized error metric and
 *  selecting the best approximation.
 * Return value:
 *  <0 when unable to find successful approximation
 *  Index of first leading edge halfwave when OK
 */
#define LC_BITS_HWBETAF 16

static int lc_hwaves_match(struct lc_halfwave_buffer *hwbuf,
			   const struct lc_refedge *refedge, int32_t *acc_ret,
			   int verbose)
{
	struct lc_halfwave *hwcur, *hwend, *hwmax, *hwb;
	int32_t accmax = 0;

	assert(hwbuf->hw_count >= LC_REFEDGE_SZ);

	hwb = hwbuf->hwaves;
	hwcur =	hwbuf->hwaves;
	hwend = hwcur + hwbuf->hw_count - LC_REFEDGE_SZ + 1;

	for(; hwcur < hwend; ++hwcur) {
		/*
		 * Calculate beta coefficient for best curve fit
		 *  beta = sum(Ri * Vi) / sum(Ri ^ 2)
		 *  Vi - sampled halfwave magnitude
		 *  Ri - reference halwave magnitude
		 *
		 * We arranged magnitudes in the way that beta always should lie
		 *  in the range [0 .. 1). If this condition doesn't hold, we
		 *  just abandon calculations and maybe raise error condition..
		 *
		 * As |beta| is always < 1, it should be stored as fractional
		 *  number. We've chosen fixed point format for that, with 16
		 *  bit fraction part. So that when we multiply Vi by beta, the
		 *  result fits into 32 bit register (given Vi is up to 16 bits
		 *  wide).
		 */
		int32_t beta;

		/*
		 * When calculating beta, we assume that sum(Ri^2) and
		 *  sum(Ri * Vi) both fit into 32 bit unsigned int.
		 *  When max(Ri) = max(Vi) = 4096 (as it is with 12 bit ADC:
		 *   in the worst case, when DC removal didn't work for some
		 *   reason. Most of the time Vi <= 2048) and REFEDGE_SZ is 15,
		 *  these conditions are met.
		 */
		int32_t sum_rivi = 0;
		for(size_t i = 0; i < LC_REFEDGE_SZ; ++i)
			sum_rivi+= hwcur[i].mag * refedge->amps[i];

		beta = sum_rivi / (refedge->cvar >> LC_BITS_HWBETAF);
		if(verbose) {
			printf("Beta: %d + %d/%d\t REFAMP: %4d\t",
			       (beta >> LC_BITS_HWBETAF),
			       (beta & ((1<<LC_BITS_HWBETAF) - 1)),
			       (1<<LC_BITS_HWBETAF),
			       (abs(beta) * LC_REFEDGE_MAG) >> LC_BITS_HWBETAF);
		}

		/*
		 * Residual squares sum - cumulative fitting error measure
		 * RSS = sum((beta * Ri - Vi) ^ 2)
		 * As we can fit both sum(Ri ^ 2) and sum(Vi ^ 2) into int32,
		 *  this beast definitely fits there too..
		 */
		int32_t rss = 0;

		for(size_t i = 0; i < LC_REFEDGE_SZ; ++i) {
			int32_t b_ri;
			b_ri = beta * refedge->amps[i];
			b_ri >>= LC_BITS_HWBETAF;	/* remove fraction */
			int32_t psum;
			psum = (int32_t) b_ri - (int32_t) hwcur[i].mag;
			rss+= psum * psum;
		}
		if(verbose)
			printf("RSS: %5d\t", rss);

		/*
		 * Now we need to normalize RSS in some way, so that we are
		 * able to determine how good is approximation that we obtained
		 * at current halfwave..
		 *
		 * sqrt(RSS / REFEDGE_SZ) = sqrt(sum((b*Ri - Vi)^2) / N) -
		 *  this definitely looks like standard deviation of fitting
		 *  errors measured at each point.
		 *
		 * STD can be easily compared with amplitude (or STD..) of
		 *  either beta*Ri, or Vi series
		 *
		 * Lets make following assumption: if
		 *  beta * STD(Ri) / sqrt(rss / N) > thres
		 * then we have good approximation at current point.
		 *
		 * Simplifying above inequality, we can obtain following:
		 *  beta^2 * sum(Ri ^ 2) / RSS > thres ^ 2
		 */
		int32_t comp_lh;

		//comp_lh = (LC_HWAVE_VAL_C >> LC_BITS_HWBETAF) * beta;
		//comp_lh = (comp_lh >> LC_BITS_HWBETAF) * beta;
		comp_lh = ((beta * beta) >> LC_BITS_HWBETAF);
		comp_lh*= (refedge->cvar >> LC_BITS_HWBETAF);
		if(verbose)
			printf("LH: %8d\t", comp_lh);

		int32_t accuracy;
		if(rss)
			accuracy = comp_lh / rss;
		else
			accuracy = comp_lh;

		if(!accmax || hwcur - hwmax < 5) {
			if(accuracy > 80 && accuracy > accmax) {
				accmax = accuracy;
				hwmax = hwcur;
				if(verbose)
					printf("'@' ");
			} else {
				if(verbose)
					printf("'-' ");
			}
		} else {
			if(verbose)
				printf("'x' ");
		}

		if(verbose) {
			printf("ACC: %5d\t", accuracy);
			printf("HW %3d start %3u end %3u mag %4d\n",
			       hwcur - hwb,
			       hwcur == hwb ? 0 : ((hwcur-1)->last + 1),
			       hwcur->last,
			       hwcur->mag);
		}
		/*
		 * If we just touched last halfwave in window or last point
		 *  that is served by current window, exit..
		 */
		if(hwcur->last >= LC_STA_WINSHIFT - 1)
			break;
	}

	if(accmax) {
		if(acc_ret)
			*acc_ret = accmax;
		return hwmax - hwb;
	}

	return -1;
}

/*
 * Calculate variance (signal power estimate) in window
 *  Can either sum absolute sample values (mean deviation) or
 *   sum sample squares (standard deviation square (variance)).
 *  However, taking square root is diffucult, so we leave it as is,
 *   just dividing by sample count.
 *  If window size is < 256 points, and absolute sample value is <4096,
 *   sum(Vi^2) fits inside uint32 without overflow.
 *
 *  NOTE: If window is greater than 256 points, or samples are > 12 bits,
 *   rewrite and use something clever here!
 */
static uint32_t lc_calc_var(lc_type_comb_s *win, size_t first, size_t last)
{
	uint32_t res = 0;
	assert(last >= first);

	for(size_t i = first; i <= last; ++i) {
		int32_t sv = (*win++) >> LC_BITS_COMBF;
		//res += sv*sv;
		res += abs(sv);
	}

	res /= ((last - first) >> 3) + 1;//(LC_STA_WINSZ/16);

	return res;
}

#define LC_WVAR_COMPARABLE(var, ref) ((var) < ((ref) + (ref) / 4) &&\
				      (var) > ((ref) - (ref) / 4))

/*
 * Determines PC (phase code) in sample interval [first, last].
 * If OK returns phase code index (>=0), if phase code can't be certainly
 *  deduced, returns negative number.
 */
#define LC_PCDET_INCORRECT	-1	/* Return on no correlation with PC */
#define LC_PCDET_SPURIOUS	-2	/* Return on correlation with >1 PCs */
static int lc_determine_pc(struct lc_station *sta, size_t first, size_t last)
{
	assert(last >= first);
	/* Calculate variances of signal in all windows */
	sta->win0_var = lc_calc_var(sta->win0, first, last);
	for(int pc = 0; pc < LC_PC_CNT; ++pc)
		sta->xcorr_var[pc] = lc_calc_var(sta->xcorr[pc], first, last);

	/*
	 * First step: determine if one of "normal" XCorr windows has variance
	 *  similar to that of win0
	 */
	int pc_cand_idx = 0;
	uint32_t pc_cand_var = 0;
	for(int pc = 0; pc < LC_PC_NCNT; ++pc)
		if(LC_WVAR_COMPARABLE(sta->xcorr_var[pc], sta->win0_var)) {
			pc_cand_idx = pc;
			pc_cand_var = sta->xcorr_var[pc];
			break;
		}

	if(!pc_cand_var)
		return LC_PCDET_INCORRECT;

	/*
	 * Next step: try to find PC with variance similar to that of candidate.
	 *  If we succeed, then candidate is considered false :-(
	 * NOTE: now synthetic codes' windows are considered too!
	 */
	for(int pc = 0; pc < LC_PC_CNT; ++pc) {
		if(pc == pc_cand_idx)
			continue;
		if(sta->xcorr_var[pc] > (pc_cand_var - pc_cand_var / 4))
			return LC_PCDET_SPURIOUS;
	}

	return pc_cand_idx;
}

/*******************************************************************************
 * Complex routines that do the thing in different lc_station lifecycle states
 */

/*
 * SEEK processing:
 *  When called, station windows have already accumulated enough samples.
 *  After return from routine, station will be probably moved to the next point
 *   in seek interval.
 *  The task is to analyze present data contained in windows for presence of
 *   proper Loran-C pulse group, i.e. perform:
 *  - Amplitude selection
 *  - Phase code extraction
 *  - Pulse leading edge search
 *  - Calculation of pulse shape quality normalized metric
 *
 * Upon return from routine, caller probably will write back start-of-pulse
 *  position into list of candidates along with quality metric.
 *
 * ------
 *
 * LOCKING processing:
 *  Routine is called with station positioned in a way that previously found
 *   start-of-pulse is PRESUMED to be at the center of window.
 *  Task is to ensure that pulse group signal is stable: has the same PC, always
 *   has high accuracy of pulse leading edge, and doesn't wander back and forth
 *   in time... Also routine should determine offset of SOP point from window
 *   center and require station shift if necessary.
 *  Routine is called every N FRIs (N can be ~10)
 *
 * ------
 *
 * TRACKING processing:
 *  Routine is called with station previously positioned in a way that
 *   start-of-pulse is EXACTLY at the window center (+- shift introduced by
 *    1 FRI delay).
 *  Routine is called _every_ FRI, so performance should be considered too.
 *
 *  The task is to scan part of window near its center to determine necessary
 *   shift, find samples that surround SZC point and perform precise SZC
 *   interpolation. Also signal integrity should be estimated (amplitude,
 *   pulse edge quality, PC correlation thresholds satisfied)
 */

void lc_process_seek(struct lc_chain *chain, struct lc_station *sta)
{
	/* Determine PC */
	int pc = lc_determine_pc(sta, 0, LC_STA_WINSZ - 1);

	if(pc == LC_PCDET_SPURIOUS)
		printf("\t--SPURIOUS DETECTED--\n");
	if(pc < 0)
		return;

	if(sta->state != LC_STAST_LOCKING) {
		printf("STA %1d PC: %1u, OFFSET: %5u\n",
		       sta->idx, pc, sta->offset);

		printf("\t      WIN0 level: %7u\n", sta->win0_var);
		for(size_t pc = 0; pc < LC_PC_CNT; ++pc)
			printf("\tXCORR PC %u level: %7u\n", pc,
			       sta->xcorr_var[pc]);
	}
	/*
	 * Here we can finally assume that station windows at current offset
	 *  contain correct (from phase coding point of view) pulse group.
	 * Now pulse edge detection and its correctness estimation should happen
	 */
	int32_t edge_acc;
	int sophw;

	lc_hwaves_win_gen(sta->xcorr[pc], 0, LC_STA_WINSZ - 1, &pulwinhw);
	sophw = lc_hwaves_match(&pulwinhw, chain->refedge, &edge_acc,
				sta->state == LC_STAST_SEEK);
	if(sophw < 0)
		return;

	size_t sop_offset = lc_hwave_find_peak(sophw, sta->xcorr[pc],
					       &pulwinhw);

	/* Add new candidate to chain's list */
	if(chain->seek_cand_cnt == LC_CAND_MAXCNT)
		return;

	struct lc_cand *cand = &chain->seek_cand[chain->seek_cand_cnt++];
	cand->offset = (sta->offset + sop_offset) % chain->grin;
	cand->pc = pc;
	cand->accuracy = edge_acc;
}
