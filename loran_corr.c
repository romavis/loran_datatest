#include "loran.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

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

static struct lc_halfwave pulwinhw[LC_WIN_HWAVESCNT];

/**/
void lc_hwaves_win_gen(lc_type_comb_s *wbuf, struct lc_halfwave *hwbuf,
		       size_t peakidx)
{
	/*
	 * Go from end to beginning: probably there is more correct 100khz
	 *  sinewave at the end of pulse than before its beginning :)
	 */

	size_t hwend, hwstart;
	struct lc_halfwave *hwcur;

//	hwend = LC_STA_WINSZ - 1;
//	hwcur =	hwbuf + LC_WIN_HWAVESCNT - 1;
//	for(size_t i = hwend; i > 0 && hwcur >= hwbuf; --i) {
//		if(LC_SGN(wbuf[i]) == LC_SGN(wbuf[i - 1]))
//			continue;

//		lc_type_comb_s hwamp = 0;
//		hwstart = i;
//		for(size_t j = hwstart; j <= hwend; ++j)
//			if(abs(wbuf[j]) > abs(hwamp))
//				hwamp = wbuf[j];
//		hwcur->mag = hwamp >> LC_BITS_COMBF;
//		hwcur->end = hwend;
//		hwcur->begin = hwstart;
//		//printf("HWAVE %3d: %u .. %u , amp %d\n",
//		//       hwcur - hwbuf, hwstart, hwend, hwamp);
//		hwcur--;
//		hwend = i - 1;
//	}
	hwstart = 0;
	hwcur = hwbuf;
	hwbuf = hwcur + LC_WIN_HWAVESCNT;

	for(size_t i = hwstart; i < (LC_STA_WINSZ - 1) && hwcur < hwbuf; ++i) {
		if(LC_SGN(wbuf[i]) == LC_SGN(wbuf[i + 1]))
			continue;

		lc_type_comb_s hwamp = 0;
		hwend = i;
		for(size_t j = hwstart; j <= hwend; ++j)
			if(abs(wbuf[j]) > abs(hwamp))
				hwamp = wbuf[j];
		hwcur->mag = hwamp >> LC_BITS_COMBF;
		hwcur->last = hwend;

		hwcur++;
		hwstart = i + 1;
	}
	/* Finish halfwaves list */
	hwcur->last = LC_STA_WINSZ - 1;
}

#define LC_BITS_HWBETAF 16

void lc_hwaves_match(struct lc_station *sta, struct lc_halfwave *hwbuf,
		     int verbose)
{
	struct lc_halfwave *hwcur, *hwbend;
	int accmax = 0, hwmax = 0;


	hwcur =	hwbuf;
	hwbend = hwbuf + LC_WIN_HWAVESCNT - LC_REFEDGE_SZ + 1;

	for(hwcur = hwbuf; hwcur < hwbend; ++hwcur) {
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
			sum_rivi+= hwcur[i].mag * lc_refedge[i];

		beta = sum_rivi / (LC_REFEDGE_CVAR >> LC_BITS_HWBETAF);
//		if(verbose)
//			printf("Beta: %d + %d/%d\t", (beta >> LC_BITS_HWBETAF),
//			       (beta & ((1<<LC_BITS_HWBETAF) - 1)),
//			       (1<<LC_BITS_HWBETAF));
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
			b_ri = beta * lc_refedge[i];
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
		comp_lh*= (LC_REFEDGE_CVAR >> LC_BITS_HWBETAF);
		if(verbose)
			printf("LH: %8d\t", comp_lh);

		int32_t accuracy;
		if(rss)
			accuracy = comp_lh / rss;
		else
			accuracy = comp_lh;

		if(!accmax || hwcur - hwbuf - hwmax < 5){
			if(accuracy > 80 && accuracy > accmax) {
				accmax = accuracy;
				hwmax = hwcur - hwbuf;
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
			       hwcur - hwbuf,
			       hwcur == hwbuf ? 0 : ((hwcur-1)->last + 1),
			       hwcur->last,
			       hwcur->mag);
		}
		/*
		 * If we just touched last halfwave in window or last point
		 *  that is served by current window, exit..
		 */
		if(hwcur[LC_REFEDGE_SZ].last == LC_STA_WINSZ - 1 ||
		   hwcur->last >= LC_STA_WINSHIFT - 1)
			break;
	}
	if(verbose) {
		if(!accmax)
			printf("\t--UNABLE TO INTERPOLATE--\n");
		else
			printf("\t--Selected interpolation at %d (ACC %d)--\n",
			       sta->offset + hwbuf[hwmax].last, accmax);
	}
	if(accmax && sta->state == LC_STAST_LOCKING && hwbuf[hwmax].last != sta->pbegin_last) {
		printf("WARNING %1d Begin of pulse ACC %5d offset %3d, prev %3d, HW %2d\n",
		       sta->idx, accmax, hwbuf[hwmax].last, sta->pbegin_last, hwmax);
		sta->pbegin_last = hwbuf[hwmax].last;
	}
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
static uint32_t lc_calc_var(lc_type_comb_s *win)
{
	uint32_t res = 0;
	for(size_t i = 0; i < LC_STA_WINSZ; ++i) {
		//res += abs((*win++) >> LC_BITS_COMBF);
		int32_t sv = (*win++) >> LC_BITS_COMBF;
		//res += sv*sv;
		res += abs(sv);
	}

	res /= (LC_STA_WINSZ/16);

	return res;
}

#define LC_WVAR_COMPARABLE(var, ref) ((var) < ((ref) + (ref) / 4) &&\
				      (var) > ((ref) - (ref) / 4))

void lc_corr_sta(struct lc_station *sta)
{
	/*
	 * Find maximum (peak) of signal in win0 and ensure that it is in the
	 *  middle of window, except if station processes beginning of seek
	 *  interval
	 */
	size_t maxidx = 0;
	int    maxval = 0;

	for(size_t i = 0; i < LC_STA_WINSZ; ++i)
		if(abs(sta->win0[i]) > maxval) {
			maxval = abs(sta->win0[i]);
			maxidx = i;
		}
//	if(maxidx < LC_WCEN_BEGIN || maxidx >= LC_WCEN_END)
//	if(0 && sta->state == LC_STAST_SEEK && maxidx < LC_WIN_RH &&
//	   sta->offset != sta->seek_cur_int->begin)
//		return;

	/* Calculate variances of signal in all windows */
	sta->win0_var = lc_calc_var(sta->win0);
	for(size_t pc = 0; pc < LC_PC_CNT; ++pc)
		sta->xcorr_var[pc] = lc_calc_var(sta->xcorr[pc]);

	/*
	 * First step: determine if one of "normal" XCorr windows has variance
	 *  similar to that of win0
	 */
	size_t pc_cand_idx = 0;
	uint32_t pc_cand_var = 0;
	for(size_t pc = 0; pc < LC_PC_NCNT; ++pc)
		if(LC_WVAR_COMPARABLE(sta->xcorr_var[pc], sta->win0_var)) {
			pc_cand_idx = pc;
			pc_cand_var = sta->xcorr_var[pc];
			break;
		}

	if(!pc_cand_var)
		return;

	if(sta->state != LC_STAST_LOCKING) {
		printf("STA %1d PC CANDIDATE: %1u, OFFSET: %5u, MAXVAL %4d, MAXIDX %4u, MAXPOS %5u\n",
		       sta->idx, pc_cand_idx, sta->offset,
		       maxval >> LC_BITS_COMBF, maxidx, maxidx + sta->offset);

		printf("\t      WIN0 level: %7u\n", sta->win0_var);
		for(size_t pc = 0; pc < LC_PC_CNT; ++pc)
			printf("\tXCORR PC %u level: %7u\n", pc,
			       sta->xcorr_var[pc]);
	}
	/*
	 * Next step: try to find PC with variance similar to that of candidate.
	 *  If we succeed, then candidate is considered false :-(
	 * NOTE: now synthetic codes' windows are considered too!
	 */
	for(size_t pc = 0; pc < LC_PC_CNT; ++pc) {
		if(pc == pc_cand_idx)
			continue;
		//if(LC_WVAR_COMPARABLE(sta->xcorr_var[pc], pc_cand_var)) {
		if(sta->xcorr_var[pc] > (pc_cand_var - pc_cand_var / 4)) {
			printf("\t--SPURIOUS DETECTED--\n");
			return;
		}
	}

	/*
	 * Here we can finally assume that station windows at current offset
	 *  contain correct (from phase coding point of view) pulse group.
	 * Now pulse edge detection and its correctness estimation should happen
	 */

	lc_hwaves_win_gen(sta->xcorr[pc_cand_idx], pulwinhw, maxidx);
	lc_hwaves_match(sta, pulwinhw, sta->state == LC_STAST_SEEK);
	if(!sta->lock_enabled) {
		sta->lock_point = sta->offset + maxidx;
		sta->lock_pc = pc_cand_idx;
		sta->lock_enabled = 1;
	}

//	lc_type_comb_s	wmax = 0;	/* Maximum value in windows */

//	for(size_t j = 0; j < LC_STA_WINSZ; ++j)
//		if(abs(sta->win0[j]) > wmax)
//			wmax = abs(sta->win0[j]);

//	wmax = wmax >> LC_BITS_COMBF;

//	if(wmax < 90)
//		return;

//	int xer = 0;

//	/* Find maximum among PC correlation windows */
//	for(size_t pc = 0; pc < LC_PC_CNT-1; ++pc) {
//		int max = abs(sta->xcorr[pc][0]);
//		size_t maxidx = 0;
//		for(size_t i = 0; i < LC_STA_WINSZ; ++i) {
//			if(abs(sta->xcorr[pc][i]) > max) {
//				maxidx = i;
//				max = abs(sta->xcorr[pc][i]);
//			}
//		}
//		/*
//		 * Ensure that our peak is within central part of window
//		 *  This reduces ambiguity in pulse detection (no more than one
//		 *   window contains given pulse in each seek interval) and also
//		 *   guarantees that there is enough space before detected peek
//		 *   and after it (LC_WCEN_PRE and LC_WCEN_POST)
//		 */
//		if(maxidx < LC_WCEN_BEGIN || maxidx >= LC_WCEN_END)
//			continue;

//		max = max >> LC_BITS_COMBF;
//		int cnorm = (max << LC_BITS_XCNORMF) / (int) wmax;

//		if(cnorm > LC_XCORR_THRES) {
//			if(sta->state == LC_STAST_SEEK) {
//				xer = 1;
//				printf("STA %1d Max for PC %u is at index %u, "
//				       "val %u, norm %u AMP %d GRI offset %u "
//				       "(sta offset %u)\n",
//				       sta->idx, pc, maxidx, max, cnorm, wmax,
//				       sta->offset + maxidx, sta->offset);
//			}

//			lc_hwaves_win_gen(sta->win0, pulwinhw, maxidx);
//			lc_hwaves_match(sta, pulwinhw, sta->state == LC_STAST_SEEK);
//			if(!sta->lock_enabled) {
//				sta->lock_point = sta->offset + maxidx;
//				sta->lock_pc = pc;
//				sta->lock_enabled = 1;
//			}
//		}
//	}

//	if(xer) {
//		printf("\t      WIN0 DEVIATION: %5d\n", lc_calc_dev_comb(sta->win0));
//		for(size_t pc = 0; pc < LC_PC_CNT; ++pc)
//			printf("\tXCORR PC %u DEVIATION: %5d\n", pc, lc_calc_dev_comb(sta->xcorr[pc]));
//	}

	/* Output windows content to files... */
//	if(sta->offset == 33356) {
//		char	fname[128] = "data_win_33356.txt";

//		FILE* dfile = fopen(fname, "wb");
//		if(!dfile) {
//			perror("fopen() while writing windows contents");
//			return;
//		}
//		for(size_t j = 0; j < LC_STA_WINSZ; ++j)
//			fprintf(dfile, "%d ", (sta->win0[j] >> LC_BITS_COMBF));
//		fputs("\n", dfile);

//		for(size_t pc = 0; pc < LC_PC_CNT; ++pc) {
//			for(size_t j = 0; j < LC_STA_WINSZ; ++j)
//				fprintf(dfile, "%d ", sta->xcorr_pc[pc][j] >> LC_BITS_COMBF);
//			fputs("\n", dfile);
//		}
//		fclose(dfile);
//	}
}
