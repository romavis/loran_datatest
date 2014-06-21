#include "loran.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#define LC_HWAVE_SZT	5U
#define LC_HWAVE_SZ	(LC_HWAVE_SZT * LC_FS / 1000000)
#define LC_WIN_HWAVESCNT (LC_STA_WINSZ * 2 * 100000 / LC_FS)

#define LC_SGN(x) ((x) < 0 ? -1 : 1)

struct lc_halfwave
{
	/* First and last sample indexes (in some window) */
	uint16_t	begin;
	uint16_t	end;
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

	hwend = LC_STA_WINSZ - 1;
	hwcur =	hwbuf + LC_WIN_HWAVESCNT - 1;
	for(size_t i = hwend; i > 0 && hwcur >= hwbuf; --i) {
		if(LC_SGN(wbuf[i]) == LC_SGN(wbuf[i - 1]))
			continue;

		lc_type_comb_s hwamp = 0;
		hwstart = i;
		for(size_t j = hwstart; j <= hwend; ++j)
			if(abs(wbuf[j]) > abs(hwamp))
				hwamp = wbuf[j];
		hwcur->mag = hwamp >> LC_BITS_COMBF;
		hwcur->end = hwend;
		hwcur->begin = hwstart;
		//printf("HWAVE %3d: %u .. %u , amp %d\n",
		//       hwcur - hwbuf, hwstart, hwend, hwamp);
		hwcur--;
		hwend = i - 1;
	}
}

#define LC_BITS_HWBETAF 16

void lc_hwaves_match(struct lc_halfwave *hwbuf)
{
	struct lc_halfwave *hwcur;

	hwcur =	hwbuf + LC_WIN_HWAVESCNT - LC_REFEDGE_SZ;

	for(; hwcur >= hwbuf; --hwcur) {
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
		printf("Beta: %d.%d\t", (beta >> LC_BITS_HWBETAF),
		       (beta & ((1<<LC_BITS_HWBETAF) - 1)));

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
		printf("RSS: %d\t", rss);

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
		 *  beta^2 * sum(Ri ^ 2) / (thres ^ 2) > RSS
		 */
		int32_t comp_lh;

		//comp_lh = (LC_HWAVE_VAL_C >> LC_BITS_HWBETAF) * beta;
		//comp_lh = (comp_lh >> LC_BITS_HWBETAF) * beta;
		comp_lh = ((beta * beta) >> LC_BITS_HWBETAF);
		comp_lh*= (LC_REFEDGE_CVAR >> LC_BITS_HWBETAF);
		printf("LH: %d\t", comp_lh);

		int32_t accuracy;
		if(rss)
			accuracy = comp_lh / rss;
		else
			accuracy = comp_lh;

		printf("ACC: %d\t", accuracy);

		printf("HW %d start %u end %u mag %d\n",
		       hwcur - hwbuf, hwcur->begin, hwcur->end, hwcur->mag);
	}
}

void lc_corr_sta(struct lc_station *sta)
{
	lc_type_comb_s	wmax = 0;	/* Maximum value in windows */

	for(size_t j = 0; j < LC_STA_WINSZ; ++j)
		if(abs(sta->win0[j]) > wmax)
			wmax = abs(sta->win0[j]);

	wmax = wmax >> LC_BITS_COMBF;

	if(wmax < 90)
		return;

	/* Find maximum among PC correlation windows */
	for(size_t pc = 0; pc < LC_PC_CNT; ++pc) {
		int max = abs(sta->xcorr_pc[pc][0]);
		size_t maxidx = 0;
		for(size_t i = 0; i < LC_STA_WINSZ; ++i) {
			if(abs(sta->xcorr_pc[pc][i]) > max) {
				maxidx = i;
				max = abs(sta->xcorr_pc[pc][i]);
			}
		}
		/*
		 * Ensure that our peak is within central part of window
		 *  This reduces ambiguity in pulse detection (no more than one
		 *   window contains given pulse in each seek interval) and also
		 *   guarantees that there is enough space before detected peek
		 *   and after it (LC_WCEN_PRE and LC_WCEN_POST)
		 */
		if(maxidx < LC_WCEN_BEGIN || maxidx >= LC_WCEN_END)
			continue;

		max = max >> LC_BITS_COMBF;
		int cnorm = (max << LC_BITS_XCNORMF) / (int) wmax;

		if(cnorm > LC_XCORR_THRES) {
			printf("Max for PC %u is at index %u, val %u, norm %u"
			       " AMP %d GRI offset %u (sta offset %u)\n",
			       pc, maxidx, max, cnorm, wmax,
			       sta->offset + maxidx, sta->offset);
			lc_hwaves_win_gen(sta->win0, pulwinhw, maxidx);
			lc_hwaves_match(pulwinhw);
			if(!sta->lock_enabled) {
				sta->lock_point = sta->offset + maxidx;
				sta->lock_pc = pc;
				sta->lock_enabled = 1;
			}
		}
	}

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
