#ifndef LORAN_H
#define LORAN_H

#include <stdint.h>
#include <stddef.h>

#include "loran_reference.h"

#define LC_MAX(x,y)	((x) > (y) ? (x) : (y))
#define LC_MIN(x,y)	((x) < (y) ? (x) : (y))

/*******************************************************************************
 * Basic info about Loran-C signal processing stages
 *
 * FRI - frame repetition interval
 * GRI - group repetition interval
 *
 * FRIn, GRIn - number of samples in one FRI/GRI
 *
 * Loran-C signal fully repeats every FRI - that is, each FRIn samples
 * Therefore, all time-related information can be represented in terms of offset
 *  inside FRI. Moreover, because of the fact that FRI consists of 2 GRIs, and
 *  only difference between them are used phase codes (all timings are the same)
 *  we can express all timings in terms of offsets inside GRI.
 *
 * The most natural way to filter Loran-C signal from noise and interfering
 *  signals is to integrate it over FRI - that is, use integrating (peaking)
 *  comb filters.
 * However, we are limited in available memory, and integrating ALL the signal
 *  over FRI is very memory-consuming activity (we have to store complete FRIn
 *  samples for comb filter delay line). Therefore, we are using integration in
 *  windows. We store only that parts of FRI (and integrate signal over them)
 *  that we are interested in.
 *  These parts are Loran-C pulse leading edges (even not complete pulse
 *  groups!).
 * How many pulses we have on FRI? For example, lets consider a chain that has
 *  5 stations: 1 master and 4 slaves. Each station emits pulse group each GRI -
 *  that is 8 pulses. We need to store 2 GRIs (1 FRI), because each GRI has
 *  its own pulse group phase code. That is, 16 pulses per station.
 *  All in all, we need to store 5*16 = 80 pulse windows.
 *
 * Also, we need some mechanism that will distribute signal samples among these
 *  windows - we need offsets (coordinates!) of these windows inside FRI and
 *  their lengths. 16 windows that correspond to each station have constant
 *  offsets relative to first station's pulse in FRI - 1 ms, 2 ms, ... GRI + 0,
 *  GRI + 1ms, .. and so on. Therefore, the key moment is offset of first
 *  station's window inside GRI - that is called "station offset".
 *  Naturally speaking, this is offset of station pulse group inside GRI.
 *
 * Therefore, main tasks of signal processing are:
 *  1) receive a buffer with signal samples
 *  2) distribute samples to correct windows (some samples are discarded)
 *  3) for each window: perform DC component removal and integrate with comb
 *      filter
 *  4) thats all.
 *
 * All these things should be performed fast and probably in real-time.
 * What to do with integrated and filtered parts of signal we have in windows?
 * Naturally, try to find Loran pulses there!
 *
 * This task (pulse finding and phase code decoding) is done by finding cross-
 *  correlation between signal in windows and reference pulse edge that is
 *  hard-coded in table.
 *
 * TODO!!
 *
 * NOTES on data types:
 *  -ADC samples come in biased, with unsigned type lc_type_sample. There are
 *   only lower LC_BITS_D bits meaningful, all other bits are expected to be set
 *   to zeros. For example, with LC_BITS_D = 12, sample values can range from
 *   0 to 4095 inclusive.
 *  -DC removal filter removes DC bias from samples, so that in _perfect_ case
 *   their values would lie in range [-2048, 2048]. However, we shouldn't expect
 *   that to be true all the time, and SHOULDN'T rely on that - there can be
 *   sudden glitches on ADC input caused by some unknown reasons that would
 *   cause values on DC filter output to jump up to (approximately) 4095 in
 *   magnitude. Therefore, DC filter output is guaranteed to be in range
 *    [-4095, 4095], NOT [-2048, 2048]!
 *  -
 */

/*******************************************************************************
 * Data types:
 *
 * lc_type_sample - ADC sample type, usually unsigned integer :-)
 * lc_type_data_s - data type on output from DC blocker (usually signed counter-
 *	part of lc_type_sample)
 * lc_type_comb_s - data type inside comb filter buffers, usually - fixed point
 *	with LC_BITS_COMBF fractional part length (signed, of course)
 *	We need fixed point here because just integral type introduces high
 *	quantization error when using high P values
 * lc_type_var - data type for storing variance estimate
 *	should be 2x larger than sample type
 *	It is the type used for cross-correlation window data
 */

typedef uint16_t	lc_type_sample;
typedef int16_t		lc_type_data_s;
typedef int32_t		lc_type_comb_s;

typedef int32_t		lc_type_var;
/**/

typedef uint16_t	lc_type_data_u;	/* Filtered data has this type */
typedef int16_t		lc_type_data_sf;/* Fixed-point data type
					   Also see LC_BITS_FRACT */
typedef	uint32_t	lc_type_var_u;	/* Signed variance data type */
typedef int32_t		lc_type_var_s;	/* Unsigned variance data type */


/*
 * Significant bits in data sample
 *  This value propagates thoughout entire filtering "pipeline"..
 */

#define LC_BITS_D	12

/*******************************************************************************
 * Cross-correlation related
 *
 * Cross-correlation buffer size (between signal window and reference edge)
 */
#define LC_XCORR_SZ	((size_t) (LC_STA_WINSZ - LC_REFPULSE_SZ + 1))

/*
 * Some bit counting magic
 *
 * We calculate variance (and cross-correlation) as follows:
 *	X(f,g,N) = 1/N * sum_i(f[i] * g[i]), where i = [0..N-1]
 * If we sum all partial products first and then divide by N, we are running
 * into overflow. If we first divide each partial product by N and then sum all
 * them up, we are running into underflow :-(
 *
 * So, we are going to sum each M products and divide these partial sums by N
 * What is M? M is the number obtained from equation:
 *	M * 2^(2*D) < 2^L
 *	D is count of significant bits in sample type:
 *		for whole-ranged u16 it is 16, for s16 it is 15
 *	L is size of lc_type_var: for u32 it is 32, for s32 it is 31, so
 *		2^L - 1 is maximum representable by lc_type_var value
 *
 * If we are using 12-bit ADC sample with DC removed, then D is 12, L is 31,
 *	partial product can have value up to 2^24. lc_type_var can store value
 *	up to 2^31 - 1. Then we can sum up to 2^(31 - 24) - 1 = 127 partial sums
 *	and guarantee to have no overflow, supposing that sample value is in
 *	range -4096 to 4096 (-2^12 to +2^12)
 */

#define LC_BITS_STYPE(type) (sizeof(type) * 8U - 1)
#define LC_BITS_UTYPE(type) (sizeof(type) * 8U)

#define LC_BITS_L	LC_BITS_STYPE(lc_type_var)

#define LC_CORR_M	((1U << (LC_BITS_L - 2 * LC_BITS_D)) - 1)

/*
 * Normalized cross-correlation value square maximum should be represented
 * using fixed point number, because its magnitude is usually low (0..4)
 * Length of fractional part is specified by LC_BITS_XCNORMF
 */
#define LC_BITS_XCNORMF	8

/*
 * Threshold (normalized xcorr square value) for station to become
 *  candidate for locking procedure
 *
 * NOTE: 1. This is default value
 *       2. This is fixed point number value, so take fraction part into account
 */

#define LC_XCORR_THRES	((lc_type_var) (0.7 * (1<<LC_BITS_XCNORMF)))

/*
 * Threshold for estimating whether halfwaves amplitudes fitting is 'good' or
 *  'bad'
 */
#define LC_HWAVE_THRES	15

/*******************************************************************************
 * Comb filter related
 *
 *  Equation: y_i = k * x_i + (1-k) * y_(i-d), where
 *	k(0..1) specifies frequency response sharpness (lower k means more sharp
 *       response)
 *	d is filter period in samples (in Loran-C case this is FRI duration)
 *
 * To greatly simplify calculations we use k = 2^(-P)
 *  P = 0 means no filtering
 *  Minimal reasonable P is 2 (k = 0.25)
 *  Maximal reasonable P is 10 (k ~= 0.001)
 */

#define LC_COMB_KP_MIN	2
#define LC_COMB_KP_MAX	10

/*
 * Comb filter step response timings
 *
 * Calculated as N_periods > ln(1-gamma) / ln(1-k)
 *	gamma is required relative signal level output by filter (i.e. 0.6)
 *	k is filter sharpness coefficient (i.e. 0.05)
 *	N_periods is count of repetition intervals passed through filter
 *		that is required to reach specified signal level at output
 * These values are pre-calculated (gamma = 0.7) and stored in table
 *  Access: lc_comb_delay[p - 1]
 */

extern size_t	lc_comb_delay[LC_COMB_KP_MAX+1];

/*
 * Comb filter delay line (window data) actually contains fixed-point data type,
 * not integral. This is needed to overcome underflow caused by truncating P
 * least significant bits of samples when calculating (sample >> P)
 *
 * We set it to P_MAX, because in this case there can't be any underflow..
 *
 * Note! LC_BITS_COMBF+LC_BITS_D should fit inside lc_type_comb_s type!
 */
#define LC_BITS_COMBF	LC_COMB_KP_MAX

/*******************************************************************************
 * DC removing filter related
 *
 * We implement DC removing filter as 1-zero 1-pole notch with f0 = 0
 * y_i = b * x_i - b * x_i-1 + (1-k) * y_i-1
 *	k adjusts notch sharpness, b is calculated from k to normalize response
 *	b = 1 - k/2
 *
 * Again, we use k = 2^-P to simplify computations
 * Then, b = 1 - 2^(-P-1)
 *
 * P is fixed in this case
 */

#define LC_DCREM_KP	3

/*******************************************************************************
 * Essential values that have global meaning
 */

/* Loran-C data sampling frequency */
#define LC_FS		400000U

/* Group repetitions in one code repetition interval (FRI) */
#define LC_GRCNT	2U
/* Window count in one GRI for 1 station */
#define LC_STA_GRIWCNT	8U
/* Overall station window count */
#define LC_STA_WINCNT	(LC_GRCNT * LC_STA_GRIWCNT)
/* Interval between two consequtive windows */
#define LC_STA_WININTT	1000U					/* in uS */
#define LC_STA_WININT	(LC_FS * LC_STA_WININTT / 1000000)	/* in samples */
/* Single window size */
#define LC_STA_WINSZT	500U					/* in uS */
#define LC_STA_WINSZ	(LC_FS * LC_STA_WINSZT / 1000000)	/* in samples */

/*
 * Approximate pulse leading edge length
 * Used to ensure that leading edge fits in window in detection phase
 */
#define LC_PLE_SZT	100U					/* in uS */
#define LC_PLE_SZ	(LC_FS * LC_PLE_SZT / 1000000)		/* in samples */

/* Window (station) shift step used in seek process */
#define LC_STA_WINSHIFT	(LC_STA_WINSZ - LC_PLE_SZ)

/* Station single GRI windows "pack" duration (in samples) */
#define LC_STA_SPAN	(LC_STA_WININT * (LC_STA_GRIWCNT - 1) + LC_STA_WINSZ)
/* 9.9ms interval (in samples) used to apply timing constraints */
/*
 * LOWERED: Chayka GRI 8000 doesn't satisfy USCG 9.9ms constraint!
 *  Instead, TD difference is not smaller than ~8.5 ms in RSDN-3/10..
 */
//#define LC_INT_99MS	(LC_FS * 99 / 10000)
#define LC_INT_99MS	(LC_FS * 85 / 10000)

/*******************************************************************************
 * Phase codes
 *
 * There are 4 8-bit codes defined in Loran-C standard: Master A,B and Slave A,B
 * These codes are used for corresponding stations' (master/slave) pulse groups
 *  switching between A and B codes in every next GRI. FRI contains two
 *  consequtive GRIs, during which, obviously, every station have used both A
 *  and B codes. We construct 4 16-bit codes that represent each of possible
 *  cases of phase code sequences for every station on FRI:
 *  0 - Master A then Master B
 *  1 - Master B then Master A
 *  2 - Slave A then Slave B
 *  3 - Slave B then Slave A
 *
 * One another code (4) was handcrafted to help solving pulse group search
 *  ambiguity in one very special but yet too common case.
 */

#define LC_PC_MASTER_A	0
#define LC_PC_MASTER_B	1
#define LC_PC_SLAVE_A	2
#define LC_PC_SLAVE_B	3
#define LC_PC_SYNTH	4

#define LC_PC_NCNT	4	/* "Normal" PC count */
#define LC_PC_CNT	5	/* Overall PC count, including synthetic ones */

/*
 * Phase codes themself
 * Values meaning: >0 - coefficient '+1', <0 - '-1', 0 - '0'
 */
extern int		lc_pc[LC_PC_CNT][LC_GRCNT * LC_STA_GRIWCNT];

/*******************************************************************************
 * Data structures
 */

/*
 * Maximal station count in one chain. Refers to count of lc_station structures,
 *  intervals, seek intervals for each station actually allocated.
 * Allocation is completely static and compile-time here
 */
#define LC_STA_MAXCNT	8

/*
 * Maximal candidates (Loran pulse groups that are found during seek process)
 *  count
 */
#define LC_CAND_MAXCNT	(LC_STA_MAXCNT * 2)

/* Represents interval inside GRI: [begin, end). NOTE: end >= begin (always) */
struct lc_int {
	size_t		begin;
	size_t		end;
};

/* Station states */
#define LC_STAST_IDLE	0
#define LC_STAST_SEEK	1
#define LC_STAST_BUSY	2
#define LC_STAST_LOCKING 3

#define LC_STA_IS_BUSY(x) ((x)->state == LC_STAST_BUSY || (x)->state == LC_STAST_LOCKING)
#define LC_STA_IS_IDLE(x) ((x)->state == LC_STAST_IDLE)

/*
 * Spanzone: interval on FRI (strictly in FRI coordinates: 0..FRIn-1)
 *  that station "occupies".
 * At max there are 3 spanzones, minimal count is 2.
 * Spanzones follow ordering:
 *  zone 0 - first zone (begin = offset)
 *  zone 1 - next zone (begin = GRIn+offset)
 *  ...
 *  zone N - exists when zone 1 "overlaps" over end of FRI and represents
 *   overlapping part of zone 1 (begin = 0)
 */

/* Struct represents spanzone - interval of FRI that station occupies */
struct lc_spanzone {
	/* In FRI "coordinates" (i.e. 0..FRIn-1 for begin, ..FRIn for end) */
	size_t		begin;
	size_t		end;

	size_t		win_idx;
	size_t		win_offset;
};

/* Struct represents single candidate (position in GRI) for locking */
struct lc_cand {
	size_t		offset;		/* Start of pulse offset in GRI */
	uint32_t	accuracy;	/* Pulse leading edge accuracy */
	int		pc;		/* Phase code index */
	_Bool		passed;		/* Selection process pass flag */
};

/* Struct represents single station */
struct lc_station {
	/*** GENERIC ***/
	/* Station index */
	size_t		idx;
	/* Offset inside GRI: from 0 to GRIn-1 */
	size_t		offset;

	struct lc_spanzone spans[LC_GRCNT + 1];
	uint8_t		spans_cnt;
	/* Next spanzone to be processed by filling routine */
	uint8_t		span_next;

	/*
	 * Sets to true when spanzone 0 begins being filled - passed FRIs are
	 *  counted from this moment
	 */
	_Bool		span0_open;

	/* Number of FRIs passed "through" station (internal use) */
	size_t		fris_passed;
	/* Number of last passed-through GRI (internal use) */
	size_t		last_gri_idx;
	/* Period (FRIs count) to call station service routine with */
	size_t		ssr_call_period;


	/* P to use for comb filter */
	int		comb_kp;

	/* DC removal filter delay line */
	lc_type_sample	dcrem_buf;

	/*
	 * Window that stores comb-filtered representation of first pulse (first
	 *  bit of PC)
	 */
	lc_type_comb_s	win0[LC_STA_WINSZ];

	/*
	 * PC correlator temporary result (over one FRI) accumulation windows
	 * These accumulate sums of samples directly from DC removal filter,
	 *  not yet filtered with comb filter..
	 */
	lc_type_data_s	xcorr_tmp[LC_PC_CNT][LC_STA_WINSZ];

	/*
	 * Cross-correlation results - for each phase code
	 * These are result of passing xcorr_tmp windows through a comb filter
	 * during several FRIs
	 */
	lc_type_comb_s	xcorr[LC_PC_CNT][LC_STA_WINSZ];

	/*
	 * Win0 and PC cross-correlation windows' signal variance.
	 *  Can be thought of as a measure of in-window signal power
	 */
	uint32_t	win0_var;
	uint32_t	xcorr_var[LC_PC_CNT];

	/* Station state */
	int		state;

	/*** SEEK related ***/
	/* Intervals to seek through assigned to this station - sorted (ASC) */
	size_t		seek_ints_cnt;
	struct lc_int	seek_ints[LC_STA_MAXCNT + 1];
	/* Current seeking interval */
	struct lc_int	*seek_cur_int;
	/* Flag - set when all seek intervals are processed */
	_Bool		seek_complete;

	/*** LOCKING related ***/
	size_t		lock_point;
	size_t		lock_pc;
	_Bool		lock_enabled;
	int32_t		pbegin_last;
};

/* Struct represents single Loran-C chain */
struct lc_chain {
	uint16_t	gri;
	size_t		station_cnt;

	uint32_t	frin;
	uint32_t	grin;

	size_t		next_sample;	/* from 0 to frin-1 */

	size_t		free_int_cnt;
	struct lc_int	free_int[LC_STA_MAXCNT + 1];

	const struct lc_refedge *refedge;

	/* SEEK related */
	size_t		seek_sta_cnt;

	size_t		seek_cand_cnt;
	struct lc_cand	seek_cand[LC_CAND_MAXCNT];

	/* Chain stations themself */
	struct lc_station sta[LC_STA_MAXCNT];
};

/*
 * Process new samples
 */
void lc_init(struct lc_chain *chain);
void lc_new_samples(struct lc_chain *chain, lc_type_sample *buf, size_t count);

/*
 * Processing functions for different station states
 */
void lc_process_seek(struct lc_chain *chain, struct lc_station *sta);

#endif // LORAN_H
