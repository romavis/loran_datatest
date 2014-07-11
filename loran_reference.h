#ifndef LORAN_REFERENCE_H
#define LORAN_REFERENCE_H

//#include "loran.h"

#include <stdint.h>
#include <stddef.h>

/*
 * Here are reference pulse edge definitions
 *
 * Reference pulse edge is represented as an array of first pulse's N sine
 *  halfwaves amplitudes
 */
typedef	int16_t	lc_type_refedge;

#define LC_REFEDGE_MAG	4096U

#define LC_REFEDGE_SZ	14

/*
 * struct lc_refedge
 * Represents reference pulse envelope rising edge model.
 * As it models envelope, it is given in magnitudes of carrier halfwaves, not
 *  RF sample values directly.
 *
 * Fields:
 *  -amps: model halfwaves magnitudes,
 *  -cvar: sum of squares of amps
 *  -szc_offset: number of halfwave that precedes standard zero crossing (SZC),
 *   i.e. if it is 5, then SZC is on the edge of 5 and 6 halfwaves
 */

struct lc_refedge
{
	lc_type_refedge	amps[LC_REFEDGE_SZ];	/* Ri values */
	int32_t		cvar;			/* Precalculated sum(Ri^2) */
	uint8_t		szc_offset;
};

extern const struct lc_refedge lc_refedge_rsdn310;
extern const struct lc_refedge lc_refedge_rsdn5bm;
extern const struct lc_refedge lc_refedge_uscg;

#define lc_refedge_def lc_refedge_rsdn310


#endif // LORAN_REFERENCE_H
