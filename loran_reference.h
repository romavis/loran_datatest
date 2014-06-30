#ifndef LORAN_REFERENCE_H
#define LORAN_REFERENCE_H

//#include "loran.h"

#include <stdint.h>
#include <stddef.h>

/*
 * These are reference pulse edge definitions
 *
 * Reference pulse edge is represented as an array of first pulse's N sine
 *  halfwaves amplitudes
 */
typedef	int16_t	lc_type_refedge;

#define LC_REFEDGE_MAG	4096U

#define LC_REFEDGE_SZ	14

/*
 * Pre-calculated sum(Ri ^ 2) - can be thought of as some sort of variance..
 */
//#define LC_REFEDGE_CVARFLT 7.3179	/* USCG */
#define LC_REFEDGE_CVARFLT 6.4986	/* RSDN-3/10 */
#define LC_REFEDGE_CVAR (int32_t) (LC_REFEDGE_CVARFLT *\
	(LC_REFEDGE_MAG * LC_REFEDGE_MAG))

extern lc_type_refedge	lc_refedge_uscg[LC_REFEDGE_SZ];
extern lc_type_refedge	lc_refedge_rsdn310[LC_REFEDGE_SZ];
extern lc_type_refedge	lc_refedge_rsdn5bm[LC_REFEDGE_SZ];

#define lc_refedge lc_refedge_rsdn310


#endif // LORAN_REFERENCE_H
