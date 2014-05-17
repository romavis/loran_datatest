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
typedef	uint16_t	lc_type_refedge;

#define LC_REFEDGE_MAG	4096U

#define LC_REFEDGE_SZ	15

/*
 * Pre-calculated sum(Ri ^ 2) - can be thought of as some sort of variance..
 */
#define LC_REFEDGE_CVARFLT 5.816949
#define LC_REFEDGE_CVAR (uint32_t) (LC_REFEDGE_CVARFLT *\
	(LC_REFEDGE_MAG * LC_REFEDGE_MAG))

extern lc_type_refedge	lc_refedge[LC_REFEDGE_SZ];

#endif // LORAN_REFERENCE_H
