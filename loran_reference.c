#include "loran.h"
#include "loran_reference.h"

/*
 * NOTE: Chayka edge from GOST R 53168-2008
 *	 Loran-C from USCG Loran-C 1994 Signal Specification
 */

/*
 * Chayka RSDN-3/10 reference pulse leading edge halfwave amplitudes,
 * distorted by RF circuits in receiver (taken from SPICE simulation)
 */

const struct lc_refedge lc_refedge_rsdn310 = {
	.amps = {
		(0.0878	* LC_REFEDGE_MAG),
		(-0.1385* LC_REFEDGE_MAG),
		(0.2066	* LC_REFEDGE_MAG),
		(-0.2911* LC_REFEDGE_MAG),
		(0.3873	* LC_REFEDGE_MAG),
		(-0.4894* LC_REFEDGE_MAG),
		(0.5951	* LC_REFEDGE_MAG),
		(-0.6925* LC_REFEDGE_MAG),
		(0.7817	* LC_REFEDGE_MAG),
		(-0.8603* LC_REFEDGE_MAG),
		(0.9237	* LC_REFEDGE_MAG),
		(-0.9683* LC_REFEDGE_MAG),
		(0.9894	* LC_REFEDGE_MAG),
		(-1.0	* LC_REFEDGE_MAG)
	},
	.cvar = (6.4986 * (LC_REFEDGE_MAG * LC_REFEDGE_MAG)),
	.szc_offset = 6
};


/*
 * Chayka RSDN-5BM reference pulse leading edge halfwave amplitudes,
 * distorted by RF circuits in receiver (taken from SPICE simulation)
 */

const struct lc_refedge lc_refedge_rsdn5bm = {
	.amps = {
		(0.1915	* LC_REFEDGE_MAG),
		(-0.2603* LC_REFEDGE_MAG),
		(0.3330	* LC_REFEDGE_MAG),
		(-0.4180* LC_REFEDGE_MAG),
		(0.5079	* LC_REFEDGE_MAG),
		(-0.5989* LC_REFEDGE_MAG),
		(0.6794	* LC_REFEDGE_MAG),
		(-0.7651* LC_REFEDGE_MAG),
		(0.8423	* LC_REFEDGE_MAG),
		(-0.8995* LC_REFEDGE_MAG),
		(0.9471	* LC_REFEDGE_MAG),
		(-0.9735* LC_REFEDGE_MAG),
		(0.9937	* LC_REFEDGE_MAG),
		(-1.0	* LC_REFEDGE_MAG)
	},
	.cvar = (7.3179 * (LC_REFEDGE_MAG * LC_REFEDGE_MAG)),
	.szc_offset = 5
};

/*
 * US Coast Guard reference pulse leading edge halfwave amplitudes, distorted by
 * RF circuits in receiver (taken from SPICE simulation)
 */

const struct lc_refedge lc_refedge_uscg = {
	.amps = {
		(0.1901	* LC_REFEDGE_MAG),
		(-0.2598* LC_REFEDGE_MAG),
		(0.3411	* LC_REFEDGE_MAG),
		(-0.4234* LC_REFEDGE_MAG),
		(0.5111	* LC_REFEDGE_MAG),
		(-0.5924* LC_REFEDGE_MAG),
		(0.6843	* LC_REFEDGE_MAG),
		(-0.7635* LC_REFEDGE_MAG),
		(0.8321	* LC_REFEDGE_MAG),
		(-0.8902* LC_REFEDGE_MAG),
		(0.9219	* LC_REFEDGE_MAG),
		(-0.9694* LC_REFEDGE_MAG),
		(0.9905	* LC_REFEDGE_MAG),
		(-1.0	* LC_REFEDGE_MAG)
	},
	.cvar = (7.3179 * (LC_REFEDGE_MAG * LC_REFEDGE_MAG)),
	.szc_offset = 5
};

