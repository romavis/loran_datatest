#include "loran.h"

//lc_type_refpulse	lc_refpulse[LC_REFPULSE_SZ] = {
//	0	,61	,0	,-108	,0	,205	,0	,-375	,
//	0	,545	,0	,-784	,0	,1015	,0	,-1261	,
//	0	,1474	,0	,-1670	,0	,1842	,0	,-2012	,
//	0	,2012	,0	,-2047
//};

/*
 * Values taken from:
 * GOST R 53168-2008: Chayka Navigation System. Transmitted signal specification
 *                    par. 3.1.1.4, page 5
 */

//lc_type_data_u		lc_refedge[LC_REFEDGE_SZ] = {
//	(0.0	* LC_REFEDGE_MAG),
//	(0.030	* LC_REFEDGE_MAG),
//	(0.053	* LC_REFEDGE_MAG),
//	(0.100	* LC_REFEDGE_MAG),
//	(0.183	* LC_REFEDGE_MAG),
//	(0.266	* LC_REFEDGE_MAG),
//	(0.383	* LC_REFEDGE_MAG),
//	(0.496	* LC_REFEDGE_MAG),
//	(0.616	* LC_REFEDGE_MAG),
//	(0.720	* LC_REFEDGE_MAG),
//	(0.816	* LC_REFEDGE_MAG),
//	(0.900	* LC_REFEDGE_MAG),
//	(0.983	* LC_REFEDGE_MAG),
//	(0.983	* LC_REFEDGE_MAG),
//	(1.0	* LC_REFEDGE_MAG)
//};

/*
 * US Coast Guard reference pulse leading edge halfwave amplitudes, distorted by
 * RF circuits in receiver (taken from SPICE simulation)
 */

lc_type_data_s			lc_refedge_uscg[LC_REFEDGE_SZ] = {
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
	(-1.0	* LC_REFEDGE_MAG),

};

