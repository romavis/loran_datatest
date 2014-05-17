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

lc_type_data_u		lc_refedge[LC_REFEDGE_SZ] = {
	(0.0	* LC_REFEDGE_MAG),
	(0.030	* LC_REFEDGE_MAG),
	(0.053	* LC_REFEDGE_MAG),
	(0.100	* LC_REFEDGE_MAG),
	(0.183	* LC_REFEDGE_MAG),
	(0.266	* LC_REFEDGE_MAG),
	(0.383	* LC_REFEDGE_MAG),
	(0.496	* LC_REFEDGE_MAG),
	(0.616	* LC_REFEDGE_MAG),
	(0.720	* LC_REFEDGE_MAG),
	(0.816	* LC_REFEDGE_MAG),
	(0.900	* LC_REFEDGE_MAG),
	(0.983	* LC_REFEDGE_MAG),
	(0.983	* LC_REFEDGE_MAG),
	(1.0	* LC_REFEDGE_MAG)
};

