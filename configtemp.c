#ifdef CONFIGT
/* ======================================================================== */
/* configtemp.cpp                                                           */
/*                                                                          */
/*		This subroutine calculates the configurational temperature of the . */
/* system. The numerator and the denominator are accumulated separately     */
/* in the structure config and then the ratio computed to give the inst.    */
/* configurational temperature                                              */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						k   :		The box number                          */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

void configtemp(int k){
#ifndef LJ
	config[k].hesr = config[k].hesx + config[k].hesy + config[k].hesz;
#endif // if LJ is defined hesr is already computed in cnonbond
	double force2_sum = 0.0;
	for (int i=0; i<box[k].boxns;i++){
		force2_sum += ff[k][i].x * ff[k][i].x + ff[k][i].y*ff[k][i].y + ff[k][i].z*ff[k][i].z;
	}

/* --------------------------------------------------------	*/
/* If not CTDOS then accumulate num and denominators and	*/
/* compute accumulated Configurational Temperature			*/
/* -------------------------------------------------------- */
#ifndef CTDOS
	config[k].num += force2_sum;
	config[k].den += config[k].hesr;
	config[k].T    = config[k].num/config[k].den/RG/0.001;
}
#endif

/* --------------------------------------------------------	*/
/* If CTDOS is defined then do not accumulate num and den.	*/
/* Just compute the current values and save in the variable */
/* Later assign these values to the correspoding DOS vars.	*/
/* -------------------------------------------------------- */
#ifdef CTDOS
	config[k].num  = force2_sum;
	config[k].den  = config[k].hesr;
	config[k].T    = config[k].num/config[k].den/RG/0.001;
}
#endif
#endif

