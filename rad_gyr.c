/* ======================================================================== */
/* gyration.cpp		                                                        */
/*                                                                          */
/* This subroutine calculates the radius of gyration based on the specific  */
/* sites mentioned in the gyration.inp file                                 */
/* ======================================================================== */

#include "defines.h"

void rad_gyr(int k){

	double r_gy = 0.0;
	if(gyration[k].n != 0){	
		
		double xc = 0.0;
		double yc = 0.0;
		double zc = 0.0;
		double mc = 0.0;

		for(int i=0; i<gyration[k].n; i++) {
			int a = gyration[k].sites[i];
			xc += atnopbc[k][a].x * pott[k][a].mas;
			yc += atnopbc[k][a].y * pott[k][a].mas;
			zc += atnopbc[k][a].z * pott[k][a].mas;
			mc   += pott[k][a].mas;
		}
		xc /= mc;
		yc /= mc;
		zc /= mc;
		

		for(int i=0; i<gyration[k].n; i++){
			int a = gyration[k].sites[i];
			double mass = pott[k][a].mas;
			double dmx  = atnopbc[k][a].x - xc;
			double dmy  = atnopbc[k][a].y - yc;
			double dmz  = atnopbc[k][a].z - zc;

			r_gy += (dmx*dmx + dmy*dmy + dmz*dmz)*mass;
		}
		r_gy = r_gy/mc;
		ordparam[k].gyr = sqrt(r_gy);
	}
}
