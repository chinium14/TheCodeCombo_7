/* ======================================================================== */
/* smd_work.cpp		                                                        */
/*                                                                          */
/*		This subroutine calculates the differential work done per unit time */
/* when one site is fixed and other is pulled by a harmonic spring.			*/ 
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"
#ifdef SMD
/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

void smd_work(int ibox){
	int k = ibox;

	double r1x = steermd[k].r1x;
	double r1y = steermd[k].r1y;
	double r1z = steermd[k].r1z;
	double r2x = steermd[k].r2x;
	double r2y = steermd[k].r2y;
	double r2z = steermd[k].r2z;
	double ex  = steermd[k].ex;
	double ey  = steermd[k].ey;
	double ez  = steermd[k].ez;
	
	int ia = steermd[k].site1;
	int ib = steermd[k].site2;


	steermd[k].effr1 = ex*r1x + ey*r1y + ez*r1z;

	double r2v2	= (r2x*vv[k][ib].x + r2y*vv[k][ib].y + r2z*vv[k][ib].z)*sim.dt;
	steermd[k].dw1	 += 2.0*steermd[k].k2*r2v2;			// work done on peptide being pulled

	/* I am not sure if the following is correct. ideally the work done on the system should	*/ 
	/* use the smd input velocity as against the velocity of the last site						*/

	double effv	= (ex*vv[k][ib].x + ey*vv[k][ib].y + ez*vv[k][ib].z);
	steermd[k].effr2 = ex*r2x + ey*r2y + ez*r2z;
	steermd[k].dw2	 += 2.0*effv * steermd[k].k2* steermd[k].effr2 *sim.dt; // total work done on system 

	steermd[k].effr1_av +=	steermd[k].effr1;
	steermd[k].effr2_av +=	steermd[k].effr2;
	steermd[k].count ++;
}

void init_smd(int k){
	if(steermd[k].count==0){
		steermd[k].effr1_av =	0.0;
		steermd[k].effr2_av	=	0.0;
	}
	else{
		steermd[k].effr1_av /=	steermd[k].count;
		steermd[k].effr2_av	/=	steermd[k].count;
		steermd[k].count =1;
	}
}

#endif

#ifdef NSMD

void smd_work(int ibox){
	int k = ibox;
	steermd[k].effr1_av +=	ordparam[k].d_nc;
	steermd[k].f2_av	+=	ordparam[k].force_2;
	steermd[k].count ++;

}

void init_smd(int k){
	if(steermd[k].count==0){
		steermd[k].effr1_av =	0.0;
		steermd[k].f2_av	=	0.0;
	}
	else{
		steermd[k].effr1_av /=	steermd[k].count;
		steermd[k].f2_av	/=	steermd[k].count;
		steermd[k].count =1;
	}
}

#endif

