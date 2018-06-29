/* ======================================================================== */
/* csmd.cpp		                                                            */
/*                                                                          */
/*		This subroutine calculates energies and forces associated with the  */
/* one site being fixed and other being pulled by a harmonic spring.        */ 
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

void smd_setup(int ibox){
	int k = ibox;

	int ia = steermd[k].site1;
	int ib = steermd[k].site2;

	steermd[k].x1 = atnopbc[k][ia].x;
	steermd[k].y1 = atnopbc[k][ia].y;
	steermd[k].z1 = atnopbc[k][ia].z;

	steermd[k].x2 = atnopbc[k][ib].x;
	steermd[k].y2 = atnopbc[k][ib].y;
	steermd[k].z2 = atnopbc[k][ib].z;

	steermd[k].r2x = 0.0;
	steermd[k].r2y = 0.0;
	steermd[k].r2z = 0.0;
	steermd[k].dw1 = 0.0;
	steermd[k].dw2 = 0.0;

  /* ================================================================== */
  /* Calculate the distance (r0) between the two restrain points        */
  /* ================================================================== */

	double dx	  =	steermd[k].x2 - steermd[k].x1;
	double dy	  =	steermd[k].y2 - steermd[k].y1;
	double dz	  =	steermd[k].z2 - steermd[k].z1;

	steermd[k].r0 = sqrt(dx * dx + dy * dy + dz * dz);
	steermd[k].ex = dx/steermd[k].r0;
	steermd[k].ey = dy/steermd[k].r0;
	steermd[k].ez = dz/steermd[k].r0;

	steermd[k].count	= 0;
	steermd[k].effr1_av	= 0.0;
	steermd[k].effr2_av	= 0.0;

}

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
double csmd (int ibox)
{

  int k = ibox;
  double fx2, fy2, fz2;
  /* ================================================================== */
  /* Declare the energy accumulator                                     */
  /* ================================================================== */
  double esmd = 0.0;

  /* ================================================================== */
  /* Zero out the virial accumulators                                   */
  /* ================================================================== */
#ifdef PRESSURE
  double wxx = 0.0;
  double wyx = 0.0;
  double wyy = 0.0;
  double wzx = 0.0;
  double wzy = 0.0;
  double wzz = 0.0;
#endif
	  /* ================================================================== */
	  /* Calculate the energy and force associated with each spring.        */
	  /* ================================================================== */

	  /* ----------------------------------------- */
	  /* Assign the atoms involved in the bonding. */
	  /* ----------------------------------------- */
		int ia = steermd[k].site1;
		int ib = steermd[k].site2;

	  /* ----------------------------------------- */
	  /* Calculate the distances between the atoms */
	  /* and the spring.                           */
	  /* ----------------------------------------- */
		double dx2  = atnopbc[k][ib].x - steermd[k].x2;
		double dy2  = atnopbc[k][ib].y - steermd[k].y2;
		double dz2  = atnopbc[k][ib].z - steermd[k].z2;

		steermd[k].r2x = -dx2;
		steermd[k].r2y = -dy2;
		steermd[k].r2z = -dz2;


		double drr2 = dx2*dx2 + dy2*dy2 + dz2*dz2;

		double dr2  = sqrt(drr2);

		double fr2 = -steermd[k].k2 * (dr2) * 2.0;
    
		esmd += -((fr2 * dr2)) * 0.5;

	  /* ----------------------------------------- */
	  /* Calculate the forces.                     */
	  /* ----------------------------------------- */    
		if(dr2 == 0.0){
			fx2 = 0.0;
			fy2 = 0.0;
			fz2 = 0.0;
		}
		else{
			fx2 = fr2 * dx2 / dr2;
			fy2 = fr2 * dy2 / dr2;
			fz2 = fr2 * dz2 / dr2;
		}
		
		ff[k][ib].x += fx2;
		ff[k][ib].y += fy2;
		ff[k][ib].z += fz2;

	  /* ----------------------------------------- */
	  /* If PRESSURE is defined, accumulate the    */
	  /* components of the virial tensor.          */
	  /* ----------------------------------------- */
	#ifdef PRESSURE
		double sxx =  fx2 * dx2;
		double syx =  fy2 * dx2;
		double syy =  fy2 * dy2;
		double szx =  fz2 * dx2;
		double szy =  fz2 * dy2;
		double szz =  fz2 * dz2;
		wxx += sxx;
		wyx += syx;
		wyy += syy;
		wzx += szx;
		wzy += szy;
		wzz += szz;
	#ifdef ZHOU
		stresses[k][ia][0][0] += sxx * 0.5;
		stresses[k][ia][1][0] += syx * 0.5;
		stresses[k][ia][1][1] += syy * 0.5;
		stresses[k][ia][2][0] += szx * 0.5;
		stresses[k][ia][2][1] += szy * 0.5;
		stresses[k][ia][2][2] += szz * 0.5;

		stresses[k][ib][0][0] += sxx * 0.5;
		stresses[k][ib][1][0] += syx * 0.5;
		stresses[k][ib][1][1] += syy * 0.5;
		stresses[k][ib][2][0] += szx * 0.5;
		stresses[k][ib][2][1] += szy * 0.5;
		stresses[k][ib][2][2] += szz * 0.5;
	#endif
	#endif

	  /* ================================================================== */
	  /*                                                                    */
	  /* After the loop, if PRESSURE is defined, store the bonding          */
	  /* contribution to the virial tensor to pvir[k].bond[i].              */
	  /*                                                                    */
	  /* ================================================================== */
	#ifdef PRESSURE
	  pvir[k].bond[0] = wxx;
	  pvir[k].bond[1] = wyx;
	  pvir[k].bond[2] = wyy;
	  pvir[k].bond[3] = wzx;
	  pvir[k].bond[4] = wzy;
	  pvir[k].bond[5] = wzz;
	#endif
  
	 steermd[k].r2 = dr2;

  return (esmd);

}
#endif
#ifdef NSMD
/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

void smd_setup(int ibox){
	int k = ibox;

	int ia = steermd[k].site1;
	int ib = steermd[k].site2;

	steermd[k].x1 = atnopbc[k][ia].x;
	steermd[k].y1 = atnopbc[k][ia].y;
	steermd[k].z1 = atnopbc[k][ia].z;

	steermd[k].x2 = atnopbc[k][ib].x;
	steermd[k].y2 = atnopbc[k][ib].y;
	steermd[k].z2 = atnopbc[k][ib].z;


  /* ================================================================== */
  /* Calculate the distance (r0) between the two restrain points        */
  /* ================================================================== */

	double dx	  =	steermd[k].x2 - steermd[k].x1;
	double dy	  =	steermd[k].y2 - steermd[k].y1;
	double dz	  =	steermd[k].z2 - steermd[k].z1;

	steermd[k].r0 = sqrt(dx * dx + dy * dy + dz * dz);
	steermd[k].ex = dx/steermd[k].r0;
	steermd[k].ey = dy/steermd[k].r0;
	steermd[k].ez = dz/steermd[k].r0;
//	steermd[k].ex = 0.0;
//	steermd[k].ey = 0.0;
//	steermd[k].ez = -1.0;
	
	steermd[k].count	= 0;
	steermd[k].effr1_av	= 0.0;
	steermd[k].f2_av	= 0.0;

}

#endif

