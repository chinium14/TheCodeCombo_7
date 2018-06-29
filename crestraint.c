/* ======================================================================== */
/* crestraint.cpp                                                           */
/*                                                                          */
/*		This subroutine calculates energies and forces associated with the  */
/* sites being restrained by springs.                                       */ 
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                                                                          */
/* ======================================================================== */
#ifdef REST
#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
double crestraint (int ibox)
{

  int k = ibox;
  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;

  /* ================================================================== */
  /* Declare the energy accumulator                                     */
  /* ================================================================== */
  double erestraint = 0.0;

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
	  /* Calculate the energy and force associated with each restrained     */
      /* single atom (rest1).                                               */
	  /* ================================================================== */
  for(int i=0; i<rest1[k].n; i++){
	  /* ----------------------------------------- */
	  /* Assign the atom involved in the restraint */
	  /* ----------------------------------------- */
		int ia = rest1[k].site[i];

	  /* ----------------------------------------- */
	  /* Calculate the distances between the atoms */
	  /* and their original position.              */
	  /* ----------------------------------------- */
		double dx = atom[k][ia].x - rest1[k].x0[i];
		double dy = atom[k][ia].y - rest1[k].y0[i];
		double dz = atom[k][ia].z - rest1[k].z0[i];

  /* ----------------------------------------- */
  /* Apply minimum image convention.           */
  /* ----------------------------------------- */
#ifndef PR_NPT
      if(dx >  hx) dx -= lx;
      else if(dx < -hx) dx += lx;
      if(dy >  hy) dy -= ly;
      else if(dy < -hy) dy += ly;
      if(dz >  hz) dz -= lz;
      else if(dz < -hz) dz += lz;
#endif
#ifdef PR_NPT
	min_image_npt_full(k, &dx, &dy, &dz);
#endif

	  /* ----------------------------------------- */
	  /* Calculate the lengths the each spring is  */
	  /* removed from its equilibrium length of 0. */
	  /* ----------------------------------------- */

	#ifndef UMBP
		double drr = sqrt(dx*dx + dy*dy + dz*dz);
		double dr  = drr - rest1[k].req[i];
	#endif

	  /* ------------------------------------------	*/
	  /* If umbrella potential is defined, the the	*/
	  /* the restraint is 1-D in Z dir for all the	*/
	  /* restraints that have site1==site2			*/
	  /* ------------------------------------------	*/
	#ifdef UMBP
		int S[16]		=	{1,18,25,35,54,64,80,96,112,136,155,179,198,222,241,265};
	    double rcm	=	0.0;
		int rcm_count	=	0;
	  		for(int ii=0; ii<16; ii++){
				double rcmi	 =	atom[k][S[ii]-1].z + atom[k][S[ii]+275].z;
				rcm			+=	rcmi;
				rcm_count ++;	
			}
		rcm		/= rcm_count;
		dz	 = dz - rcm;
		double drr = fabs(dz);
		double dr  = drr - rest1[k].req[i];
	#endif

		double fr = -rest1[k].k[i] * (dr) * 2.0;
		erestraint += -(fr * dr * 0.5);

	  /* ----------------------------------------- */
	  /* Calculate the forces.                     */
	  /* ----------------------------------------- */    
		double fx, fy, fz;
		if(dr == 0.0){
			fx = 0.0;
			fy = 0.0;
			fz = 0.0;
		}
		else{
			fx = fr * dx / dr;
			fy = fr * dy / dr;
			fz = fr * dz / dr;
		}

		ff[k][ia].x += fx;
		ff[k][ia].y += fy;
		ff[k][ia].z += fz;

	  /* ----------------------------------------- */
	  /* If PRESSURE is defined, accumulate the    */
	  /* components of the virial tensor.          */
	  /* ----------------------------------------- */
	#ifdef PRESSURE
		double sxx = fx * dx;
		double syx = fy * dx;
		double syy = fy * dy;
		double szx = fz * dx;
		double szy = fz * dy;
		double szz = fz * dz;
		wxx += sxx;
		wyx += syx;
		wyy += syy;
		wzx += szx;
		wzy += szy;
		wzz += szz;
	#ifdef ZHOU
		stresses[k][ia][0][0] += sxx;
		stresses[k][ia][1][0] += syx;
		stresses[k][ia][1][1] += syy;
		stresses[k][ia][2][0] += szx;
		stresses[k][ia][2][1] += szy;
		stresses[k][ia][2][2] += szz;

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
  
  }

	  /* ================================================================== */
	  /* Calculate the energy and force associated with each restrained     */
      /* atom pair (rest2).                                                 */
	  /* ================================================================== */
  for(int i=0; i<rest2[k].n; i++){
	  /* ----------------------------------------- */
	  /* Assign the atoms involved in the bonding. */
	  /* ----------------------------------------- */
		int ia = rest2[k].site1[i];
		int ib = rest2[k].site2[i];

	  /* ----------------------------------------- */
	  /* Calculate the distances between the atoms */
	  /* and the springs.                          */
	  /* ----------------------------------------- */
		double dx = atom[k][ib].x - atom[k][ia].x;
		double dy = atom[k][ib].y - atom[k][ia].y;
		double dz = atom[k][ib].z - atom[k][ia].z;

  /* ----------------------------------------- */
  /* Apply minimum image convention.           */
  /* ----------------------------------------- */
#ifndef PR_NPT
      if(dx >  hx) dx -= lx;
      else if(dx < -hx) dx += lx;
      if(dy >  hy) dy -= ly;
      else if(dy < -hy) dy += ly;
      if(dz >  hz) dz -= lz;
      else if(dz < -hz) dz += lz;
#endif
#ifdef PR_NPT
	min_image_npt_full(k, &dx, &dy, &dz);
#endif
	  /* ----------------------------------------- */
	  /* Calculate the lengths the each spring is  */
	  /* removed from its equilibrium length of 0. */
	  /* ----------------------------------------- */
		double drr = sqrt(dx*dx + dy*dy + dz*dz);
		double dr  = drr - rest2[k].req[i];
		double fr = -rest2[k].k[i] * (dr) * 2.0;
		erestraint += -(fr * dr * 0.5);

	  /* ----------------------------------------- */
	  /* Calculate the forces.                     */
	  /* ----------------------------------------- */    
		double fx, fy, fz;
		if(dr == 0.0){
			fx = 0.0;
			fy = 0.0;
			fz = 0.0;
		}
		else{
			fx = fr * dx / drr;
			fy = fr * dy / drr;
			fz = fr * dz / drr;
		}

		ff[k][ia].x -= fx;
		ff[k][ia].y -= fy;
		ff[k][ia].z -= fz;
		ff[k][ib].x += fx;
		ff[k][ib].y += fy;
		ff[k][ib].z += fz;

	  /* ----------------------------------------- */
	  /* If PRESSURE is defined, accumulate the    */
	  /* components of the virial tensor.          */
	  /* ----------------------------------------- */
	#ifdef PRESSURE
		double sxx = fx * dx; 
		double syx = fy * dx;
		double syy = fy * dy;
		double szx = fz * dx;
		double szy = fz * dy;
		double szz = fz * dz;
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

  }
  return (erestraint);

}
#endif
