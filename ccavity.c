/* ======================================================================== */
/* ccavity.cpp                                                              */
/*            This potential function restricts the chain within a spherical*/
/* cavity. The cavity has its center at the origin of the coordinate system */
/* Based on the potential used by Klimov and Thirumalai, PNAS, 99, no 12,   */
/* p8019, 2002.																*/		
/* Written by Nitin Rathore , January 19, 2005.								*/
/* ======================================================================== */
#ifdef CAVITY
#include "defines.h"

  /* ====================================================== */
  /* This calculates a purely repulsive wall potential.     */
  /* ====================================================== */


double ccavity(int ibox){
  double ecavity= 0.0;
  int k			= ibox;

  double sig	= golik[k].req;
  double eps	= golik[k].eps * 1.0;
  double Rs		= box[k].lx*0.1;			/* The cavity size is 20% the size of box	*/
  double Rs2	= Rs * Rs;
  double factor = 0.8*PI*eps;
  double alpha	= 0.01;

  for(int b=0; b<box[k].boxns; b++){
		double dR2 = atom[k][b].x*atom[k][b].x + atom[k][b].y*atom[k][b].y + atom[k][b].z*atom[k][b].z;
		
		double dR	= sqrt(dR2);
//		if (dR+3.0>= Rs) {
//			fprintf(stdout,"Amino acids %d is within %lf Angstrom of the wall, %d\n",b,Rs-dR, n_iter);
//		}

#ifdef ACAVITY
		double dra  = Rs-dR;
		double drb	= Rs+dR;
		double rt	= Rs/dR;
		double sra  = (sig) / dra;
		double srb  = (sig) / drb;
		double sra2 = sra*sra;
		double srb2 = srb*srb;
		double sra3 = sra*sra2;
		double srb3 = srb*srb2;
		double sra4 = sra2*sra2;
		double srb4 = srb2*srb2;
		double sra6 = sra3*sra3;
		double srb6 = srb3*srb3;
		double sra9 = sra6*sra3;
		double srb9 = srb6*srb3;
		double sra10= sra6*sra4;
		double srb10= srb6*srb4;
		double force;
		ecavity		+= factor*(sra10 - srb10)*rt;
 		force		=  factor/sig*rt*((rt - 11.0)*sra10*sra -(rt + 11.0)*srb10*srb);
				
		  double fx          = force * atom[k][b].x/dR;
		  ff[k][b].x += fx;

		  double fy          = force * atom[k][b].y/dR;
		  ff[k][b].y += fy;

		  double fz          = force * atom[k][b].z/dR;
		  ff[k][b].z += fz;

#endif
#ifdef BCAVITY
		double dra  = Rs-dR;
		double force;
		ecavity		+= (alpha/Rs)*((exp(-dra))*(dR-1.0)-(dR*dR/2.0));
 		force		= -(alpha/Rs)*((exp(-dra)-1.0));
		  double fx          = force * atom[k][b].x;
		  ff[k][b].x += fx;

		  double fy          = force * atom[k][b].y;
		  ff[k][b].y += fy;

		  double fz          = force * atom[k][b].z;
		  ff[k][b].z += fz;
#endif
  }

  return (ecavity);
}
#endif
