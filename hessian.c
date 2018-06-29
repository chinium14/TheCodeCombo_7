/* ======================================================================== */
/* Hessian.cpp	    														*/
/*                                                                          */
/*		This subroutine calculates the hessian matrix for the current 		*/
/* set of coordinates. It is valid only when sim.NB=1. Does not support     */
/* more than one box. Therefore k set tom "0" in the beginning.				*/
/*                                                                          */
/* ======================================================================== */
#include "defines.h"
#ifdef HESSIAN

void forces(int);

double DEL = 1E-7 ;

void hess(){
	int k=0;
  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;
    /*---------------------------------------------	*/
	/* First back up the original forces and crds	*/
	/*---------------------------------------------	*/
	for(int i =0; i< box[k].boxns; i++){
		atom_temp[k][i]		= atom[k][i];			
		ff_temp[k][i]		= ff[k][i];
	}

	for (int a =0; a< box[k].boxns; a++){

		/*---------------------------------------------	*/
		/* Now perturb the particle "a" in x direction	*/
		/*---------------------------------------------	*/
		atom[k][a].x += DEL;
		
		if (atom[k][a].x >  hx) atom[k][a].x -= lx;

		forces(k);
		for( int b=a; b< box[k].boxns; b++){
			hessian[3*a+0][3*b+0]= hessian[3*b+0][3*a+0] = (ff_temp[k][b].x - ff[k][b].x)/DEL;
			hessian[3*a+0][3*b+1]= hessian[3*b+1][3*a+0] = (ff_temp[k][b].y - ff[k][b].y)/DEL;
			hessian[3*a+0][3*b+2]= hessian[3*b+2][3*a+0] = (ff_temp[k][b].z - ff[k][b].z)/DEL;
		}
		
		/*---------------------------------------------	*/
		/* Now restore the particle "a" in x direction	*/
		/* and perturb the particle "a" in y direction	*/
		/*---------------------------------------------	*/
		atom[k][a].x	= atom_temp[k][a].x;
		atom[k][a].y   += DEL;
		if (atom[k][a].y >  hy) atom[k][a].y -= ly;
		
		forces(k);
		for( int b=a; b< box[k].boxns; b++){
			
			hessian[3*a+1][3*b+0]= hessian[3*b+0][3*a+1] = (ff_temp[k][b].x - ff[k][b].x)/DEL;
			hessian[3*a+1][3*b+1]= hessian[3*b+1][3*a+1] = (ff_temp[k][b].y - ff[k][b].y)/DEL;
			hessian[3*a+1][3*b+2]= hessian[3*b+2][3*a+1] = (ff_temp[k][b].z - ff[k][b].z)/DEL;
		}
	
		/*---------------------------------------------	*/
		/* Now restore the particle "a" in y direction	*/
		/* and perturb the particle "a" in z direction	*/
		/*---------------------------------------------	*/
		atom[k][a].y	= atom_temp[k][a].y;
		atom[k][a].z   += DEL;
		if (atom[k][a].z >  hz) atom[k][a].z -= lz;

		forces(k);
		for( int b=a; b< box[k].boxns; b++){
			
			hessian[3*a+2][3*b+0]= hessian[3*b+0][3*a+2] = (ff_temp[k][b].x - ff[k][b].x)/DEL;
			hessian[3*a+2][3*b+1]= hessian[3*b+1][3*a+2] = (ff_temp[k][b].y - ff[k][b].y)/DEL;
			hessian[3*a+2][3*b+2]= hessian[3*b+2][3*a+2] = (ff_temp[k][b].z - ff[k][b].z)/DEL;
		}
		/*---------------------------------------------	*/
		/* Now restore the particle "a" in z direction	*/
		/*---------------------------------------------	*/
		atom[k][a].z	= atom_temp[k][a].z;
	}
	for(int i =0; i< box[k].boxns; i++){
		ff[k][i]		= ff_temp[k][i];
	}
}
#endif
