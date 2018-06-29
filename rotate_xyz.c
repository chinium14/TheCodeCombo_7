/* ======================================================================== */
/* rotate.cpp	                                                            */
/* This suboutine rotates the system so that vector between SITE1 and SITE2 */
/* is aligned to the axis passed as a parameter (0=x axis, 1=y axis, 2=z    */
/* axis).  ROTN must be defined.                                            */
/*																		    */
/* Written by Thomas Knotts 22 Oct 03                                       */
/* ======================================================================== */
#ifdef ROTN
#include "defines.h"
void pbc_all(int);


void rotate_xyz(int ibox, int axis)
{
	int k	  =	ibox;

	/* ---------------------------------------------------- */
	/* Determine the angles.                                */
	/* ---------------------------------------------------- */

	double lx	= atnopbc[k][SITE2].x -atnopbc[k][SITE1].x;
	double ly	= atnopbc[k][SITE2].y -atnopbc[k][SITE1].y;
	double lz	= atnopbc[k][SITE2].z -atnopbc[k][SITE1].z;
	
	double costheta, costheta2;
	double sintheta, sintheta2;
	double cosphi, cosphi2;
	double sinphi;
	
	if(axis == 2){
		costheta		=	lz/sqrt(lx*lx + ly*ly + lz*lz);
		costheta2	=   costheta*costheta;
		sintheta2	=	1.0- costheta2;
		sintheta		=	sqrt (sintheta2);
		cosphi		=	lx/sqrt(lx*lx + ly*ly);
		sinphi		=	ly/sqrt(lx*lx + ly*ly);
	}
	else if(axis == 1){
		costheta		=	ly/sqrt(lx*lx + ly*ly + lz*lz);
		costheta2	=   costheta*costheta;
		sintheta2	=	1.0- costheta2;
		sintheta		=	sqrt (sintheta2);
		cosphi		=	lx/sqrt(lx*lx + lz*lz);
		sinphi		=	lz/sqrt(lx*lx + lz*lz);
	}
	else{
		costheta		=	lx/sqrt(lx*lx + ly*ly + lz*lz);
		costheta2	=   costheta*costheta;
		sintheta2	=	1.0- costheta2;
		sintheta		=	sqrt (sintheta2);
		cosphi		=	lz/sqrt(lz*lz + ly*ly);
		sinphi		=	ly/sqrt(lz*lz + ly*ly);
	}


	/* ---------------------------------------------------- */
	/* Determine the rotation matrix.                       */
	/* ---------------------------------------------------- */
	double R11	=	costheta*cosphi;
	double R12	=	costheta*sinphi;
	double R13	=	-sintheta;
	double R21	=	-sinphi;
	double R22	=	cosphi;
	double R23	=	0.0;
	double R31	=	cosphi*sintheta;
	double R32	=	sinphi*sintheta;
	double R33	=	costheta;	

	/* ---------------------------------------------------- */
	/* Apply the rotation matrix to all coordinates of the  */
	/* system.                                              */
	/* ---------------------------------------------------- */
	for (int i=0; i<box[k].boxns; i++){
		double pos_x	=	R11*atnopbc[k][i].x + R12*atnopbc[k][i].y+ R13*atnopbc[k][i].z;
		double pos_y	=	R21*atnopbc[k][i].x + R22*atnopbc[k][i].y+ R23*atnopbc[k][i].z;
		double pos_z	=	R31*atnopbc[k][i].x + R32*atnopbc[k][i].y+ R33*atnopbc[k][i].z;
			
		atnopbc[k][i].x	=	pos_x;
		atnopbc[k][i].y	=	pos_y;
		atnopbc[k][i].z	=	pos_z;
		atom[k][i].x	=	atnopbc[k][i].x;
		atom[k][i].y	=	atnopbc[k][i].y;
		atom[k][i].z	=	atnopbc[k][i].z;
	}
	/* ---------------------------------------------------- */
	/* Apply pbc's.                                         */
	/* ---------------------------------------------------- */
		pbc_all(k);




}

#endif

		
