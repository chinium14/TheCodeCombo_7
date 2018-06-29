/* ======================================================================== */
/* pbc_all.cpp                                                              */
/*                                                                          */
/*		This subroutine applys periodic boundary conditions to all the      */
/* site/atom coordinates.  It is before forces are calculated.              */
/*                                                                          */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"
#ifdef STATUS
void curr_status		  (int,int);
#endif

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void pbc_all (int ibox)
{
  int k = ibox;
#ifdef STATUS
	  curr_status(k,1);
#endif
  /* ================================================================== */
  /*                                                                    */
  /* Box information                                                    */
  /*                                                                    */
  /* ================================================================== */
  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;
//  fprintf(stdout, "pbc_all: line 40");
  /* ================================================================== */
  /*                                                                    */
  /* Determine if each atom's coordinates lie inside the box. If they   */
  /* do not, add or subtract a boxlength to their coordinates.          */
  /*                                                                    */
  /* ================================================================== */
  for(int i=0; i<box[k].boxns; i++) {

    double dx = atom[k][i].x;
    double dy = atom[k][i].y;
    double dz = atom[k][i].z;
    
    double pbcx = 0.0;
    double pbcy = 0.0;
    double pbcz = 0.0;


	while (dx<-hx || dx>hx){
		if (dx< -hx) pbcx = +lx;
		if (dx> hx)  pbcx = -lx;
		dx += pbcx;
		#ifdef STATUS
			  curr_status(k,2);
		#endif
	}
	while (dy<-hy || dy>hy){
		if (dy< -hy) pbcy = +ly;
		if (dy> hy)  pbcy = -ly;
		dy += pbcy;
		#ifdef STATUS
			  curr_status(k,3);
		#endif
	}
	while (dz<-hz || dz>hz){
		if (dz< -hz) pbcz = +lz;
		if (dz> hz)  pbcz = -lz;
		dz += pbcz;
		#ifdef STATUS
			curr_status(k,4);
		#endif
	}
    atom[k][i].x = dx;
    atom[k][i].y = dy;
    atom[k][i].z = dz;

  }
//  fprintf(stdout, "pbc_all: line 87");
}  

