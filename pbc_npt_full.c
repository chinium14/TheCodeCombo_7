/* ======================================================================== */
/* pbc_npt_full.cpp                                                         */
/*                                                                          */
/*		This subroutine applys periodic boundary conditions to all the      */
/* site/atom coordinates when using Parinello Rahman NPT dynamics.          */
/*                                                                          */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                                                                          */
/* ======================================================================== */
#ifdef PR_NPT
#include "defines.h"


/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void pbc_npt_full(int ibox)
{
	int k = ibox;
	double sx, sy, sz;

	for(int i=0; i<box[k].boxns; i++){
  /* ================================================================== */
  /*                                                                    */
  /* Obtain the scaled coordinates for each particle.                   */
  /*                                                                    */
  /* ================================================================== */
		sx = axesi[k][0] * atom[k][i].x + axesi[k][1] * atom[k][i].y + axesi[k][2] * atom[k][i].z;
		sy = axesi[k][3] * atom[k][i].x + axesi[k][4] * atom[k][i].y + axesi[k][5] * atom[k][i].z;
		sz = axesi[k][6] * atom[k][i].x + axesi[k][7] * atom[k][i].y + axesi[k][8] * atom[k][i].z;

  /* ================================================================== */
  /*                                                                    */
  /* Apply pbc for length between -0.5 and 0.5.                         */
  /*                                                                    */
  /* ================================================================== */
/*		while (sx< 0.0 || sx>1.0){
			if (sx < 0.0) sx = sx + 1.0;
			if (sx >  1.0) sx = sx - 1.0;
		}
		while (sy< 0.0 || sy>1.0){
			if (sy < 0.0) sy = sy + 1.0;
			if (sy >  1.0) sy = sy - 1.0;
		}
		while (sz<0.0 || sz>1.0){
			if (sz < 0.0) sz = sz + 1.0;
			if (sz >  1.0) sz = sz - 1.0;
		}
*/
		while (sx< -.5 || sx>.5){
			if (sx < -.5) sx = sx + 1.0;
			if (sx >  .5) sx = sx - 1.0;
		}
		while (sy< -.5 || sy>.5){
			if (sy < -.5) sy = sy + 1.0;
			if (sy >  .5) sy = sy - 1.0;
		}
		while (sz<-.5 || sz>.5){
			if (sz < -.5) sz = sz + 1.0;
			if (sz >  .5) sz = sz - 1.0;
		}

  /* ================================================================== */
  /*                                                                    */
  /* Obtain the real coordinates.                                       */
  /*                                                                    */
  /* ================================================================== */
		atom[k][i].x = axes[k][0] * sx + axes[k][1] * sy + axes[k][2] * sz;
		atom[k][i].y = axes[k][3] * sx + axes[k][4] * sy + axes[k][5] * sz;
		atom[k][i].z = axes[k][6] * sx + axes[k][7] * sy + axes[k][8] * sz;
	}



}
#endif
