/* ======================================================================== */
/* min_image_npt_full.cpp                                                   */
/*                                                                          */
/*		This subroutine calculates the minimum image between two particles. */
/*                                                                          */
/* Passed Parameters:                                                       */
/*					  ibox:		box number
/*						dx:		x distance between two particles            */
/*                      dy:		y distance between two particles			*/
/*						dz:		z distance between two particles            */
/*                                                                          */
/* ======================================================================== */
#ifdef PR_NPT
#include "defines.h"


/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void min_image_npt_full(int ibox, double *dx, double *dy, double *dz)
{
	int k=ibox;	
	double sx, sy, sz;

  /* ================================================================== */
  /*                                                                    */
  /* Obtain the scaled coordinates for each particle.                   */
  /*                                                                    */
  /* ================================================================== */
	sx = axesi[k][0] * *dx + axesi[k][1] * *dy + axesi[k][2] * *dz;
	sy = axesi[k][3] * *dx + axesi[k][4] * *dy + axesi[k][5] * *dz;
	sz = axesi[k][6] * *dx + axesi[k][7] * *dy + axesi[k][8] * *dz;

  /* ================================================================== */
  /*                                                                    */
  /* Apply minimum image for length between -0.5 and 0.5.               */
  /*                                                                    */
  /* ================================================================== */
	if (sx < -0.5) sx = sx + 1.0;
	if (sx >  0.5) sx = sx - 1.0;
	if (sy < -0.5) sy = sy + 1.0;
	if (sy >  0.5) sy = sy - 1.0;
	if (sz < -0.5) sz = sz + 1.0;
	if (sz >  0.5) sz = sz - 1.0;

  /* ================================================================== */
  /*                                                                    */
  /* Obtain the real coordinates.                                       */
  /*                                                                    */
  /* ================================================================== */
	*dx = axes[k][0] * sx + axes[k][1] * sy + axes[k][2] * sz;
	*dy = axes[k][3] * sx + axes[k][4] * sy + axes[k][5] * sz;
	*dz = axes[k][6] * sx + axes[k][7] * sy + axes[k][8] * sz;
	

}//end subroudine
#endif


