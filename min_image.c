/* ======================================================================== */
/* min_image.cpp                                                            */
/*                                                                          */
/*    This subroutine calculates the minimum image distances.               */
/* Written by Thomas Knotts 30 Aug 2005.                                    */
/* ======================================================================== */
#include "defines.h"
void min_image(int ibox, double *dx_ptr, double *dy_ptr, double *dz_ptr)
{
  int k = ibox;
  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;
  
  if(*dx_ptr >  hx) *dx_ptr -= lx;
  else if(*dx_ptr < -hx) *dx_ptr += lx;
  if(*dy_ptr >  hy) *dy_ptr -= ly;
  else if(*dy_ptr < -hy) *dy_ptr += ly;
  if(*dz_ptr >  hz) *dz_ptr -= lz;
  else if(*dz_ptr < -hz) *dz_ptr += lz;
}


/* ======================================================================== */
/* min_image_npt_full.cpp                                                   */
/*                                                                          */
/*    This subroutine calculates the minimum image between two particles if */
/* the box is not square or rectangular.                                    */
/*                                                                          */
/* Passed Parameters:                                                       */
/*    ibox:  box number                                                     */
/*    dx:    x distance between two particles                               */
/*    dy:    y distance between two particles                               */
/*    dz:    z distance between two particles                               */
/*                                                                          */
/* ======================================================================== */
#ifdef PR_NPT
void min_image_npt_full(int ibox, double *dx, double *dy, double *dz)
{
	int k=ibox;	
	double sx, sy, sz;

  /* ================================================================== */
  /* Obtain the scaled coordinates for each particle.                   */
  /* ================================================================== */
	sx = axesi[k][0] * *dx + axesi[k][1] * *dy + axesi[k][2] * *dz;
	sy = axesi[k][3] * *dx + axesi[k][4] * *dy + axesi[k][5] * *dz;
	sz = axesi[k][6] * *dx + axesi[k][7] * *dy + axesi[k][8] * *dz;

  /* ================================================================== */
  /* Apply minimum image for length between -0.5 and 0.5.               */
  /* ================================================================== */
	if (sx < -0.5) sx = sx + 1.0;
	if (sx >  0.5) sx = sx - 1.0;
	if (sy < -0.5) sy = sy + 1.0;
	if (sy >  0.5) sy = sy - 1.0;
	if (sz < -0.5) sz = sz + 1.0;
	if (sz >  0.5) sz = sz - 1.0;

  /* ================================================================== */
  /* Obtain the real coordinates.                                       */
  /* ================================================================== */
	*dx = axes[k][0] * sx + axes[k][1] * sy + axes[k][2] * sz;
	*dy = axes[k][3] * sx + axes[k][4] * sy + axes[k][5] * sz;
	*dz = axes[k][6] * sx + axes[k][7] * sy + axes[k][8] * sz;
	

}//end subroudine
#endif


