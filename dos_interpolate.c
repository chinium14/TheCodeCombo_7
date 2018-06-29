#ifdef TDXEDOS
/* ======================================================================== */
/* dos_interpolate.cpp			                                            */
/* This file contains subroutines to interpolate the value of the density   */
/* of states between two bins.					                            */
/*                                                                          */
/* Written by Thomas Knotts           8 June 2004                           */
/*																			*/		
/* ======================================================================== */

#include "defines.h"

/* ======================================================================== */
/* dos_interp_2																*/
/* This suboutine does bilinear interpolation.  See Numerical Recipies      */
/* section 3.6.                                                             */
/* ======================================================================== */

double dos_interp_2(int ibox,double x1_val,double x2_val)
{

	int k = ibox;
	double x1 = x1_val;
	double x2 = x2_val;
	double x1_width = sim_dos[k].x1_width;
	double x2_width = sim_dos[k].x2_width;
	double x1_width_h = x1_width/2.0;
	double x2_width_h = x2_width/2.0;
	
	double g_val;

	/* -------------------------------------------------------- */
	/* Determine the bins.  Note: we cannot simply use the bins */
	/* determined in the calling subroutines because the        */
	/* interpolation grid begins at sim_dos[k].x1_begin +       */
	/* sim_dos[k].x1_width*0.5 and not at sim_dos[k].x1_begin.  */
	/* -------------------------------------------------------- */
	double x1_begin = sim_dos[k].x1_begin + x1_width_h; 
	double x2_begin = sim_dos[k].x2_begin + x2_width_h; 

	int bin_1 = (int)((x1 - x1_begin)/x1_width);
	int bin_2 = (int)((x2 - x2_begin)/x2_width);

	if(bin_1 < 0) bin_1 = 0;
	else if(bin_1 == sim_dos[k].x1_bins-1) bin_1 = sim_dos[k].x1_bins-2;

	if(bin_2 < 0) bin_2 = 0;
	else if(bin_2 == sim_dos[k].x2_bins-1) bin_2 = sim_dos[k].x2_bins-2;

	/* -------------------------------------------------------- */
	/* Calculate the interpolated/extrapolated value of the     */
	/* density of states using bilinear interpolation.          */
	/* -------------------------------------------------------- */
	double y1 = dos_hist[k][bin_1][bin_2].g_of_l;
	double y2 = dos_hist[k][bin_1+1][bin_2].g_of_l;
	double y3 = dos_hist[k][bin_1+1][bin_2+1].g_of_l;
	double y4 = dos_hist[k][bin_1][bin_2+1].g_of_l;


	double x1a_begin = dos_hist[k][bin_1][bin_2].x1_mid;
	double x1a_end   = dos_hist[k][bin_1+1][bin_2].x1_mid;
	double x2a_begin = dos_hist[k][bin_1][bin_2].x2_mid;
	double x2a_end   = dos_hist[k][bin_1][bin_2+1].x2_mid;

	double t = (x1-x1a_begin)/(x1a_end-x1a_begin);
	double u = (x2-x2a_begin)/(x2a_end-x2a_begin);

	g_val = (1.0-t)*(1.0-u)*y1 + t*(1.0-u)*y2 + t*u*y3 + (1.0-t)*u*y4;


	return (g_val);

}
#endif

#ifdef XEDOS
/* ======================================================================== */
/* dos_interp															                                 	*/
/* This suboutine does linear interpolation.                                */
/* ======================================================================== */
#include "defines.h"
double dos_interp(int ibox, int bin, double val)
{
  double g;
  int k = ibox;
  if(bin ==0){
	  g= dos_hist[k][bin].g_of_l + (dos_hist[k][bin+1].g_of_l-dos_hist[k][bin].g_of_l)*
			   (val-dos_hist[k][bin].l_mid)/sim_dos[k].l_width;
  }
  else if(bin ==sim_dos[k].l_bins-1){
	  g= dos_hist[k][bin].g_of_l + (dos_hist[k][bin].g_of_l-dos_hist[k][bin-1].g_of_l)*
			   (val-dos_hist[k][bin].l_mid)/sim_dos[k].l_width;
  }
  else {
	  g= dos_hist[k][bin].g_of_l + (dos_hist[k][bin+1].g_of_l-dos_hist[k][bin-1].g_of_l)*
			   (val-dos_hist[k][bin].l_mid)*0.5/sim_dos[k].l_width;
  }

  return(g);
}
#endif

#if defined(DOS) || defined(CTDOS)
/* ======================================================================== */
/* dos_interp															                                 	*/
/* This suboutine does linear interpolation.                                */
/* ======================================================================== */
#include "defines.h"
double dos_interp(int ibox, int bin, double val)
{
  double g;
  int k = ibox;
  if(bin ==0){
	  g= dos_hist[k][bin].g_of_e + (dos_hist[k][bin+1].g_of_e-dos_hist[k][bin].g_of_e)*
			   (val-dos_hist[k][bin].e_mid)/sim_dos[k].e_width;
  }
  else if(bin ==sim_dos[k].e_bins-1){
	  g= dos_hist[k][bin].g_of_e + (dos_hist[k][bin].g_of_e-dos_hist[k][bin-1].g_of_e)*
			   (val-dos_hist[k][bin].e_mid)/sim_dos[k].e_width;
  }
  else {
	  g= dos_hist[k][bin].g_of_e + (dos_hist[k][bin+1].g_of_e-dos_hist[k][bin-1].g_of_e)*
			   (val-dos_hist[k][bin].e_mid)*0.5/sim_dos[k].e_width;
  }
  return(g);
}
#endif
