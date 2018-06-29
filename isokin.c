/* ======================================================================== */
/* isokin.cpp                                                               */
/*                                                                          */
/*		This subroutine scales the velocities of the particles to the       */
/* desired temperature producing the isokinetic ensemble.                   */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                      count:		The number of temperature values        */
/*                                  accumulated in temp_accum               */
/*                      temp_accum:	The accumulated temperature             */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void isokin (int ibox)
{
  int k = ibox;
  double scale;
  scale = sqrt(sim.T[k]/box[k].temp);
  for(int i=0; i<box[k].boxns; i++){
	vv[k][i].x *= scale;
  vv[k][i].y *= scale;
  vv[k][i].z *= scale;
	uu[k][i].x *= scale;
	uu[k][i].y *= scale;
	uu[k][i].z *= scale;
  }
}
