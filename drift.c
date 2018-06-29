/* ======================================================================== */
/* drift.cpp                                                                */
/*                                                                          */
/*		This subroutine adjusts the velocites to remove drift so that the   */
/* net momentum rematains zero.                                             */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void zero_angular_momentum(int);

void drift (int ibox)
{
  int k = ibox;
#ifndef ZEROAM
  /* ================================================================== */
  /*                                                                    */
  /* Zero out momentum accumulators.                                    */
  /*                                                                    */
  /* ================================================================== */
  double vvx = 0.0;
  double vvy = 0.0;
  double vvz = 0.0;

  /* ================================================================== */
  /*                                                                    */
  /* Calculate net momentum in each direction.                          */
  /*                                                                    */
  /* ================================================================== */
  for(int i=0; i<box[k].boxns; i++) {
    vvx += pott[k][i].mas * vv[k][i].x;
    vvy += pott[k][i].mas * vv[k][i].y;
    vvz += pott[k][i].mas * vv[k][i].z;
  }

  /* ================================================================== */
  /*                                                                    */
  /* Calculate the correction to the momentum.                          */
  /*                                                                    */
  /* ================================================================== */
#ifndef KONS
  vvx /= (1.0*box[k].boxns);
  vvy /= (1.0*box[k].boxns);
  vvz /= (1.0*box[k].boxns);
#endif
#ifdef KONS
  vvx /= (1.0*(box[k].boxns-cons[k].n));
  vvy /= (1.0*(box[k].boxns-cons[k].n));
  vvz /= (1.0*(box[k].boxns-cons[k].n));
#endif
  //double vvv = sqrt(vvx*vvx + vvy*vvy + vvz*vvz);

  /* ================================================================== */
  /*                                                                    */
  /* Subtract the correction form the atomic velocities.                */
  /*                                                                    */
  /* ================================================================== */
  for(int i=0; i<box[k].boxns; i++) {
    double imass = 1.0 / pott[k][i].mas;
    vv[k][i].x -= vvx * imass;
    vv[k][i].y -= vvy * imass;
    vv[k][i].z -= vvz * imass;
  }
#ifdef KONS
  for(int i=0; i<cons[k].n; i++){
	int j=cons[k].site[i];
	vv[k][j].x = 0.0;
	vv[k][j].y = 0.0;
	vv[k][j].z = 0.0;
  }
#endif

#else //Do this if ZEROAM is defined
	zero_angular_momentum(k);
        //Note: zero_angular_momentum removes the linear momemtum.
#endif
  /* ================================================================== */
  /*                                                                    */
  /* Check for net linear momentum.                                     */
  /* (Un-comment the code below if needed for debugging or if you want  */
  /* to add error handling.)                                            */
  /*                                                                    */
  /* ================================================================== */
/*
  vvx = 0.0;
  vvy = 0.0;
  vvz = 0.0;
  for(int i=0; i<box[k].boxns; i++) {
    vvx += pott[k][i].mas * vv[k][i].x;
    vvy += pott[k][i].mas * vv[k][i].y;
    vvz += pott[k][i].mas * vv[k][i].z;
  }
  vvx /= (1.0*box[k].boxns);
  vvy /= (1.0*box[k].boxns);
  vvz /= (1.0*box[k].boxns);
  vvv = sqrt(vvx*vvx + vvy*vvy + vvz*vvz);
*/
}

