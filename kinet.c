/* ======================================================================== */
/* kinet.cpp                                                                */
/*                                                                          */
/*		This subroutine calculates the instantaneous kinetic energy.        */
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
void kinet (int ibox)
{
  int k = ibox;
  /* ================================================================== */
  /*                                                                    */
  /* Zero out kinetic energy tensor.                                    */
  /*                                                                    */
  /* ================================================================== */

  for(int i=0; i<6; i++)
	{
		en[k].ken[i] = 0.0;
	}


  /* ================================================================== */
  /*                                                                    */
  /* Accumulate the kinectic energy tensor.                             */
  /*                                                                    */
  /* ================================================================== */
  for(int i=0; i<box[k].boxns; i++) {
    double fac = pott[k][i].mas;
    double vvx = vv[k][i].x;
    double vvy = vv[k][i].y;
    double vvz = vv[k][i].z;
    double uux = uu[k][i].x;
    double uuy = uu[k][i].y;
    double uuz = uu[k][i].z;
    en[k].ken[0] += fac * (vvx + uux)*(vvx + uux);
    en[k].ken[1] += fac * (vvy + uuy)*(vvx + uux);
    en[k].ken[2] += fac * (vvy + uuy)*(vvy + uuy);
    en[k].ken[3] += fac * (vvz + uuz)*(vvx + uux);
    en[k].ken[4] += fac * (vvz + uuz)*(vvy + uuy);
    en[k].ken[5] += fac * (vvz + uuz)*(vvz + uuz);
  }
  /* ================================================================== */
  /*                                                                    */
  /* Convert kinetic energy to kJ/mol                                   */
  /*                                                                    */
  /* ================================================================== */
  for(int i=0; i<6; i++) en[k].ken[i] *= (0.125) * (NA*AU*10.0);
  /* ================================================================== */
  /*                                                                    */
  /* Calculate the scalar kinetic energy.                               */
  /*                                                                    */
  /* ================================================================== */
  en[k].kinet = en[k].ken[0] + en[k].ken[2] + en[k].ken[5];



}
