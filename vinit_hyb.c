/* ======================================================================== */
/* vinit_hyb.cpp                                                            */
/*                                                                          */
/*		This subroutine generates initial velocites from a gaussian         */
/* distribution whose mean corresponds to the simulation temperature. It    */
/* is called if no simul.vel# is found in the INPUT folder. The net linear  */
/* momentum is not set to zero.												*/
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

  /* ================================================================== */
  /*                                                                    */
  /* Function Prototypes                                                */
  /*                                                                    */
  /* ================================================================== */
double dgauss(double,double);
double ran2();

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void vinit_hyb (int k) 
{


  /* ================================================================== */
  /*                                                                    */
  /* Initialize velocites from a Gaussing destribution.                 */
  /*                                                                    */
  /* ================================================================== */
// for(int k=0; k<sim.NB; k++) {
    for(int i=0; i<box[k].boxns; i++) {
      double imass	= 1.0/pott[k][i].mas;
	  double ave	= 0.0;
      double var	= sqrt(RG*sim.T[k]*0.1*imass);
      vv[k][i].x	= 1.0*dgauss(ave,var);
      vv[k][i].y	= 1.0*dgauss(ave,var);
      vv[k][i].z	= 1.0*dgauss(ave,var);
    }
 // }


//  for(int k=0; k<sim.NB; k++) {

  
  /* ================================================================== */
  /*                                                                    */
  /* Determine the current kinetic energy                               */
  /*                                                                    */
  /* ================================================================== */
    double kinet = 0.0;
    for(int i=0; i<box[k].boxns; i++) {
  //    double imass = 1.0 / pott[k][i].mas;
      kinet += pott[k][i].mas * (vv[k][i].x * vv[k][i].x +
				 vv[k][i].y * vv[k][i].y +
				 vv[k][i].z * vv[k][i].z);
    }

    kinet *= 0.5;

    double nfree = (3.0 * box[k].boxns);

  /* ================================================================== */
  /*                                                                    */
  /* Determine the scale for the velocites to get the desired average   */
  /* temperature and scale the temperatures.                            */
  /*                                                                    */
  /* ================================================================== */
    double kin0 = sim.T[k] * nfree * 0.5 * RG * 0.10;
	double scal = sqrt(kin0 / kinet);

    for(int i=0; i<box[k].boxns; i++) {
      vv[k][i].x *= scal;
      vv[k][i].y *= scal;
      vv[k][i].z *= scal;
    }

  /* ================================================================== */
  /*                                                                    */
  /* Set both velocity arrays to the new values.                        */
  /*                                                                    */
  /* ================================================================== */
    for(int i=0; i<box[k].boxns; i++) {
      uu[k][i].x = vv[k][i].x;
      uu[k][i].y = vv[k][i].y;
      uu[k][i].z = vv[k][i].z;
    }

//  }

}
