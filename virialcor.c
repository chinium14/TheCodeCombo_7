/* ======================================================================== */
/* virialcor.cpp                                                            */
/*                                                                          */
/*		This subroutine calculates the long-range correction for the        */
/* virial.  It is called from main() and returns a double.                  */         
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		Simulation box number                   */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
double virialcor (int ibox)
{
  /* ================================================================== */
  /*                                                                    */
  /* Variables needed in the subroutine.                                */
  /*                                                                    */
  /* ================================================================== */
  int k = ibox;
  double rc    = sqrt(box[k].rc2);
  double sum6  = 0.0;
  double sum12 = 0.0;

  /* ================================================================== */
  /*                                                                    */
  /* Loop over the sites and accumulate the correction.                 */
  /*                                                                    */
  /* ================================================================== */
  for(int i=0; i<box[k].boxns-1; i++) {
    double epsi = pott[k][i].eps;
    double sigi = pott[k][i].sig;
    for(int j=i+1; j<box[k].boxns; j++) {
      double epsj = pott[k][j].eps;
      double sigj = pott[k][j].sig;
      double ss   = 0.5 * (sigi + sigj);
      double sa   = ss * ss * ss;
      double sb   = sa * sa;
      double ee   = sqrt(epsi * epsj) * sb;

      sum6  += ee;
      sum12 += ee * sb;
    }
  }
#ifdef LJ
  double aa = 32.0 * PI / (rc*rc*rc);
  double bb = 2.0 / (3.0 * pow(rc,6.0));
#endif
#ifndef LJ
  double aa = 64.0 * PI / (rc*rc*rc);
  double bb = 1.0 / (3.0 * pow(rc,6.0));
#endif

  double cc = 1.0;

  /* ----------------------------------------------- */
  /* Calculate the correction.  0.5 factor avoids    */
  /* double counting.                                */
  /* ----------------------------------------------- */    
  double vircor =  aa * (bb * sum12 - cc * sum6);
//  double vircor2 = 16.0*PI/3.0*(box[k].boxns)*(box[k].boxns-1)*pow(pott[k][0].sig, 3.0)*pott[k][0].eps*(2*pow(rc/pott[k][0].sig,-9.0)-3*pow(rc/pott[k][0].sig, -3.0));
  
  return(vircor);

}

