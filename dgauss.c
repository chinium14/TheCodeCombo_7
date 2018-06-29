/* ======================================================================== */
/* dgauss.cpp                                                               */
/*                                                                          */
/*	Generates random numbers from a Gaussian distribution.  It returns  */
/* a double.                                                                */
/*                                                                          */
/* See: Allen, M. P., and D. J. Tildesley, Computer Simulations of Liquids, */
/*      p347 (1987).							    */	
/*                                                                          */
/* Passed Parameters:                                                       */
/*			ave:		The average of the distribution     */
/*                      var:		The varience of the distribution.   */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

  /* ================================================================== */
  /*                                                                    */
  /* Function Prototypes                                                */
  /*                                                                    */
  /* ================================================================== */
double ran2 (void);

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

double dgauss (double ave, double var)
{
 
  /* ================================================================== */
  /*                                                                    */
  /* Generate two uniform random numbers on the interval (0,1).         */
  /*                                                                    */
  /* ================================================================== */
  double x = ran2();
  double y = ran2();
  
  /* ================================================================== */
  /*                                                                    */
  /* Transform the uniform random numbers generated to a gaussian       */
  /* distribution.                                                      */
  /*                                                                    */
  /* ================================================================== */
  double dgauss = ave + var * sqrt(-2.0*log(x)) * cos(2.0*PI*y);
  
  return (dgauss);
}
