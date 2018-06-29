/* ======================================================================== */
/* velscale.cpp                                                             */
/*                                                                          */
/*		This subroutine calculates the scaling factor for the velocity      */
/* using the Berendsen thermostat(weak coupling) scheme.  This subroutine   */
/* is used if an NVT ensemble is selected and a Nose Hoover thermostat is   */
/* not selected in simul.input.  The subroutine returns the scaling factor  */
/* as a double and is called from main();                                   */
/*                                                                          */
/* See: Berensen et al.  1984                                               */
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
double velscale (int ibox)
{
  /* ================================================================== */
  /*                                                                    */
  /* Variables needed in the subroutine.                                */
  /*                                                                    */
  /* ================================================================== */
  int k = ibox;
  double temp  = box[k].temp;
  double temp0 = sim.T[k];

  /* ================================================================== */
  /*                                                                    */
  /* Calculate the scaling factor using the Berendsen scheme.           */
  /*                                                                    */
  /* ================================================================== */
  double tmp = 1.0 + ((temp0/temp) - 1.0) * (sim.dt / sim.Ttau[k]);
  if(tmp < 0.0) {
    fprintf(stdout,"===> problem with temperature control (negative argument)\n");
    exit(1);
  }

  double scale = sqrt(tmp);

  
  return(scale);
}
