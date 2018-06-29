/* ======================================================================== */
/* force_short.cpp                                                          */
/*                                                                          */
/*		This subroutine is part of the driver for the force calculations    */
/* if the multiple time step (RESPA) integration is selected in             */
/* simul.input. It calls the subroutines that calculate the fast-changing   */
/* forces.  It zero's out the forces on the atoms (e.g. ff[k][i].x) so that */
/* the contributions to the different energy modes can be added in their    */
/* respective subroutines (e.g. cbond()).  It also zero's out the virial    */
/* accumulators.                                                            */
/*                                                                          */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

  /* ================================================================== */
  /*                                                                    */
  /* Function Prototypes                                                */
  /*                                                                    */
  /* ================================================================== */
double cbond (int);
double cbend (int);
#ifndef NEUTRAL
double curey(int);
#endif
double ctorsion (int);


#ifdef SASA
double csasa(int);
#endif
#ifdef REST
double crestraint(int);
#endif

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

void force_short (int ibox)
{
  int k = ibox;

  /* ================================================================== */
  /*                                                                    */
  /* Zero out atomic forces arrays.                                     */
  /*                                                                    */
  /* ================================================================== */
  for(int i=0; i<box[k].boxns; i++) {
    ff[k][i].x = 0.0;
    ff[k][i].y = 0.0;
    ff[k][i].z = 0.0;
  }

  /* ================================================================== */
  /*                                                                    */
  /* Zero out the virial accumulators                                   */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  for(int i=0; i<6; i++) {
    pvir[k].bond[i]  = 0.0;
    pvir[k].bend[i]  = 0.0;
#ifdef SASA
	pvir[k].sasa[i]  = 0.0;
#endif

  }
#endif

  /* ================================================================== */
  /*                                                                    */
  /* Call subroutines that calculate energy and forces from the         */
  /* different contributions to the force fields that are fast          */
  /* changing.															*/
  /*                                                                    */
  /* ================================================================== */
  
  en[k].bond = cbond(k);
#ifdef REST
  en[k].bond += crestraint(k);
#endif
  en[k].bend = cbend(k);
#ifndef NEUTRAL
  en[k].urey = curey(k);
#endif
#ifdef SASA
  en[k].esasa= csasa(k);
#endif



  /* ----------------------------------------------------------	*/
  /* If FLIM is defined, check to make sure the forces are not  */
  /* too great.  If they are, reduce them to FORCE_LIMIT.       */
  /* ----------------------------------------------------------	*/
#ifdef FLIM
  for(int i=0; i<box[k].boxns; i++){
	if(ff[k][i].x>FORCE_LIMIT) ff[k][i].x = FORCE_LIMIT;
	else if (ff[k][i].x<-FORCE_LIMIT) ff[k][i].x = -FORCE_LIMIT;
	if(ff[k][i].y>FORCE_LIMIT) ff[k][i].y = FORCE_LIMIT;
	else if (ff[k][i].y<-FORCE_LIMIT) ff[k][i].y = -FORCE_LIMIT;
	if(ff[k][i].z>FORCE_LIMIT) ff[k][i].z = FORCE_LIMIT;
	else if (ff[k][i].z<-FORCE_LIMIT) ff[k][i].z = -FORCE_LIMIT;
  }
#endif//FLIM

  /* ----------------------------------------------------------	*/
  /* If KONS is defined, set the forces and velocities of the  */
  /* constrained atoms to zero.                                 */
  /* ----------------------------------------------------------	*/
#ifdef KONS
  for(int i=0; i<cons[k].n; i++){
	int j = cons[k].site[i];
	ff[k][j].x = 0.0;
	ff[k][j].y = 0.0;
	ff[k][j].z = 0.0;
  }
#endif

  /* ================================================================== */
  /*                                                                    */
  /* Assign the value for ff_short[k][i].  The slowly changing forces.  */
  /*                                                                    */
  /* ================================================================== */
  for(int i=0; i<box[k].boxns; i++) {
    ff_short[k][i].x = ff[k][i].x;
    ff_short[k][i].y = ff[k][i].y;
    ff_short[k][i].z = ff[k][i].z;
  }

}
