/* ======================================================================== */
/* force_long.c                                                             */
/*                                                                          */
/*		This subroutine is part of the driver for the force calculations    */
/* if the multiple time step (RESPA) integration is selected in             */
/* simul.input. It calls the subroutines that calculate the slowly-changing */
/* forces.  It zero's out the forces on the atoms (e.g. ff[k][i].x) so that */
/* the contributions to the different energy modes can be added in their    */
/* respective subroutines (e.g. cnonbond()).  It also zero's out the virial */
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
double ctorsion (int);
double cimproper (int);
#ifdef SASA
double csasa(int);
#endif
#ifndef NLIST
void   cnonbond (int,double*);
#endif
#ifdef NLIST
void   cnonbondnl (int,double*);
#endif

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

void force_long (int ibox)
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
    pvir[k].nbond[i] = 0.0;
    pvir[k].tors[i]  = 0.0;
	pvir[k].impr[i]  = 0.0;
  }
#endif

  /* ================================================================== */
  /*                                                                    */
  /* Call subroutines that calculate energy and forces from the         */
  /* different contributions to the force fields that are slowly        */
  /* changing.															*/
  /*                                                                    */
  /* ================================================================== */
  
  en[k].tors = ctorsion(k);
  en[k].impr = cimproper(k);
  double interen[8];
#ifndef NLIST
  cnonbond(k,interen);
#endif
#ifdef NLIST
  cnonbondnl(k,interen);
#endif
  en[k].nbond  = interen[0];
  en[k].nbonds = interen[1];

#ifdef COULOMB
  en[k].coulomb = interen[2];
#endif

#ifdef EWALD
  en[k].coulomb = interen[2];
  en[k].tewald  = interen[2];
  en[k].rewald  = interen[7];
  en[k].rewalds = interen[3];
  en[k].kewald  = interen[4];
  en[k].sewald  = interen[5];
  en[k].iewald  = interen[6];
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
  /* Assign the value for ff_long[k][i].  The slowly changing forces.   */
  /*                                                                    */
  /* ================================================================== */
  for(int i=0; i<box[k].boxns; i++) {
    ff_long[k][i].x = ff[k][i].x;
    ff_long[k][i].y = ff[k][i].y;
    ff_long[k][i].z = ff[k][i].z;
  }
  

}
