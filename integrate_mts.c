/* ======================================================================== */
/* integrate_mts.cpp                                                        */
/*                                                  						*/
/*		This subroutine is one of several that integrate the equations of   */
/* motion.  The subroutine called is controled by sim.ID2 which is          */
/* specified in simul.input.  This subroutine integrates the equations of   */
/* motion using a reversible, multiple time step, Velocity Verlet           */
/* algorithm(RESPA).                                                        */
/*                                                                          */
/* See: Tuckerman, M., et. al., "Reversible multiple time scale molecular   */
/* dynamics," J. Chem. Phys., 97 (3) 1990-2001 1 August 1992.               */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                      vscale:		The scaling value for the velocities    */
/*                                  needed for NVT simulations.  It is      */
/*                                  calculated in vscale().  For NVE,       */
/*                                  vscal = 1.0.                            */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

  /* ================================================================== */
  /*                                                                    */
  /* Function Prototypes                                                */
  /*                                                                    */
  /* ================================================================== */
void force_short(int);
void force_long (int);
void pbc_all(int);

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void integrate_mts (int ibox, double vscal)
{
  int k = ibox;

  /* ================================================================== */
  /*                                                                    */
  /* Assign the time steps needed for integration.                      */
  /*                                                                    */
  /* ================================================================== */
  double dt_short	=	sim.dt;
  double dt_shorth	=	0.5* sim.dt;
//  double dt_long	=	sim.dtlong;
  double dt_longh	=	0.5*sim.dtlong;

  /* ================================================================== */
  /*                                                                    */
  /* Apply the RESPA mutiple time step reversible algorithm.            */
  /*                                                                    */
  /* ================================================================== */
  for(int i=0; i<box[k].boxns; i++) {

	  double imass=100.0/pott[k][i].mas;

  /* ----------------------------------------------- */
  /* Update the velocities to half a long time step  */
  /* using the long(slowly changing) forces.         */
  /* ----------------------------------------------- */
   	  vv[k][i].x = vscal *(vv[k][i].x + dt_longh * imass*  ff_long[k][i].x);
      vv[k][i].y = vscal *(vv[k][i].y + dt_longh * imass*  ff_long[k][i].y);
      vv[k][i].z = vscal *(vv[k][i].z + dt_longh * imass*  ff_long[k][i].z);
  }

  /* ----------------------------------------------- */
  /* Loop over n short time steps.                   */
  /* ----------------------------------------------- */
  for ( int j=0; j< sim.nsteps; j++){
	  for (int i=0; i< box[k].boxns; i++){
	  	  double imass=100.0/pott[k][i].mas;

       /* *********************************************** */
       /* Update the velocities to half a short time step */
       /* using the short(fast changing) forces.          */
       /* *********************************************** */
		  vv[k][i].x = (vv[k][i].x + vscal *dt_shorth * imass*  ff_short[k][i].x);
		  vv[k][i].y = (vv[k][i].y + vscal *dt_shorth * imass*  ff_short[k][i].y);
		  vv[k][i].z = (vv[k][i].z + vscal *dt_shorth * imass*  ff_short[k][i].z);

       /* *********************************************** */
       /* Update the positions to a full short time step. */
       /* *********************************************** */
		  double ddx = dt_short * vv[k][i].x;
		  double ddy = dt_short * vv[k][i].y;
		  double ddz = dt_short * vv[k][i].z;
    
		  atnopbc[k][i].x += ddx;
		  atnopbc[k][i].y += ddy;
		  atnopbc[k][i].z += ddz;
		  atom[k][i].x += ddx;
	 	  atom[k][i].y += ddy;
		  atom[k][i].z += ddz;
	  }
       /* *********************************************** */
       /* Apply periodic boundary conditions.             */
       /* *********************************************** */
	  pbc_all(k);

       /* *********************************************** */
       /* Update the short range forces.                  */
       /* *********************************************** */
	  force_short(k);

       /* *********************************************** */
       /* Update the velocities to half a short time step */
       /* using the short(fast changing) forces.          */
       /* *********************************************** */
	  for ( int i=0; i< box[k].boxns; i++){
		  double imass=100.0/pott[k][i].mas;
		  vv[k][i].x = (vv[k][i].x + vscal *dt_shorth * imass*  ff_short[k][i].x);
		  vv[k][i].y = (vv[k][i].y + vscal *dt_shorth * imass*  ff_short[k][i].y);
		  vv[k][i].z = (vv[k][i].z + vscal *dt_shorth * imass*  ff_short[k][i].z);
	  }
  }

  /* ----------------------------------------------- */
  /* Update the long range forces.                   */
  /* ----------------------------------------------- */
  force_long(k);

  /* ----------------------------------------------- */
  /* Update the velocities to half a long time step  */
  /* using the long(slowly changing) forces.         */
  /* ----------------------------------------------- */
  for(int i=0; i<box[k].boxns; i++) {
  	  double imass=100.0/pott[k][i].mas;
      vv[k][i].x = (vv[k][i].x + vscal *dt_longh * imass*  ff_long[k][i].x);
      vv[k][i].y = (vv[k][i].y + vscal *dt_longh * imass*  ff_long[k][i].y);
      vv[k][i].z = (vv[k][i].z + vscal *dt_longh * imass*  ff_long[k][i].z);
	  uu[k][i].x = vv[k][i].x;
      uu[k][i].y = vv[k][i].y;
      uu[k][i].z = vv[k][i].z;
  }

}
