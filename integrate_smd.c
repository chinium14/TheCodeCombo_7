/* ======================================================================== */
/* integrate_smd.cpp                                                    */
/*                                                                          */
/*		This subroutine is one of several that integrate the equations of   */
/* motion.  The subroutine called is controled by sim.ID2 which is          */
/* specified in simul.input.  This subroutine integrates the equations of   */
/* motion using the Velocity Verlet algorithm.                              */
/*                                                                          */
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
#ifdef SMD
  /* ================================================================== */
  /*                                                                    */
  /* Function Prototypes                                                */
  /*                                                                    */
  /* ================================================================== */
void forces(int);
void force_short(int);
void force_long (int);
void pbc_all(int);

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void integrate_smd (int ibox, double vscal)
{
  int k = ibox;

  /* ================================================================== */
  /*                                                                    */
  /* Assign the time steps needed for integration.                      */
  /*                                                                    */
  /* ================================================================== */
  double dt	=	sim.dt;
  double dt_h	=	0.5* sim.dt;

  /* ================================================================== */
  /*                                                                    */
  /* Apply the Velocity Verlet algorithm to advance the momenta and     */
  /* postions.                                                          */
  /*                                                                    */
  /* ================================================================== */
  
  double vel	= steermd[k].v;
  steermd[k].x2	+= vel*dt*steermd[k].ex;
  steermd[k].y2	+= vel*dt*steermd[k].ey;
  steermd[k].z2	+= vel*dt*steermd[k].ez;


  for(int i=0; i<box[k].boxns; i++) {
	  if(i!= steermd[k].site1){
		  double imass=100.0/pott[k][i].mas;

	  /* ----------------------------------------------- */
	  /* Update the velocities to half a time step.      */
	  /* ----------------------------------------------- */
   		  vv[k][i].x = vscal *(vv[k][i].x + dt_h * imass*  ff[k][i].x);
		  vv[k][i].y = vscal *(vv[k][i].y + dt_h * imass*  ff[k][i].y);
		  vv[k][i].z = vscal *(vv[k][i].z + dt_h * imass*  ff[k][i].z);

	  /* ----------------------------------------------- */
	  /* Update the postions to a whole time step.       */
	  /* ----------------------------------------------- */
		  double ddx = dt* vv[k][i].x;
		  double ddy = dt* vv[k][i].y;
		  double ddz = dt* vv[k][i].z;
    
		  atnopbc[k][i].x += ddx;
		  atnopbc[k][i].y += ddy;
		  atnopbc[k][i].z += ddz;
		  atom[k][i].x += ddx;
		  atom[k][i].y += ddy;
		  atom[k][i].z += ddz;
	  }
  }
	  /* ----------------------------------------------- */
	  /* Apply periodic boundary conditions.             */
	  /* ----------------------------------------------- */

	  pbc_all(k);

	  /* ----------------------------------------------- */
	  /* Update the forces.                              */
	  /* ----------------------------------------------- */
	  forces(k);

	  /* ----------------------------------------------- */
	  /* Update the velocities to half a time step.      */
	  /* ----------------------------------------------- */
	  for ( int i=0; i< box[k].boxns; i++){
		  if(i== steermd[k].site1){
			  vv[k][i].x = 0.0;
			  vv[k][i].y = 0.0;
			  vv[k][i].z = 0.0;
			  uu[k][i].x = vv[k][i].x;
			  uu[k][i].y = vv[k][i].y;
			  uu[k][i].z = vv[k][i].z;
		  }
		  else{
			  double imass=100.0/pott[k][i].mas;
			  vv[k][i].x = (vv[k][i].x + vscal *dt_h * imass*  ff[k][i].x);
			  vv[k][i].y = (vv[k][i].y + vscal *dt_h * imass*  ff[k][i].y);
			  vv[k][i].z = (vv[k][i].z + vscal *dt_h * imass*  ff[k][i].z);
			  uu[k][i].x = vv[k][i].x;
			  uu[k][i].y = vv[k][i].y;
			  uu[k][i].z = vv[k][i].z;
		  }

	  }
//  }
}
 
#endif// SMD

#ifdef NSMD
  /* ================================================================== */
  /*                                                                    */
  /* Function Prototypes                                                */
  /*                                                                    */
  /* ================================================================== */
void forces(int);
void force_short(int);
void force_long (int);
void pbc_all(int);

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void integrate_smd (int ibox, double vscal)
{
  int k = ibox;

  /* ================================================================== */
  /*                                                                    */
  /* Assign the time steps needed for integration.                      */
  /*                                                                    */
  /* ================================================================== */
  double dt	=	sim.dt;
  double dt_h	=	0.5* sim.dt;

  /* ================================================================== */
  /*                                                                    */
  /* Apply the Velocity Verlet algorithm to advance the momenta and     */
  /* postions.                                                          */
  /*                                                                    */
  /* ================================================================== */
  
  double vel	= steermd[k].v;
  atnopbc[k][steermd[k].site2].x	+= vel*dt*steermd[k].ex;
  atnopbc[k][steermd[k].site2].y	+= vel*dt*steermd[k].ey;
  atnopbc[k][steermd[k].site2].z	+= vel*dt*steermd[k].ez;

  atom[k][steermd[k].site2].x	+= vel*dt*steermd[k].ex;
  atom[k][steermd[k].site2].y	+= vel*dt*steermd[k].ey;
  atom[k][steermd[k].site2].z	+= vel*dt*steermd[k].ez;


  for(int i=0; i<box[k].boxns; i++) {
	  if(i!= steermd[k].site1 && i!=steermd[k].site2){
		  double imass=100.0/pott[k][i].mas;

	  /* ----------------------------------------------- */
	  /* Update the velocities to half a time step.      */
	  /* ----------------------------------------------- */
   		  vv[k][i].x = vscal *(vv[k][i].x + dt_h * imass*  ff[k][i].x);
		  vv[k][i].y = vscal *(vv[k][i].y + dt_h * imass*  ff[k][i].y);
		  vv[k][i].z = vscal *(vv[k][i].z + dt_h * imass*  ff[k][i].z);

	  /* ----------------------------------------------- */
	  /* Update the postions to a whole time step.       */
	  /* ----------------------------------------------- */
		  double ddx = dt* vv[k][i].x;
		  double ddy = dt* vv[k][i].y;
		  double ddz = dt* vv[k][i].z;
    
		  atnopbc[k][i].x += ddx;
		  atnopbc[k][i].y += ddy;
		  atnopbc[k][i].z += ddz;
		  atom[k][i].x += ddx;
		  atom[k][i].y += ddy;
		  atom[k][i].z += ddz;
	  }
  }
	  /* ----------------------------------------------- */
	  /* Apply periodic boundary conditions.             */
	  /* ----------------------------------------------- */

	  pbc_all(k);

	  /* ----------------------------------------------- */
	  /* Update the forces.                              */
	  /* ----------------------------------------------- */
	  forces(k);

	  /* ----------------------------------------------- */
	  /* Update the velocities to half a time step.      */
	  /* ----------------------------------------------- */
	  for ( int i=0; i< box[k].boxns; i++){
		 if(i!= steermd[k].site1 && i!=steermd[k].site2){
			  double imass=100.0/pott[k][i].mas;
			  vv[k][i].x = (vv[k][i].x + vscal *dt_h * imass*  ff[k][i].x);
			  vv[k][i].y = (vv[k][i].y + vscal *dt_h * imass*  ff[k][i].y);
			  vv[k][i].z = (vv[k][i].z + vscal *dt_h * imass*  ff[k][i].z);
			  uu[k][i].x = vv[k][i].x;
			  uu[k][i].y = vv[k][i].y;
			  uu[k][i].z = vv[k][i].z;
		  }

	  }
//  }
}
 
#endif//NSMD
