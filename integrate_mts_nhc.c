/* ======================================================================== */
/* integrate_mts_nhc.cpp                                                    */
/*                                                                          */
/*		This subroutine is one of several that integrate the equations of   */
/* motion.  The subroutine called is controled by sim.ID2 which is          */
/* specified in simul.input.  This subroutine employs the XI and XO RESPA   */
/* NVT multiple time step algorithms (This mode(XI/XO) is also specified in */
/* simul.input).                                                            */
/*                                                                          */
/* See: Tuckerman, M., et. al., "Reversible multiple time scale molecular   */
/* dynamics," J. Chem. Phys., 97 (3) 1990-2001 1 August 1992.               */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                      vscale:		The scaling value for the velocities    */
/*                                  needed for NVT simulations.  For this   */
/*                                  subroutine, it is always 1.0.           */
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
void integrate_nhc(double, int);

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void integrate_mts_nhc (int ibox, double vscal)
{
  int k = ibox;

  /* ================================================================== */
  /*                                                                    */
  /* Assign the time steps needed for integration.                      */
  /*                                                                    */
  /* ================================================================== */
  double dt_short	=	sim.dt;
  double dt_shorth	=	0.5* sim.dt;
  //double dt_long	=	sim.dtlong;
  double dt_longh	=	0.5*sim.dtlong;
 
  /* ================================================================== */
  /*                                                                    */
  /* The NVT-XO-RESPA/NVT-XI-RESPA algorithm.                           */
  /*                                                                    */
  /* ================================================================== */

  /* -------------------------------------------------------- */
  /* Update the particle velocites, thermostat velocites      */
  /* and thermostat positions to half a long time step for    */
  /* XO-RESPA and half a short time step for XI-RESPA.        */
  /* -------------------------------------------------------- */
  if(sim.xRESPA[k] == 0) integrate_nhc(dt_longh,k); 
  for (int j=0; j<sim.nsteps; j++){
	if(sim.xRESPA[k] == 1) integrate_nhc(dt_shorth,k);

  /* -------------------------------------------------------- */
  /* Update the particle velocities with the long time force  */
  /* to half a long time step if this is the first time       */
  /* through the loop.                                        */
  /* -------------------------------------------------------- */
	if(j==0)
	{
		for(int i=0; i<box[k].boxns; i++) {

			double imass=100.0/pott[k][i].mas;

   			vv[k][i].x = (vv[k][i].x + dt_longh * imass*  ff_long[k][i].x);
			vv[k][i].y = (vv[k][i].y + dt_longh * imass*  ff_long[k][i].y);
			vv[k][i].z = (vv[k][i].z + dt_longh * imass*  ff_long[k][i].z);
		}//for i
	}//if(j==0) 
	
  /* -------------------------------------------------------- */
  /* Update the particle velocites with the short force for   */
  /* half a short time step.                                  */
  /* -------------------------------------------------------- */
	for (int i=0; i< box[k].boxns; i++){
	  	double imass=100.0/pott[k][i].mas;

		vv[k][i].x = (vv[k][i].x + dt_shorth * imass*  ff_short[k][i].x);
		vv[k][i].y = (vv[k][i].y + dt_shorth * imass*  ff_short[k][i].y);
		vv[k][i].z = (vv[k][i].z + dt_shorth * imass*  ff_short[k][i].z);

  /* -------------------------------------------------------- */
  /* Update the particle positions to a whole short time      */
  /* step.                                                    */
  /* -------------------------------------------------------- */
		double ddx = dt_short * vv[k][i].x;
		double ddy = dt_short * vv[k][i].y;
		double ddz = dt_short * vv[k][i].z;
    
		atnopbc[k][i].x += ddx;
		atnopbc[k][i].y += ddy;
		atnopbc[k][i].z += ddz;
		atom[k][i].x += ddx;
	 	atom[k][i].y += ddy;
		atom[k][i].z += ddz;
	}//for i
		
  /* -------------------------------------------------------- */
  /* Apply periodic boundary conditions.                      */
  /* -------------------------------------------------------- */
	pbc_all(k);

  /* -------------------------------------------------------- */
  /* Update the short (fast changing) forces.                 */
  /* -------------------------------------------------------- */
	force_short(k);
	
  /* -------------------------------------------------------- */
  /* Update the particle velocites with the short force for   */
  /* half a short time step.                                  */
  /* -------------------------------------------------------- */
	for ( int i=0; i< box[k].boxns; i++){
		double imass=100.0/pott[k][i].mas;
		vv[k][i].x = (vv[k][i].x + dt_shorth * imass*  ff_short[k][i].x);
		vv[k][i].y = (vv[k][i].y + dt_shorth * imass*  ff_short[k][i].y);
		vv[k][i].z = (vv[k][i].z + dt_shorth * imass*  ff_short[k][i].z);
	}//for i
  
  /* -------------------------------------------------------- */
  /* Update the particle velocities with the long time force  */
  /* to half a long time step if this is the last time        */
  /* through the loop.                                        */
  /* -------------------------------------------------------- */
	if(j==sim.nsteps-1)	
	{
       /* *********************************************** */
       /* Update the long (slowly changing) forces.       */
       /* *********************************************** */
		force_long(k);
		for(int i=0; i<box[k].boxns; i++) {
  			double imass=100.0/pott[k][i].mas;
			vv[k][i].x = (vv[k][i].x + dt_longh * imass*  ff_long[k][i].x);
			vv[k][i].y = (vv[k][i].y + dt_longh * imass*  ff_long[k][i].y);
			vv[k][i].z = (vv[k][i].z + dt_longh * imass*  ff_long[k][i].z);
			uu[k][i].x = vv[k][i].x;
			uu[k][i].y = vv[k][i].y;
			uu[k][i].z = vv[k][i].z;
		}//for i
	}//if(j==sim.nsteps-1)	

  /* -------------------------------------------------------- */
  /* Update the particle velocites, thermostat velocites      */
  /* and thermostat positions to half a long time step for    */
  /* XO-RESPA and half a short time step for XI-RESPA.        */
  /* -------------------------------------------------------- */	
	if(sim.xRESPA[k] == 1) integrate_nhc(dt_shorth,k);
  }//for j
  if(sim.xRESPA[k]== 0) integrate_nhc(dt_longh,k);
}
