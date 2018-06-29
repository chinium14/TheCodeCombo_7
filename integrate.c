/* ======================================================================== */
/* integrate.cpp                                                            */
/*                                                                          */
/*		This subroutine is one of several that integrate the equations of   */
/* motion.  The subroutine called is controled by sim.ID2 which is          */
/* specified in simul.input.  This subroutine integrates the equations of   */
/* motion using the Leaq-Frog algorithm.                                    */
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

  /* ================================================================== */
  /*                                                                    */
  /* Function Prototypes                                                */
  /*                                                                    */
  /* ================================================================== */
void forces(int);
void pbc_all(int);
void force_short(int);
void force_long(int);

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

void integrate (int ibox, double vscal)
{
  int k = ibox;

  /* ================================================================== */
  /*                                                                    */
  /* Assign the time step needed for integration.                       */
  /*                                                                    */
  /* ================================================================== */
  double dt = sim.dt;
  double rdt = sim.dt * 100.0;    // 100.0 factor for correct units
 
  /* ================================================================== */
  /*                                                                    */
  /* Apply periodic boundary conditions.                                */
  /*                                                                    */
  /* ================================================================== */
 pbc_all(k);

  /* ================================================================== */
  /*                                                                    */
  /* Calculate the forces on the atoms/sites.                           */
  /*                                                                    */
  /* ================================================================== */
 forces(k);

  /* ================================================================== */
  /*                                                                    */
  /* Apply Leap-Frog algorithm to advance the momenta and postions.     */
  /*                                                                    */
  /* ================================================================== */
  for(int i=0; i<box[k].boxns; i++) {
    double imass = 1.0 / pott[k][i].mas;
  /* ----------------------------------------------- */
  /* Assign old velocity to holder uu. (Needed to    */
  /* calculate velocity at time t rather than -t/2   */
  /* or t/2, the velocities used in the Leap Frog    */
  /* algorithm).                                     */
  /* ----------------------------------------------- */
    uu[k][i].x = vv[k][i].x;
    uu[k][i].y = vv[k][i].y;
    uu[k][i].z = vv[k][i].z;

  /* ----------------------------------------------- */
  /* Update velocities.                              */
  /* ----------------------------------------------- */
    vv[k][i].x = vscal * (vv[k][i].x + rdt * imass * ff[k][i].x);
    vv[k][i].y = vscal * (vv[k][i].y + rdt * imass * ff[k][i].y);
    vv[k][i].z = vscal * (vv[k][i].z + rdt * imass * ff[k][i].z);

  /* ----------------------------------------------- */
  /* Update positions.                               */
  /* ----------------------------------------------- */
    double ddx = dt * vv[k][i].x;
    double ddy = dt * vv[k][i].y;
    double ddz = dt * vv[k][i].z;
    atnopbc[k][i].x += ddx;
    atnopbc[k][i].y += ddy;
    atnopbc[k][i].z += ddz;
    atom[k][i].x += ddx;
    atom[k][i].y += ddy;
    atom[k][i].z += ddz;
  }

}

