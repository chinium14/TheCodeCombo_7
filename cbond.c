/* ======================================================================== */
/* cbond.cpp                                                                */
/*                                                                          */
/*		This subroutine calculates the bond stretching energies and forces. */
/* It returns the value of the potential energy due to the bonds to be      */
/* stored in en[k].ebond. It is called from either forces() or              */
/* force_short().                                                           */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

#ifdef PR_NPT
void min_image_npt_full(int, double*, double*, double*);
#endif


/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

double cbond (int ibox)
{

  int k = ibox;
  /* ================================================================== */
  /*                                                                    */
  /* Zero out the energy accumulator                                    */
  /*                                                                    */
  /* ================================================================== */
  double ebond = 0.0;

  /* ================================================================== */
  /*                                                                    */
  /* Zero out the virial accumulators                                   */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  double wxx = 0.0;
  double wyx = 0.0;
  double wyy = 0.0;
  double wzx = 0.0;
  double wzy = 0.0;
  double wzz = 0.0;
#endif

  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;
  
  /* ================================================================== */
  /*                                                                    */
  /* Calculate the energy and force associated with each bond.  The     */
  /* loop is over the number of bonds.                                  */
  /*                                                                    */
  /* ================================================================== */
  for(int i=0; i<bondN[k]; i++) {
  /* ----------------------------------------- */
  /* Assign the atoms involved in the bonding. */
  /* ----------------------------------------- */
	int ia = bond[k][i].a;
    int ib = bond[k][i].b;

  /* ----------------------------------------- */
  /* Calculate the distances between the atoms */
  /* involved in the bonding.                  */
  /* ----------------------------------------- */
    double dx = atom[k][ib].x - atom[k][ia].x;
    double dy = atom[k][ib].y - atom[k][ia].y;
    double dz = atom[k][ib].z - atom[k][ia].z;

  /* ----------------------------------------- */
  /* Apply minimum image convention.           */
  /* ----------------------------------------- */
#ifndef PR_NPT
      if(dx >  hx) dx -= lx;
      else if(dx < -hx) dx += lx;
      if(dy >  hy) dy -= ly;
      else if(dy < -hy) dy += ly;
      if(dz >  hz) dz -= lz;
      else if(dz < -hz) dz += lz;
#endif
#ifdef PR_NPT
	min_image_npt_full(k, &dx, &dy, &dz);
#endif
  /* ----------------------------------------- */
  /* Calculate the length of the bond and then */
  /* the energy using the CHARMM force field.  */
  /* ----------------------------------------- */

    double drr = dx*dx + dy*dy + dz*dz;
    double dr  = sqrt(drr);
    double dri = 1.0 / dr;

    double rr = dr - bond[k][i].req;
    double fr = -bond[k][i].krbond * rr*2.0;
    ebond += -(fr * rr * 0.5);

  /* ----------------------------------------- */
  /* Calculate the forces.                     */
  /* ----------------------------------------- */
    double fx = fr * dri * dx;
    double fy = fr * dri * dy;
    double fz = fr * dri * dz;

    ff[k][ia].x -= fx;
    ff[k][ia].y -= fy;
    ff[k][ia].z -= fz;
    ff[k][ib].x += fx;
    ff[k][ib].y += fy;
    ff[k][ib].z += fz;

  /* ----------------------------------------- */
  /* If PRESSURE is defined, accumulate the    */
  /* components of the virial tensor.          */
  /* ----------------------------------------- */
#ifdef PRESSURE
    double sxx = fx * dx;
    double syx = fy * dx;
    double syy = fy * dy;
    double szx = fz * dx;
    double szy = fz * dy;
    double szz = fz * dz;
    wxx += sxx;
    wyx += syx;
    wyy += syy;
    wzx += szx;
    wzy += szy;
    wzz += szz;

/* ----------------------------------------- */
/* If ZHOU is defined, accumulate the        */
/* components of the stress tensor.          */
/* ----------------------------------------- */
#ifdef ZHOU
        stresses[k][ia][0][0] += 0.5*sxx;
        stresses[k][ia][1][0] += 0.5*syx;
        stresses[k][ia][1][1] += 0.5*syy;
        stresses[k][ia][2][0] += 0.5*szx;
        stresses[k][ia][2][1] += 0.5*szy;
        stresses[k][ia][2][2] += 0.5*szz;

        stresses[k][ib][0][0] += 0.5*sxx;
        stresses[k][ib][1][0] += 0.5*syx;
        stresses[k][ib][1][1] += 0.5*syy;
        stresses[k][ib][2][0] += 0.5*szx;
        stresses[k][ib][2][1] += 0.5*szy;
        stresses[k][ib][2][2] += 0.5*szz;
        
#endif  

#endif
	}

  /* ================================================================== */
  /*                                                                    */
  /* After the loop, if PRESSURE is defined, store the bonding          */
  /* contribution to the virial tensor to pvir[k].bond[i].              */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  pvir[k].bond[0] = wxx;
  pvir[k].bond[1] = wyx;
  pvir[k].bond[2] = wyy;
  pvir[k].bond[3] = wzx;
  pvir[k].bond[4] = wzy;
  pvir[k].bond[5] = wzz;
#endif

  return (ebond);

}
