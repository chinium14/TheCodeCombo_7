/* ======================================================================== */
/* curey2.cpp                                                               */
/*                                                                          */
/*	This subroutine calculates the Urey-Bradley 1-3 energies and forces.    */
/* It returns the value of the potential energy due to the ureys to be      */
/* stored in en[k].eurey. It is called from either forces() or              */
/* force_short().                                                           */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                                                                          */
/* Written by Thomas A. Knotts IV		29 Aug 03							*/
/* ======================================================================== */
#ifndef NEUTRAL
#include "defines.h"
#ifdef CONFIGT
#ifdef PR_NPT
void min_image_npt_full(int, double*, double*, double*);
#endif


/* ************************************************************************ */
/* ======================== Begin Subroutine ============================== */
/* ************************************************************************ */

double curey2 (int ibox)
{

  int k = ibox;
  /* ================================================================== */
  /* Zero out the energy accumulator                                    */
  /* ================================================================== */
  double eurey = 0.0;

  /* ================================================================== */
  /* Zero out the virial accumulators                                   */
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
  /* Calculate the energy and force associated with each bond.  The     */
  /* loop is over the number of bonds.                                  */
  /* ================================================================== */
  for(int i=0; i<ureyN[k]; i++) {
  /* ----------------------------------------- */
  /* Assign the atoms involved in the bonding. */
  /* ----------------------------------------- */
	int ia = urey[k][i].a;
    int ib = urey[k][i].b;

  /* ----------------------------------------- */
  /* Calculate the distances between the atoms */
  /* involved in the bonding.                  */
  /* ----------------------------------------- */
    double xab = atom[k][ib].x - atom[k][ia].x;
    double yab = atom[k][ib].y - atom[k][ia].y;
    double zab = atom[k][ib].z - atom[k][ia].z;

  /* ----------------------------------------- */
  /* Apply minimum image convention.           */
  /* ----------------------------------------- */
#ifndef PR_NPT
    if(xab >  hx) xab -= lx;
    else if(xab < -hx) xab += lx;
    if(yab >  hy) yab -= ly;
    else if(yab < -hy) yab += ly;
    if(zab >  hz) zab -= lz;
    else if(zab < -hz) zab += lz;
#endif
#ifdef PR_NPT
	min_image_npt_full(k, &xab, &yab, &zab);
#endif
  /* ----------------------------------------- */
  /* Calculate the length of the bond and then */
  /* the energy using the CHARMM force field.  */
  /* ----------------------------------------- */

    double rab2 = xab*xab + yab*yab + zab*zab;
    double rab  = sqrt(rab2);
    double dri  = 1.0 / rab;
	double dt   = rab -urey[k][i].Seq;
    double fr = -urey[k][i].k * dt*2.0;
    eurey += -(fr * dt * 0.5);
  /*------------------------------------------ */
  /* Calculate the contribution of bonds to	   */
  /* to diagonal of hessian					   */
  /*------------------------------------------ */
	double deddt	= -fr;
	double d2eddt2	= 2.0*urey[k][i].k;
	double de		= deddt/rab;
	double term		= (d2eddt2-de) / rab2;
	double termx	= term * xab;
    double termy	= term * yab;
	double termz	= term * zab;
    double d2e11	= termx*xab + de;
	double d2e22	= termy*yab + de;
	double d2e33	= termz*zab + de;

	config[k].hesx	   += 2.0* d2e11;
	config[k].hesy	   += 2.0* d2e22;
	config[k].hesz	   += 2.0* d2e33;


  /* ----------------------------------------- */
  /* Calculate the forces.                     */
  /* ----------------------------------------- */
    double fx = fr * dri * xab;
    double fy = fr * dri * yab;
    double fz = fr * dri * zab;

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
    double sxx = fx * xab;
    double syx = fy * xab;
    double syy = fy * yab;
    double szx = fz * xab;
    double szy = fz * yab;
    double szz = fz * zab;
	wxx += sxx;
	wyx += syx;
	wyy += syy;
	wzx += szx;
	wzy += szy;
	wzz += szz;
#ifdef ZHOU
	stresses[k][ia][0][0] += sxx * 0.5;
	stresses[k][ia][1][0] += syx * 0.5;
	stresses[k][ia][1][1] += syy * 0.5;
	stresses[k][ia][2][0] += szx * 0.5;
	stresses[k][ia][2][1] += szy * 0.5;
	stresses[k][ia][2][2] += szz * 0.5;

	stresses[k][ib][0][0] += sxx * 0.5;
	stresses[k][ib][1][0] += syx * 0.5;
	stresses[k][ib][1][1] += syy * 0.5;
	stresses[k][ib][2][0] += szx * 0.5;
	stresses[k][ib][2][1] += szy * 0.5;
	stresses[k][ib][2][2] += szz * 0.5;

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
  pvir[k].urey[0] = wxx;
  pvir[k].urey[1] = wyx;
  pvir[k].urey[2] = wyy;
  pvir[k].urey[3] = wzx;
  pvir[k].urey[4] = wzy;
  pvir[k].urey[5] = wzz;
#endif

  return (eurey);

}
#endif//CONFIGT
#endif//NEUTRAL

