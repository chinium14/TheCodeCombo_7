#ifdef MC
#ifndef NEUTRAL
#include "defines.h"



/* ************************************************************************ */
/* ======================== Begin Subroutine ============================== */
/* ************************************************************************ */

double eurey_mc (int ibox)
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
    double dx = atom[k][ib].x - atom[k][ia].x;
    double dy = atom[k][ib].y - atom[k][ia].y;
    double dz = atom[k][ib].z - atom[k][ia].z;

  /* ----------------------------------------- */
  /* Apply minimum image convention.           */
  /* ----------------------------------------- */
    if(dx >  hx) dx -= lx;
    else if(dx < -hx) dx += lx;
    if(dy >  hy) dy -= ly;
    else if(dy < -hy) dy += ly;
    if(dz >  hz) dz -= lz;
    else if(dz < -hz) dz += lz;
  /* ----------------------------------------- */
  /* Calculate the length of the bond and then */
  /* the energy using the CHARMM force field.  */
  /* ----------------------------------------- */

    double drr = dx*dx + dy*dy + dz*dz;
    double dr  = sqrt(drr);
    double dri = 1.0 / dr;

    double rr = dr - urey[k][i].Seq;
    double fr = -urey[k][i].k * rr*2.0;
    eurey += -(fr * rr * 0.5);

  /* ----------------------------------------- */
  /* Calculate the forces.                     */
  /* ----------------------------------------- */
    double fx = fr * dri * dx;
    double fy = fr * dri * dy;
    double fz = fr * dri * dz;

    ffox[k][ia] -= fx;
    ffoy[k][ia] -= fy;
    ffoz[k][ia] -= fz;
    ffox[k][ib] += fx;
    ffoy[k][ib] += fy;
    ffoz[k][ib] += fz;

  /* ----------------------------------------- */
  /* If PRESSURE is defined, accumulate the    */
  /* components of the virial tensor.          */
  /* ----------------------------------------- */
#ifdef PRESSURE
    wxx += fx * dx;
    wyx += fy * dx;
    wyy += fy * dy;
    wzx += fz * dx;
    wzy += fz * dy;
    wzz += fz * dz;
#endif

	}

  /* ================================================================== */
  /*                                                                    */
  /* After the loop, if PRESSURE is defined, store the bonding          */
  /* contribution to the virial tensor to pvir[k].bond[i].              */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  pviro[k].urey[0] = wxx;
  pviro[k].urey[1] = wyx;
  pviro[k].urey[2] = wyy;
  pviro[k].urey[3] = wzx;
  pviro[k].urey[4] = wzy;
  pviro[k].urey[5] = wzz;
#endif

  return (eurey);

}
#endif //NEUTRAL
#endif//MC
