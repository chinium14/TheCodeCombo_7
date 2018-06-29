/* ======================================================================== */
/* pbc.c                                                                    */
/*            Apply periodic boundary conditions                            */
/* ======================================================================== */

#include "defines.h"

void pbc (int ibox, double *da, double *db)
{

  int k = ibox;
  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;

  double pbcx = 0.0;
  double pbcy = 0.0;
  double pbcz = 0.0;

  if(da[0] >  hx) pbcx =- lx;
  if(da[0] < -hx) pbcx =+ lx;
  if(da[1] >  hy) pbcy =- ly;
  if(da[1] < -hy) pbcy =+ ly;
  if(da[2] >  hz) pbcz =- lz;
  if(da[2] < -hz) pbcz =+ lz;
  
  db[0] = da[0] + pbcx;
  db[1] = da[1] + pbcy;
  db[2] = da[2] + pbcz;
  
}
