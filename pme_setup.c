/* ======================================================================== */
/* pme_setup.cpp                                                            */
/*            Setup routine for PME                                         */
/* ======================================================================== */

#include "defines.h"

#ifdef SPME

void fillspline (double,double*,double*);
void pme_dftmod (double*,double*,int);
void pme_bspline (void);
void pme_init (int);

void pme_chkbox (int k, int swap)
{


  int bx = (int)(box[k].lx);
  int by = (int)(box[k].ly);
  int bz = (int)(box[k].lz);

  if(oddeven(bx)) bx += 1;
  if(oddeven(by)) by += 1;
  if(oddeven(bz)) bz += 1;

  int chx = 0;
  int chy = 0;
  int chz = 0;

  if(sim.grid[0] != bx) {sim.grid[0] = bx; chx = 1;};
  if(sim.grid[1] != by) {sim.grid[1] = by; chy = 1;};
  if(sim.grid[2] != bz) {sim.grid[2] = bz; chz = 1;};
  
  if(chx || chy || chz || swap) {

    if( (bx >= gridmx) || (by >= gridmx) || (bz >= gridmx) ) 
      fprintf(stderr,"WARNING: box dimensions (%i*%i*%i) exceeds PME grid (%i) !!!\n", bx, by, bz, gridmx); 

    int ngrid = bx*by*bz;
    for(int i=0; i<ngrid; i++) {
      qgrida[i][0] = 0.0;   qgrida[i][1] = 0.0;
    }

    fftw_destroy_plan(qgrid_forward);
    fftw_destroy_plan(qgrid_back);
    qgrid_forward = fftw_plan_dft_3d(bx,by,bz,qgrida,qgrida,FFTW_FORWARD,FFTW_PATIENT);
    qgrid_back    = fftw_plan_dft_3d(bx,by,bz,qgrida,qgrida,FFTW_BACKWARD,FFTW_PATIENT);

    pme_bspline();

  }


}



void pme_init (int k)
{

  int ORDER = sim.order;
  int GRIDX = sim.grid[0];
  int GRIDY = sim.grid[1];
  int GRIDZ = sim.grid[2];

  qgrid_forward = fftw_plan_dft_3d(GRIDX,GRIDY,GRIDZ,qgrida,qgrida,FFTW_FORWARD,FFTW_PATIENT);
  qgrid_back    = fftw_plan_dft_3d(GRIDX,GRIDY,GRIDZ,qgrida,qgrida,FFTW_BACKWARD,FFTW_PATIENT);

  pme_bspline();

}


void pme_clean (void)
{

  fftw_destroy_plan(qgrid_forward);
  fftw_destroy_plan(qgrid_back);

  fftw_free(qgrida);

}


void pme_bspline (void)
{

  int ORDER = sim.order;
  int GRIDX = sim.grid[0];
  int GRIDY = sim.grid[1];
  int GRIDZ = sim.grid[2];

  double w = 0.0;
  fillspline(w,mn,dmn);

  int ngrid = greater(GRIDX,greater(GRIDY,GRIDZ))+1;

  for(int i=1; i<=ngrid; i++) bsp_arr[i] = 0.0; 

  for(int i=2; i<=(ORDER+1); i++) bsp_arr[i] = mn[i-1];

  pme_dftmod(bsp_modx,bsp_arr,GRIDX);
  pme_dftmod(bsp_mody,bsp_arr,GRIDY);
  pme_dftmod(bsp_modz,bsp_arr,GRIDZ);

}


void pme_dftmod (double *bmod, double *barr, int grid)
{

  double twopi = 2.0*PI;
  double tiny = 1.0E-7;

  for(int i=1; i<=grid; i++) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    for(int j=1; j<=grid; j++) {
      double arg = twopi*(i-1)*(j-1)/grid;
      sum1 += barr[j] * cos(arg);
      sum2 += barr[j] * sin(arg);
    }
    bmod[i] = sum1*sum1 + sum2*sum2;
  }

  //for(int i=1; i<=grid; i++) 
  //  if(bmod[i] < tiny) bmod[i] = 0.5*(bmod[i-1] + bmod[i+1]);

}

#endif
