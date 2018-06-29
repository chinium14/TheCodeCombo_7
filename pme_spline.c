/* ======================================================================== */
/* pme_spline.cpp                                                           */
/*            Calculate B-spline function and differential for PME          */
/* ======================================================================== */

#include "defines.h"

#ifdef SPME

void fillspline (double w, double *array, double *darray)
{

  int ORDER = sim.order;

  for(int i=3; i<=ORDER; i++) array[i] = 0.0;

  array[2] = w;
  array[1] = 1.0 - w;
  array[0] = 0.0;

  // generate ORDER spline
  for(int i=3; i<=ORDER-1; i++) {
    double div = (double)(1.0 / (i-1));
    array[i] = div * w * array[i-1];
    double di = (double)(i);
    for(int j=1; j<=(i-2); j++) {
      double dj = (double)(j);
      array[i-j] = div * ((w+dj)*array[i-j-1] + (di-dj-w)*array[i-j]);
    }
    array[1] *= div*(1.0-w);
  }

  // calculate differential
  darray[0] = 0.0;
  darray[1] = - array[1];
  for(int i=2; i<=ORDER; i++) darray[i] = array[i-1] - array[i];

  // generate order spline
  int kk = ORDER;
  double div = (double)(1.0 / (kk-1));
  array[kk] = div * w * array[kk-1];
  double dkk = (double)(kk);
  for(int i=1; i<=(kk-2); i++) {
    double di = (double)(i);
    array[kk-i] = div * ((w+di)*array[kk-i-1] + (dkk-di-w)*array[kk-i]);
  }
  array[1] *= div*(1.0-w);

}

#endif
