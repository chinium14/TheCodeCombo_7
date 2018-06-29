#ifdef WALL
/* ======================================================================== */
/* wall_angle.cpp		                                                    */
/* This suboutine calculates the angle between the a site sent as a         */
/* parameter and a wall at z=0.                                             */
/*                                                                          */
/* Written by Thomas Knotts 7 June 04                                       */	
/* ======================================================================== */
#include "defines.h"
double wall_angle(int ibox, int isite){
  int k=ibox;
  int site = isite;
  double x = atom[k][site].x;
  double y = atom[k][site].y;
  double z = atom[k][site].z;
  return(atan(z/sqrt(x*x + y*y)));
}
#endif

