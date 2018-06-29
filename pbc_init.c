/* ======================================================================== */
/* pbc_init.c                                                               */
/*                                                                          */
/*		This subroutine applys periodic boundary conditions to all the      */
/* site/atom coordinates of the initial configuration. It is called only    */
/* once in read_config().                                                   */
/*                                                                          */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void pbc_init (int ibox)
{
  int k = ibox;

  /* ================================================================== */
  /*                                                                    */
  /* Box information                                                    */
  /*                                                                    */
  /* ================================================================== */
  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  //double hx = box[k].hx;
  //double hy = box[k].hy;
  //double hz = box[k].hz;

  double ilx = 1.0/lx;
  double ily = 1.0/ly;
  double ilz = 1.0/lz;


  /* ================================================================== */
  /*                                                                    */
  /* The following operations make the center of the box (0,0,0), so    */
  /* the coordinates of the sites/atoms goes from -(box length)/2.0 to  */
  /* +(box length)/2.0.													*/
  /*                                                                    */
  /* ================================================================== */
  for(int i=0; i<box[k].boxns; i++) {

    double dx = atnopbc[k][i].x;
    double dy = atnopbc[k][i].y;
    double dz = atnopbc[k][i].z;


    double qx = dx*ilx;
    int nx = (int)(qx < 0.0 ? qx-0.5 : qx+0.5);
    double qy = dy*ily;
    int ny = (int)(qy < 0.0 ? qy-0.5 : qy+0.5);
    double qz = dz*ilz;
    int nz = (int)(qz < 0.0 ? qz-0.5 : qz+0.5);

    dx -= (nx * lx);
    dy -= (ny * ly);
    dz -= (nz * lz);

    

    atom[k][i].x = dx;
    atom[k][i].y = dy;
    atom[k][i].z = dz;

  }


}  

