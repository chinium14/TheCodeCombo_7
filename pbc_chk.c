/* ======================================================================== */
/* pbc_chk.cpp                                                              */
/*                                                                          */
/*		This subroutine applys periodic boundary conditions to all the      */
/* site/atom coordinates.  It is before forces are calculated.              */
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
void pbc_chk (int iter)
{

  /* ================================================================== */
  /*                                                                    */
  /* Box information                                                    */
  /*                                                                    */
  /* ================================================================== */
  for(int k=0; k<sim.NB; k++) {
  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;

	double er = 1E-4;
  
	  FILE *err1;
	  char namepbc[40];
#ifdef MPI
	  sprintf(namepbc, "./OUTPUT/BOX%d/error%d.log",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	  sprintf(namepbc, "./OUTPUT/BOX%d/error%d.log",k,k);
#endif
	  err1= fopen(namepbc,"a");

	  for(int i=0; i<box[k].boxns; i++) {

		double dx = fabs(atnopbc[k][i].x - atom[k][i].x) + er/5.0;
		double dy = fabs(atnopbc[k][i].y - atom[k][i].y) + er/5.0;
		double dz = fabs(atnopbc[k][i].z - atom[k][i].z) + er/5.0;
    
		double fdx = fmodf((float)dx,(float)lx);
		double fdy = fmodf((float)dy,(float)ly);
		double fdz = fmodf((float)dz,(float)lz);
      
		if ((fabs(fdx) > er) || (fabs(fdy) > er) ||(fabs(fdz) > er)){
			fprintf(err1,"%lf	%lf	%lf\n",dx,dy,dz);
			fprintf(err1,"%lf	%lf	%lf\n",lx,ly,lz);
			fprintf(err1,"At iteration %d, there is error in coords of site %d:	(%lf),(%lf),(%lf)\n",iter,i,fabs(fdx),fabs(fdy),fabs(fdz));
			fprintf(err1,"pbc coords are:	x=%lf	y=%lf	z=%lf\n",atom[k][i].x,atom[k][i].y,atom[k][i].z);
			fprintf(err1,"nopbc coords are:	x=%lf	y=%lf	z=%lf\n",atnopbc[k][i].x,atnopbc[k][i].y,atnopbc[k][i].z);
			fprintf(err1,"\n");

		}
	  }
	  fclose(err1);
  }
}  

