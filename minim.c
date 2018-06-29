/*********************************************/
/*====   Potential Energy Minimization   ====*/
/*===========     Copyright 2003    =========*/
/*=====  Tushar Jain and Nitin Rathore ======*/
/*=== Adapted from Numerical Recipes in C ===*/
/*
	frprmn.cpp - Conjugate gradient method
	linmin.cpp - Line Minimization along 
				 a direction
	mnbrak.cpp - Bracket a minima along the 
				 direction
	 brent.cpp - 
*/
/*********************************************/

#include "defines.h"

#ifdef NMA
void force_nma(double *pos_nma,double *for_nma);
double energy_nma(double *pos_nma);
void frprmn(double p[], int n, double ftol, int *iter, double *fret,double (*func)(double []), void (*dfunc)(double [], double []),int itera);
void output(int,int);
void xyz_config(int,unsigned long);
void movie (int,unsigned long);
double minim()
{
 double enermin;
 int itera,iter,i;
 double *x,*f;

 f=(double*) calloc(3*box[0].boxns+1,sizeof(double));
 x=(double*) calloc(3*box[0].boxns+1,sizeof(double));

 for(i=0;i<box[0].boxns;i++){
	x[i*3+1]=(double)(atom[0][i].x);
	x[i*3+2]=(double)(atom[0][i].y);
	x[i*3+3]=(double)(atom[0][i].z);
  }

 force_nma(x,f);
 frprmn(x,3*box[0].boxns,0,&iter,&enermin,energy_nma,force_nma,10);
 for(i=0;i<box[0].boxns;i++){
	atnopbc[0][i].x=atom[0][i].x;
	atnopbc[0][i].y=atom[0][i].y;
	atnopbc[0][i].z=atom[0][i].z;
  }
 output(iter,1);
 xyz_config(1,iter);
 movie(1,iter);
 free(f);
 free(x);
 return enermin;
}
#endif

