/* ======================================================================== */
/* quatfit.cpp                                                              */
/*            Superpose two molecules and find the best fit RMSD			*/
/* ======================================================================== */
#ifdef RMSD
#include "defines.h"


/*
 The program to superimpose atoms of two molecules by quaternion method
 This version is adapted from the original source code by

 David J. Heisterberg
 The Ohio Supercomputer Center
 1224 Kinnear Rd.
 Columbus, OH  43212-1163
 (614)292-6036
 djh@osc.edu    djh@ohstpy.bitnet    ohstpy::djh

 Translated to C from fitest.f program and interfaced with Xmol program
 by Jan Labanowski,  jkl@osc.edu   jkl@ohstpy.bitnet   ohstpy::jkl

 Copyright: Ohio Supercomputer Center, David J. Heisterberg, 1990.
 The program can be copied and distributed freely, provided that
 this copyright in not removed. You may acknowledge the use of the
 program in published material as:
 David J. Heisterberg, 1990, unpublished results.

*/

#include <ctype.h>

#define MAXLINELEN    250
  

/*====================================================================
CENTER
 center or translate a molecule. 
 n - number of atoms
 x - on input  - original xyz coordinates of a molecule
     on output - moved xyz coordinates (see io for modes).

 w - if io=1, weights of atoms
     if io=2 or 3, unused

 io - 1 weighted geometric center of the molecule will be at (0,0,0)
      2 molecule will be moved by a vector -o (i.e., components of a vector o
        will be subtracted from atom coordinates). 
      3 molecule will be moved by a vector +o (i.e., components of a vector o
        will be added atom coordinates). 

 o - if io=1, output, center of original coordinates
     if io=2, input, vector o will be subtracted from atomic coordinates
     if io=3, input, vector o will be added to atomic coordinates

=====================================================================*/

void center(int n, double x[4][MAXPOINTS], double w[MAXPOINTS], int io, double o[4])
{
 double wnorm, modif;
 int i;

 if (io == 2) {
   modif = -1.0;
   }
 else if (io == 3) {
   modif = 1.0;
   }
 else {
   modif = -1.0;
   o[1] = 0.0;
   o[2] = 0.0;
   o[3] = 0.0;
   wnorm = 0.0;
   for (i = 1; i <= n; i++) {
     o[1] = o[1] + x[1][i] * sqrt(w[i]);
     o[2] = o[2] + x[2][i] * sqrt(w[i]);
     o[3] = o[3] + x[3][i] * sqrt(w[i]);
     wnorm = wnorm + sqrt(w[i]);
     }
   o[1] = o[1] / wnorm;
   o[2] = o[2] / wnorm;
   o[3] = o[3] / wnorm;
   }


 for (i = 1; i <= n; i++) {
   x[1][i] = x[1][i] + modif*o[1];
   x[2][i] = x[2][i] + modif*o[2];
   x[3][i] = x[3][i] + modif*o[3];
   }

}

/*================================
 ROTMOL
 rotate a molecule
 n - number of atoms
 x - input coordinates
 y - rotated coordinates y = u * x
 u - left rotation matrix
==================================*/

void rotmol (int n, double x[4][MAXPOINTS], double y[4][MAXPOINTS], double u[4][4])
{
 double yx, yy, yz;
 int i;

 for (i = 1; i <= n; i++) {
   yx = u[1][1] * x[1][i] + u[1][2] * x[2][i] + u[1][3] * x[3][i];
   yy = u[2][1] * x[1][i] + u[2][2] * x[2][i] + u[2][3] * x[3][i];
   yz = u[3][1] * x[1][i] + u[3][2] * x[2][i] + u[3][3] * x[3][i];

   y[1][i] = yx;
   y[2][i] = yy;
   y[3][i] = yz;
   }
}

/*=======================================================
 JACOBI
 Jacobi diagonalizer with sorted output. It is only good for 4x4 matrices.
 a - input: matrix to diagonalize
 v - output: eigenvectors
 d - output: eigenvalues
 nrot - input: maximum number of sweeps
=========================================================*/

void jacobi (double a[4][4], double d[4], double v[4][4], int nrot)
{
 double onorm, dnorm;
 double b, dma, q, t, c, s;
 double atemp, vtemp, dtemp;
 int i, j, k, l;

 for (j = 0; j <= 3; j++) {
   for (i = 0; i <= 3; i++) {
     v[i][j] = 0.0;
     }
   v[j][j] = 1.0;
   d[j] = a[j][j];
   }

 for (l = 1; l <= nrot; l++) {
   dnorm = 0.0;
   onorm = 0.0;
   for (j = 0; j <= 3; j++) {
     dnorm = dnorm + fabs(d[j]);
     for (i = 0; i <= j - 1; i++) {
       onorm = onorm + fabs(a[i][j]);
       }
     }
   if((onorm/dnorm) <= 1.0e-12) goto Exit_now;
   for (j = 1; j <= 3; j++) {
     for (i = 0; i <= j - 1; i++) {
       b = a[i][j];
       if(fabs(b) > 0.0) {
         dma = d[j] - d[i];
         if((fabs(dma) + fabs(b)) <=  fabs(dma)) {
           t = b / dma;
           }
         else {
           q = 0.5 * dma / b;
           t = 1.0/(fabs(q) + sqrt(1.0+q*q));
           if(q < 0.0) {
             t = -t;
             }
           }
         c = 1.0/sqrt(t * t + 1.0);
         s = t * c;
         a[i][j] = 0.0;
         for (k = 0; k <= i-1; k++) {
           atemp = c * a[k][i] - s * a[k][j];
           a[k][j] = s * a[k][i] + c * a[k][j];
           a[k][i] = atemp;
           }
         for (k = i+1; k <= j-1; k++) {
           atemp = c * a[i][k] - s * a[k][j];
           a[k][j] = s * a[i][k] + c * a[k][j];
           a[i][k] = atemp;
           }
         for (k = j+1; k <= 3; k++) {
           atemp = c * a[i][k] - s * a[j][k];
           a[j][k] = s * a[i][k] + c * a[j][k];
           a[i][k] = atemp;
           }
         for (k = 0; k <= 3; k++) {
           vtemp = c * v[k][i] - s * v[k][j];
           v[k][j] = s * v[k][i] + c * v[k][j];
           v[k][i] = vtemp;
           }
         dtemp = c*c*d[i] + s*s*d[j] - 2.0*c*s*b;
         d[j] = s*s*d[i] + c*c*d[j] +  2.0*c*s*b;
         d[i] = dtemp;
         }  /* end if */
       } /* end for i */
     } /* end for j */
   } /* end for l */
 
Exit_now:

 nrot = l;

 for (j = 0; j <= 2; j++) {
   k = j;
   dtemp = d[k];
   for (i = j+1; i <= 3; i++) {
     if(d[i] < dtemp) {
       k = i;
       dtemp = d[k];
       }
     }

   if(k > j) {
     d[k] = d[j];
     d[j] = dtemp;
     for (i = 0; i <= 3; i++) {
       dtemp = v[i][k];
       v[i][k] = v[i][j];
       v[i][j] = dtemp;
       }
     }
   }
}



/*==========================================
 Q2MAT
 Generate a left rotation matrix from a normalized quaternion

 INPUT
   q      - normalized quaternion

 OUTPUT
   u      - the rotation matrix
===========================================*/

void q2mat (double q[4], double u[4][4])
{
 u[1][1] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
 u[2][1] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
 u[3][1] = 2.0 * (q[1] * q[3] + q[0] * q[2]);

 u[1][2] = 2.0 * (q[2] * q[1] + q[0] * q[3]);
 u[2][2] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
 u[3][2] = 2.0 * (q[2] * q[3] - q[0] * q[1]);

 u[1][3] = 2.0 *(q[3] * q[1] - q[0] * q[2]);
 u[2][3] = 2.0 * (q[3] * q[2] + q[0] * q[1]);
 u[3][3] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
}


/*==========================================
 QTRFIT
 Find the quaternion, q,[and left rotation matrix, u] that minimizes

   |qTXq - Y| ^ 2  [|uX - Y| ^ 2]

 This is equivalent to maximizing Re (qTXTqY).

 This is equivalent to finding the largest eigenvalue and corresponding
 eigenvector of the matrix

 [A2   AUx  AUy  AUz ]
 [AUx  Ux2  UxUy UzUx]
 [AUy  UxUy Uy2  UyUz]
 [AUz  UzUx UyUz Uz2 ]

 where

   A2   = Xx Yx + Xy Yy + Xz Yz
   Ux2  = Xx Yx - Xy Yy - Xz Yz
   Uy2  = Xy Yy - Xz Yz - Xx Yx
   Uz2  = Xz Yz - Xx Yx - Xy Yy
   AUx  = Xz Yy - Xy Yz
   AUy  = Xx Yz - Xz Yx
   AUz  = Xy Yx - Xx Yy
   UxUy = Xx Yy + Xy Yx
   UyUz = Xy Yz + Xz Yy
   UzUx = Xz Yx + Xx Yz

 The left rotation matrix, u, is obtained from q by

   u = qT1q

 INPUT
   n      - number of points
   x      - fitted molecule coordinates
   y      - reference molecule coordinates
   w      - weights

 OUTPUT
   q      - the best-fit quaternion
   u      - the best-fit left rotation matrix
   nr     - max number of jacobi sweeps

=====================================*/

void qtrfit (int n, double x[4][MAXPOINTS], double y[4][MAXPOINTS], double w[MAXPOINTS], double q[4], double u[4][4], int nr)
{
 double xxyx, xxyy, xxyz;
 double xyyx, xyyy, xyyz;
 double xzyx, xzyy, xzyz;
 double c[4][4], v[4][4];
 double d[4];
 int i, j;


/* generate the upper triangle of the quadratic form matrix */

 xxyx = 0.0;
 xxyy = 0.0;
 xxyz = 0.0;
 xyyx = 0.0;
 xyyy = 0.0;
 xyyz = 0.0;
 xzyx = 0.0;
 xzyy = 0.0;
 xzyz = 0.0;
 
 for (i = 1; i <= n; i++) {
   xxyx = xxyx + x[1][i] * y[1][i] * w[i];
   xxyy = xxyy + x[1][i] * y[2][i] * w[i];
   xxyz = xxyz + x[1][i] * y[3][i] * w[i];
   xyyx = xyyx + x[2][i] * y[1][i] * w[i];
   xyyy = xyyy + x[2][i] * y[2][i] * w[i];
   xyyz = xyyz + x[2][i] * y[3][i] * w[i];
   xzyx = xzyx + x[3][i] * y[1][i] * w[i];
   xzyy = xzyy + x[3][i] * y[2][i] * w[i];
   xzyz = xzyz + x[3][i] * y[3][i] * w[i];
   }
 
 for(i = 0; i <= 3; i++) {
   for(j = 0; j <= 3; j++) {
      c[i][j] = 0.0;
      }
   }

 c[0][0] = xxyx + xyyy + xzyz;

 c[0][1] = xzyy - xyyz;
 c[1][1] = xxyx - xyyy - xzyz;

 c[0][2] = xxyz - xzyx;
 c[1][2] = xxyy + xyyx;
 c[2][2] = xyyy - xzyz - xxyx;

 c[0][3] = xyyx - xxyy;
 c[1][3] = xzyx + xxyz;
 c[2][3] = xyyz + xzyy;
 c[3][3] = xzyz - xxyx - xyyy;

/* diagonalize c */

 jacobi (c, d, v, nr);

/* extract the desired quaternion */

 q[0] = v[0][3];
 q[1] = v[1][3];
 q[2] = v[2][3];
 q[3] = v[3][3];

/* generate the rotation matrix */

 q2mat (q, u);

}

/*=======================================================*/


/* FITEST
 rigid fit test driver
 reads in data, fits, and writes out
 k	= box number
 mode_file	=	0	if rmsd is to be computed on the fly
 mode_file	=	1	if rmsd is to be computed between given files ref.pdb and fit.pdb
 mode_file	=	2	if rmsd is to be computed from./INPUT/ref.pdb & current coords 
 mode_rmsd	=	0	if rmsd is based over all atoms
 mode_rmsd	=	1	if rmsd is based over all non hydrogen atoms
 mode_rmsd	=	2	if rmsd is based over all C-alpha atoms only
 mode_rmsd	=	3	if rmsd is based over all backbone atoms
 mode_excl	=	0	if no residues excluded from rmsd calculations
 mode_excl	=	1	if residues are excluded from rmsd calculations
*/
/*=======================================================*/

int quatfit(int k, int mode_file, int mode_rmsd)
{
 int nat_r;                     /* number of all atoms in reference molecule */
 int nat_f;                     /* number of all atoms in fitted molecule	 */
 int npairs;                    /* no of fitted atom pairs */
 int atoms_r[MAXPOINTS];        /* atoms of ref. molecule to be superimposed */
 int atoms_f[MAXPOINTS];        /* atoms of fit. molecule to be superimposed */
// int res_excl[MAXPOINTS];		/* This will read in the residues number of the excluded residues*/
// int count_excl=0;				/* counter to keep count of number of residues to be excluded */
 int max_sweeps;                /* max number of iterations in jacobi */
 
 double xyz_r[4][MAXPOINTS];    /* coordinates for reference molecule */
 double xyz_f[4][MAXPOINTS];    /* coordinates for fitted molecule */
 double ref_xyz[4][MAXPOINTS];  /* ref. molecule atom coordinates to fit */
 double fit_xyz[4][MAXPOINTS];  /* fit. molecule atom coordinates to fit */
 double weight[MAXPOINTS];      /* fitted atom pair weights */
 double ref_center[4];          /* center of ref. molecule fitted atoms */
 double fit_center[4];          /* center of ref. molecule fitted atoms */
 double q[4];                   /* quaternion */
 double u[4][4];                /* left rotation matrix for coordinates */
 double s, d, wd, rms, wnorm;   /* aux variables */
 
// char name2[50];				/* string to read in the name of input files */
 char name3[50];				/* string to read in the name of input files */
 char name4[50];				/* string to read in the name of input files */
 char tt[80];
 char temp[5];
 char temp1[5]="END";
 char temp2[5]="REMA";
char amino_name[47][4]={"ABC","GLY","ALA","SER","CYS","VAL","THR","ILE","PRO","MET","ASP","ASN","LEU","LYS","GLU",
		"GLN","ARG","HIS","PHE","TYR","TRP","CYT","GUA","ADE","THY","HOH","SOD","CLA","LIP","CHO","EAM","POT","TRE", "GOL","F00", "F01",
		"F02","F03","F04","F05","F06","F07","F08","F09","F10","F11","MEO"};

struct atom_rmsd fit[MAXPOINTS];

 
 /* file with fit goodness values */
/*
  FILE *statfile;                
#ifdef MPI
	  sprintf(name2,"./OUTPUT/BOX%d/rmsd%d.output",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	  sprintf(name2,"./OUTPUT/BOX%d/rmsd%d.output",k,k);
#endif
	statfile= fopen(name2,"a");
*/
 /* set defaults */
 max_sweeps = 30;
 

 /*==========================================================================*/
/* read in exclusions if exclusion mode is on */
/*
 if (mode_excl==0){
	 for(int i=0; i<MAXPOINTS; i++){
		 res_excl[i]=0;
	 }
 }
 else if (mode_excl ==1){
	 FILE *exc1;

  if( NULL == (exc1=fopen("./INPUT/rmsd_excl.inp","r")) ) {
    fprintf(stdout,"input file rmsd_excl.inp does not exist!!!\n");
	
    exit(1);
  }  
     while (!feof(exc1)){
		 count_excl++;
		 fscanf(exc1, "%d\n",&res_excl[count_excl]);
	 }
	 fclose(exc1);
 }
*/
/*==========================================================================*/

 
 nat_r=0;
 if(box[k].boxns>= MAXPOINTS) {
   fprintf(stderr,  
   "Error: Molecule too big. Recompile program with larger MAXPOINTS (at least %i)\n", box[k].boxns+1);
   exit(8);
   }

 /* read coordinates of ref molecule */
if( mode_file==0){						// read from within simulations
	 FILE *das;
#ifdef MPI
	  sprintf(name3,"./OUTPUT/BOX%d/simulinitial%d.pdb",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	  sprintf(name3,"./OUTPUT/BOX%d/simulinitial%d.pdb",k,k);
#endif
	  if( NULL == (das=fopen(name3,"r")) ) {			// open reference file
		fprintf(stdout,"input file simulinitial.pdb (%s) does not exist!!!\n", name3);
		
		exit(10);
	  }  	  	
	  	  int i =0;	
	  
	  while (!feof(das)){
		  fscanf(das, "%4s",   temp);
		  if (strcmp(temp,temp1)==0 || strcmp(temp,temp2)==0 ) {
		  fgets(tt,80,das);
			  continue;
		  }
		  fscanf(das, "%d",   &i); nat_r =i;
		  if(nat_r >= MAXPOINTS) {
			fprintf(stderr,"Error: Molecule too big. Recompile program with larger MAXPOINTS (at least %i)\n", nat_r+1);
			exit(12);
		  }
		  fscanf(das, "%4s",   ref[i].name);
		  fscanf(das, "%3s",   ref[i].res);
		  fscanf(das, "%d",   &ref[i].nres);
		  fscanf(das,"%lf",  &ref[i].x);xyz_r[1][i]=ref[i].x;
		  fscanf(das,"%lf",  &ref[i].y);xyz_r[2][i]=ref[i].y;
		  fscanf(das,"%lf",  &ref[i].z);xyz_r[3][i]=ref[i].z;fgets(tt,80,das);
	  }
	  fclose(das);
	  i =0;	
	  FILE *das2;
#ifdef MPI
	  sprintf(name4,"./OUTPUT/BOX%d/simulfinal%d.pdb",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	  sprintf(name4,"./OUTPUT/BOX%d/simulfinal%d.pdb",k,k);
#endif
	  if( NULL == (das2=fopen(name4,"r")) ) {			// open fitted file
		fprintf(stdout,"input file simulfinal.pdb (%s) does not exist!!!\n", name4);
		
		exit(11);
	  } 
	  while (!feof(das2)){
		  fscanf(das2, "%4s",   temp);
		  if (strcmp(temp,temp1)==0 || strcmp(temp,temp2)==0 ) {
		  fgets(tt,80,das);
			  continue;
		  }
		  fscanf(das2, "%d",   &i); nat_f=i;
		  if(nat_f >= MAXPOINTS) {
			fprintf(stderr,"Error: Molecule too big. Recompile program with larger MAXPOINTS (at least %i)\n", nat_f+1);
			exit(13);
		  }
		  fscanf(das2, "%4s",   fit[i].name);
		  if (strcmp(fit[i].name,ref[i].name)!=0){
			  fprintf(stdout,"ERROR: atom names do not match in ref (%s) and fit (%s) pdb files at site number %d!!!\n",ref[i].name, fit[i].name, i);
			  exit(14);
		  }
		  fscanf(das2,"%3s",  fit[i].res);
		  fscanf(das2,"%d",   &fit[i].nres);
		  fscanf(das2,"%lf",  &fit[i].x);	xyz_f[1][i]=fit[i].x;
		  fscanf(das2,"%lf",  &fit[i].y);	xyz_f[2][i]=fit[i].y;
		  fscanf(das2,"%lf",  &fit[i].z);	xyz_f[3][i]=fit[i].z;
		  fgets(tt,80,das2);
	  }
	  
	  fclose(das2);
	  if (nat_r !=nat_f){
			  fprintf(stdout,"ERROR: atom numbers do not match in ref (%i) and fit (%i) pdb files !!!\n", nat_r, nat_f);
			  exit(15);
	  }
}// mode_file==0

else if (mode_file ==1){				// read from ref.pdb and fit.pdb in INPUT folder   
  
	  FILE *das;

	  if( NULL == (das=fopen("./INPUT/ref.pdb","r")) ) {			// open reference file
		fprintf(stdout,"input file ref.pdb does not exist!!!\n");
		exit(10);
	  } 
	  FILE *das2;

	  if( NULL == (das2=fopen("./INPUT/fit.pdb","r")) ) {			// open fitted file
		fprintf(stdout,"input file fit.pdb does not exist!!!\n");
		exit(11);
	  }  
	  int i =0;	
	  
	  while (!feof(das)){
		  fscanf(das, "%4s",   temp);
		  if (strcmp(temp,temp1)==0 || strcmp(temp,temp2)==0 ) {
		  fgets(tt,80,das);
			  continue;
		  }
		  fscanf(das, "%d",   &i); nat_r =i;
		  if(nat_r >= MAXPOINTS) {
			fprintf(stderr,"Error: Molecule too big. Recompile program with larger MAXPOINTS (at least %i)\n", nat_r+1);
			exit(12);
		  }
		  fscanf(das,"%4s",  ref[i].name);
		  fscanf(das,"%3s",  ref[i].res);
		  fscanf(das,"%d",   &ref[i].nres);
		  fscanf(das,"%lf",  &ref[i].x);	xyz_r[1][i]=ref[i].x;
		  fscanf(das,"%lf",  &ref[i].y);	xyz_r[2][i]=ref[i].y;
		  fscanf(das,"%lf",  &ref[i].z);	xyz_r[3][i]=ref[i].z;
		  fgets(tt,80,das);
	  }

	  i =0;	  
	  while (!feof(das2)){
		  fscanf(das2, "%4s",   temp);
		  if (strcmp(temp,temp1)==0 || strcmp(temp,temp2)==0 ) {
		  fgets(tt,80,das);
			  continue;
		  }
		  fscanf(das2, "%d",   &i); nat_f=i;
		  if(nat_f >= MAXPOINTS) {
			fprintf(stderr,"Error: Molecule too big. Recompile program with larger MAXPOINTS (at least %i)\n", nat_f+1);
			exit(13);
		  }
		  fscanf(das2, "%4s",   fit[i].name);
		  if (strcmp(fit[i].name,ref[i].name)!=0){
			  fprintf(stdout,"ERROR: atom names do not match in ref (%s) and fit (%s) pdb files at site number %d!!!\n",ref[i].name, fit[i].name, i);
			  exit(14);
		  }
		  fscanf(das2,"%3s",  fit[i].res);
		  fscanf(das2,"%d",   &fit[i].nres);
		  fscanf(das2,"%lf",  &fit[i].x);	xyz_f[1][i]=fit[i].x;
		  fscanf(das2,"%lf",  &fit[i].y);	xyz_f[2][i]=fit[i].y;
		  fscanf(das2,"%lf",  &fit[i].z);	xyz_f[3][i]=fit[i].z;
		  fgets(tt,80,das2);
	  }
	  fclose(das);
	  fclose(das2);
	  if (nat_r !=nat_f){
			  fprintf(stdout,"ERROR: atom numbers do not match in ref (%i) and fit (%i) pdb files !!!\n", nat_r, nat_f);
			  exit(15);
	  }
}//mode_file==1
else if (mode_file ==2){				// read from ref.pdb in INPUT folder and compare with current coords in atnopbc 
  

	  for(int ii=0; ii<MAXPOINTS; ii++){
		xyz_r[1][ii]=ref[ii].x;
		xyz_r[2][ii]=ref[ii].y;
		xyz_r[3][ii]=ref[ii].z;
	  }
	  nat_r = refN;
	int rescount =0;double pdb =1.00;
    int count = 0;
    for(int m=0; m<sim.NC; m++) {
      for(int i=0; i<bp[k][m].nbox; i++) {
		  for (int kk=0; kk<mol[m].Nres; kk++) {
			for(int j=0; j<residue[k][rescount].Nsite; j++) {
				nat_f=count+1;
				if(nat_f >= MAXPOINTS) {
					fprintf(stderr,"Error: Molecule too big. Recompile program with larger MAXPOINTS (at least %i)\n", nat_f+1);
					exit(13);
				}
				strcpy(fit[count+1].name,atom[k][count].name);

				if (strcmp(fit[count+1].name,ref[count+1].name)!=0){
					  fprintf(stdout,"ERROR: atom names do not match in ref (%s) and fit (%s) pdb files at site number %d!!!\n",ref[count+1].name, fit[count+1].name, count+1);
					  exit(14);
				}
				strcpy(fit[count+1].res	,amino_name[residue[k][rescount].type]);
				fit[count+1].nres	= rescount+1;
				fit[count+1].x		= atnopbc[k][count].x;	xyz_f[1][count+1]=fit[count+1].x;
				fit[count+1].y		= atnopbc[k][count].y;	xyz_f[2][count+1]=fit[count+1].y;
				fit[count+1].z		= atnopbc[k][count].z;	xyz_f[3][count+1]=fit[count+1].z;	  	  
				count++;
			}
			rescount++;
		  }
      }// loop i ends
	}// loop m ends
//	fclose(das);
	  if (nat_r !=nat_f){
			  fprintf(stdout,"ERROR: atom numbers do not match in ref (%i) and fit (%i) pdb files !!!\n", nat_r, nat_f);
			  exit(15);
	  }
}//mode_file==2

/*==========================================================================*/
/* Having read the reference and fitted molecule apply exclusions */

 
 npairs=0;
 for(int i=1; i<=nat_r; i++){
	 int flag=0;
	 for(int j=1; j<=count_excl; j++){
		 if(ref[i].nres == res_excl[j]){
			 flag=1;
			 break;
		 }
	 }// loop j ends here

	 if (flag==1) continue;
	 else if (flag==0){
			 if(mode_rmsd == 0) {  
				 npairs++;
				 atoms_r[npairs] = i;
				 atoms_f[npairs] = i;
				 weight[npairs]  = 1.0;
			 }
			 else if (mode_rmsd ==1){
				 if(ref[i].name[0]=='H') continue; 
				 npairs++;
				 atoms_r[npairs] = i;
				 atoms_f[npairs] = i;
				 weight[npairs]  = 1.0;
			 }
			 else if (mode_rmsd ==2){
				 if(strcmp(ref[i].name,"CA")!=0) continue;
				 npairs++;
				 atoms_r[npairs] = i;
				 atoms_f[npairs] = i;
				 weight[npairs]  = 1.0;
			 }
			 else if (mode_rmsd ==3){
				 if   ((strcmp(ref[i].name,"C")==0) ||(strcmp(ref[i].name,"N")==0 )||(strcmp(ref[i].name,"CA")==0)) {
				 npairs++;
				 atoms_r[npairs] = i;
				 atoms_f[npairs] = i;
				 weight[npairs]  = 1.0;
				 }
			 }
	 }// flag ==0 ends
 }// loop i ends here

/*=======================================================*/

 /* extract fitted atoms to tables */
 for (int i = 1; i <= npairs; i++) {
   for (int j = 1; j <= 3; j++) {
     ref_xyz[j][i] = xyz_r[j][atoms_r[i]];
     fit_xyz[j][i] = xyz_f[j][atoms_f[i]];
     }
   }

 /* ===  Atom coordinates are fit in both modes === */
 /* center ref molecule fitted atoms around (0,0,0) */
 center (npairs, ref_xyz, weight, 1, ref_center);

 /* center fitted molecule fitted atoms around (0,0,0) */
 center (npairs, fit_xyz, weight, 1, fit_center);

 /* fit specified atoms of fit_molecule to those of ref_molecule */
 qtrfit(npairs, fit_xyz, ref_xyz, weight, q, u, max_sweeps);

 /* subtract coordinates of the center of fitted atoms of the fitted molecule
    from all atom coordinates of the fitted molecule (note that weight is
    a dummy parameter) */
 center(nat_f, xyz_f, weight, 2, fit_center);

 /* rotate the fitted molecule by the rotation matrix u */
 rotmol(nat_f, xyz_f, xyz_f, u);
 /* same with set of fitted atoms of the fitted molecule */
 rotmol(npairs, fit_xyz, fit_xyz, u);

 
 /* translate atoms of the fitted molecule to the center
      of fitted atoms of the reference molecule */
 center(nat_f, xyz_f, weight, 3, ref_center);
 /* same with set of fitted atoms of the fitted molecule */
 center(npairs, fit_xyz, weight, 3, ref_center);
 /* translate fitted atoms of reference molecule to their orig. location */
 center(npairs, ref_xyz, weight, 3, ref_center);

// fprintf(statfile,"\nDistances and weighted distances between fitted atoms\n");
// fprintf(statfile,"Ref.At. Fit.At.  Distance  Dist*sqrt(weight)  weight\n");

 rms = 0.0;
 wnorm = 0.0;
 for (int i = 1; i <= npairs; i++) {
   d = 0.0;
   for (int j = 1; j <= 3; j++) {
     s = ref_xyz[j][i] - fit_xyz[j][i];
     d += s*s;
     }
   d = sqrt(d);
   s = sqrt(weight[i]);
   wd = s*d;
   rms += wd*wd;
   wnorm += s;
//   fprintf(statfile, "  %3d    %3d  %11.6f   %11.6f    %11.6f\n",
//  atoms_r[i], atoms_f[i], d, wd, weight[i]);
   }
  
 rms = sqrt(rms/wnorm);
 ordparam[k].rmsd = rms;
/*
 fprintf(statfile, "\n\nWeighted root mean square=%10.6f\n\n", rms);
 fprintf(statfile, "\n\nCenter of reference molecule fitted atoms\n");
 fprintf(statfile, "Xc = %11.6f Yc = %11.6f Zc = %11.6f\n",
         ref_center[1], ref_center[2], ref_center[3]);
  
 fprintf(statfile, "\n\nCenter of fitted molecule fitted atoms\n");
 fprintf(statfile, "Xc = %11.6f Yc = %11.6f Zc = %11.6f\n",
         fit_center[1], fit_center[2], fit_center[3]);
   
 fprintf(statfile,"\n\nLeft rotation matrix\n");
 for (int i = 1; i <= 3; i++) {
   fprintf(statfile, " %11.6f  %11.6f  %11.6f\n", 
           u[1][i], u[2][i], u[3][i]);

  }
*/
// fclose(statfile);
 return(0);
 }
#endif
