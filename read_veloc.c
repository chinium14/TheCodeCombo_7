/* ======================================================================== */
/* read_veloc.cpp                                                           */
/*                                                                          */
/*		This subroutine reads in the velocities of all the sites if the     */
/* file simul.vel# is found in the INPUT folder.  If it is not, the global  */
/* variable vtest will be 1 and main() will call vinit to assign random,    */
/* gaussian velocites to the sites.                                         */
/*                                                                          */
/* Note: There must be a separate simul.vel# file for each box.             */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void read_veloc (int *vtest) 
{

  /* ================================================================== */
  /*                                                                    */
  /* Variables and pointers needed in the subroutine.                   */
  /*                                                                    */
  /* ================================================================== */
  int qw;
  char t2[5], tt[80], nameb[20];
  FILE *das;
  
  *vtest = 1;

  for(int k=0; k<sim.NB; k++) {
  /* ================================================================== */
  /*                                                                    */
  /* Check to see if the file exists. If it does, read in the           */
  /* velocities. If it doesn't, return to main with vtest = 1.          */
  /*                                                                    */
  /* ================================================================== */
#ifdef MPI
	sprintf(nameb,"./INPUT/simul.vel%d",mpi.my_rank);
#endif
#ifndef MPI
	sprintf(nameb,"./INPUT/simul.vel%d",k);
#endif
    if( (das=fopen(nameb,"r")) ) {
      fgets(tt,80,das);
      fgets(tt,80,das);
      fgets(tt,80,das);
      fgets(tt,80,das); 
      int count = 0;
      for(int m=0; m<sim.NC; m++) {
		for(int i=0; i<bp[k][m].nbox; i++) {
			for(int j=0; j<mol[m].Nsite; j++) {
				fscanf(das,"%d",   &qw);
				fscanf(das,"%4s",  t2);
				fscanf(das,"%lf",  &vv[k][count].x);
				fscanf(das,"%lf",  &vv[k][count].y);
				fscanf(das,"%lf\n",&vv[k][count].z);
				count++;
			}//for j
		}//for i
      }// for m
      fclose(das);    
      *vtest = 0;
    }//das=fopen(nameb,"r")

#ifdef PR_NPT
	for(int i=0; i<box[k].boxns; i++){
		double mass = pott[k][i].mas;
		vvab[k][0] = vvab[k][0] + vv[k][i].x * vv[k][i].x * mass;
		vvab[k][1] = vvab[k][1] + vv[k][i].x * vv[k][i].y * mass;
		vvab[k][2] = vvab[k][2] + vv[k][i].x * vv[k][i].z * mass;
		vvab[k][4] = vvab[k][4] + vv[k][i].y * vv[k][i].y * mass;
		vvab[k][5] = vvab[k][5] + vv[k][i].y * vv[k][i].z * mass;
		vvab[k][8] = vvab[k][8] + vv[k][i].z * vv[k][i].z * mass;
	}
  /* ----------------------------------------------- */
  /* Convert gm/mol to kg.                           */
  /* ----------------------------------------------- */
	    vvab[k][0] = vvab[k][0] / 1000.0 / NA;
		vvab[k][1] = vvab[k][1] / 1000.0 / NA;
		vvab[k][2] = vvab[k][2] / 1000.0 / NA;
		vvab[k][4] = vvab[k][4] / 1000.0 / NA;
		vvab[k][5] = vvab[k][5] / 1000.0 / NA;
		vvab[k][8] = vvab[k][8] / 1000.0 / NA;


	//to have a symmetric vvab matrix;
	vvab[k][3] = vvab[k][1];
	vvab[k][6] = vvab[k][2];
	vvab[k][7] = vvab[k][5];

#endif
  /* ================================================================== */
  /*                                                                    */
  /* Set both velocity arrays to the new values.                        */
  /*                                                                    */
  /* ================================================================== */
    for(int i=0; i<box[k].boxns; i++) {
      uu[k][i].x = vv[k][i].x;
      uu[k][i].y = vv[k][i].y;
      uu[k][i].z = vv[k][i].z;
    }
  
  }//for k

}
