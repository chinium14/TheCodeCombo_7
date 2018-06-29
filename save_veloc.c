/* ======================================================================== */
/* save_veloc.cpp                                                           */
/*                                                                          */
/*		This subroutine saves the velocities for each site. It is called    */
/* from main() with a frequency set by sim.blockc.  This frequency is       */
/* set in the simul.input file.  The file output is the same format as the  */
/* simul.vel# file so a simulation can be started from a previous           */
/* simulation.  The file saves the velocities in all of the boxes to        */
/* different files.  The output file is named simul.vel#.sav. To use it to  */
/* start a simulation, simply delete .sav from the name.                    */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						block:		The iteration number                    */
/*                      cycle:		The number of iterations in the         */
/*                                  simulation.                             */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void save_veloc (unsigned long block, unsigned long cycle) 
{
  /* ================================================================== */
  /*                                                                    */
  /* Variables and pointers needed in the subroutine.                   */
  /*                                                                    */
  /* ================================================================== */
  char name[50];
//  char da[6] = "ABCDE";
  time_t now;
  FILE *fas;
  
  for(int k=0; k<sim.NB; k++) {
    time(&now);
  /* ================================================================== */
  /*                                                                    */
  /* Open the file and write the output.                                */
  /*                                                                    */
  /* ================================================================== */
#ifdef MPI
	  sprintf(name,"./OUTPUT/BOX%d/simul.vel%d.sav",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	  sprintf(name,"./OUTPUT/BOX%d/simul.vel%d.sav",k,k);
#endif
    fas = fopen(name,"w");
  /* ----------------------------------------------- */
  /* Write heading and simulation parameters.        */
  /* ----------------------------------------------- */    
    fprintf(fas,"Last saved: %lu blocks %lu of cycles @ %.24s\n",
	    block,cycle,ctime(&now));
    fprintf(fas,"\n");
    fprintf(fas,"velocities of sites [A/ps]\n");
    fprintf(fas,"\n");
//    int count  = 0;
  /* ----------------------------------------------- */
  /* Loop around the number of simulation            */
  /* components, molecule m in box k, and number of  */
  /* sites in molecule m and write velocities.       */
  /* ----------------------------------------------- */    
	for(int i=0; i<box[k].boxns; i++){
			fprintf(fas,"%-4d ",    i);
			fprintf(fas," %2s   ",  atom[k][i].t);
			fprintf(fas," %16.8e ", vv[k][i].x);
			fprintf(fas," %16.8e ", vv[k][i].y);
			fprintf(fas," %16.8e\n",vv[k][i].z);
	}
  /*    for(int m=0; m<sim.NC; m++) {
      for(int i=0; i<bp[k][m].nbox; i++) {
		for(int j=0; j<mol[m].Nsite; j++) {
			fprintf(fas,"%-4d ",    count);
			fprintf(fas," %2s   ",  atom[k][count].t);
			fprintf(fas," %16.8e ", vv[k][count].x);
			fprintf(fas," %16.8e ", vv[k][count].y);
			fprintf(fas," %16.8e\n",vv[k][count].z);
			count++;
		}// for j
		fprintf(fas,"\n");
      }// for i
    }// for m
 */
  fclose(fas);

   }// for k
}
