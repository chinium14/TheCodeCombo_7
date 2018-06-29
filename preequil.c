/* ======================================================================== */
/* preequil.cpp                                                             */
/*                                                                          */
/*		This subroutine is the main loop for the preequilibration of the      */
/* system.  If no velocity input file is found, the system is ramped from   */
/* T=0.001 K to the input temperature in 10,000 steps.  It is done using    */
/* the isokinetice ensemble.                                                */
/*                                                                          */
/* Written by Thomas A. Knotts on 28 Apr 05.                                */
/* ======================================================================== */

#include "defines.h"
  /* ================================================================== */
  /* Function prototypes                                                */
  /* ================================================================== */
void pbc_all(int);
void forces(int);
void force_long(int);
void force_short(int);
void isokin(int);
void kinet(int);
void drift(int);
void nblist(int);
void vinit();
void integrate_vverlet(int, double);
void calcvalue(int);
double velscale(int);
int ran_int(int, int);
#ifdef TLATE
void translate(int);
#endif
#ifdef DLIST
void nl_check(int, int, int);
#endif

void preequil(){
  int Nsteps = 10000;
  /* ----------------------------------------------- */
  /* Back up the simulation temperatures and         */
  /* calculate the correct temperature increment for */
  /* each time step.                                 */
  /* ----------------------------------------------- */
  double *simT,*dT;
  simT = (double*)calloc(sim.NB,sizeof(double));
  dT = (double*)calloc(sim.NB,sizeof(double));

  for(int k=0; k<sim.NB; k++){
    simT[k] = sim.T[k];
    dT[k] = (sim.T[k]-0.001)/(double)Nsteps;
    sim.T[k] = 0.001;
  }

  vinit();

  /* ================================================================== */
  /* Pre-equilibration uses the isokinetic ensemble so vscal[k] is      */
  /* always 1.0 and all scaling is handled in the subroutine isokin().  */
  /* ================================================================== */

  /* ----------------------------------------------- */
  /* Do the pre-equilibration steps.                 */
  /* ----------------------------------------------- */
	for(int i=0; i<Nsteps; i++){
		n_iter = i;
    for(int k=0; k<sim.NB; k++){
  /* ----------------------------------------------- */
  /* Update the neighborlist with frequency          */
  /* sim.nlist.                                      */
  /* ----------------------------------------------- */
			#ifdef NLIST
			#ifndef DLIST
				if(i%sim.nlist==0) nblist(k);
			#endif
			#ifdef DLIST
				nl_check(0,box[k].boxns,k);
				if(nl_flag[k] == 1) nblist(k);
			#endif
			#endif	

			/* ----------------------------------------------- */
			/* Itegration of the Equations of Motion using     */
      /* velocity verlet with a scaling of 1.0.          */
      /* ----------------------------------------------- */
       integrate_vverlet(k, 1.0);
			
			/* ----------------------------------------------- */
			/* Calculate the properties of the system.         */
			/* ----------------------------------------------- */
			kinet(k);    	            
			calcvalue(k);         	     

			/* ----------------------------------------------- */
			/* Correct for momentum drift and translate        */
			/* ----------------------------------------------- */
      if((i%1000) == 0){
        drift(k); 
			  #ifdef TLATE
			  translate(k);
			  #endif
      }

			/* ----------------------------------------------- */
			/* Increment tempererature by sim.dtemp.           */
			/* ----------------------------------------------- */
			sim.T[k] += dT[k]; 
			isokin(k); 

		}//end for k
    

  }//end i



  /* ================================================================== */
  /* Apply pbc's, recalculate forces and properties of the system to    */
  /* prepare for equilibration run.                                     */
  /* ================================================================== */
  for(int k=0; k<sim.NB; k++){
    sim.T[k] = simT[k];
		pbc_all(k);
#ifndef GOLIK
		force_short(k); force_long(k); 
#endif
		forces(k);   
		kinet(k);
		calcvalue(k);  
  }
  fprintf(stdout,"Pre-equilibration completed!!!\n");
}

