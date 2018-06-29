/* ======================================================================== */
/* mcmove_hmc.cpp                                                           */
/*     This move is a hybrid molecular dynmaics/Monte Carlo move.           */
/* Written by Thomas Knotts, 1 Sep 2005.                                    */ 
/* ======================================================================== */
	
#include "defines.h"
void   vinit_hyb          (int);
void   pbc_all            (int);
void   integrate_vverlet  (int, double);
void   integrate_dos_mts  (int, double);
void   kinet              (int);
void   calcvalue          (int);
void   nblist             (int);
double ran2               (void);
double velscale           (int);
void   mcmove_bkp         (int,int,int);
void   mcmove_rstr        (int,int,int);
void   forces		          (int);
void   force_short        (int);
void   force_long         (int);

#ifdef DLIST
void   nl_check           (int,int,int);
#endif


void mcmove_hmc(int ibox)
{

  int k = ibox;
  int boxns = box[k].boxns;

  /* ====================================================== */
  /* Initialize the velocities.                             */
  /* ====================================================== */
  vinit_hyb(k);
  kinet(k);
  calcvalue(k);

  /* ====================================================== */
  /* Back up the necessary information so that the system   */
  /* can be restored if the move is rejected.               */
  /* ====================================================== */
  mcmove_bkp(0,boxns,k);

  /* ====================================================== */
  /* Set the pre-move variables.                            */
  /* ====================================================== */
  double U_initial;
  double U_final;
  double K_initial;
  double K_final;
  double E_initial;
  double E_final;

  double vscal = 1.0;
  double  beta	= 1.0/sim.kT[k];
  
	U_initial		= en[k].potens;
	K_initial		= en[k].kinet;
  E_initial   = U_initial + K_initial;

  /* ====================================================== */
  /* Now propose the move.                                  */
  /* ====================================================== */

	/* ----------------------------------------------- */
	/* Propagate the system with verlet integrator.    */
	/* ----------------------------------------------- */	
  for (int j=0;j < sim_hyb.cyc_hybmd; j++){

    if(sim.nsteps ==1){
      integrate_vverlet(k, vscal);
      #ifdef DLIST
	      nl_check(0,box[k].boxns,k);
	      if(nl_flag[k] == 1) nblist(k);
      #endif
    }

    else {
      if(j==0){
        force_short(k);
        for(int l=0; l< boxns; l++){
          ff_long[k][l].x		= ff_temp[k][l].x-ff_short[k][l].x;
          ff_long[k][l].y		= ff_temp[k][l].y-ff_short[k][l].y;
          ff_long[k][l].z		= ff_temp[k][l].z-ff_short[k][l].z;
	      }
	    }
      integrate_dos_mts(k,vscal);
	    #ifdef DLIST
        nl_check(0,boxns,k);
        if(nl_flag[k] == 1) nblist(k);
	    #endif

    }
  }
  
  /* ----------------------------------------------- */
  /* Calculate the properties of the system.         */
  /* ----------------------------------------------- */
	kinet(k);    											
	calcvalue(k);
	U_final	= en[k].potens;
	K_final	= en[k].kinet;
	E_final = U_final+K_final;

  /* ====================================================== */
  /* Check to see if the move is accepted.                  */
  /* ====================================================== */
  double arg = (E_initial-E_final)*beta;
  if (arg>0) mc_hmc_acc[k]++; //accepted
	else if(exp(arg) > ran2()) mc_hmc_acc[k]++; //accepted
  else mcmove_rstr(0,boxns,k); //rejected

}
