/* ======================================================================== */
/* mcmove_disp.cpp                                                          */
/*     This move is a typical atom displacement move.  It is called with    */
/* frequency mc_rand[k].  It is similiar to the xedos_linear move and the   */
/* (ct)dos_trans move.  In the xedos case, the move is also called with     */
/* frequency mc_rand[k] but in the (ct)dos case, it is called with          */
/* mc_trans[k].                                                             */
/* Written by Thomas Knotts, 1 Sep 2005.                                    */ 
/* ======================================================================== */

#include "defines.h"

void   svalues        (int);
void   pbc_all        (int);
void   calcvalue      (int);
void   nblist_pivot   (int,int);
double ran2           (void);
int	   ran_int	      (int,int);     
void   mcmove_bkp     (int,int,int);
void   mcmove_rstr    (int,int,int);
#ifdef MC
void   forces_mc      (int,int,int,int);
#else
void   forces		      (int);
#endif

#ifdef DLIST
void   nl_check       (int,int,int);
#endif



void mcmove_displ(int ibox)
{
	int k	=	ibox;

	/* ---------------------------------------------------- */
	/* First, pick a random solute molecule to move along   */
	/* the reaction coordinate.  Then determine the upper   */
	/* and lower bounds for the sites to move.              */
	/* ---------------------------------------------------- */
	int lb;
	int ub;
	int rand_site =0;
	if(molec.Nsolute == 1){		
		lb = 0;
		ub = molec.lsite[0];
		rand_site		=   ran_int(lb,ub);
	}
	else{
		int molecule = ran_int(0, molec.Nsolute);
		lb = molec.fsite[molecule];//first site of molecule
		ub = molec.lsite[molecule];//first site of next molecule
		rand_site		=   ran_int(lb,ub);
	}

  /* ====================================================== */
  /* Back up the necessary information so that the system   */
  /* can be restored if the move is rejected.               */
  /* ====================================================== */
  mcmove_bkp(rand_site,rand_site+1,k);

  /* ====================================================== */
  /* Set the pre-move variables.                            */
  /* ====================================================== */
	double U_initial;
	double U_final;

  #ifdef MC
	forces_mc(rand_site,rand_site+1,0,k);
  #endif
	calcvalue(k);
	
	double beta	= 1.0/sim.kT[k];			

  #ifdef MC		
	U_initial		= enmc[k].o_potens;
  #else		
	U_initial		= en[k].potens;
  #endif

  /* ====================================================== */
  /* Now propose the move.                                  */
  /* ====================================================== */

	/* ------------------------------------------------------ */
	/* Randomly displace the atom along each direction.       */
	/* ------------------------------------------------------ */
  double dx, dy,dz;
  dx = (ran2()-0.5)*mc_rand.delta[k]; 
  dy = (ran2()-0.5)*mc_rand.delta[k];
  dz = (ran2()-0.5)*mc_rand.delta[k];

  atnopbc[k][rand_site].x += dx;
  atnopbc[k][rand_site].y += dy;
  atnopbc[k][rand_site].z += dz;
  atom[k][rand_site].x += dx;
  atom[k][rand_site].y += dy;
  atom[k][rand_site].z += dz;


	/* ----------------------------------------------- */
	/* Apply periodic boundary conditions.             */
	/* ----------------------------------------------- */
	pbc_all(k);
	#ifndef DLIST
	nblist_pivot(k,rand_site+1); // call Nblist 
	#else
	nl_check(lb,ub,k);
	if(nl_flag[k] == 1) nblist_pivot(k,rand_site+1);
	#endif

	/* ----------------------------------------------- */
	/* Update the forces.                              */
	/* ----------------------------------------------- */
  #ifdef MC
	forces_mc(rand_site, rand_site+1,1,k);
  #else
	forces(k);
  #endif

	calcvalue(k);

  #ifdef MC
	U_final	= enmc[k].n_potens;
  #else
	U_final	= en[k].potens;
  #endif


  /* ====================================================== */
  /* Check to see if the move is accepted.                  */
  /* ====================================================== */
  double arg = (U_initial-U_final)*beta;
  if (arg>0) mc_rand.rand_acc[k]++; //accepted
	else if(exp(arg) > ran2()) mc_rand.rand_acc[k]++; //accepted
  else mcmove_rstr(rand_site,rand_site+1,k); //rejected
  
}

