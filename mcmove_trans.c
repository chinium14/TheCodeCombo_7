/* ======================================================================== */
/* mcmove_trans.cpp                                                         */
/*     This move is is specifically used to move two separate               */
/* molecules, such as ssDNA together or apart.  It does so along the the    */
/* vector defined by SITE1 and SITE2.  This move could be generalized so    */
/* any two solute molecules of interest could be randomly picked and moved  */
/* along the vector defined by their centers of mass.  It is called with    */
/* probability mc_trans[k] from nvt_mc.cpp.  mc_trans[k] is read in from    */
/* mc0.input as the second  move in of the five possibilies.  This          */
/* move is essentialy the same thing as xedos_trans/(ct)dos_assoc_trans.    */
/* A key note: xedos_trans is called with the probability mc_trans[k] and   */
/* xedos_linear which is a single atom displacement, is called with         */
/* mc_rand.  Thus, there is some incongruity in the two types of density of */
/* states.                                                                  */
/* Written by Thomas Knotts, 31 Aug 2005.                                   */ 
/*                                                                          */
/* 6 Jul 09 Made this subroutine so that it moves a random molecule in a    */
/* random direction and not along the reaction coordinate given by SITE1    */
/* and SITE2.  This is more general for multiple molecules.  You can go     */
/* and uncomment the "old" code if you want.                                */
/* ======================================================================== */

#include "defines.h"

void   svalues        (int);
void   pbc_all        (int);
void   calcvalue      (int);
void   nblist_pivot   (int,int);
double ran2           (void);
int	   ran_int	      (int,int);
double calc_dist      (int,int,int);     
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



void mcmove_trans(int ibox)
{
	int k	=	ibox;

	/* ---------------------------------------------------- */
	/* First, pick a random solute molecule to move along   */
	/* the reaction coordinate.  Then determine the upper   */
	/* and lower bounds for the sites to move.              */
	/* ---------------------------------------------------- */
  int lb;
  int ub;
  if(molec.Nsolute == 1){		
    lb = 0;
    ub = molec.lsite[0];
  }
  else{
    int molecule = ran_int(0, molec.Nsolute);
    lb = molec.fsite[molecule];//first site of molecule
    ub = molec.lsite[molecule];//first site of next molecule
  }

  /* ====================================================== */
  /* Back up the necessary information so that the system   */
  /* can be restored if the move is rejected.               */
  /* ====================================================== */
  mcmove_bkp(lb,ub,k);

  /* ====================================================== */
  /* Set the pre-move variables.                            */
  /* ====================================================== */
	double U_initial;
	double U_final;
  double L_initial;
  double L_final;

  #ifdef MC
	forces_mc(lb,ub,0,k);
  #endif
	calcvalue(k);
	
	double beta	= 1.0/sim.kT[k];			

  #ifdef MC		
	U_initial		= enmc[k].o_potens;
  #else		
	U_initial		= en[k].potens;
  #endif

	L_initial = calc_dist(k,SITE1,SITE2); 

  /* ====================================================== */
  /* Now propose the move.                                  */
  /* ====================================================== */

	/* ------------------------------------------------------ */
	/* Determine the distance between the sites on the        */
	/* reaction coordinate.                                   */
	/* ------------------------------------------------------ */
/*	double dx	= atnopbc[k][SITE2].x -atnopbc[k][SITE1].x;
	double dy	= atnopbc[k][SITE2].y -atnopbc[k][SITE1].y;
	double dz	= atnopbc[k][SITE2].z -atnopbc[k][SITE1].z;

	double r = sqrt(dx*dx+dy*dy+dz*dz);
*/
	/* ------------------------------------------------------ */
	/* Randomly increase or decrease the reaction coordinate  */
	/* some fraction of the strain value input in dos#.input  */
	/* First the  distance is determined.  Then, one strand   */
	/* is chosen at random to move while the other is chosen  */
	/* to remain where it is.                                 */
	/* ------------------------------------------------------ */
/*	double ex_dist = (ran2()-0.5)*2.0 * mc_rand.delta[k];

	double scalar = ex_dist/r;

	dx *= scalar;
	dy *= scalar;
	dz *= scalar;
*/

        double dx, dy,dz;
        dx = (ran2()-0.5)*mc_trans.delta[k];
        dy = (ran2()-0.5)*mc_trans.delta[k];
        dz = (ran2()-0.5)*mc_trans.delta[k];

	for (int i=lb; i<ub; i++){
		atnopbc[k][i].x	+=	dx;
		atnopbc[k][i].y	+=	dy;
		atnopbc[k][i].z	+=	dz;
		atom[k][i].x	=	atnopbc[k][i].x;
		atom[k][i].y	=	atnopbc[k][i].y;
		atom[k][i].z	=	atnopbc[k][i].z;
	}

	/* ----------------------------------------------- */
	/* Apply periodic boundary conditions.             */
	/* ----------------------------------------------- */
	pbc_all(k);
	#ifdef NLIST
	#ifndef DLIST
	nblist_pivot(k,ub); // call Nblist 
	#endif
	#ifdef DLIST
	nl_check(lb,ub,k);
	if(nl_flag[k] == 1) nblist_pivot(k,ub);
  #endif
	#endif

	/* ----------------------------------------------- */
	/* Update the forces.                              */
	/* ----------------------------------------------- */
  #ifdef MC
	forces_mc(lb, ub, 1,k);
  #else
	forces(k);
  #endif

	calcvalue(k);

  #ifdef MC
	U_final	= enmc[k].n_potens;
  #else
	U_final	= en[k].potens;
  #endif

	L_final = calc_dist(k,SITE1,SITE2);


  /* ====================================================== */
  /* Check to see if the move is accepted.                  */
  /* ====================================================== */
  double arg = (U_initial-U_final)*beta+2.0 * log(L_final/L_initial);
  if (arg>0) mc_trans.trans_acc[k]++; //accepted
	else if(exp(arg) > ran2()) mc_trans.trans_acc[k]++; //accepted
  else mcmove_rstr(lb,ub,k); //rejected

}
