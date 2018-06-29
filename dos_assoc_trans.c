#ifdef DOS
/* ======================================================================== */
/* dos_assoc_trans.cpp			                                                */
/* This move is for WLDOS.  It is specifically used to move two separate    */
/* molecules, such as ssDNA together or apart.  It does so along the the    */
/* vector defined by SITE1 and SITE2.  This move could be generalized so    */
/* any two solute molecules of interest could be randomly picked and moved  */
/* along the vector defined by their centers of mass.  It is called with    */
/* probability mc_rand[k] from dos.cpp.  mc_rand[k] is read in from         */
/* dos0.input as the second to last move in of the five possibilies.  This  */
/* move is essentialy the same thing as xedos_trans.  A key note:           */
/* xedos_trans is called with the probability mc_trans[k] and xedos_linear  */
/* which is a single atom displacement, is called with mc_rand.  Thus,      */
/* there is some incongruity in the two types of density of states.         */
/* Written by Thomas Knotts, 30 Aug 2005.                                   */ 
/* ======================================================================== */

#include "defines.h"

void   dos_svalues    (int);
void   pbc_all        (int);
void   calcvalue      (int);
void   nblist_pivot   (int,int);
void   nblist         (int);
double ran2           (void);
int	   ran_int	      (int,int);
void   flat_histogram (int);
double dos_interp     (int,int,double);
double calc_dist      (int,int,int);     

#ifdef MC
void   forces_mc      (int,int,int,int);
#else
void   forces		      (int);
#endif

#ifdef STATUS
void   curr_status    (int,int);
#endif

#ifdef DLIST
void   nl_check       (int,int,int);
#endif



void dos_assoc_trans(int ibox)
{
	int k	  =	ibox;

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
	
	int boxns = box[k].boxns;

	/* ====================================================== */
	/* Determine the distance between the sites on the        */
	/* reaction coordinate.                                   */
	/* ====================================================== */
	double dx	= atnopbc[k][SITE2].x -atnopbc[k][SITE1].x;
	double dy	= atnopbc[k][SITE2].y -atnopbc[k][SITE1].y;
	double dz	= atnopbc[k][SITE2].z -atnopbc[k][SITE1].z;

	double r = sqrt(dx*dx+dy*dy+dz*dz);

	/* ====================================================== */
	/* Randomly increase or decrease the reaction coordinate  */
	/* some fraction of the strain value input in dos#.input  */
	/* First the  distance is determined.  Then, one strand   */
	/* is chosen at random to move while the other is chosen  */
	/* to remain where it is.                                 */
	/* ====================================================== */
	double ex_dist = (ran2()-0.5)*2.0 * mc_rand.delta[k];

	double scalar = ex_dist/r;

	dx *= scalar;
	dy *= scalar;
	dz *= scalar;

#ifdef MC
	forces_mc(lb,ub,0,k);
#endif

	double U_initial;
	double U_final;
  double L_initial;
  double L_final;
	
	int old_bin;
	int new_bin;
	double beta_old	= 1.0/sim.kT[k];			
	double beta_new	= beta_old;

	calcvalue(k);
	for(int l=lb; l<ub; l++){
		atom_temp[k][l]	= atom[k][l];				    /* Back up coordinates with pdb			*/
		atnopbc_temp[k][l]	= atnopbc[k][l];		/* Back up coordinates without pdb		*/
	}
	for(int l=0; l<boxns; l++){
		ff_temp[k][l]		= ff[k][l];				      /* Back up all the force  components	*/
	}
#ifdef PRESSURE
		pvir_temp[k]	= pvir[k];					      /* Back up all the virial components	*/
#endif
		en_temp[k]		= en[k];					        /* Back up all the energy components	*/
		
//#ifdef MC		
//		U_initial		= enmc[k].o_potens;
//#else		
		U_initial		= en[k].potens;
//#endif
		L_initial = calc_dist(k,SITE1,SITE2);   /* Distance betweeb SITE1 and SITE2		*/

	  old_bin			= (int) ((U_initial - sim_dos[k].e_begin)/sim_dos[k].e_width);	

	/* Now propose a rescaling move	*/

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
#endif
#ifndef MC
		forces(k);
#endif

		calcvalue(k);
#ifdef MC
		U_final	= en_temp[k].potens - enmc[k].o_potens + enmc[k].n_potens;
#endif
#ifndef MC
		U_final	= en[k].potens;
#endif
		L_final = calc_dist(k,SITE1,SITE2);


		/* ----------------------------------------------- */
		/* Check the acceptance of the move. First check   */
    /* to see if the system moved out of range, then   */
    /* check to see if it is accepted.                 */
		/* ----------------------------------------------- */

		if(U_final >= sim_dos[k].e_begin && U_final < sim_dos[k].e_end){
      new_bin	=  (int) ((U_final - sim_dos[k].e_begin)/sim_dos[k].e_width);
      double g_old = dos_interp(k,old_bin,U_initial);
      double g_new = dos_interp(k,new_bin,U_final);
      double arg = g_old - g_new + 2.0 * log(L_final/L_initial);

		  if (arg>0){
				  dos_hist[k][new_bin].h_of_e ++;
				  dos_hist[k][new_bin].g_of_e += sim_dos[k].mod_f;
				  flat_histogram(k);
				  mc_rand.rand_acc[k] ++;
		  }//accepted
		  else if (exp(arg) > ran2()){
			  dos_hist[k][new_bin].h_of_e ++;
			  dos_hist[k][new_bin].g_of_e += sim_dos[k].mod_f;
			  flat_histogram(k);
			  mc_rand.rand_acc[k] ++;
		  }//accepted
		  else{
			  for(int l=lb; l<ub; l++){
				  atom[k][l]		= atom_temp[k][l];	 /* Back up coordinates with pdb		*/
				  atnopbc[k][l]	= atnopbc_temp[k][l];/* Back up coordinates without pdb		*/
			  }
			  for(int l=0; l< box[k].boxns; l++){
				  ff[k][l]		= ff_temp[k][l];		/* Back up all the force  components	*/
			  }
			  #ifdef PRESSURE
			  pvir[k]		= pvir_temp[k];				/* Back to old  virial components		*/
			  #endif
			  en[k]		= en_temp[k];				/* Back to old  energy components		*/
			  #ifdef NLIST
			  #ifndef DLIST
				  nblist_pivot(k,ub); // call Nblist 
			  #endif
			  #ifdef DLIST
				  nl_check(lb,ub,k);
				  if(nl_flag[k] == 1) nblist_pivot(k,ub);
			  #endif
			  #endif
			  dos_hist[k][old_bin].h_of_e ++;
			  dos_hist[k][old_bin].g_of_e += sim_dos[k].mod_f;
			  flat_histogram(k);
		  }//rejected		
    }
    else{
			  for(int l=lb; l<ub; l++){
				  atom[k][l]		= atom_temp[k][l];	 /* Back up coordinates with pdb		*/
				  atnopbc[k][l]	= atnopbc_temp[k][l];/* Back up coordinates without pdb		*/
			  }
			  for(int l=0; l< box[k].boxns; l++){
				  ff[k][l]		= ff_temp[k][l];		/* Back up all the force  components	*/
			  }
			  #ifdef PRESSURE
			  pvir[k]		= pvir_temp[k];				/* Back to old  virial components		*/
			  #endif
			  en[k]		= en_temp[k];				/* Back to old  energy components		*/
			  #ifdef NLIST
			  #ifndef DLIST
				  nblist_pivot(k,ub); // call Nblist 
			  #endif
			  #ifdef DLIST
				  nl_check(lb,ub,k);
				  if(nl_flag[k] == 1) nblist_pivot(k,ub);
			  #endif
			  #endif
			  dos_hist[k][old_bin].h_of_e ++;
			  dos_hist[k][old_bin].g_of_e += sim_dos[k].mod_f;
			  flat_histogram(k);
    }
    calcvalue(k);
		dos_svalues(k);
}



#endif

