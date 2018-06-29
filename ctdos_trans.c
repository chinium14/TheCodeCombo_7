#ifdef CTDOS
#include "defines.h"

void dos_svalues    (int);
void pbc_all        (int);
#ifdef MC
void forces_mc      (int,int,int,int);
#endif
#ifndef MC
void forces         (int);
#endif
void calcvalue      (int);
void nblist_pivot   (int, int);
double ran2         (void);
int ran_int         (int,int);
void flat_histogram (int);
double dos_interp   (int,int,double);
#ifdef STATUS
void curr_status    (int,int);
#endif
#ifdef DLIST
void nl_check       (int, int, int);
#endif

void ctdos_trans(int k){

	/* ---------------------------------------------------- */
	/* First, pick a random solute molecule and then a      */
	/* random residue of that molecule to move.             */
	/* ---------------------------------------------------- */
	int res_lb;
	int res_ub;
	if(molec.Nsolute == 1){		
		res_lb = 0;
		res_ub = molec.lres[0];
	}
	else{
		int molecule = ran_int(0, molec.Nsolute);
		res_lb = molec.fres[molecule];//first site of molecule
		res_ub = molec.lres[molecule];//first site of next molecule
	}
	
	int res = ran_int(res_lb, res_ub);	
	
	double dx = (ran2()-0.5)*mc_trans.delta[k];; 
	double dy = (ran2()-0.5)*mc_trans.delta[k];
	double dz = (ran2()-0.5)*mc_trans.delta[k];

	int lb=0;
	for (int i =0; i< res; i++){
		 lb+=residue[k][i].Nsite;
	}
	int ub = lb + residue[k][res].Nsite;

#ifdef MC
	forces_mc(lb,ub,0,k);
#endif
	double U_initial;
	double U_final;
	int old_bin;
	int new_bin;
	
	calcvalue(k);
	for(int l =lb; l<ub; l++){
		atom_temp[k][l]		= atom[k][l];				/* Back up coordinates with pdb			*/
		atnopbc_temp[k][l]	= atnopbc[k][l];		/* Back up coordinates without pdb		*/
	}
	for(int l =0; l< box[k].boxns; l++){
		ff_temp[k][l]		= ff[k][l];				/* Back up all the force  components	*/
	}
	#ifdef PRESSURE
		pvir_temp[k]		= pvir[k];					/* Back up all the virial components	*/
	#endif
	en_temp[k]		= en[k];					/* Back up all the energy components	*/
	config_temp[k]  = config[k];
//#ifdef MC		
//		U_initial		= enmc[k].o_potens;
//#endif
//#ifndef MC		
		U_initial		= en[k].potens;
//#endif

	old_bin			=  (int) ((U_initial - sim_dos[k].e_begin)/sim_dos[k].e_width);

	for (int i =lb; i<ub; i++){
	  atnopbc[k][i].x += dx;
	  atnopbc[k][i].y += dy;
	  atnopbc[k][i].z += dz;
	  atom[k][i].x += dx;
	  atom[k][i].y += dy;
	  atom[k][i].z += dz;
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

		/* ----------------------------------------------- */
		/* Check the acceptance of the move. First check   */
    /* to see if the system moved out of range, then   */
    /* check to see if it is accepted.                 */
		/* ----------------------------------------------- */

		if(U_final >= sim_dos[k].e_begin && U_final < sim_dos[k].e_end){
      new_bin	=  (int) ((U_final - sim_dos[k].e_begin)/sim_dos[k].e_width);
      double g_old = dos_interp(k,old_bin,U_initial);
      double g_new = dos_interp(k,new_bin,U_final);
      double arg = g_old - g_new;

		  if (arg>0){
				  dos_hist[k][new_bin].h_of_e ++;
				  dos_hist[k][new_bin].g_of_e += sim_dos[k].mod_f;
				  flat_histogram(k);
				  mc_trans.trans_acc[k] ++;
		  }//accepted
		  else if (exp(arg) > ran2()){
			  dos_hist[k][new_bin].h_of_e ++;
			  dos_hist[k][new_bin].g_of_e += sim_dos[k].mod_f;
			  flat_histogram(k);
			  mc_trans.trans_acc[k] ++; 
		  }//accepted
		  else{
			  for(int l =lb; l<ub; l++){
				  atom[k][l]		= atom_temp[k][l];	 /* Back up coordinates with pdb		*/
				  atnopbc[k][l]	= atnopbc_temp[k][l];/* Back up coordinates without pdb		*/
			  }
			  for(int l =0; l< box[k].boxns; l++){
				  ff[k][l]		= ff_temp[k][l];		/* Back up all the force  components	*/
			  }
			  #ifdef PRESSURE
			  pvir[k]		= pvir_temp[k];				/* Back to old  virial components		*/
			  #endif
			  en[k]		= en_temp[k];				/* Back to old  energy components		*/
			  config[k]   = config_temp[k];
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
			  for(int l =lb; l<ub; l++){
				  atom[k][l]		= atom_temp[k][l];	 /* Back up coordinates with pdb		*/
				  atnopbc[k][l]	= atnopbc_temp[k][l];/* Back up coordinates without pdb		*/
			  }
			  for(int l =0; l< box[k].boxns; l++){
				  ff[k][l]		= ff_temp[k][l];		/* Back up all the force  components	*/
			  }
			  #ifdef PRESSURE
			  pvir[k]		= pvir_temp[k];				/* Back to old  virial components		*/
			  #endif
			  en[k]		= en_temp[k];				/* Back to old  energy components		*/
			  config[k]   = config_temp[k];
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
    calcvalue(k);
		dos_svalues(k);
}
#endif

	
