#ifdef CTDOS
/* ======================================================================== */
/* ctdos_hmc.cpp			                                                */
/* This suboutine performs density of states simulation. In the beginning a */
/* random walk is performed in U by assigning initial veolocities based on  */
/* a beta of the concerned box. In additon to the pivot moves, swap moves   */
/* are implemented to facilitate sampling and also to achieve faster		*/
/* convergence.																*/		
/* ======================================================================== */
	
#include "defines.h"

double ran2           (void);
void vinit_hyb        (int);
void dos_svalues      (int);
void pbc_all          (int);
void pbc_chk          (int);
void kinet            (int);
void calcvalue        (int);
void nblist           (int);
void flat_histogram   (int);
void integrate_vverlet(int, double);
void integrate_dos_mts(int,double);
void force_short      (int);
double dos_interp     (int, int, double);

#ifdef CONFIGT
  void configtemp     (int);
#endif
#ifdef DLIST
  void nl_check       (int, int, int);
#endif

void ctdos_hmc(int k)
{
	int		old_bin;
	int		new_bin;
	double	U_initial;
	double	U_final;
	double	K_initial;
	double	K_final;
	double	vscal[NBOXES];
	double	beta_old;
	double	beta_new;
	vscal[k] = 1.0;					// initial scale for velocity


	sim.T[k]	= sim_dos[k].T_begin;
	sim.kT[k]	= RG * sim.T[k]*.001; // the units are in kJ/mol , same as energy
	beta_old	= 1.0/sim.kT[k];			
	beta_new	= beta_old;
	vinit_hyb(k);
	kinet(k);
	calcvalue(k);
	for(int l =0; l< box[k].boxns; l++){
		atom_temp[k][l]		= atom[k][l];		/* Back up coordinates with pdb			*/
		atnopbc_temp[k][l]	= atnopbc[k][l];	/* Back up coordinates without pdb		*/
		ff_temp[k][l]		= ff[k][l];			/* Back up all the force  components	*/
		vv_temp[k][l]		= vv[k][l];
	}

	#ifdef PRESSURE
	pvir_temp[k]	= pvir[k];					/* Back up all the virial components	*/
	#endif
	en_temp[k]		= en[k];					/* Back up all the energy components	*/
	U_initial		= en[k].potens;
	K_initial		= en[k].kinet;
	config_temp[k]  = config[k];
	old_bin			=  (int)((U_initial - sim_dos[k].e_begin)/sim_dos[k].e_width);			

	/* ----------------------------------------------- */
	/* Propagate the system with verlet integrator.    */
	/* ----------------------------------------------- */	
	for (int j=0;j < sim_hyb.cyc_hybmd; j++){

		if(sim.nsteps ==1){
			integrate_vverlet(k, vscal[k]);
			#ifdef DLIST
				nl_check(0,box[k].boxns,k);
				if(nl_flag[k] == 1) nblist(k);
			#endif
		}

		else {
			if(j==0){
			  force_short(k);
				for(int l=0; l< box[k].boxns; l++){
					ff_long[k][l].x		= ff_temp[k][l].x-ff_short[k][l].x;
					ff_long[k][l].y		= ff_temp[k][l].y-ff_short[k][l].y;
					ff_long[k][l].z		= ff_temp[k][l].z-ff_short[k][l].z;
				}
			}
			integrate_dos_mts(k,vscal[k]);
			#ifdef DLIST
				nl_check(0,box[k].boxns,k);
				if(nl_flag[k] == 1) nblist(k);
			#endif

		}
		#ifdef STATUS
		curr_status(k,7);
		#endif
	}
  
  /* ----------------------------------------------- */
  /* Calculate the properties of the system.         */
  /* ----------------------------------------------- */
	kinet(k);    											
	calcvalue(k);
	U_final	= en[k].potens;
	K_final	= en[k].kinet;

	/* ----------------------------------------------- */
	/* Check the acceptance of the move. First check   */
  /* to see if the system moved out of range, then   */
  /* check to see if it is accepted.                 */
	/* ----------------------------------------------- */

	if(U_final >= sim_dos[k].e_begin && U_final < sim_dos[k].e_end){
    new_bin	=  (int) ((U_final - sim_dos[k].e_begin)/sim_dos[k].e_width);
    double g_old = dos_interp(k,old_bin,U_initial);
    double g_new = dos_interp(k,new_bin,U_final);
    double arg = g_old - g_new + K_initial*beta_old - K_final*beta_new;

	  if (arg>0){
			  dos_hist[k][new_bin].h_of_e ++;
			  dos_hist[k][new_bin].g_of_e += sim_dos[k].mod_f;
			  flat_histogram(k);
			  #ifdef NLIST
			  #ifndef DLIST
				  if(sim_dos[k].dos_acc %sim.nlist==0) nblist(k); // call Nblist every nblist accepted moves
			  #endif
			  #ifdef DLIST
				  nl_check(0,box[k].boxns,k);
				  if(nl_flag[k] == 1) nblist(k);
			  #endif
			  #endif
			  sim_dos[k].dos_acc ++;
	  }//accepted
	  
	  else if (exp(arg) > ran2()){
			  dos_hist[k][new_bin].h_of_e ++;
			  dos_hist[k][new_bin].g_of_e += sim_dos[k].mod_f;
			  flat_histogram(k);
			  #ifdef NLIST
			  #ifndef DLIST
				  if(sim_dos[k].dos_acc %sim.nlist==0) nblist(k); // call Nblist every nblist accepted moves
			  #endif
			  #ifdef DLIST
				  nl_check(0,box[k].boxns,k);
				  if(nl_flag[k] == 1) nblist(k);
			  #endif
			  #endif
			  sim_dos[k].dos_acc ++;
	  }//accepted
	  
	  else{
			  for(int l =0; l< box[k].boxns; l++){
				  atom[k][l]		= atom_temp[k][l];		/* Back up coordinates with pdb			*/
				  atnopbc[k][l]	= atnopbc_temp[k][l];	/* Back up coordinates without pdb		*/
				  ff[k][l]		= ff_temp[k][l];		/* Back to old  force  components		*/
				  vv[k][l]		= vv_temp[k][l];
				  uu[k][l]		= vv[k][l];
			  }
			  #ifdef PRESSURE
			  pvir[k]		= pvir_temp[k];					/* Back to old  virial components		*/
			  #endif
			  en[k]		= en_temp[k];					/* Back to old  energy components		*/
			  config[k]  = config_temp[k];
			  dos_hist[k][old_bin].h_of_e ++;
			  dos_hist[k][old_bin].g_of_e += sim_dos[k].mod_f;
			  flat_histogram(k);
	  }//rejected
  }
  else{
		for(int l =0; l< box[k].boxns; l++){
			atom[k][l]		= atom_temp[k][l];		/* Back up coordinates with pdb			*/
			atnopbc[k][l]	= atnopbc_temp[k][l];	/* Back up coordinates without pdb		*/
			ff[k][l]		= ff_temp[k][l];		/* Back to old  force  components		*/
			vv[k][l]		= vv_temp[k][l];
			uu[k][l]		= vv[k][l];
		}
		#ifdef PRESSURE
		pvir[k]		= pvir_temp[k];					/* Back to old  virial components		*/
		#endif
		en[k]		= en_temp[k];					/* Back to old  energy components		*/
		config[k]  = config_temp[k];
		dos_hist[k][old_bin].h_of_e ++;
		dos_hist[k][old_bin].g_of_e += sim_dos[k].mod_f;
		flat_histogram(k);
  }

	calcvalue(k); // needed to update box.variable correctly
	dos_svalues(k);
  }
#endif

