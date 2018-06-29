#ifdef XEDOS
/* ======================================================================== */
/* xedos_hmc.cpp			                                                */
/* This suboutine performs density of states simulation. In the beginning a */
/* random walk is performed in U by assigning initial veolocities based on  */
/* a beta of the concerned box. In additon to the pivot moves, swap moves   */
/* are implemented to facilitate sampling and also to achieve faster		*/
/* convergence.																*/
/* ======================================================================== */
	
#include "defines.h"


double ran2				(void);
void vinit_hyb			(int);
void dos_svalues		(int);
void pbc_all			(int);
void pbc_chk			(int);
void kinet				(int);
void calcvalue			(int);
void nblist				(int);
void flat_histogram		(int);
void integrate_vverlet	(int, double);
void integrate_dos_mts	(int, double);
void end_to_end			(int,int,int);
void force_short(int);
double dos_interp (int,int,double);
#ifdef PRESSURE
  double virialcor (int);
#endif
#ifdef DLIST
void nl_check(int,int,int);
#endif


void xedos_hmc(int k)
{
	int		old_bin;
	int		new_bin;
	int		boxns = box[k].boxns;
	double  L_initial;
	double  L_final;
	double	U_initial;
	double	U_final;
	double	K_initial;
	double	K_final;
	double  E_initial;
	double  E_final;
	double	vscal[NBOXES];
	double	beta_old;
	double	beta_new;
	vscal[k] = 1.0;					// initial scale for velocity

	beta_old	= 1.0/sim.kT[k];			
	beta_new	= beta_old;
	vinit_hyb(k);
	kinet(k);
	calcvalue(k);

	/* ------------------------------------	*/
	/* Back up coords, force, virial and	*/
	/* energy components.					*/
	/* ------------------------------------	*/

	for(int l =0; l< boxns; l++){
		atom_temp[k][l]		= atom[k][l];		
		atnopbc_temp[k][l]	= atnopbc[k][l];	
		ff_temp[k][l]		= ff[k][l];			
		vv_temp[k][l]		= vv[k][l];
	}
	#ifdef PRESSURE
	pvir_temp[k]	= pvir[k];					
	#endif
	en_temp[k]		= en[k];					
	U_initial		= en[k].potens;
	K_initial		= en[k].kinet;
	E_initial		= U_initial+K_initial;
	end_to_end(k,SITE1,SITE2);
	L_initial		= ordparam[k].d_nc;
	old_bin			= (int) ((L_initial - sim_dos[k].l_begin)/sim_dos[k].l_width);			
	
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
	E_final = U_final+K_final;
	end_to_end(k,SITE1,SITE2);
	L_final = ordparam[k].d_nc;

#ifdef WALL
		int beyond_wall_flag = 0; //flag to see if particle moved past wall
		//for(int i=0; i<wall[k].n; i++){
			for(int j=0; j<box[k].boxns; j++){  
				if(atom[k][j].z < wall[k].z[0]){
					beyond_wall_flag = 1;
					break;
				}
			}
		//}
	/* ----------------------------------------------- */
	/* Check the acceptance of the move. First check   */
  /* to see if the system moved out of range, then   */
  /* check to see if it is accepted.                 */
	/* ----------------------------------------------- */

		if (L_final >= sim_dos[k].l_begin && L_final < sim_dos[k].l_end && beyond_wall_flag == 0){
#else 
		if (L_final >= sim_dos[k].l_begin && L_final < sim_dos[k].l_end){	
#endif

	  new_bin	= (int) ((L_final - sim_dos[k].l_begin)/sim_dos[k].l_width);
		double g_old = dos_interp(k,old_bin,L_initial); 
    double g_new = dos_interp(k,new_bin,L_final);
		double arg = g_old - g_new + E_initial*beta_old-E_final*beta_new;

		if (arg>0.0){
			dos_hist[k][new_bin].h_of_l ++;
			dos_hist[k][new_bin].g_of_l += sim_dos[k].mod_f;
			flat_histogram(k);
			#ifdef NLIST
			#ifndef DLIST
				if(sim_dos[k].dos_acc %sim.nlist==0) nblist(k); // call Nblist every nblist accepted moves
			#endif
			#ifdef DLIST
				nl_check(0, boxns, k);
				if(nl_flag[k] == 1) nblist(k);
			#endif
			#endif
			sim_dos[k].dos_acc ++;
		}//accepted
	
		else if (exp(arg) > ran2()){
			dos_hist[k][new_bin].h_of_l ++;
			dos_hist[k][new_bin].g_of_l += sim_dos[k].mod_f;
			flat_histogram(k);
			#ifdef NLIST
			#ifndef DLIST
				if(sim_dos[k].dos_acc %sim.nlist==0) nblist(k); // call Nblist every nblist accepted moves
			#endif
			#ifdef DLIST
				nl_check(0, boxns, k);
				if(nl_flag[k] == 1) nblist(k);
			#endif
			#endif
			sim_dos[k].dos_acc ++;
		}//accepted
	
		else{
			for(int l =0; l< boxns; l++){
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
			dos_hist[k][old_bin].h_of_l ++;
			dos_hist[k][old_bin].g_of_l += sim_dos[k].mod_f;
			flat_histogram(k);
			end_to_end(k,SITE1,SITE2);				// recompute to get back to the old

		}//rejected
	}
	else{
		for(int l =0; l< boxns; l++){
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
		dos_hist[k][old_bin].h_of_l ++;
		dos_hist[k][old_bin].g_of_l += sim_dos[k].mod_f;
		flat_histogram(k);
		end_to_end(k,SITE1,SITE2);				// recompute to get back to the old

	}//rejected

	calcvalue(k); // needed to update box.variable correctly
	dos_svalues(k);
  }
#endif

