#ifdef TDXEDOS
/* ======================================================================== */
/* tdxedos_hmc.cpp			                                                */
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
void kinet				(int);
void calcvalue			(int);
void nblist				(int);
void flat_histogram_td		(int);
void integrate_vverlet	(int, double);
void integrate_dos_mts	(int, double);
void end_to_end			(int,int,int);
void force_short(int);
#ifdef DLIST
void nl_check(int,int,int);
#endif
#ifdef WALL
double wall_angle	(int, int);
#endif
double dos_interp_2	(int,double,double);

void tdxedos_hmc(int k)
{
	double U_initial;
	double U_final;
	double K_initial;
	double K_final;
	double E_initial;
	double E_final;
	double x1_initial;
	double x1_final;
	double x2_initial;
	double x2_final;

	int old_bin_1;
	int new_bin_1;
	int old_bin_2;
	int new_bin_2;

	int	boxns = box[k].boxns;

	/*------------------------------------*/
	/* Call subroutines that calculate    */
	/* the properties defining the rxn-   */
	/* coordinate and assign x1 and x2    */
	/* initial values.                    */
	/*------------------------------------*/
	end_to_end(k,SITE1,SITE2);
	x1_initial = ordparam[k].d_nc;
	x2_initial = wall_angle(k,wall[k].angle_site);
	old_bin_1  = (int)((x1_initial - sim_dos[k].x1_begin)/sim_dos[k].x1_width);
	old_bin_2  = (int)((x2_initial - sim_dos[k].x2_begin)/sim_dos[k].x2_width);

	double beta_old	= 1.0/sim.kT[k];			
	double beta_new	= beta_old;

	double	vscal[NBOXES];
	vscal[k] = 1.0;					// initial scale for velocity
	vinit_hyb(k);

	kinet(k);
	calcvalue(k);

	/* ---------------------------------------------------- */
	/* Back up the necessary information, and calculate the */
	/* system parameters before the move.                   */
	/* ---------------------------------------------------- */

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
	
	/* ----------------------------------------------- */
	/* Propagate the system with verlet integrator.    */
	/* ----------------------------------------------- */	
	for (int j=0; j<sim_hyb.cyc_hybmd; j++){
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

	/*------------------------------------*/
	/* Call subroutines that calculate    */
	/* the properties defining the rxn-   */
	/* coordinate and assign x1 and x2    */
	/* final values.                      */
	/*------------------------------------*/
	x1_final = ordparam[k].d_nc;
	x2_final = wall_angle(k,wall[k].angle_site);
	new_bin_1  = (int)((x1_final - sim_dos[k].x1_begin)/sim_dos[k].x1_width);
	new_bin_2  = (int)((x2_final - sim_dos[k].x2_begin)/sim_dos[k].x2_width);		

	#ifdef WALL
	int beyond_wall_flag = 0; //flag to see if particle moved past wall
		for(int j=0; j<box[k].boxns; j++){  
			if(atom[k][j].z < wall[k].z[0]){
				beyond_wall_flag = 1;
				break;
			}
		}

	if (x1_final >= sim_dos[k].x1_begin && x1_final < sim_dos[k].x1_end &&
		x2_final >= sim_dos[k].x2_begin && x2_final < sim_dos[k].x2_end && beyond_wall_flag == 0){
	#else 
	if (x1_final >= sim_dos[k].x1_begin && x1_final < sim_dos[k].x1_end &&
		x2_final >= sim_dos[k].x2_begin && x2_final < sim_dos[k].x2_end){
	#endif
		/* ------------------------------------------------ */
		/* Determine the new and old density of states by   */
		/* 2D interpolation.                                */
		/* ------------------------------------------------ */
		double g_old; double g_new;
		
		g_old = dos_interp_2(k,x1_initial,x2_initial);
		g_new = dos_interp_2(k,x1_final,x2_final);

		double arg = g_old - g_new + E_initial*beta_old-E_final*beta_new;


		/* ------------------------------------------------ */
		/* Accept or reject the move and perform the        */
		/* necessary updates.                               */
		/* ------------------------------------------------ */

		if (arg>0.0){
			dos_hist[k][new_bin_1][new_bin_2].h_of_l ++;
			dos_hist[k][new_bin_1][new_bin_2].g_of_l += sim_dos[k].mod_f;
			flat_histogram_td(k);
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
			dos_hist[k][new_bin_1][new_bin_2].h_of_l ++;
			dos_hist[k][new_bin_1][new_bin_2].g_of_l += sim_dos[k].mod_f;
			flat_histogram_td(k);
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
			for(int l=0; l<boxns; l++){
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
			dos_hist[k][old_bin_1][old_bin_2].h_of_l ++;
			dos_hist[k][old_bin_1][old_bin_2].g_of_l += sim_dos[k].mod_f;
			flat_histogram_td(k);
			end_to_end(k,SITE1,SITE2);				// recompute to get back to the old
			ordparam[k].x1 = ordparam[k].d_nc;
			ordparam[k].x2 = wall_angle(k,wall[k].angle_site);
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
		dos_hist[k][old_bin_1][old_bin_2].h_of_l ++;
		dos_hist[k][old_bin_1][old_bin_2].g_of_l += sim_dos[k].mod_f;
		flat_histogram_td(k);
		end_to_end(k,SITE1,SITE2);				// recompute to get back to the old
		ordparam[k].x1 = ordparam[k].d_nc;
		ordparam[k].x2 = wall_angle(k,wall[k].angle_site);
	}//rejected

	calcvalue(k); // needed to update box.variable correctly
	dos_svalues(k);
}
#endif

