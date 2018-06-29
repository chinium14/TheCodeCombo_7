#ifdef TDXEDOS	
/* ======================================================================== */
/* tdxedos_linear.cpp		                                                */
/* This suboutine performs a move where a random site is picked and moved   */
/* by some random displacement. The move is accepted according to the		*/
/* tdxedos.                                                                 */
/* Written by Thomas Knotts 14 June 2004.                                   */
/* ======================================================================== */

#include "defines.h"

void dos_svalues  (int);
void pbc_all	  (int);
#ifdef MC
void forces_mc(int,int,int,int);
#else
void forces		  (int);
#endif
void calcvalue	  (int);
void nblist_pivot (int,int);
void nblist		  (int);
double ran2		  (void);
int	   ran_int	  (int,int);
void flat_histogram_td(int);
#ifdef STATUS
void curr_status  (int,int);
#endif
void end_to_end	  (int,int,int);
#ifdef DLIST
void nl_check	  (int,int,int);
#endif
#ifdef WALL
double wall_angle	(int, int);
#endif
double dos_interp_2	(int,double,double);

void tdxedos_linear(int ibox){
	int k	  =	ibox;

	/* ---------------------------------------------------- */
	/* First, pick a random solute molecule and then a      */
	/* random site on that molecule to move.                */
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

	
	/* ---------------------------------------------------- */
	/* Calculate the values needed for the move.            */
	/* ---------------------------------------------------- */
	double U_initial;
	double U_final;
	double x1_initial;
	double x1_final;
	double x2_initial;
	double x2_final;

	int old_bin_1;
	int new_bin_1;
	int old_bin_2;
	int new_bin_2;

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

	#ifdef MC
		forces_mc(rand_site,rand_site+1,0,k);
	#endif	

	double beta_old	= 1.0/sim.kT[k];			
	double beta_new	= beta_old;

	calcvalue(k);

	/* ---------------------------------------------------- */
	/* Back up the necessary information, and calculate the */
	/* system parameters before the move.                   */
	/* ---------------------------------------------------- */
	atom_temp[k][rand_site]		= atom[k][rand_site];			/* Back up coordinates with pdb			*/
	atnopbc_temp[k][rand_site]	= atnopbc[k][rand_site];		/* Back up coordinates without pdb		*/

	for(int l=0; l< box[k].boxns; l++){
		ff_temp[k][l]		= ff[k][l];				/* Back up all the force  components	*/
	}
	#ifdef PRESSURE
	pvir_temp[k]	= pvir[k];					/* Back up all the virial components	*/
	#endif

	en_temp[k]		= en[k];					/* Back up all the energy components	*/
	
	#ifdef MC		
	U_initial		= enmc[k].o_potens;
	#else		
	U_initial		= en[k].potens;
	#endif

	/* ---------------------------------------------------- */
	/* Move the site by some random displacement.           */
	/* ---------------------------------------------------- */
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
	#ifdef NLIST
	#ifndef DLIST
	nblist_pivot(k,ub); // call Nblist 
	#else
	nl_check(lb, ub, k);
	if(nl_flag[k] == 1) nblist_pivot(k,ub);
	#endif
	#endif

	/* ----------------------------------------------- */
	/* Update the forces and system parameters.        */
	/* ----------------------------------------------- */
	#ifdef MC
	forces_mc(rand_site, rand_site+1, 1,k);
	#else
	forces(k);
	#endif

	calcvalue(k);

	#ifdef MC
	U_final	= enmc[k].n_potens;
	#else
	U_final	= en[k].potens;
	#endif
	/*------------------------------------*/
	/* Call subroutines that calculate    */
	/* the properties defining the rxn-   */
	/* coordinate and assign x1 and x2    */
	/* final values.                      */
	/*------------------------------------*/
	end_to_end(k,SITE1,SITE2);
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

		double arg = g_old- g_new + (U_initial-U_final)*beta_old; 

		/* ------------------------------------------------ */
		/* Accept or reject the move and perform the        */
		/* necessary updates.                               */
		/* ------------------------------------------------ */
		if (arg>0.0){
			dos_hist[k][new_bin_1][new_bin_2].h_of_l ++;
			dos_hist[k][new_bin_1][new_bin_2].g_of_l += sim_dos[k].mod_f;
			flat_histogram_td(k);
			mc_rand.rand_acc[k] ++;
		}//accepted
		else if (exp(arg) > ran2()){
			dos_hist[k][new_bin_1][new_bin_2].h_of_l ++;
			dos_hist[k][new_bin_1][new_bin_2].g_of_l += sim_dos[k].mod_f;
			flat_histogram_td(k);
			mc_rand.rand_acc[k] ++;	
		}//accepted
		else{
			atom[k][rand_site]		= atom_temp[k][rand_site];	 /* Back up coordinates with pdb		*/
			atnopbc[k][rand_site]	= atnopbc_temp[k][rand_site];/* Back up coordinates without pdb		*/
			for(int l=0; l< box[k].boxns; l++){
				ff[k][l]		= ff_temp[k][l];	 /* Back up all the force  components	*/
			}
			#ifdef PRESSURE
			pvir[k]		= pvir_temp[k];				 /* Back to old  virial components		*/
			#endif

			en[k]		= en_temp[k];				 /* Back to old  energy components		*/

			#ifdef NLIST
			#ifndef DLIST
			nblist_pivot(k,ub); // call Nblist 
			#else
			nl_check(lb, ub, k);
			if(nl_flag[k] == 1) nblist_pivot(k,ub);
			#endif
			#endif

			dos_hist[k][old_bin_1][old_bin_2].h_of_l ++;
			dos_hist[k][old_bin_1][old_bin_2].g_of_l += sim_dos[k].mod_f;
			flat_histogram_td(k);
			end_to_end(k,SITE1,SITE2);				// recompute to get back to the old
			ordparam[k].x1 = ordparam[k].d_nc;
			ordparam[k].x2 = wall_angle(k,wall[k].angle_site);
		}//rejected (due to arg)
	}
	else{
		atom[k][rand_site]		= atom_temp[k][rand_site];	 /* Back up coordinates with pdb		*/
		atnopbc[k][rand_site]	= atnopbc_temp[k][rand_site];/* Back up coordinates without pdb		*/
		for(int l=0; l< box[k].boxns; l++){
			ff[k][l]		= ff_temp[k][l];	 /* Back up all the force  components	*/
		}
		#ifdef PRESSURE
		pvir[k]		= pvir_temp[k];				 /* Back to old  virial components		*/
		#endif

		en[k]		= en_temp[k];				 /* Back to old  energy components		*/

		#ifdef NLIST
		#ifndef DLIST
		nblist_pivot(k,ub); // call Nblist 
		#else
		nl_check(lb, ub, k);
		if(nl_flag[k] == 1) nblist_pivot(k,ub);
		#endif
		#endif

		dos_hist[k][old_bin_1][old_bin_2].h_of_l ++;
		dos_hist[k][old_bin_1][old_bin_2].g_of_l += sim_dos[k].mod_f;
		flat_histogram_td(k);
		end_to_end(k,SITE1,SITE2);				// recompute to get back to the old
		ordparam[k].x1 = ordparam[k].d_nc;
		ordparam[k].x2 = wall_angle(k,wall[k].angle_site);
	}//rejected (outside range)
	calcvalue(k);
	dos_svalues(k);

}

#endif

		
