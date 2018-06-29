#ifdef TDXEDOS
/* ======================================================================== */
/* tdxedos_trans.cpp	                                                    */
/* This suboutine performs a move and evaluates its acceptance for tdxedos.	*/
/* It is designed to be used when investigating such things as the          */
/* hybrization of DNA or protein ligand binding.  It takes two separate     */
/* molecules and moves them apart from or closer to each other along the    */
/* the reaction coordinate line defined by SITE1 and SITE2 in mol.input.    */
/*																		    */
/* Written by Thomas Knotts 14 June 2004                                    */
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
void nl_check(int,int,int);
#endif
#ifdef WALL
double wall_angle	(int, int);
#endif
double dos_interp_2	(int,double,double);

void tdxedos_trans(int ibox){
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

	/* ---------------------------------------------------- */
	/* Calculate the values needed for the move.            */
	/* ---------------------------------------------------- */
	double U_initial;
	double U_final;
	double L_final;
	double L_initial;
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
	forces_mc(lb,ub,0,k);
	#endif

	double beta_old	= 1.0/sim.kT[k];			
	double beta_new	= beta_old;

	calcvalue(k);

	/* ---------------------------------------------------- */
	/* Back up the necessary information, and calculate the */
	/* system parameters before the move.                   */
	/* ---------------------------------------------------- */
	for(int l=lb; l<ub; l++){
		atom_temp[k][l]	= atom[k][l];				/* Back up coordinates with pdb			*/
		atnopbc_temp[k][l]	= atnopbc[k][l];		/* Back up coordinates without pdb		*/
	}
	for(int l=0; l<boxns; l++){
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

	L_initial		= ordparam[k].d_nc;

	/* ---------------------------------------------------- */
	/* Perform a scaling move.                              */
	/* ---------------------------------------------------- */

	double dx	= atnopbc[k][SITE2].x -atnopbc[k][SITE1].x;
	double dy	= atnopbc[k][SITE2].y -atnopbc[k][SITE1].y;
	double dz	= atnopbc[k][SITE2].z -atnopbc[k][SITE1].z;

	double r = sqrt(dx*dx+dy*dy+dz*dz);

	double ex_dist = (ran2()-0.5)*2.0 * mc_trans.delta[k];

	double scalar = ex_dist/r;

	dx *= scalar;
	dy *= scalar;
	dz *= scalar;

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
	#else
	nl_check(lb,ub,k);
	if(nl_flag[k] == 1) nblist_pivot(k,ub);
	#endif
	#endif

	/* ----------------------------------------------- */
	/* Update the forces and system parameters.        */
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

	/*------------------------------------*/
	/* Call subroutines that calculate    */
	/* the properties defining the rxn-   */
	/* coordinate and assign x1 and x2    */
	/* final values.                      */
	/*------------------------------------*/
	end_to_end(k,SITE1,SITE2);
	L_final = ordparam[k].d_nc;
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

		double arg = g_old- g_new + (U_initial-U_final)*beta_old + 2.0*log(L_final/L_initial);

		if (arg>0.0){
			dos_hist[k][new_bin_1][new_bin_2].h_of_l ++;
			dos_hist[k][new_bin_1][new_bin_2].g_of_l += sim_dos[k].mod_f;
			flat_histogram_td(k);
			mc_trans.trans_acc[k] ++;
		}//accepted
		else if (exp(arg) > ran2()){
			dos_hist[k][new_bin_1][new_bin_2].h_of_l ++;
			dos_hist[k][new_bin_1][new_bin_2].g_of_l += sim_dos[k].mod_f;
			flat_histogram_td(k);
			mc_trans.trans_acc[k] ++;
		}//accepted
		else{
			for(int l=lb; l<ub; l++){
				atom[k][l]		= atom_temp[k][l];	 /* Back up coordinates with pdb		*/
				atnopbc[k][l]	= atnopbc_temp[k][l];/* Back up coordinates without pdb		*/
			}
			for(int l=0; l<boxns; l++){
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
			nl_check(lb,ub,k);
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
		for(int l=lb; l<ub; l++){
			atom[k][l]		= atom_temp[k][l];	 /* Back up coordinates with pdb		*/
			atnopbc[k][l]	= atnopbc_temp[k][l];/* Back up coordinates without pdb		*/
		}
		for(int l=0; l<boxns; l++){
			ff[k][l]		= ff_temp[k][l];	 /* Back up all the force  components	*/
		}
		#ifdef PRESSURE
		pvir[k]		= pvir_temp[k];				 /* Back to old  virial components		*/
		#endif

		en[k]		= en_temp[k];				 /* Back to old  energy components		*/

		#ifdef NLIST
		#ifndef DLIST
		nblist_pivot(k,ub); // call Nblist 
		#endif
		#ifdef DLIST
		nl_check(lb,ub,k);
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
		
