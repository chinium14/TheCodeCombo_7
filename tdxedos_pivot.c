#ifdef TDXEDOS
/* ======================================================================== */
/* tdxedos_pivot.cpp                                                        */
/* This subroutine attempts a pivot move on one of the molecules defined as */
/* a solute in the mol.input file.  It then accepts the move based upon the */
/* tdxedos criterion.                                                       */
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
void nblist_pivot (int, int);
void nblist (int);
double ran2		  (void);
int	 ran_int	  (int,int);
void flat_histogram_td(int);
#ifdef STATUS
void curr_status		  (int,int);
#endif
void end_to_end(int, int, int);
#ifdef DLIST
void nl_check(int,int,int);
#endif
#ifdef WALL
double wall_angle	(int, int);
#endif
double dos_interp_2	(int,double,double);

void tdxedos_pivot(int ibox)
{
	int k			=	ibox;
	/* ---------------------------------------------------- */
	/* Here are some variable needed for a pivot move.      */
	/* ---------------------------------------------------- */
	double theta	=	ran2()*PI*mc_pivot.theta[k];
	double costheta	=	cos(theta);
	double sintheta	=	sin(theta);
	
	/*--------------------------------------*/
	/* INITIALIZATION OF PIVOT MATRICES		*/
	/*--------------------------------------*/

	int ii,jj, kk,ll;	
	for(ii=0; ii<3; ii++){
		for(jj=0; jj<2; jj++){
			for(kk=0; kk<3; kk++){
				for(ll=0; ll<3; ll++) PIVOT[ii][jj][kk][ll]=0;
			}
		}
	}
	// intialize all 54 enteries

	
	PIVOT[0][0][0][0]= 1.0;				// Rx (clockwise)
	PIVOT[0][0][1][1]= costheta;
	PIVOT[0][0][1][2]=-sintheta;
	PIVOT[0][0][2][1]= sintheta;
	PIVOT[0][0][2][2]= costheta;

	PIVOT[0][1][0][0]= 1.0;				// Rx (counter clockwise)
	PIVOT[0][1][1][1]= costheta;
	PIVOT[0][1][1][2]= sintheta;
	PIVOT[0][1][2][1]=-sintheta;
	PIVOT[0][1][2][2]= costheta;
	
	PIVOT[1][0][0][0]= costheta;
	PIVOT[1][0][0][2]= sintheta;		// Ry(clockwise)
	PIVOT[1][0][1][1]= 1.0;
	PIVOT[1][0][2][0]=-sintheta;
	PIVOT[1][0][2][2]= costheta;

	PIVOT[1][1][0][0]= costheta;
	PIVOT[1][1][0][2]=-sintheta;		// Ry(counter clockwise)
	PIVOT[1][1][1][1]= 1.0;
	PIVOT[1][1][2][0]= sintheta;
	PIVOT[1][1][2][2]= costheta;
	
	PIVOT[2][0][0][0]= costheta;
	PIVOT[2][0][0][1]= sintheta;		// Rz(clockwise)
	PIVOT[2][0][1][0]=-sintheta;
	PIVOT[2][0][1][1]= costheta;
	PIVOT[2][0][2][2]= 1.0;

	PIVOT[2][1][0][0]= costheta;
	PIVOT[2][1][0][1]=-sintheta;		// Rz(counter clockwise)
	PIVOT[2][1][1][0]= sintheta;
	PIVOT[2][1][1][1]= costheta;
	PIVOT[2][1][2][2]= 1.0;

	/* ---------------------------------------------------- */
	/* First, pick a random solute molecule to apply the    */
	/* move to.  Then pick a random site on that molecule.  */
	/* Then, pick the axis, side, and angle to rotate.      */
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
		if((ub-lb) < 5){
			lb = 0;
			ub = molec.lsite[0];
		}
	}
	int p;
	p = ran_int(lb,ub);		// pick the pivot site : pick Calpha or main_chain N atom only
	int pivot_site = strcmp(atom[k][p].name,"CA")*strcmp(atom[k][p].name,"N")*strcmp(atom[k][p].name,"P")*
		strcmp(atom[k][p].name,"C5'")*strcmp(atom[k][p].name,"O5'")*strcmp(atom[k][p].name,"S");
//	while (strcmp(atom[k][p].name,"CA")!=0 && strcmp(atom[k][p].name,"N")!=0){
	while(pivot_site != 0){
		p = ran_int(lb,ub);
		#ifdef STATUS
			curr_status(k,0);
		#endif
		pivot_site = strcmp(atom[k][p].name,"CA")*strcmp(atom[k][p].name,"N")*strcmp(atom[k][p].name,"P")*
			strcmp(atom[k][p].name,"C5'")*strcmp(atom[k][p].name,"O5'")*strcmp(atom[k][p].name,"S");
	}

	int AXIS	=	ran_int(0,3);				// rotation about x, y or z axis

	#ifdef WALL
	int SIDE = 0;
	#else
	int SIDE	=	ran_int(0,2);				// amino or carboxy side with amino =0 and carboxy =1 side
	#endif	

	int ANGLE	=	ran_int(0,2);				// 0 is 90 1 is -90

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

	/* ---------------------------------------------------- */
	/* Variable to evaluate if move is accepted.            */
        /* ---------------------------------------------------- */
	#ifdef MC
	forces_mc(lb,ub,0,k);
	#endif

	double beta_old	= 1.0/sim.kT[k];

	/* ---------------------------------------------------- */
	/* Back up the necessary information, and calculate the */
	/* system parameters before the move.                   */
	/* ---------------------------------------------------- */
	calcvalue(k);
	for(int l=lb; l<ub; l++){
		atom_temp[k][l]	= atom[k][l];				/* Back up coordinates with pdb			*/
		atnopbc_temp[k][l]	= atnopbc[k][l];		/* Back up coordinates without pdb		*/
	}
	for(int l =0; l< box[k].boxns; l++){
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
	/* Now, propose a pivot move.                           */
	/* ---------------------------------------------------- */	
	if (SIDE==0){
		for (int i=lb; i<p; i++){
			double pos_x	=	PIVOT[AXIS][ANGLE][0][0]*(atnopbc[k][i].x-atnopbc[k][p].x)+ 
								PIVOT[AXIS][ANGLE][0][1]*(atnopbc[k][i].y-atnopbc[k][p].y)+
								PIVOT[AXIS][ANGLE][0][2]*(atnopbc[k][i].z-atnopbc[k][p].z)+ atnopbc[k][p].x;
			double pos_y	=	PIVOT[AXIS][ANGLE][1][0]*(atnopbc[k][i].x-atnopbc[k][p].x)+
								PIVOT[AXIS][ANGLE][1][1]*(atnopbc[k][i].y-atnopbc[k][p].y)+
								PIVOT[AXIS][ANGLE][1][2]*(atnopbc[k][i].z-atnopbc[k][p].z)+ atnopbc[k][p].y;
			double pos_z	=	PIVOT[AXIS][ANGLE][2][0]*(atnopbc[k][i].x-atnopbc[k][p].x)+
								PIVOT[AXIS][ANGLE][2][1]*(atnopbc[k][i].y-atnopbc[k][p].y)+
								PIVOT[AXIS][ANGLE][2][2]*(atnopbc[k][i].z-atnopbc[k][p].z)+ atnopbc[k][p].z;
			atnopbc[k][i].x	=	pos_x;
			atnopbc[k][i].y	=	pos_y;
			atnopbc[k][i].z	=	pos_z;
			atom[k][i].x	=	atnopbc[k][i].x;
			atom[k][i].y	=	atnopbc[k][i].y;
			atom[k][i].z	=	atnopbc[k][i].z;
		}
	}
	else{	//if SIDE ==1
		for (int i=p+1; i<ub; i++){
			double pos_x	=	PIVOT[AXIS][ANGLE][0][0]*(atnopbc[k][i].x-atnopbc[k][p].x)+ 
								PIVOT[AXIS][ANGLE][0][1]*(atnopbc[k][i].y-atnopbc[k][p].y)+
								PIVOT[AXIS][ANGLE][0][2]*(atnopbc[k][i].z-atnopbc[k][p].z)+ atnopbc[k][p].x;
			double pos_y	=	PIVOT[AXIS][ANGLE][1][0]*(atnopbc[k][i].x-atnopbc[k][p].x)+
								PIVOT[AXIS][ANGLE][1][1]*(atnopbc[k][i].y-atnopbc[k][p].y)+
								PIVOT[AXIS][ANGLE][1][2]*(atnopbc[k][i].z-atnopbc[k][p].z)+ atnopbc[k][p].y;
			double pos_z	=	PIVOT[AXIS][ANGLE][2][0]*(atnopbc[k][i].x-atnopbc[k][p].x)+
								PIVOT[AXIS][ANGLE][2][1]*(atnopbc[k][i].y-atnopbc[k][p].y)+
								PIVOT[AXIS][ANGLE][2][2]*(atnopbc[k][i].z-atnopbc[k][p].z)+ atnopbc[k][p].z;
			atnopbc[k][i].x	=	pos_x;
			atnopbc[k][i].y	=	pos_y;
			atnopbc[k][i].z	=	pos_z;
			atom[k][i].x	=	atnopbc[k][i].x;
			atom[k][i].y	=	atnopbc[k][i].y;
			atom[k][i].z	=	atnopbc[k][i].z;
		}
	}

	/* ----------------------------------------------- */
	/* Apply periodic boundary conditions.             */
	/* ----------------------------------------------- */
	pbc_all(k);
	#ifdef NLIST
	#ifndef DLIST
	nblist_pivot(k,ub); // call Nblist
	#else
	if (SIDE==0) nl_check(lb,p,k);
	else nl_check(p+1,ub,k);
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
			mc_pivot.pivot_acc[k] ++;
		}//accepted
		else if (exp(arg) > ran2()){
			dos_hist[k][new_bin_1][new_bin_2].h_of_l ++;
			dos_hist[k][new_bin_1][new_bin_2].g_of_l += sim_dos[k].mod_f;
			flat_histogram_td(k);
			mc_pivot.pivot_acc[k] ++;
		}//accepted
		else{
			for(int l=lb; l<ub; l++){
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

			#ifdef NLIST
			#ifndef DLIST
			nblist_pivot(k,ub); // call Nblist
			#endif
			#ifdef DLIST
			nl_check(lb,p,k);
			if(nl_flag[k] == 1) nblist_pivot(k,ub);
			#endif
			#endif

			dos_hist[k][old_bin_1][old_bin_2].h_of_l ++;
			dos_hist[k][old_bin_1][old_bin_2].g_of_l += sim_dos[k].mod_f;
			flat_histogram_td(k);
			end_to_end(k,SITE1,SITE2);
			ordparam[k].x1 = ordparam[k].d_nc;
			ordparam[k].x2 = wall_angle(k,wall[k].angle_site);
		}//rejected
	}
	else{
		for(int l=lb; l<ub; l++){
			atom[k][l]		= atom_temp[k][l];	 /* Back up coordinates with pdb		*/
			atnopbc[k][l]	= atnopbc_temp[k][l];/* Back up coordinates without pdb		*/
		}
		for(int l =0; l< box[k].boxns; l++){
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
		nl_check(lb,p,k);
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

		
