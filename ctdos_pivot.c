#ifdef CTDOS
/* ======================================================================== */
/* ctdos_pivot.cpp															*/
/* ======================================================================== */
	
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
void nblist_pivot   (int,int);
void nblist         (int);
double ran2         (void);
int ran_int	        (int,int);
void flat_histogram (int);
double dos_interp   (int, int, double);
#ifdef STATUS
void curr_status    (int,int);
#endif
#ifdef DLIST
void nl_check       (int, int, int);
#endif

void ctdos_pivot(int ibox)
{
	int k			=	ibox;
	double theta	=	ran2()*PI*mc_pivot.theta[k];
	double costheta	=	cos(theta);
	double sintheta	=	sin(theta);
	
	/*--------------------------------------*/
	/* INITIALIZATION OF PIVOT MATRICES		*/
	/*--------------------------------------*/

	int ii,jj, kk,ll;	
	for(ii=0; ii<3; ii++)
	{
		for(jj=0; jj<2; jj++)
		{
			for(kk=0; kk<3; kk++)
			{
				for(ll=0; ll<3; ll++)
					PIVOT[ii][jj][kk][ll]=0;
			}
		}
	}
	
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
		if((ub-lb) < 5){
			lb = 0;
			ub = molec.lsite[0];
		}
	}
	int p;
	p = ran_int(lb,ub);		// pick the pivot site : pick Calpha or main_chain N atom only
	while (strcmp(atom[k][p].name,"CA")!=0 && strcmp(atom[k][p].name,"N")!=0 && 
		   strcmp(atom[k][p].name,"P") !=0 && strcmp(atom[k][p].name,"S")!=0){
		p = ran_int(lb,ub);
		#ifdef STATUS
			curr_status(k,0);
		#endif
	}

	int AXIS	=	ran_int(0,3);				// rotation about x, y or z axis
#ifdef WALL
	int SIDE = 0;
#else
	int SIDE	=	ran_int(0,2);				// amino or carboxy side with amino =0 and carboxy =1 side
#endif	
	int ANGLE	=	ran_int(0,2);				// 0 is 90 1 is -90

#ifdef MC
	forces_mc(lb,ub,0,k);
#endif

	double U_initial;
	double U_final;
	int old_bin;
	int new_bin;

	calcvalue(k);
	for(int l =lb; l<ub; l++){
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
		config_temp[k]  = config[k];
//#ifdef MC		
//		U_initial		= enmc[k].o_potens;
//#endif
//#ifndef MC		
		U_initial		= en[k].potens;
//#endif

		old_bin			=  (int) ((U_initial - sim_dos[k].e_begin)/sim_dos[k].e_width);


	/* Now propose a pivot move	*/


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
	#endif
	#ifdef DLIST
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
			  mc_pivot.pivot_acc[k] ++;
		  }//accepted
		  else if (exp(arg) > ran2()){
			  dos_hist[k][new_bin].h_of_e ++;
			  dos_hist[k][new_bin].g_of_e += sim_dos[k].mod_f;
			  flat_histogram(k);
			  mc_pivot.pivot_acc[k] ++;
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
    }
	  calcvalue(k); // needed to update box.variable correctly
	  dos_svalues(k);

}

#endif

		
