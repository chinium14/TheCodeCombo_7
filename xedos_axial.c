
#ifdef XEDOS	
#include "defines.h"

void dos_svalues  (int);
void pbc_all	  (int);
#ifdef MC
void forces_mc(int,int,int,int);
#endif
#ifndef MC
void forces		  (int);
#endif
void calcvalue	  (int);
void nblist_pivot (int, int);
void nblist		  (int);
double ran2		  (void);
int	   ran_int	  (int,int);
void flat_histogram(int);
double dos_interp(int,int,double);
#ifdef STATUS
void curr_status  (int,int);
#endif
void end_to_end	  (int,int,int);
#ifdef DLIST
void nl_check(int,int,int);
#endif



void xedos_axial(int ibox)
{
	int k	  =	ibox;

/* ============= Note =================== */
/* This move only makes sense if you only */
/* have one solute molecule.  The         */
/* selction of lb and ub need to be       */
/* checked for multiple solute molecules  */
/* ====================================== */
	/* ---------------------------------------------------- */
	/* Set the bounds on which to do the move from 0 to the */
	/* last site of the last solute molecule.               */
	/* ---------------------------------------------------- */
	int lb;
	int ub;
	if(molec.Nsolute == 1){
		lb = 0;
		ub = molec.lsite[0];
	}
	else{
		lb = 0;//first site of first solute molecule
		ub = molec.lsite[molec.Nsolute-1];//last site of last solute molecule
	}

	int N = ub-lb;

	double lx	= atnopbc[k][SITE2].x -atnopbc[k][SITE1].x;
	double ly	= atnopbc[k][SITE2].y -atnopbc[k][SITE1].y;
	double lz	= atnopbc[k][SITE2].z -atnopbc[k][SITE1].z;

	double costheta		=	lz/sqrt(lx*lx + ly*ly + lz*lz);
	double costheta2	=   costheta*costheta;
	double sintheta2	=	1.0- costheta2;
	double sintheta		=	sqrt (sintheta2);
	double cos2theta	=	2.0*costheta2 - 1.0;
	double cosphi		=	lx/sqrt(lx*lx + ly*ly);
	double cosphi2		=   cosphi*cosphi;
	double sinphi2		=	1.0-cosphi2;
//	double sinphi		=	sqrt(sinphi2);					//why not this??
	double sinphi		=	ly/sqrt(lx*lx + ly*ly);

//	double s		=   1.0 + (ran2()-0.5)*mc_axial.strain[k]*2.0;
//	double s = 1.1;

	double L_initial;
	end_to_end(k,SITE1,SITE2);				/* Distance betweeb SITE1 and SITE2		*/
	L_initial		= ordparam[k].d_nc;
	double s		=   1.0 + (ran2()-0.5)*mc_axial.strain[k]*2.0/L_initial;
//		double s		=   1.0 + mc_axial.strain[k]*2.0/L_initial;

	double R11	=	sinphi2+cosphi2*(costheta2+s*sintheta2);
	double R12	=	(s-1.0)*cosphi*sinphi*sintheta2;
	double R13	=	(s-1.0)*cosphi*costheta*sintheta;
	double R21	=	R12;
	double R22	=	cosphi2+sinphi2*(costheta2+sintheta2*s);
	double R23	=	(s-1.0)*costheta*sinphi*sintheta;
	double R31	=	R13;
	double R32	=	R23;
	double R33	=	0.5*(1.0+s+(s-1.0)*cos2theta);

/*	This is the basis vector transformation Matrix
	double R11	=	costheta*cosphi;
	double R12	=	costheta*sinphi;
	double R13	=	-sintheta;
	double R21	=	-sinphi;
	double R22	=	cosphi;
	double R23	=	0.0;
	double R31	=	cosphi*sintheta;
	double R32	=	sinphi*sintheta;
	double R33	=	costheta;	
*/
	

#ifdef MC
	forces_mc(lb,ub,0,k);
#endif
	double U_initial;
	double U_final;
	double L_final;
	
	int old_bin;
	int new_bin;
	double beta_old	= 1.0/sim.kT[k];			
	double beta_new	= beta_old;

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
#endif
#ifndef MC		
		U_initial		= en[k].potens;
#endif
		end_to_end(k,SITE1,SITE2);				/* Distance betweeb SITE1 and SITE2		*/
		L_initial		= ordparam[k].d_nc;
		old_bin			= (int) ((L_initial - sim_dos[k].l_begin)/sim_dos[k].l_width);			

	/* Now propose a rescaling move	*/

		for (int i=lb; i<ub; i++){
			double pos_x	=	R11*atnopbc[k][i].x + R12*atnopbc[k][i].y+ R13*atnopbc[k][i].z;
			double pos_y	=	R21*atnopbc[k][i].x + R22*atnopbc[k][i].y+ R23*atnopbc[k][i].z;
			double pos_z	=	R31*atnopbc[k][i].x + R32*atnopbc[k][i].y+ R33*atnopbc[k][i].z;
			
			atnopbc[k][i].x	=	pos_x;
			atnopbc[k][i].y	=	pos_y;
			atnopbc[k][i].z	=	pos_z;
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
		U_final	= enmc[k].n_potens;
#endif
#ifndef MC
		U_final	= en[k].potens;
#endif

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
			double arg = g_old- g_new + (U_initial-U_final)*beta_old + N*log(L_final/L_initial);

			if (arg>0.0){
				dos_hist[k][new_bin].h_of_l ++;
				dos_hist[k][new_bin].g_of_l += sim_dos[k].mod_f;
				flat_histogram(k);
				mc_axial.axial_acc[k] ++;
			}//accepted
			else if (exp(arg) > ran2()){
				dos_hist[k][new_bin].h_of_l ++;
				dos_hist[k][new_bin].g_of_l += sim_dos[k].mod_f;
				flat_histogram(k);
				mc_axial.axial_acc[k] ++;
			}//accepted
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
				#endif
				#ifdef DLIST
					nl_check(lb,ub,k);
					if(nl_flag[k] == 1) nblist_pivot(k,ub);
				#endif
				#endif
				dos_hist[k][old_bin].h_of_l ++;
				dos_hist[k][old_bin].g_of_l += sim_dos[k].mod_f;
				flat_histogram(k);
				end_to_end(k,SITE1,SITE2);				// recompute to get back to the old
			}//rejected (due to arg)
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
			#endif
			#ifdef DLIST
				nl_check(lb,ub,k);
				if(nl_flag[k] == 1) nblist_pivot(k,ub);
			#endif
			#endif
			dos_hist[k][old_bin].h_of_l ++;
			dos_hist[k][old_bin].g_of_l += sim_dos[k].mod_f;
			flat_histogram(k);
			end_to_end(k,SITE1,SITE2);				// recompute to get back to the old
		}//rejected (outside range)
		calcvalue(k);
		dos_svalues(k);
}

#endif

		
