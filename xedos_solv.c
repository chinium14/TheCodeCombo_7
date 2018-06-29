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
void nblist_pivot (int,int);
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
void nl_check	  (int,int,int);
#endif


void xedos_solv(int ibox)
{
#define PMS_tr		0.5					// Probability of translation vs Rotation

#ifndef SLIST
	fprintf(stdout,"You cannot use solvent move without SLIST defined.\n");
	exit(111);
#endif
#ifdef SLIST
		int k	  =	ibox;
	/* ---------------------------------------------------- */
	/* First, pick a random solute molecule and then a      */
	/* random site on that molecule to move.                */
	/* ---------------------------------------------------- */
	int lb;
	int ub;
	int rand_sol =0;

	if(molec.Nsolute == 1){		
		lb = 0;
		ub = molec.lsite[0];
		rand_sol		=   ran_int(lb,ub);
	}
	else{
		int molecule = ran_int(0, molec.Nsolute);
		lb = molec.fsite[molecule];//first site of molecule
		ub = molec.lsite[molecule];//first site of next molecule
		rand_sol		=   ran_int(lb,ub);
	}
	/* ---------------------------------------------------	*/
	/* Now select a water molecule in the vicinity of the	*/
	/* selected solute site.								*/
	/* ----------------------------------------------------	*/

	int solvent		=	ran_int(0,slist[k].count[rand_sol]);
	int rand_site	=	slist[k].list[rand_sol][solvent];
	lb	=	rand_site;
	ub	=	rand_site+3;

	/* ---------------------------------------------------- */
	/* Now select whether a translation or rotational move	*/
	/* to perform based on probability listed on top		*/
	/* ----------------------------------------------------	*/
	int tr_flag =0;
	if (ran2()<PMS_tr) tr_flag =1;
	else tr_flag=0;
	double pos_x, pos_y, pos_z; double dr2=0.0;
	double dx, dy,dz;
	if(tr_flag ==1){								/* Move the selected site randomly in x,y,z	*/
	
		dx = (ran2()-0.5)*mc_solv.delta[k]; 
		dy = (ran2()-0.5)*mc_solv.delta[k];
		dz = (ran2()-0.5)*mc_solv.delta[k];

		pos_x = atnopbc[k][rand_site].x+dx;
		pos_y = atnopbc[k][rand_site].y+dy;
		pos_z = atnopbc[k][rand_site].z+dz;
		dr2	=	(atnopbc[k][rand_sol].x-pos_x)*(atnopbc[k][rand_sol].x-pos_x)+
					(atnopbc[k][rand_sol].y-pos_y)*(atnopbc[k][rand_sol].y-pos_y)+
					(atnopbc[k][rand_sol].z-pos_z)*(atnopbc[k][rand_sol].z-pos_z);
	}
	else{											/* Perform a rotation move on the water about oxygen	*/
		dx = atnopbc[k][rand_site].x;
		dy = atnopbc[k][rand_site].y;
		dz = atnopbc[k][rand_site].z;
	}


	double rc	=	sim.rc + (sim.rclist - sim.rc)*0.5;
	double rc2	= rc*rc;

	
	double L_initial;
	end_to_end(k,SITE1,SITE2);						/* Distance betweeb SITE1 and SITE2		*/
	L_initial		= ordparam[k].d_nc;
	int old_bin			= (int) ((L_initial - sim_dos[k].l_begin)/sim_dos[k].l_width);
	
	if(dr2<rc2){
		
	#ifdef MC
		forces_mc(lb,ub,0,k);
	#endif
		
		double U_initial;
		double U_final;
		double L_final;
		
		
		int new_bin;
		double beta_old	= 1.0/sim.kT[k];			
		double beta_new	= beta_old;
		calcvalue(k);

		for(int l=lb; l<ub; l++){
			atom_temp[k][l]	= atom[k][l];				/* Back up coordinates with pdb			*/
			atnopbc_temp[k][l]	= atnopbc[k][l];		/* Back up coordinates without pdb		*/
		}

		for(int l=0; l< box[k].boxns; l++){
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
						
			if(tr_flag ==1){
				for(int l=lb; l<ub; l++){
					atnopbc[k][l].x += dx;
					atnopbc[k][l].y += dy;
					atnopbc[k][l].z += dz;
					atom[k][l].x	+= dx;
					atom[k][l].y	+= dy;
					atom[k][l].z	+= dz;
				}
			}
			else{
				/* ---------------------------------------------------- */
				/* Here are some variable needed for a pivot move.      */
				/* ---------------------------------------------------- */
				double theta	=	(ran2()-0.5)*PI*2.0*0.2;
				double costheta	=	cos(theta);
				double sintheta	=	sin(theta);
				
				/*--------------------------------------*/
				/* INITIALIZATION OF PIVOT MATRICES		*/
				/*--------------------------------------*/

				int ii,jj, kk,ll;	
				for(ii=0; ii<3; ii++){
					for(jj=0; jj<2; jj++){
						for(kk=0; kk<3; kk++){
								for(ll=0; ll<3; ll++)
									PIVOT[ii][jj][kk][ll]=0;
							}
					}
				}
				// intialize all 54 enteries

				
				PIVOT[0][0][0][0]= 1.0;			
				PIVOT[0][0][1][1]= costheta;
				PIVOT[0][0][1][2]=-sintheta;
				PIVOT[0][0][2][1]= sintheta;
				PIVOT[0][0][2][2]= costheta;

				PIVOT[1][0][0][0]= costheta;
				PIVOT[1][0][0][2]= sintheta;	
				PIVOT[1][0][1][1]= 1.0;
				PIVOT[1][0][2][0]=-sintheta;
				PIVOT[1][0][2][2]= costheta;
				
				PIVOT[2][0][0][0]= costheta;
				PIVOT[2][0][0][1]= sintheta;
				PIVOT[2][0][1][0]=-sintheta;
				PIVOT[2][0][1][1]= costheta;
				PIVOT[2][0][2][2]= 1.0;

				int AXIS	=	ran_int(0,3);
				int ANGLE	=	0;

				pos_x	=	PIVOT[AXIS][ANGLE][0][0]*dx + PIVOT[AXIS][ANGLE][0][1]*dy + PIVOT[AXIS][ANGLE][0][2]*dz;
				pos_y	=	PIVOT[AXIS][ANGLE][1][0]*dx + PIVOT[AXIS][ANGLE][1][1]*dy + PIVOT[AXIS][ANGLE][1][2]*dz;
				pos_z	=	PIVOT[AXIS][ANGLE][2][0]*dx + PIVOT[AXIS][ANGLE][2][1]*dy + PIVOT[AXIS][ANGLE][2][2]*dz;

				for(int i=lb+1; i<ub; i++){

					pos_x	=	PIVOT[AXIS][ANGLE][0][0]*(atnopbc[k][i].x- dx)+ 
											PIVOT[AXIS][ANGLE][0][1]*(atnopbc[k][i].y-dy)+
											PIVOT[AXIS][ANGLE][0][2]*(atnopbc[k][i].z-dz)+ dx;
					pos_y	=	PIVOT[AXIS][ANGLE][1][0]*(atnopbc[k][i].x-dx)+
											PIVOT[AXIS][ANGLE][1][1]*(atnopbc[k][i].y-dy)+
											PIVOT[AXIS][ANGLE][1][2]*(atnopbc[k][i].z-dz)+ dy;
					pos_z	=	PIVOT[AXIS][ANGLE][2][0]*(atnopbc[k][i].x-dx)+
											PIVOT[AXIS][ANGLE][2][1]*(atnopbc[k][i].y-dy)+
											PIVOT[AXIS][ANGLE][2][2]*(atnopbc[k][i].z-dz)+ dz;

					atnopbc[k][i].x	=	pos_x;
					atnopbc[k][i].y	=	pos_y;
					atnopbc[k][i].z	=	pos_z;
					atom[k][i].x	=	pos_x;
					atom[k][i].y	=	pos_y;
					atom[k][i].z	=	pos_z;
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
			nl_check(lb, ub, k);
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
		  double arg = g_old- g_new + (U_initial-U_final)*beta_old; 

			if (arg>0.0){
				dos_hist[k][new_bin].h_of_l ++;
				dos_hist[k][new_bin].g_of_l += sim_dos[k].mod_f;
				flat_histogram(k);
				mc_solv.solv_acc[k] ++;
			}//accepted
			else if (exp(arg) > ran2()){
				dos_hist[k][new_bin].h_of_l ++;
				dos_hist[k][new_bin].g_of_l += sim_dos[k].mod_f;
				flat_histogram(k);
				mc_solv.solv_acc[k] ++;
			}//accepted
			else{

				for(int l=lb; l<ub; l++){
					atom[k][l]		= atom_temp[k][l];	 /* Back up coordinates with pdb		*/
					atnopbc[k][l]	= atnopbc_temp[k][l];/* Back up coordinates without pdb		*/
				}
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
				#endif
				#ifdef DLIST
					nl_check(lb, ub, k);
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
			#endif
			#ifdef DLIST
				nl_check(lb, ub, k);
				if(nl_flag[k] == 1) nblist_pivot(k,ub);
			#endif
			#endif
			dos_hist[k][old_bin].h_of_l ++;
			dos_hist[k][old_bin].g_of_l += sim_dos[k].mod_f;
			flat_histogram(k);
			end_to_end(k,SITE1,SITE2);				// recompute to get back to the old
		}//rejected (outside range)
	}
	else{
		dos_hist[k][old_bin].h_of_l ++;
		dos_hist[k][old_bin].g_of_l += sim_dos[k].mod_f;
		flat_histogram(k);
		end_to_end(k,SITE1,SITE2);				// recompute to get back to the old
	}// rejected because selected solvent moved out of the defined radius
 	calcvalue(k);
	dos_svalues(k);
#endif//SLIST
}

#endif

		
