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
double dos_interp (int,int,double);
#ifdef STATUS
void curr_status  (int,int);
#endif
void end_to_end	  (int,int,int);
#ifdef DLIST
void nl_check	  (int,int,int);
#endif


void xedos_linear(int ibox)
{
	int k	  =	ibox;

	

	/* ---------------------------------------------------- */
	/* First, pick a random solute molecule and then a      */
	/* random site on that molecule to move.                */
	/* ---------------------------------------------------- */
	int lb;
	int ub;
	int rand_site =0;

#ifdef IONC
	/* ----------------------------------------------------	*/
	/* If IONC is defined, the random move is either done	*/
	/* on protein-Cwater system or the relevant ion (SITE2)	*/
	/* ----------------------------------------------------	*/
	double ran_num	= ran2();
	if (ran_num<0.8){
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
	}
	else{
		rand_site	=	SITE2;
		lb = rand_site;
		ub = rand_site+1;
	}
#endif//IONC


#ifndef IONC
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

#endif //no IONC	

	/* Move the selected site randomly in x,y,z*/

	double dx, dy,dz;
	dx = (ran2()-0.5)*mc_rand.delta[k]; 
	dy = (ran2()-0.5)*mc_rand.delta[k];
	dz = (ran2()-0.5)*mc_rand.delta[k];
	

#ifdef MC
	forces_mc(rand_site,rand_site+1,0,k);
#endif
	double L_initial;
	end_to_end(k,SITE1,SITE2);						/* Distance betweeb SITE1 and SITE2		*/
	
	L_initial		= ordparam[k].d_nc;
	double U_initial;
	double U_final;
	double L_final;
	
	int old_bin;
	int new_bin;
	double beta_old	= 1.0/sim.kT[k];			
	double beta_new	= beta_old;
	calcvalue(k);


	atom_temp[k][rand_site]		= atom[k][rand_site];				/* Back up coordinates with pdb			*/
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
#endif
#ifndef MC		
		U_initial		= en[k].potens;
#endif
		old_bin			= (int) ((L_initial - sim_dos[k].l_begin)/sim_dos[k].l_width);			

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
	    nblist_pivot(k,rand_site+1); // call Nblist 
	#endif
	#ifdef DLIST
		nl_check(lb, ub, k);
		if(nl_flag[k] == 1) nblist_pivot(k,rand_site+1);
	#endif
	#endif
		/* ----------------------------------------------- */
		/* Update the forces.                              */
		/* ----------------------------------------------- */
#ifdef MC
		forces_mc(rand_site, rand_site+1, 1,k);
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
				mc_rand.rand_acc[k] ++;
			}//accepted
			else if (exp(arg) > ran2()){
				dos_hist[k][new_bin].h_of_l ++;
				dos_hist[k][new_bin].g_of_l += sim_dos[k].mod_f;
				flat_histogram(k);
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
					nblist_pivot(k,rand_site+1); // call Nblist 
				#endif
				#ifdef DLIST
					nl_check(lb, ub, k);
					if(nl_flag[k] == 1) nblist_pivot(k,rand_site+1);
				#endif
				#endif
				dos_hist[k][old_bin].h_of_l ++;
				dos_hist[k][old_bin].g_of_l += sim_dos[k].mod_f;
				flat_histogram(k);
				end_to_end(k,SITE1,SITE2);				// recompute to get back to the old
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
				nblist_pivot(k,rand_site+1); // call Nblist 
			#endif
			#ifdef DLIST
				nl_check(lb, ub, k);
				if(nl_flag[k] == 1) nblist_pivot(k,rand_site+1);
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

		
