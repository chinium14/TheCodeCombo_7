#ifdef FX_EDOS
	
#include "defines.h"

void init_block			(int);
void ofile				(void);
double ran2				(void);
void vinit_hyb			(int);
void xyz_config			(int,int);
void movie				(int,int);
void dos_svalues		(int);
void dos_vblock			(unsigned long);
void pbc_all			(int);
void pbc_chk			(int);
void kinet				(int);
void calcvalue			(int);
void nblist				(int);
void flat_histogram		(int);
void pivot				(int);
void trans				(int);
void save_config		(unsigned long,unsigned long);
void save_veloc			(unsigned long,unsigned long);
void dos_output			(int,int);
#ifdef TLATE
void translate			(int);
#endif
void integrate_vverlet	(int, double);
#ifdef MPI
int	 swap_mpi			(int,int,int);
#endif
#ifdef STATUS
void curr_status		  (int,int);
#endif
#ifdef CONFIGT
void configtemp(int);
#endif
#ifdef ROTN
void rotate_xyz(int,int);
#endif
#ifdef DLIST
void nl_check(int,int,int);
#endif

void fx_edos_hmc(int k)
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
	vscal[k] = 1.0;// initial scale for velocity


				sim.T[k]	= sim_dos[k].T_begin;
				sim.kT[k]	= RG * sim.T[k]*.001;
				beta_old	= 1.0/sim.kT[k];			
				beta_new	= beta_old;
				vinit_hyb(k);
				kinet(k);
				calcvalue(k);

				/* ------------------------------------	*/
				/* Back up coords, force, virial and	*/
				/* energy components.					*/
				/* ------------------------------------	*/

				for(int l =0; l< box[k].boxns; l++){
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
				old_bin			=  (int) ((U_initial - sim_dos[k].e_begin)/sim_dos[k].e_width);			
		
				for (int j=0;j < sim_hyb.cyc_hybmd; j++){
					integrate_vverlet(k, vscal[k]);
				#ifdef DLIST
					nl_check(0,box[k].boxns,k);
					if(nl_flag[k] == 1) nblist(k);
				#endif

				#ifdef STATUS
					curr_status(k,7);
				#endif
				}// j loop ends here
			  
			  /* ----------------------------------------------- */
			  /* Calculate the new properties of the system		 */
			  /* ----------------------------------------------- */
				kinet(k);    											
				calcvalue(k);
				U_final	= en[k].potens;
				K_final	= en[k].kinet;
				new_bin	=  (int) ((U_final - sim_dos[k].e_begin)/sim_dos[k].e_width);

			  /* ----------------------------------------------- */
			  /* Propose a MC move and check for acceptance. If  */
			  /* If the moves is accepted then accumulate the	 */
			  /* calculated values for block averages and error  */
			  /* bars else first return to original values and   */
			  /* then accumulate statistics						 */
			  /* ----------------------------------------------- */

				if (U_final >= sim_dos[k].e_begin && U_final < sim_dos[k].e_end && (dos_hist[k][new_bin].w_of_e + 
						K_initial*beta_old)> (dos_hist[k][old_bin].w_of_e + K_final*beta_new)){

						dos_hist[k][new_bin].nw_of_e ++;

						/* Check if the sign of walker is changed and update number of trips	*/
						
						if(new_bin==sim_dos[k].e_bins-1){
							if (sim_dos[k].sign_walker == -1) sim_dos[k].trips_num ++;
							sim_dos[k].sign_walker=1;
						}
						else if(new_bin==0){
							if (sim_dos[k].sign_walker == 1) sim_dos[k].trips_num ++;
							sim_dos[k].sign_walker=-1;
						}

						if(sim_dos[k].sign_walker ==1) dos_hist[k][new_bin].np_of_e ++;
						else dos_hist[k][new_bin].nm_of_e ++;

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
				
				else if ( U_final >= sim_dos[k].e_begin && U_final < sim_dos[k].e_end &&
						exp((dos_hist[k][new_bin].w_of_e + K_initial*beta_old)-(dos_hist[k][old_bin].w_of_e + K_final*beta_new)) > ran2()){
						dos_hist[k][new_bin].nw_of_e ++;

						/* Check if the sign of walker is changed and update number of trips	*/
						
						if(new_bin==sim_dos[k].e_bins-1){
							if (sim_dos[k].sign_walker == -1) sim_dos[k].trips_num ++;
							sim_dos[k].sign_walker=1;
						}
						else if(new_bin==0){
							if (sim_dos[k].sign_walker == 1) sim_dos[k].trips_num ++;
							sim_dos[k].sign_walker=-1;
						}

						if(sim_dos[k].sign_walker ==1) dos_hist[k][new_bin].np_of_e ++;
						else dos_hist[k][new_bin].nm_of_e ++;
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
					/* ------------------------------------	*/
					/* Go back to the original backed up	*/
					/* quantities as move is rejected.		*/
					/* ------------------------------------	*/
					for(int l =0; l< box[k].boxns; l++){
						atom[k][l]		= atom_temp[k][l];		
						atnopbc[k][l]	= atnopbc_temp[k][l];	
						ff[k][l]		= ff_temp[k][l];		
						vv[k][l]		= vv_temp[k][l];
						uu[k][l]		= vv[k][l];
					}
					#ifdef PRESSURE
					pvir[k]		= pvir_temp[k];					
					#endif
					en[k]		= en_temp[k];					
					dos_hist[k][old_bin].nw_of_e ++;

					/* Check if the sign of walker is changed and update number of trips	*/
					
					if(old_bin==sim_dos[k].e_bins-1){
						if (sim_dos[k].sign_walker == -1) sim_dos[k].trips_num ++;
						sim_dos[k].sign_walker=1;
					}
					else if(old_bin==0){
						if (sim_dos[k].sign_walker == 1) sim_dos[k].trips_num ++;
						sim_dos[k].sign_walker=-1;
					}

					if(sim_dos[k].sign_walker ==1) dos_hist[k][old_bin].np_of_e ++;
					else dos_hist[k][old_bin].nm_of_e ++;
				}//rejected

				calcvalue(k); // needed to update box.variable correctly
				dos_svalues(k);
}
#endif

