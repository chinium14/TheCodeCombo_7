#ifdef DOS
/* ======================================================================== */
/* dos.cpp			                                                        */
/* This suboutine performs density of states simulation. In the beginning a */
/* random walk is performed in U by assigning initial veolocities based on  */
/* a beta of the concerned box. In additon to the pivot moves, swap moves   */
/* are implemented to facilitate sampling and also to achieve faster		*/
/* convergence.																*/		
/* ======================================================================== */
	
#include "defines.h"

void init_block			(int);
void ofile				(void);
double ran2				(void);
void vinit_hyb			(int);
void xyz_config			(int,unsigned long);
void movie				(int,unsigned long);
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
void dos_assoc_trans (int);
int  swap_dos			(unsigned long,int);
void save_config		(unsigned long, unsigned long);
void save_veloc			(unsigned long,unsigned long);
void dos_awrite			(unsigned long,int);
void dos_wwrite			(unsigned long,int);
void ioflush			(int);
void force_short                (int);
double dos_interp (int, int, double);
#ifdef TLATE
void translate			(int);
#endif
void integrate_vverlet	(int, double);
void integrate_dos_mts(int,double);
#ifdef MPI
int	 swap_mpi			(unsigned long,int,int);
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
#ifdef STATS
void dos_stats			(int);
#endif
#ifdef TRR
void write_trr(int, unsigned long, int);
#endif
void dos()
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
	char	name5[100];

	sim_hyb.flag=0;

	int status1=0;
#ifdef STATS
	for(int k=0; k<sim.NB; k++) sim_dos[k].data_count =0;
#endif

	/*------------------------------------*/
	/* If the inital energy is outside	  */
	/* the defined range; exit right away */
	/*------------------------------------*/

	for(int k=0; k<sim.NB; k++) if (en[k].potens <sim_dos[k].e_begin || 
		en[k].potens > sim_dos[k].e_end){
		fprintf(stdout,"intial energy (%f) in box %i is outside the energy window (%f, %f) !!!\n", en[k].potens, k, sim_dos[k].e_begin, sim_dos[k].e_end);
		exit(1);
	}
		  
	for(int k=0; k<sim.NB; k++) vscal[k] = 1.0;// initial scale for velocity


  /* ================================================================== */
  /* Begin production iterations.                                       */
  /* ================================================================== */
	for(unsigned long i=1; i<=sim.cyc_pr; i++) { // DOS counter
		n_iter=i;
		for(int k=0; k<sim.NB; k++) {
			#ifdef STATUS
				curr_status(k,6);
			#endif
			double random = ran2();
      if(i==129){
        int junk = 1;
      }
			/*------------------------------------------*/
			/* If trans move is picked call trans		*/
			/*------------------------------------------*/
			
			if (random< mc_trans.PMS[k]) trans(k);
				
			/*------------------------------------------*/
			/* If pivot move is picked call pivot		*/
			/*------------------------------------------*/

			else if (random < mc_pivot.PMS[k]+mc_trans.PMS[k]) pivot(k);

      else if (random < mc_rand.PMS[k]+mc_pivot.PMS[k]+mc_trans.PMS[k]) dos_assoc_trans(k);
			
			/*------------------------------------------*/
			/* otherwise call a hybrid MD-MC move		*/
			/*------------------------------------------*/
			else{
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
			  /* Calculate the new properties of the system		 */
			  /* ----------------------------------------------- */
				kinet(k);    											
				calcvalue(k);
				U_final	= en[k].potens;
				K_final	= en[k].kinet;

			  /* ----------------------------------------------- */
			  /* Check the acceptance of the move. First check   */
        /* to see if the system moved out of range, then   */
        /* check to see if it is accepted.                 */
			  /* ----------------------------------------------- */

		    if(U_final >= sim_dos[k].e_begin && U_final < sim_dos[k].e_end){
          new_bin	=  (int) ((U_final - sim_dos[k].e_begin)/sim_dos[k].e_width);
          double g_old = dos_interp(k,old_bin,U_initial);
          double g_new = dos_interp(k,new_bin,U_final);
          double arg = g_old - g_new + K_initial*beta_old - K_final*beta_new;

				  if (arg>0){
						  dos_hist[k][new_bin].h_of_e ++;
						  dos_hist[k][new_bin].g_of_e += sim_dos[k].mod_f;
						  flat_histogram(k);
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
				  else if (exp(arg) > ran2()){
						  dos_hist[k][new_bin].h_of_e ++;
						  dos_hist[k][new_bin].g_of_e += sim_dos[k].mod_f;
						  flat_histogram(k);
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
					  dos_hist[k][old_bin].h_of_e ++;
					  dos_hist[k][old_bin].g_of_e += sim_dos[k].mod_f;
					  flat_histogram(k);
				  }//rejected
        }
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
					  dos_hist[k][old_bin].h_of_e ++;
					  dos_hist[k][old_bin].g_of_e += sim_dos[k].mod_f;
					  flat_histogram(k);
        }//rejected

				calcvalue(k); // needed to update box.variable correctly
				dos_svalues(k);
			} // dos with md move ends here corresponds to else
		}//k loop ends here
  /* ----------------------------------------------- */
  /* Propose a swap move and check acceptance.		 */
  /* This when MPI is not defined					 */
  /* ----------------------------------------------- */		
		#ifndef MPI
		if((i%sim_hyb.cyc_swap)==0 && sim.NB>=1){
			int n_success=0;
			struct{
				int box1;
				int box2;
			}swapped_boxes[NBOXES];

			if(sim_hyb.flag==0){
				int k=0;
				while ((k+1)<sim.NB){
					int success =swap_dos(i,k);
					if (success ==1) {
						swapped_boxes[n_success].box1 = k;
						swapped_boxes[n_success].box2 = k+1;
						n_success++;
					}
					k=k+2;
				}
				sim_hyb.flag=1;
			}
			else{
				int k=1;
				while ((k+1)<sim.NB){
					int success =swap_dos(i,k);
					if (success ==1) {
						swapped_boxes[n_success].box1 = k;
						swapped_boxes[n_success].box2 = k+1;
						n_success++;
					}
					k=k+2;
				}
				sim_hyb.flag=0;
			}

			if(n_success != 0){
				for(int j=0; j<n_success; j++){
					swap_ptr[swapped_boxes[j].box1] += sprintf(swap_ptr[swapped_boxes[j].box1],
						"%d iterations    flag=%d box%d &box%d swapped \n",
						i,sim_hyb.flag,swapped_boxes[j].box1,swapped_boxes[j].box2);
				}
				ioflush(0);
			}

		}// end swap counter here
		#else
		/*------------------------------------------------------*/
		/* If MPI is defined then go to swap_mpi and later		*/
		/* modify the value of flag from 0 to 1 or 1 to 0		*/
		/* If flag = 0 the process sends the message if flag is */
		/* one the process receive a message from previous rank */
		/*------------------------------------------------------*/
		if((i%sim_hyb.cyc_swap)==0 ){
			int n_success = 0;
			for(int k=0; k<sim.NB; k++){
				int success=swap_mpi(i,mpi.flag,0);
				if (mpi.flag==0) mpi.flag =1; else mpi.flag=0;	
				if (success ==1){
					swap_ptr[k] += sprintf(swap_ptr[k],
					" At iteration %d	box with rank %d swapped \n",i,mpi.my_rank);		
					n_success++;
				}			}
			if(n_success > 0) ioflush(0);
		}
		#endif //ifndef MPI

		/* ----------------------------------------------- */
		/* Check to see if stop.txt is present. If it is   */
		/* set the counter to exit the next iteration.     */
		/* ----------------------------------------------- */
		if((i%sim.blockd) == 0 || i==1) {
		  #ifdef WIN
		  struct _stat buf;
		  if(_stat("./INPUT/stop.txt",&buf) == 0) i=sim.cyc_pr+1;
		  #else
		  struct stat buf;
		  if(stat("./INPUT/stop.txt",&buf) == 0) i=sim.cyc_pr+1;
		  #endif
		}


		/* ----------------------------------------------- */
		/* Write the calculated data to the I/O buffers.   */
		/* ----------------------------------------------- */
		if((i%sim.blockc) == 0 || i==1) {
			dos_vblock(i);
			dos_awrite(i,1);
			ioflush(0);
		}
		
		/* ----------------------------------------------- */
		/* Write the gromacs .trr file and coordinate and  */
		/* velocity files with fequency blockt.            */
		/* ----------------------------------------------- */


		if((i%sim.blockt) == 0) {
			#ifdef TLATE
			for (int k =0; k<sim.NB; k++) translate(k);
			#endif

			#ifdef ROTN
			for (int k =0; k<sim.NB; k++) rotate_xyz(k,_AXIS);	
			#endif

			#ifdef TRR	
			for (int k =0; k<sim.NB; k++) write_trr(k,i,0);
			#endif

			dos_wwrite(i,1);
			save_config(i,sim.cyc_pr);								
			save_veloc(i,sim.cyc_pr);								
			xyz_config(1,i);
			movie(1,i);
			pbc_chk(i);	
		}

 
		#ifndef MPI		
		status1 =1;
		for(int k=0; k<sim.NB; k++) {
			if (sim_dos[k].mod_f>STOP_F) status1=0;
		}
		if (status1 ==1){

			for(int k=0; k<sim.NB; k++) {
				FILE *fptr_hist;


				sprintf(name5,"./OUTPUT/DOS/flathist%d.txt",k);  

				fptr_hist= fopen(name5,"w");
				for (int i =0; i<sim_dos[k].e_bins; i++){
					fprintf(fptr_hist, "%lu,%lf, %lf, %d\n", i, dos_hist[k][i].e_mid, dos_hist[k][i].g_of_e, dos_hist[k][i].h_of_e);
				}
				fclose(fptr_hist);
			}
			break;
		}//if status1==1
		#endif
/*----------------------------------------------------------------- */
/* I need to  add some message passing to make sure that the		*/	
/* simulation terminates only when all the processes have converged	*/
/*----------------------------------------------------------------- */

    }// end production run	MC counter i														
  /* ----------------------------------------------- */
  /* Call output on the last iteration everytime.    */
  /* ----------------------------------------------- */
	dos_awrite(sim.cyc_pr,1);
	dos_wwrite(sim.cyc_pr,1);
	ioflush(0);
  }

  
void flat_histogram(int k)
{
	double sum=0; double av_energy; int flag1=1; // int flag2=1;
	for (int i =0; i<sim_dos[k].e_bins; i++){
		sum = sum + dos_hist[k][i].h_of_e;
	}
	av_energy = sum/sim_dos[k].e_bins;
	for (int i =0; i<sim_dos[k].e_bins; i++){
		 if (dos_hist[k][i].h_of_e < (sim_dos[k].flat* av_energy) )		flag1=0;
	}
	if (flag1==1){
		if sim_dos[k].mod_f > 10E-10
		  sim_dos[k].mod_f = sim_dos[k].mod_f/2.0;
#ifdef STATS
        if(sim_dos[k].data_count < 50) dos_stats(k);
#endif

		if(sim_dos[k].mod_f >=STOP_F){
			for (int i =0; i<sim_dos[k].e_bins; i++) dos_hist[k][i].h_of_e=	0;
		}
		if(sim_dos[k].mod_f >=RESET_F)	init_block(k);
		 
	}
}
#endif

