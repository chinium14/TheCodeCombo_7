#ifdef CTDOS
/* ======================================================================== */
/* ctdos.cpp			                                                    */
/* This suboutine performs density of states simulation. In the beginning a */
/* random walk is performed in U by assigning initial veolocities based on  */
/* a beta of the concerned box. In additon to the pivot moves, swap moves   */
/* are implemented to facilitate sampling and also to achieve faster		*/
/* convergence.																*/		
/* Revised by Thomas Knotts 9 Sep 04 for I/O buffering and clean up.        */
/* ======================================================================== */
	
#include "defines.h"

void init_block			(int);
double ran2				(void);
void xyz_config			(int,unsigned long);
void movie				(int,unsigned long);
void dos_vblock			(unsigned long);
void pbc_chk			(int);
void flat_histogram		(int);
int  swap_ct_dos		(unsigned long,int);
void save_config		(unsigned long,unsigned long);
void save_veloc			(unsigned long,unsigned long);
void dos_awrite			(unsigned long,int);
void dos_wwrite			(unsigned long,int);
void ioflush			(int);
void calc_gT_of_e		(int);
void ctdos_pivot		(int);
void ctdos_trans		(int);
void ctdos_hmc			(int);
void ctdos_assoc_trans (int);

#ifdef TLATE
void translate			(int);
#endif

#ifdef MPI
int	 swap_ct_mpi		(unsigned long,int,int);
#endif

#ifdef STATUS
void curr_status		(int,int);
#endif

#ifdef STATS
void dos_stats			(int);
#endif

#ifdef ROTN
void rotate_xyz(int,int);
#endif

#ifdef TRR
void write_trr(int, unsigned long, int);
#endif

void ctdos()
{
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

	for(int k=0; k<sim.NB; k++) if (en[k].potens <sim_dos[k].e_begin || en[k].potens > sim_dos[k].e_end){
		fprintf(stdout,"intial energy (%f) in box %i is outside the energy window (%f, %f) !!!\n", en[k].potens, k, sim_dos[k].e_begin, sim_dos[k].e_end);
		exit(1);
	}
		  
  /* ================================================================== */
  /* Begin production iterations.                                       */
  /* ================================================================== */
	for(unsigned long i=1; i<=sim.cyc_pr; i++) { // CTDOS counter
		n_iter=i;
    for(int k=0; k<sim.NB; k++) {
			#ifdef STATUS
			curr_status(k,6);
			#endif
			double random = ran2();

			/*------------------------------------------*/
			/* If trans move is picked call pivot		*/
			/*------------------------------------------*/
			
			if (i>sim.cyc_eq && random< mc_trans.PMS[k])			ctdos_trans(k);
				
			/*------------------------------------------*/
			/* If pivot move is picked call pivot		*/
			/*------------------------------------------*/

			else if (i>sim.cyc_eq && random < mc_pivot.PMS[k]+mc_trans.PMS[k])	ctdos_pivot(k);

			/*------------------------------------------*/
			/* assoc_trans move                         */
			/*------------------------------------------*/
      else if (i>sim.cyc_eq && random < mc_rand.PMS[k]+mc_pivot.PMS[k]+mc_trans.PMS[k])	ctdos_assoc_trans(k);
			
			/*------------------------------------------*/
			/* otherwise call a hybrid MD-MC move		*/
			/*------------------------------------------*/
			else ctdos_hmc(k);
		}


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
					int success =swap_ct_dos(i,k);
					if (success ==1) {
						swapped_boxes[n_success].box1 = k;
						swapped_boxes[n_success].box2 = k+1;
						n_success++;
					}
					k=k+2;
				}
				sim_hyb.flag=1;
			}
			else {
			int k=1;
				while ((k+1)<sim.NB){
					int success =swap_ct_dos(i,k);
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
		#else //ifndef MPI
		/*------------------------------------------------------*/
		/* If MPI is defined then go to swap_mpi and later		*/
		/* modify the value of flag from 0 to 1 or 1 to 0		*/
		/* If flag = 0 the process sends the message if flag is */
		/* one the process receive a message from previous rank */
		/*------------------------------------------------------*/
		if((i%sim_hyb.cyc_swap)==0 ){
			int n_success = 0;
			for(int k=0; k<sim.NB; k++){
				int success=swap_ct_mpi(i,mpi.flag,0);
				if (mpi.flag==0) mpi.flag =1; else mpi.flag=0;	
				if (success ==1){
					swap_ptr[k] += sprintf(swap_ptr[k],
					" At iteration %d	box with rank %d swapped \n",i,mpi.my_rank);		
					n_success++;
				}
			}
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
		}// blockd ends here


		/* ----------------------------------------------- */
		/* Write the calculated data, configurations, and  */
		/* velocities to file with frequency sim.blockc.   */
		/* ----------------------------------------------- */
		if((i%sim.blockc) == 0 || i==1) {
			for(int k=0; k< sim.NB; k++) calc_gT_of_e(k);
			dos_vblock(i);
			dos_awrite(i,1);
			ioflush(0);
		}

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
					fprintf(fptr_hist, "%d,%lf, %lf, %d\n", i, dos_hist[k][i].e_mid, dos_hist[k][i].g_of_e, dos_hist[k][i].h_of_e);
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
		sim_dos[k].mod_f = sim_dos[k].mod_f/2.0;
#ifdef STATS
	if(sim_dos[k].data_count < 50) dos_stats(k);
#endif
		if(sim_dos[k].mod_f >=STOP_F){
			for (int i =0; i<sim_dos[k].e_bins; i++) {
				dos_hist[k][i].h_of_e=	0;
				if(sim_dos[k].mod_f>=RESET_F){		// reset config temp num and den
					dos_hist[k][i].ct_num = 0.0;
					dos_hist[k][i].ct_den = 0.0;
				}
			}
			if(sim_dos[k].mod_f >=RESET_F)	init_block(k);
		}		 
	}
}
void calc_gT_of_e(int k)
{
	dos_hist[k][1].gT_of_e = dos_hist[k][0].gT_of_e + sim_dos[k].e_width/(RG*0.001*dos_hist[k][1].config_T);

	for (int i =1; i<sim_dos[k].e_bins-1; i++){
		dos_hist[k][i+1].gT_of_e = dos_hist[k][i-1].gT_of_e + sim_dos[k].e_width*2.0/(RG*0.001*dos_hist[k][i].config_T);
	}
	for (int i =0; i<sim_dos[k].e_bins; i++){
		dos_hist_temp[k][i].gT_of_e = dos_hist[k][i].gT_of_e;
	}

}
#endif

