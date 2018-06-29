#ifdef MMDOS
/* ======================================================================== */
/* mmdos.cpp			                                                    */
/* This suboutine performs density of states simulation. In the beginning a */
/* random walk is performed in U by assigning initial veolocities based on  */
/* a beta of the concerned box. In additon to the pivot moves, swap moves   */
/* are implemented to facilitate sampling and also to achieve faster		*/
/* convergence.																*/		
/* ======================================================================== */
	
#include "defines.h"


void init_block			(int);
void ofile				(void);
void outend				(int);
double ran2				(void);
void xyz_config			(int,unsigned long);
void svalues			(int);
void vblock				(int);
void calcvalue			(int);
void calc_g_of_e   		(int);
void mmdos_pivot		(int);
void save_config		(unsigned long,unsigned long);
void save_veloc			(unsigned long,unsigned long);
void output				(int,int);
void mmdos_hmc			(int);
void integrate_vverlet	(int, double);
void mmdos_dke			(int);
void translate			(int);
int swap_dos			(unsigned long,int);
#ifdef MPI
int	 swap_mpi			(unsigned long,int,int);
#endif

#ifdef RMSD
int  quatfit (int, int, int);
#endif


void mmdos(){

    char	name4[50];
	char	name5[100];

    /*----------------------------------------------------------*/
	/* open the swap_dos.out file to save the list of accepted	*/
	/* swaps and also to save the currect status of convergence */
	/* factors in different boxes.								*/
	/*----------------------------------------------------------*/
	
	FILE *swap2;  FILE *ter;
#ifndef MPI
	sprintf(name4,"./OUTPUT/DOS/swap_dos.out");              
	swap2= fopen(name4,"w");
	fclose(swap2);
#endif
#ifdef MPI
	sprintf(name4,"./OUTPUT/DOS/swap_dos%d.out",mpi.my_rank);              
	swap2= fopen(name4,"w");
	fclose(swap2);
#endif

	int status1=0;
	sim_hyb.flag=0;

#ifdef STATS
	for(int k=0; k<sim.NB; k++) sim_dos[k].data_count =0;
#endif

	/*------------------------------------*/
	/* If the inital energy is outside	  */
	/* the defined range; exit right away */
	/*------------------------------------*/

	for(int k=0; k<sim.NB; k++) if (en[k].totals <sim_dos[k].e_begin || en[k].totals > sim_dos[k].e_end){
		fprintf(stdout,"intial energy (%f) in box %i is outside the energy window (%f, %f) !!!\n", en[k].totals, k, sim_dos[k].e_begin, sim_dos[k].e_end);
		exit(1);
	}
		  
	for(int k=0; k<sim.NB; k++)	{
		for ( int i=0; i<sim_dos[k].e_bins; i++){
			dos_hist[k][i].k_of_e =	1.5 * box[k].boxns * RG * sim_dos[k].T_begin*.001; //3/2NkT
			dos_hist_temp[k][i].k_of_e = dos_hist[k][i].k_of_e ;

		}
	}


  /* ================================================================== */
  /*                                                                    */
  /* Begin production iterations.                                       */
  /*                                                                    */
  /* ================================================================== */
	for(unsigned long i=1; i<=sim.cyc_pr; i++) { // MMDOS counter
		n_iter=i;
		
		for(int k=0; k<sim.NB; k++) {
			double random = ran2();
			/*------------------------------------------*/
			/* If deltaKE move is picked call mmdos_dke */
			/*------------------------------------------*/
			if (random < sim_dos[k].p_dke) mmdos_dke(k);
			/*------------------------------------------*/
			/* If pivot move is picked call mmdos_pivot	*/
			/*------------------------------------------*/
			else if  (random < mc_pivot.PMS[k]) mmdos_pivot(k);
			/*------------------------------------------*/
			/* otherwise call the NVE HMC move			*/
			/*------------------------------------------*/
			else mmdos_hmc(k); 
		}//k loop ends here
  /* ----------------------------------------------- */
  /* Calculate the block averages accumulated in     */
  /* svalues() and reset the accumulators.           */
  /* ----------------------------------------------- */
		if((i%sim.blockd) == 0 || i==1) {								
			for(int k=0; k<sim.NB; k++) calc_g_of_e(k);
			vblock(i);
			for(int k=0; k<sim.NB; k++) init_block(k);
			outend(i);
			if( NULL != (ter=fopen("./INPUT/stop.txt","r"))){ 
				fclose(ter);
				(i= sim.cyc_pr+1);
			}
		}

  /* ----------------------------------------------- */
  /* Propose a swap move and check acceptance.		 */
  /* This when MPI is not defined					 */
  /* ----------------------------------------------- */		
#ifndef MPI
		if((i%sim_hyb.cyc_swap)==0 && sim.NB>=1){
			swap2= fopen(name4,"a");
			if(sim_hyb.flag==0){
				int k=0;
				while ((k+1)<sim.NB){
					int success =swap_dos(i,k);
					if (success ==1) {
						fprintf(swap2,"%d iterations	flag=%d	box%d &box%d swapped \n",i,sim_hyb.flag,k,k+1);
					}
					k=k+2;
				}
				sim_hyb.flag=1;
			}
			else {
			int k=1;
				while ((k+1)<sim.NB){
					int success =swap_dos(i,k);
					if (success ==1) {
						fprintf(swap2,"%d iterations	flag=%d	box%d &box%d swapped \n",i,sim_hyb.flag,k,k+1);
					}
					k=k+2;
				}
				sim_hyb.flag=0;
			}
			fprintf(swap2,"mod_f:	");
			for(int k =0; k<sim.NB; k++){
				fprintf(swap2,"%.8lf(%d)	",exp(sim_dos[k].mod_f),k);
			}
			fprintf(swap2,"\n");
			fclose(swap2);

		}// end swap counter here
#endif
#ifdef MPI
	/*------------------------------------------------------*/
	/* If MPI is defined then go to swap_mpi and later		*/
	/* modify the value of flag from 0 to 1 or 1 to 0		*/
	/* If flag = 0 the process sends the message if flag is */
	/* one the process receive a message from previous rank */
	/*------------------------------------------------------*/
		if((i%sim_hyb.cyc_swap)==0 ){
			swap2= fopen(name4,"a");
			for(int k=0; k<sim.NB; k++){
				int success=swap_mpi(i,mpi.flag,0);
				if (mpi.flag==0) mpi.flag =1; else mpi.flag=0;	
				if (success ==1) fprintf(swap2," At iteration %d	box with rank %d swapped \n",i,mpi.my_rank);				
			}
			fclose(swap2);
		}
#endif
    /* ----------------------------------------------- */
    /* Write the calculated data, configurations, and  */
    /* velocities to file with frequency sim.blockc.   */
    /* ----------------------------------------------- */
		if((i%sim.blockc) == 0 || i==1) {
			for(int k=0; k< sim.NB; k++) translate(k);
			output(i,1);											
			save_config(i,sim.cyc_pr);								
			save_veloc(i,sim.cyc_pr);								
			xyz_config(1,i);
#ifdef RMSD
		for(int k=0; k<sim.NB; k++) {
			quatfit(k,rmsd_modes.file,rmsd_modes.atoms); 
		} // loop k ends
#endif			
		}// blockc ends here

#ifndef MPI		
		status1 =1;
		for(int k=0; k<sim.NB; k++) {
			if (sim_dos[k].mod_f>1E-20) status1=0;
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
	output(sim.cyc_pr,1);
  }

  
void calc_g_of_e(int k)
{
//	dos_hist[k][0].g_of_e = 0.0;
	dos_hist[k][1].g_of_e = dos_hist[k][0].g_of_e + sim_dos[k].e_width*1.5*box[k].boxns
							/(dos_hist[k][1].k_of_e/dos_hist[k][1].h_of_e);

	for (int i =1; i<sim_dos[k].e_bins-1; i++){
		dos_hist[k][i+1].g_of_e = dos_hist[k][i-1].g_of_e + sim_dos[k].e_width*3.0
								* box[k].boxns/(dos_hist[k][i].k_of_e/dos_hist[k][i].h_of_e);
	}
	// compute mod_f based on current estimate of g_of_e
	double mod_f =0.0;
	for (int i =0; i<sim_dos[k].e_bins; i++){
		mod_f += ((dos_hist[k][i].g_of_e - dos_hist[k][0].g_of_e)- 
				 (dos_hist_temp[k][i].g_of_e - dos_hist_temp[k][0].g_of_e))*
				 ((dos_hist[k][i].g_of_e - dos_hist[k][0].g_of_e)- 
				 (dos_hist_temp[k][i].g_of_e - dos_hist_temp[k][0].g_of_e));
	}
	sim_dos[k].mod_f = mod_f/sim_dos[k].e_bins;
	for (int i =0; i<sim_dos[k].e_bins; i++){
		dos_hist_temp[k][i].g_of_e = dos_hist[k][i].g_of_e;
	}

}
#endif

