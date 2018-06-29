#ifdef XEDOS
/* ======================================================================== */
/* xedos.cpp			                                                    */
/* This suboutine performs density of states simulation. In the beginning a */
/* random walk is performed in U by assigning initial veolocities based on  */
/* a beta of the concerned box. In additon to the pivot moves, swap moves   */
/* are implemented to facilitate sampling and also to achieve faster		*/
/* convergence.																*/
/* Revised by Thomas Knotts 9 Sep 04 for I/O buffereing and clean up.       */		
/* ======================================================================== */
	
#include "defines.h"

void init_block			(int);
double ran2				(void);
void xyz_config			(int,unsigned long);
void movie				(int,unsigned long);
void dos_vblock			(unsigned long);
void pbc_chk			(int);
void flat_histogram		(int);
int  swap_xe_dos		(unsigned long,int);
void save_config		(unsigned long,unsigned long);
void save_veloc			(unsigned long,unsigned long);
void dos_awrite			(unsigned long,int);
void dos_wwrite			(unsigned long,int);
void ioflush			(int);
void xedos_hmc			(int);
void xedos_axial		(int);
void xedos_scale		(int);
void xedos_trans		(int);
void xedos_linear		(int);
void xedos_pivot		(int);
void xedos_solv			(int);
void xedos_cosol		(int);
void end_to_end			(int,int,int);

#ifdef TLATE
void translate			(int);
#endif

#ifdef MPI
int	 swap_xe_mpi		(unsigned long,int,int);
#endif

#ifdef STATUS
void curr_status		(int,int);
#endif

#ifdef STATS
void dos_stats			(int);
#endif

#ifdef IONC
int avrg_xe_mpi			(int, int, int);
#endif

#ifdef ROTN
void rotate_xyz(int,int);
#endif

#ifdef TRR
void write_trr(int, unsigned long, int);
#endif

void xedos()
{
	char	name5[100];
	sim_hyb.flag=0;

	int status1=0;
	#ifdef STATS
	for(int k=0; k<sim.NB; k++) sim_dos[k].data_count =0;
	#endif

	/*------------------------------------*/
	/* If the length energy is outside	  */
	/* the defined range; exit right away */
	/*------------------------------------*/
	
	for(int k=0; k<sim.NB; k++){
		end_to_end(k,SITE1,SITE2);
		if (ordparam[k].d_nc <sim_dos[k].l_begin || ordparam[k].d_nc > sim_dos[k].l_end){
			fprintf(stdout,"intial length (%lf) is outside the window for box %d (0 to N-1)!!!\n",ordparam[k].d_nc,k);
			exit(1);
		}
		sim.T[k]	= sim_dos[k].T_begin;
		sim.kT[k]	= RG * sim.T[k]*.001;
	}

  /* ================================================================== */
  /* Begin production iterations.                                       */
  /* ================================================================== */
	for(unsigned long i=1; i<=sim.cyc_pr; i++) { // XEDOS counter
		n_iter=i;
		for(int k=0; k<sim.NB; k++) {
			#ifdef STATUS
			curr_status(k,6);
			#endif
			double random = ran2();

			/*------------------------------------------*/
			/* If axial move is picked call xedos_axial	*/
			/*------------------------------------------*/
			if      (random<mc_pivot.PMS[k])														xedos_pivot(k);
			else if (random<mc_trans.PMS[k] + mc_pivot.PMS[k])										xedos_trans(k);
			else if (random<mc_axial.PMS[k] + mc_trans.PMS[k] + mc_pivot.PMS[k])					xedos_scale(k);
			else if (random<mc_rand.PMS[k]  + mc_axial.PMS[k] + mc_trans.PMS[k] + mc_pivot.PMS[k])  xedos_linear(k);
			#ifdef SLIST
			else if (random<mc_rand.PMS[k]  + mc_axial.PMS[k] + mc_trans.PMS[k] + mc_pivot.PMS[k] + mc_solv.PMS[k]) xedos_solv(k);
			#else
			else if (random<mc_rand.PMS[k]  + mc_axial.PMS[k] + mc_trans.PMS[k] + mc_pivot.PMS[k] + mc_solv.PMS[k]) xedos_cosol(k);
			#endif
			else xedos_hmc(k);
		}

		#ifdef WALL
		/*------------------------------------------*/
		/* check to make sure no sites have entered */
		/* the wall.                                */
		/*------------------------------------------*/
		for(int k=0; k<sim.NB; k++){
			for(int j=0; j<box[k].boxns; j++){
				if(atnopbc[k][j].z < wall[k].z[0]){
					fprintf(stdout,"Particle %d (z coord=%f) is beyond the wall (z coord=%f) on iteration %d\n",j,atnopbc[k][j].z,wall[k].z[0],i);
                                        exit(7923);
				}
			}
		}
		#endif

		/*------------------------------------------*/
		/* Now call a swap move						*/
		/*------------------------------------------*/
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
					int success =swap_xe_dos(i,k);
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
					int success =swap_xe_dos(i,k);
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
				int success=swap_xe_mpi(i,mpi.flag,0);
				if (mpi.flag==0) mpi.flag =1; else mpi.flag=0;	
				if (success ==1){
					swap_ptr[k] += sprintf(swap_ptr[k],
					" At iteration %d	box with rank %d swapped \n",i,mpi.my_rank);		
					n_success++;
				}
			}		
#ifdef IONC
/*			if((i%(sim_hyb.cyc_swap*2))==0 ){
				for(int k=0; k<sim.NB; k++){
					int flag;
					double q = (mpi.p-1)/2.0;
					if (mpi.my_rank<q) flag =0;
					else if (mpi.my_rank > q) flag =1;
					else flag =2;
					
					// -----------------------------------------------	//
					// If flag ==0, the process will send a message		//
					// If flag ==1, the process receives a message		//
					// if flag ==2, the process does nothing		//
					// ------------------------------------------------	//
					int suc =0;
					if (flag !=2) suc=avrg_xe_mpi(i,flag,0);				
				}
			}
			if(n_success > 0) ioflush(0);
		
*/
#endif //ifdef IONC

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
		/* Write the calculated data, configurations, and  */
		/* velocities to file with frequency sim.blockc.   */
		/* ----------------------------------------------- */
		if((i%sim.blockc) == 0 || i==1) {
			dos_vblock(i);
			dos_awrite(i,1);
			ioflush(0);
		}// blockc 

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
				for (int i =0; i<sim_dos[k].l_bins; i++){
					fprintf(fptr_hist, "%d,%lf, %lf, %d\n", i, dos_hist[k][i].l_mid, dos_hist[k][i].g_of_l, dos_hist[k][i].h_of_l);
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
#ifndef IONC
	double sum=0; double av_length; int flag1=1; // int flag2=1;
	for (int i =0; i<sim_dos[k].l_bins; i++){
		sum = sum + dos_hist[k][i].h_of_l;
	}
	av_length = sum/sim_dos[k].l_bins;
	for (int i =0; i<sim_dos[k].l_bins; i++){
		 if (dos_hist[k][i].h_of_l < (sim_dos[k].flat* av_length) )		flag1=0; 
	}
	if (flag1==1){
		sim_dos[k].mod_f = sim_dos[k].mod_f/2.0;
#ifdef STATS
	if(sim_dos[k].data_count < 50) dos_stats(k);
#endif
		if(sim_dos[k].mod_f >=STOP_F){
			for (int i =0; i<sim_dos[k].l_bins; i++) {
				dos_hist[k][i].h_of_l=	0;
			}
			if(sim_dos[k].mod_f >=RESET_F)	init_block(k);
		}		 
	}
#endif//NOIONC
#ifdef IONC
	double sum=0; double av_length; int flag1=1; int min =dos_hist[k][0].h_of_l;;
	for (int i =0; i<sim_dos[k].l_bins; i++){
		sum = sum + dos_hist[k][i].h_of_l;
		if (min > dos_hist[k][i].h_of_l ) min = dos_hist[k][i].h_of_l;
	}
	av_length = sum/sim_dos[k].l_bins;
	for (int i =0; i<sim_dos[k].l_bins; i++){
		 if (dos_hist[k][i].h_of_l < (sim_dos[k].flat* av_length) )		flag1=0; 
	}
	if (min > 5000) flag1=1;

	if (flag1==1){
		sim_dos[k].mod_f = sim_dos[k].mod_f/2.0;
#ifdef STATS
	if(sim_dos[k].data_count < 50) dos_stats(k);
#endif
		if(sim_dos[k].mod_f >=STOP_F){
			for (int i =0; i<sim_dos[k].l_bins; i++) {
				dos_hist[k][i].h_of_l=	0;
			}
			if(sim_dos[k].mod_f >=RESET_F)	init_block(k);
		}		 
	}
#endif//IONC

}
#endif

