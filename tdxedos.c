#ifdef TDXEDOS
/* ======================================================================== */
/* 2Dxedos.cpp			                                                    */
/* This suboutine performs an expanded ensemble density of states           */
/* simulation in 2 dimensions. It calls the various moves including swap    */
/* moves.                                                                   */
/*                                                                          */
/* Written by Thomas Knotts 3 June 04                                       */
/* Revised by Thomas Knotts 9 Sep 04 for I/O buffereing and clean up.       */		
/* ======================================================================== */

    /*----------------------------------------------------------*/
	/* Include statements and Function Prototypes               */
	/*----------------------------------------------------------*/
	
#include "defines.h"

void init_block         (int);
double ran2             (void);
void xyz_config         (int,unsigned long);
void movie              (int,unsigned long);
void dos_vblock         (unsigned long);
void pbc_chk            (int);
void flat_histogram     (int);
void save_config        (unsigned long,unsigned long);
void save_veloc         (unsigned long,unsigned long);
void dos_awrite         (unsigned long,int);
void dos_wwrite         (unsigned long,int);
void ioflush            (int);
void tdxedos_axial      (int);
void tdxedos_hmc        (int);
void tdxedos_linear     (int);
void tdxedos_pivot      (int);
void tdxedos_scale      (int);
void tdxedos_solv       (int);
void tdxedos_trans      (int);
void end_to_end         (int,int,int);
int  ran_int			(int,int);

#ifdef MPI
int	 swap_tdxe_mpi      (int,int,int);
#else
int  swap_tdxe_dos      (int,int);
#endif

#ifdef TLATE
void translate          (int);
#endif

#ifdef STATUS
void curr_status        (int,int);
#endif

#ifdef STATS
void dos_stats          (int);
#endif

#ifdef ROTN
void rotate_xyz         (int,int);
#endif

#ifdef TRR
void write_trr          (int, unsigned long, int);
#endif

#ifdef WALL
double wall_angle         (int, int);
#endif


void tdxedos()
{
	char name5[100];  //variable for flathist.txt file

	sim_hyb.flag=0;
	
	int status1=0;
	#ifdef STATS
	for(int k=0; k<sim.NB; k++) sim_dos[k].data_count =0;
	#endif

  /* ================================================================== */
  /* If the initial reaction coordinates are outside the ranges, exit.  */
  /* ================================================================== */
	
	for(int k=0; k<sim.NB; k++){

		/*------------------------------------*/
		/* Call subroutines that calculate    */
		/* the properties defining the rxn-   */
		/* coordinate and assign ordparam[k].1*/
		/* and ordparam[k].2.                 */
		/*------------------------------------*/
		end_to_end(k,SITE1,SITE2);
		ordparam[k].x1 = ordparam[k].d_nc;
		ordparam[k].x2 = wall_angle(k,wall[k].angle_site);

		/*------------------------------------*/
		/* Make the checks and exit.          */
		/*------------------------------------*/
		if (ordparam[k].x1 <sim_dos[k].x1_begin || ordparam[k].x1 > sim_dos[k].x1_end){
		  fprintf(stdout,"Reaction Coordinate 1 (%f) is outside the window (%f, %f) in box %d!!!\n",ordparm[k].x1, sim_dos[k].x1_begin, sim_dos[k].x1_end, k);
		  exit(1);
		}

		if(ordparam[k].x2 < sim_dos[k].x2_begin || ordparam[k].x2 > sim_dos[k].x2_end){
		  fprintf(stdout,"Reaction Coordinate 2 (%f) is outside the window (%f, %f) in box %d!!!\n",ordparm[k].x2, sim_dos[k].x2_begin, sim_dos[k].x2_end, k);
		  exit(1);
		}

		sim.T[k]	= sim_dos[k].T_begin;
		sim.kT[k]	= RG * sim.T[k]*.001;
	}


  /* ================================================================== */
  /* Begin production iterations.                                       */
  /* ================================================================== */
	for(unsigned long i=1; i<=sim.cyc_pr; i++){
		n_iter=i;
		for(int k=0; k<sim.NB; k++){
			#ifdef STATUS
			curr_status(k,6);
			#endif
			double random = ran2();
			/*------------------------------------------*/
			/* Perform a move according to a specified  */
			/* probabilities.                           */
			/*------------------------------------------*/
			if      (random<mc_pivot.PMS[k])                                     tdxedos_pivot(k);
			else if (random<mc_trans.PMS[k] + mc_pivot.PMS[k])                   tdxedos_trans(k);
			else if (random<mc_axial.PMS[k] + mc_trans.PMS[k] + mc_pivot.PMS[k]) tdxedos_scale(k);
			else if (random<mc_rand.PMS[k]  + mc_axial.PMS[k] + mc_trans.PMS[k] 
					   + mc_pivot.PMS[k])                                          tdxedos_linear(k);
			else if (random<mc_rand.PMS[k]  + mc_axial.PMS[k] + mc_trans.PMS[k] 
					   + mc_pivot.PMS[k] + mc_solv.PMS[k])                         tdxedos_solv(k);
			else                                                                 tdxedos_hmc(k); 
		} //k

		#ifdef WALL
		/*------------------------------------------*/
		/* Check to make sure no sites have entered */
		/* the wall.                                */
		/*------------------------------------------*/
		for(int k=0; k<sim.NB; k++){
			for(int j=0; j<box[k].boxns; j++){
				if(atom[k][j].z < wall[k].z[0]){
					fprintf(stdout,"Particle %d (z coord = %f) is beyond the wall (z coord = %f) on iteration %d\n",j,atom[k][j].z, wall[k].z[0], i);
					exit(7923);
				}
			}
		} //k
		#endif

		/* ================================================================== */
		/* Propose a swap move.                                               */
		/* ================================================================== */
		#ifndef MPI
		if((i%sim_hyb.cyc_swap)==0 && sim.NB>1){
			int n_success=0;
			struct{
				int box1;
				int box2;
			}swapped_boxes[NBOXES];

			/*------------------------------------------*/
			/*First pick a box at random according to   */
   	 		/*correct probabilities.                    */
    		/*------------------------------------------*/
    		int box1, box2;
    		double random_box1 = ran2();

    		for(int k=0; k<sim.NB; k++){
      			if(random_box1 < swapping.pick_box_prob[k]){
        			box1 = k;
        			break;
      			}
    		}

    		/*------------------------------------------*/
    		/*Now randomly choose one the the nabors of */
    		/*box1.                                     */
    		/*------------------------------------------*/
    		int random_box2 = ran_int(0,swapping.n_nabors[box1]);
			box2 = swapping.nabor[box1][random_box2];
			int success = swap_tdxe_dos(box1, box2);
    		if (success == 1) {
				swapped_boxes[n_success].box1 = box1;
				swapped_boxes[n_success].box2 = box2;
				n_success++;
			}
			if(n_success != 0){
				for(int j=0; j<n_success; j++){
					swap_ptr[swapped_boxes[j].box1] += sprintf(swap_ptr[swapped_boxes[j].box1],
						"%d iterations    flag=%d box%d &box%d swapped \n",
						i,sim_hyb.flag,swapped_boxes[j].box1,swapped_boxes[j].box2);
				}
				ioflush(0);
			}

		}
		#else //ifndef MPI

		/*------------------------------------------------------*/
		/* If MPI is defined then go to swap_mpi and later		*/
		/* modify the value of flag from 0 to 1 or 1 to 0		*/
		/* If flag = 0 the process sends the message if flag is */
		/* one the process receive a message from previous rank */
		/*------------------------------------------------------*/
/*	    if((i%sim_hyb.cyc_swap)==0 ){
		  swap2= fopen(name4,"a");
		  for(int k=0; k<sim.NB; k++){
			int success=swap_xe_mpi(i,mpi.flag,0);
			if (mpi.flag==0) mpi.flag =1; else mpi.flag=0;	
			if (success ==1) fprintf(swap2," At iteration %d	box with rank %d swapped \n",i,mpi.my_rank);				
		  }
		  fclose(swap2);
		}
*/
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
		/* Write the gromacs .trr file with fequency       */
		/* blockt.                                         */
		/* ----------------------------------------------- */
		#ifdef TRR	
		if((i%sim.blockt) == 0) {
		  for (int k =0; k<sim.NB; k++) write_trr(k,i,0);
		}
		#endif

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

		/* ----------------------------------------------- */
		/* If all the boxes are converged, print the       */
		/* final flathist files and exit.                  */
		/* ----------------------------------------------- */
		#ifndef MPI		
		status1 =1;
		for(int k=0; k<sim.NB; k++) {
			if (sim_dos[k].mod_f>STOP_F) status1=0;
		}
		if(status1 ==1){
			for(int k=0; k<sim.NB; k++) {
				FILE *fptr_hist;
				sprintf(name5,"./OUTPUT/DOS/flathist%d.txt",k);  
				fptr_hist= fopen(name5,"w");
				for (int i=0; i<sim_dos[k].x1_bins; i++){
					for(int j=0; j<sim_dos[k].x2_bins; j++)
						fprintf(fptr_hist, "%d	%lf	%lf	%lf	%d\n", i, dos_hist[k][i][j].x1_mid,
						dos_hist[k][i][j].x2_mid, dos_hist[k][i][j].g_of_l, dos_hist[k][i][j].h_of_l);
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

  
void flat_histogram_td(int k)
{
  double sum=0; double flat_crit; int flag1=1;

  for(int i=0; i<sim_dos[k].x1_bins; i++){
    for(int j=0; j<sim_dos[k].x2_bins; j++)
      sum = sum + dos_hist[k][i][j].h_of_l;
  }

  flat_crit = sum/(sim_dos[k].x1_bins*sim_dos[k].x2_bins)*sim_dos[k].flat;


  for(int i=0; i<sim_dos[k].x1_bins; i++){
    for(int j=0; j<sim_dos[k].x2_bins; j++)
	  if (dos_hist[k][i][j].h_of_l < flat_crit) flag1=0; 
  }
	
  if (flag1==1){
    sim_dos[k].mod_f = sim_dos[k].mod_f/2.0;
    #ifdef STATS
	if(sim_dos[k].data_count < 50) dos_stats(k);
    #endif
    if(sim_dos[k].mod_f >= STOP_F){
      for(int i=0; i<sim_dos[k].x1_bins; i++){
	    for(int j=0; j<sim_dos[k].x2_bins; j++)
          dos_hist[k][i][j].h_of_l=0;
	  }
    if(sim_dos[k].mod_f >= RESET_F)	init_block(k);
	}		 
  }
}
#endif//TDXEDOS

