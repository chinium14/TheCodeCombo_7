/* ======================================================================== */
/* nvt_mc.cpp                                                               */
/*                                                                          */
/*    This subroutine is the main loop for NVT MC simulations.              */
/* Written by Thomas A. Knotts on 31 Aug 05.                                */
/* ======================================================================== */
#include "defines.h"
void    init_block          (int);
void    save_config         (unsigned long,unsigned long);
void    save_veloc          (unsigned long,unsigned long);
void    xyz_config          (int,unsigned long);
void    movie               (int,unsigned long);
void    ofile               (void);
void    output              (unsigned long,int);
void    outend              (unsigned long);
void    ioflush             (int);
void    svalues             (int);
void    vblock              (int);
void    pbc_all             (int);
void    forces              (int);
void    integrate_vverlet   (int, double);
void    integrate_nhc       (double, int);
void    kinet               (int);
void    calcvalue           (int);
void    nblist              (int);
void    configtemp          (int);
void    mcmove_pivot        (int);
void    mcmove_trans        (int);
void    mcmove_displ        (int);
void    mcmove_hmc          (int);
int     swap_box            (int);
double  ran2                (void);
double  velscale            (int);

#ifdef RMSD
int     quatfit             (int, int, int);
#endif

#ifdef TLATE
void    translate           (int);
#endif

#ifdef ROTN
void    rotate_xyz          (int,int);
#endif

#ifdef DLIST
void    nl_check            (int, int, int);
#endif

#ifdef TRR
void    write_trr           (int, unsigned long, int);
#endif

#ifdef TWHAM
void twham_init             (void);
void twham_write            (int, int);
#endif

#ifdef XWHAM
void xwham_write            (int, int);
#endif

#ifdef CHECKP
int WriteCheckpoint(const char *filename, unsigned long i, double vscal[NBOXES]);
int ReadCheckpoint(const char *filename, unsigned long *i, double vscal[NBOXES]);
#endif


void nvt_mc()
{
  /* ========================================== */
  /* Begin equilibration iterations.            */
  /* ========================================== */


	unsigned long i=1;

        //Determine if we should restart from a checkpoint.
	#ifdef CHECKP
	double vscal[NBOXES];	//Not used, but necessary for checkpointing

        if(restart_equil == 1)
        {
                ReadCheckpoint(equiname[checkptr], &i, vscal);          //Restart equilibration
		//We don't rebuild the neighborlist for monte carlo.
        }
        if(restart_equil == 2)
                i = sim.cyc_eq + 1;     //Skip equilibration
	#endif


  for( ; i<=sim.cyc_eq; i++){
 
    n_iter=i;
	#ifdef CHECKP
	if(terminating)
	{
		//Call ioflush here???
		ioflush(1);
		WriteCheckpoint(equiname[checkptr], i, vscal);
		printf("Terminating equilibration after iteration %lu\n", i-1);
		exit(0);
	}
	if(i%sim.blockc==0)
	{
		ioflush(1);
		WriteCheckpoint(equiname[checkptr], i, vscal);
		checkptr = 1 - checkptr;
	}
	#endif

    for(int k=0; k<sim.NB; k++){

      double random = ran2();

      /* ------------------------------------------ */
      /* Randomly pick and perform a move.          */
      /* ------------------------------------------ */			
      if      (random<mc_trans.PMS[k])                                mcmove_trans(k);
      else if (random<mc_pivot.PMS[k]+mc_trans.PMS[k])	              mcmove_pivot(k);
      else if (random<mc_rand.PMS[k]+mc_pivot.PMS[k]+mc_trans.PMS[k]) mcmove_displ(k);
      else                                                            mcmove_hmc(k);		
      calcvalue(k);
 			#ifdef CONFIGT
			configtemp(k);
			#endif

    }

    /* ------------------------------------------ */
    /* Attempt a swap move.                       */
    /* ------------------------------------------ */
    #ifndef MPI
    if(((i+1)%sim_hyb.cyc_swap)==0){
      int k;
      int n_success=0;
      struct{
        int box1;
        int box2;
      }swapped_boxes[NBOXES];

      if(sim_hyb.flag==0){
        k=0;
        sim_hyb.flag = 1;
      }else{
        k=1;
        sim_hyb.flag = 0;
      }


      while ((k+1)<sim.NB){
        int success =swap_box(k);
        if (success ==1) {
          swapped_boxes[n_success].box1 = k;
          swapped_boxes[n_success].box2 = k+1;
          n_success++;
          sim_hyb.swap_acc[k]++;
          sim_hyb.swap_acc[k+1]++;
        }
        k=k+2;
      }


      if(n_success != 0){
        for(int j=0; j<n_success; j++){
          swap_ptr[swapped_boxes[j].box1] += sprintf(swap_ptr[swapped_boxes[j].box1],
            "%lu iterations    flag=%d box%d &box%d swapped \n",
            i,sim_hyb.flag,swapped_boxes[j].box1,swapped_boxes[j].box2);
        }
        ioflush(0);
      }


    }// end swap counter here
    #else //ifndef MPI

    /* Need to code MPI swapping. */

    #endif //ifndef MPI


    if((i%sim.blockd) == 0){								
    /* ----------------------------------------------- */
    /* Check to see if stop.txt is present. If it is   */
    /* set the counter to exit the next iteration.     */
    /* ----------------------------------------------- */
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
    if((i%sim.blockc) == 0) {
      #ifdef RMSD
      for(int k=0; k<sim.NB; k++) quatfit(k,rmsd_modes.file,rmsd_modes.atoms);
      #endif
      output(i,0);
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

      save_config(i,sim.cyc_pr);								
      save_veloc(i,sim.cyc_pr);								
      xyz_config(1,i);
      movie(1,i);
    }

    #ifdef TWHAM
    for(int k=0; k<sim.NB; k++){
      if(en[k].potens < pe_min[k]) pe_min[k] = en[k].potens;
      if(en[k].potens > pe_max[k]) pe_max[k] = en[k].potens;
    }
    #endif

  }

i = 1;
#ifdef CHECKP
if(restart_equil != 2)          //Don't reinit production run if restarting production
#endif
{

  fprintf(stdout,"Equilibration completed!!!\n");

  /* ========================================== */
  /* Begin production iterations.               */
  /* ========================================== */
  #ifdef TWHAM
  twham_init();
  #endif
}
#ifdef CHECKP
if(restart_prod)
{
        //Restart from production checkpoint file
        ReadCheckpoint(prodname[checkptr], &i, vscal);
	//We don't rebuild the neighborlist for monte carlo.
}
#endif


  for( ; i<=sim.cyc_pr; i++){

    n_iter=i;
	#ifdef CHECKP
	if(terminating)
	{
		//Call ioflush here???
		ioflush(1);
		WriteCheckpoint(prodname[checkptr], i, vscal);
		printf("Terminating production after iteration %lu\n", i-1);
		exit(0);
	}
	if(i%sim.blockc==0)
	{
		ioflush(1);
		WriteCheckpoint(prodname[checkptr], i, vscal);
		checkptr = 1 - checkptr;
	}
	#endif

    for(int k=0; k<sim.NB; k++){

      double random = ran2();

      /* ------------------------------------------ */
      /* Randomly pick and perform a move.          */
      /* ------------------------------------------ */			
      if      (random<mc_trans.PMS[k])                                mcmove_trans(k);
      else if (random<mc_pivot.PMS[k]+mc_trans.PMS[k])                mcmove_pivot(k);
      else if (random<mc_rand.PMS[k]+mc_pivot.PMS[k]+mc_trans.PMS[k]) mcmove_displ(k);
      else                                                            mcmove_hmc(k);		
      calcvalue(k);
			#ifdef CONFIGT
			configtemp(k);
			#endif
      svalues(k);
    }

    /* ------------------------------------------ */
    /* Attempt a swap move.                       */
    /* ------------------------------------------ */
    #ifndef MPI
    if(((i+1)%sim_hyb.cyc_swap)==0){
      int k;
      int n_success=0;
      struct{
        int box1;
        int box2;
      }swapped_boxes[NBOXES];

      if(sim_hyb.flag==0){
        k=0;
        sim_hyb.flag = 1;
      }else{
        k=1;
        sim_hyb.flag = 0;
      }


      while ((k+1)<sim.NB){
        int success =swap_box(k);
        if (success ==1) {
          swapped_boxes[n_success].box1 = k;
          swapped_boxes[n_success].box2 = k+1;
          n_success++;
          sim_hyb.swap_acc[k]++;
          sim_hyb.swap_acc[k+1]++;
        }
        k=k+2;
      }


      if(n_success != 0){
        for(int j=0; j<n_success; j++){
          swap_ptr[swapped_boxes[j].box1] += sprintf(swap_ptr[swapped_boxes[j].box1],
            "%lu iterations    flag=%d box%d &box%d swapped \n",
            i,sim_hyb.flag,swapped_boxes[j].box1,swapped_boxes[j].box2);
        }
        ioflush(0);
      }


    }// end swap counter here
    #else //ifndef MPI

    /* Need to code MPI swapping. */

    #endif //ifndef MPI


    if((i%sim.blockd) == 0){								
    /* ----------------------------------------------- */
    /* Check to see if stop.txt is present. If it is   */
    /* set the counter to exit the next iteration.     */
    /* ----------------------------------------------- */
      #ifdef WIN
      struct _stat buf;
      if(_stat("./INPUT/stop.txt",&buf) == 0) i=sim.cyc_pr+1;
      #else
      struct stat buf;
      if(stat("./INPUT/stop.txt",&buf) == 0) i=sim.cyc_pr+1;
      #endif

			vblock(i);
			for(int k=0; k<sim.NB; k++) init_block(k);
			outend(i);
			ioflush(0);
    }

    /* ----------------------------------------------- */
    /* Write the calculated data, configurations, and  */
    /* velocities to file with frequency sim.blockc.   */
    /* ----------------------------------------------- */
    if((i%sim.blockc) == 0) {
      #ifdef RMSD
      for(int k=0; k<sim.NB; k++) quatfit(k,rmsd_modes.file,rmsd_modes.atoms);
      #endif
      output(i,1);
      ioflush(0);
      #ifdef TWHAM
      for(int k=0; k<sim.NB; k++) twham_write(k,0);
      #endif
      #ifdef XWHAM
      for(int k=0; k<sim.NB; k++) xwham_write(k,0);
      #endif
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

      save_config(i,sim.cyc_pr);								
      save_veloc(i,sim.cyc_pr);								
      xyz_config(1,i);
      movie(1,i);
    }
  }


    /* ----------------------------------------------- */
    /* Call output on the last iteration everytime.    */
    /* ----------------------------------------------- */
    output(sim.cyc_pr,1);
    ioflush(0);
    #ifdef TWHAM
    for(int k=0; k<sim.NB; k++) twham_write(k,1);
    #endif
    #ifdef XWHAM
    for(int k=0; k<sim.NB; k++) xwham_write(k,1);
    #endif
 
}
  
