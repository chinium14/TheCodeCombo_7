/* ======================================================================== */
/* repexch.cpp                                                              */
/* This suboutine performs parallel tempering using hybrid monte carlo      */
/* moves. After every certain bybrid MC move a swap moves is proposed and   */
/* accepted or rejected based on the standard criterion.                    */  
/* ======================================================================== */
  
#include "defines.h"
void init_all     (void);
void init_block   (int);
void outend      (unsigned long);
void xyz_config   (int,unsigned long);
void movie      (int,unsigned long);
void vblock      (int);
void pbc_all    (int);
void forces      (int);
void force_short  (int);
void force_long   (int);
void svalues    (int);
void kinet      (int);
void calcvalue    (int);
void drift      (int);
void nblist      (int);
void save_config  (unsigned long,unsigned long);
void save_veloc   (unsigned long,unsigned long);
void output      (unsigned long,int);
void ioflush    (int);
void integrate      (int, double);
void integrate_vverlet  (int, double);
void integrate_mts    (int, double);
void integrate_mts_nhc  (int, double);
void integrate_nhc    (double, int);
  
double ran2      (void);
double velscale   (int);
int     ran_int    (int, int);
#ifdef TLATE
void   translate  (int);
#endif
#ifdef SASA
double csasa (int);
#endif
#ifdef RMSD
void center(int, double , double ,int, double);
void rotmol (int, double, double, double);
void jacobi(double, double, double, int);
void q2mat(double, double);
void qtrfit(int, double,double,double,double,double,int);
int  quatfit (int, int, int);
#endif

#ifdef CONFIGT
void configtemp(int);
void hess (void);
#endif
#ifdef ROTN
void rotate_xyz(int,int);
#endif
#ifdef DLIST
void nl_check(int,int,int);
#endif
#ifdef TRR
void write_trr(int, unsigned long, int);
#endif
#ifdef MPI
int   swap_box_mpi(int);
#else
int   swap_box    (int);
#endif
#ifdef REM
void init_rem(void);
#endif
#ifdef TWHAM
void twham_init(void);
void twham_write(int,int);
#endif
#ifdef XWHAM
void xwham_write(int,int); 
#endif
#ifdef CHECKP
int WriteCheckpoint(const char *filename, unsigned long i, double vscal[NBOXES], int state);
int ReadCheckpoint(const char *filename, unsigned long *i, double vscal[NBOXES], int state);
#endif

void repexch()
{
    double vscal[NBOXES];

    sim_hyb.flag=0;

/* **********************************************************************/
/*  NVT EQUILIBRATION RUN STARTS                                        */
/* **********************************************************************/


  for(int k=0; k<sim.NB; k++) vscal[k]= 1.0;// initial scale for velocity

  /* ================================================================== */
  /*                                                                    */
  /* Begin equilibration iterations.                                    */
  /*                                                                    */
  /* ================================================================== */

unsigned long i=1;
	#ifdef CHECKP
        if(restart_equil == 1)
        {
                ReadCheckpoint(equiname[checkptr], &i, vscal, 0);          //Restart equilibration

                for(int k=0; k<sim.NB; k++)
                {
                        #ifdef NLIST
                        #ifndef SASA    /* Neighborlist will be checkpointed if SASA is defined */
                                nblist(k);
			#endif
                        #endif
                }
        }
        if(restart_equil == 2)
                i = sim.cyc_eq + 1;     //Skip equilibration
	#endif

  for( ; i<=sim.cyc_eq; i++) { // MC counter

	#ifdef CHECKP
	//Don't write checkpoints when catching TERM signals for parallel jobs
	if(i%sim.blockc==0)
	{
		//Make sure to write all checkpoints synchronously
		#ifdef MPI
		MPI_Barrier(MPI_COMM_WORLD);
		#endif
		ioflush(1);
		WriteCheckpoint(equiname[checkptr], i, vscal, 0);
		checkptr = 1 - checkptr;
	}
	#endif

    for(int k=0; k<sim.NB; k++) {
  /* ----------------------------------------------- */
  /* Update the neighborlist with frequency          */
  /* sim.nlist.                                      */
  /* ----------------------------------------------- */
      #ifdef NLIST
      #ifndef DLIST
      if(i%sim.nlist==0) nblist(k);
      #else
      nl_check(0,box[k].boxns,k);
      if(nl_flag[k] == 1){
        nblist(k);
        //printf("nblist called at iteration %d\n",i);
      }
      #endif
      #endif  
  /* ----------------------------------------------- */
  /* Itegration of the Equations of Motion.          */
  /* ----------------------------------------------- */
       /* ----------------------------------------------- */
       /* Leap Frog Algorithm                             */
       /* ----------------------------------------------- */
      if (sim.ID2 ==1 ) integrate(k,vscal[k]); 

       /* ----------------------------------------------- */
       /* Velocity Verlet Algorithm                       */
       /* ----------------------------------------------- */
      else if (sim.ID2 ==2 || sim.ID2 ==3)  integrate_vverlet(k, vscal[k]);

       /* ----------------------------------------------- */
       /* Nose Hoover Velocity Verlet Algorithm           */
       /* ----------------------------------------------- */
      else if (sim.ID2 ==4 || sim.ID2 ==5) 
      {
        integrate_nhc(sim.dt, k);
        integrate_vverlet(k, 1.0);
        integrate_nhc(sim.dt, k);
      }

      #ifdef WALL
      #ifndef SPHERE
        /*------------------------------------------*/
      /* check to make sure no sites have entered */
      /* the wall.                                */
      /*------------------------------------------*/
      for(int j=0; j<box[k].boxns; j++){
        if(atom[k][j].z < wall[k].z[0]){
          fprintf(stdout,"Particle %d (z coord = %f) is beyond the wall (z coord=%f) on iteration %lu\n",j, atom[k][j].z, wall[k].z[0], i);
	  exit(7923);
        }
      }
      #endif
      #endif
      /* ----------------------------------------------- */
      /* Calculate the properties of the system.         */
      /* ----------------------------------------------- */
      kinet(k);                  
      calcvalue(k);  
      
      /* ----------------------------------------------- */
      /* Correct for momentum drift with frequency       */
      /* sim.drift.                                      */
      /* ----------------------------------------------- */
      if((i%sim.drift) == 0) drift(k);

      /* ----------------------------------------------- */
      /* Increment tempererature by sim.dtemp.           */
      /* ----------------------------------------------- */
      sim.T[k] += sim.dtemp[k]; 
      
      /* ----------------------------------------------- */
      /* Determine the velocity scaling factor           */
      /* ----------------------------------------------- */
      if (sim.ID2!=4 && sim.ID2!=5 && sim.ID2!=8) vscal[k] = velscale(k); 
      
    }//end for k
    
    /* ----------------------------------------------- */
    /* Check to see if stop.txt is present. If it is   */
    /* set the counter to exit the next iteration.     */
    /* ----------------------------------------------- */
    if((i%sim.blockd) == 0) {              
      #ifdef WIN
      struct _stat buf;
      if(_stat("./INPUT/stop.txt",&buf) == 0) i=sim.cyc_eq+1;
      #else
      struct stat buf;
      if(stat("./INPUT/stop.txt",&buf) == 0) i=sim.cyc_eq+1;
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

    }//blockc ends

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

#ifndef MPI    
    /* ----------------------------------------------- */
    /* Propose a swap move and check acceptance.     */
    /* ----------------------------------------------- */    
    if((i%sim_hyb.cyc_swap)==0){
      int n_success=0;
      struct{
        int box1;
        int box2;
      }swapped_boxes[NBOXES];

      if(sim_hyb.flag==0){
        int k=0;
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
        sim_hyb.flag=1;
      }
      else {
        int k=1;
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
        sim_hyb.flag=0;
      }

      if(n_success != 0){
        for(int j=0; j<n_success; j++) swap_ptr[swapped_boxes[j].box1] += sprintf(swap_ptr[swapped_boxes[j].box2],"%lu\t%d\n",i,swapped_boxes[j].box1);
        ioflush(0);
      }
    }// end swap counter here
#else
                /*------------------------------------------------------*/
                /* If MPI is defined then go to swap_mpi and later              */
                /* modify the value of flag from 0 to 1 or 1 to 0               */
                /* If flag = 0 the process sends the message if flag is */
                /* one the process receive a message from previous rank */
                /*------------------------------------------------------*/
    if((i%sim_hyb.cyc_swap)==0 ){
      int n_success = 0;
      int success=swap_box_mpi(mpi.flag);
      if (success ==1){
        if (mpi.flag == 0) swap_ptr[0] += sprintf(swap_ptr[0],"%lu\t%d\n", i, mpi.my_rank-1); //only print if you are the "top" box (the one receiving the request
        n_success++;
        sim_hyb.swap_acc[0]++;
      }                      
      if (mpi.flag==0) mpi.flag =1; else mpi.flag =0;
      if(n_success > 0) ioflush(0);
    }
#endif //MPI  
  

    #ifdef TWHAM
    #ifndef SASAREX
    for(int k=0; k<sim.NB; k++){
      if(en[k].potens < pe_min[k]) pe_min[k] = en[k].potens;
      if(en[k].potens > pe_max[k]) pe_max[k] = en[k].potens;
    }
    #else
    for(int k=0; k<sim.NB; k++){
      if(en[k].esasa < pe_min[k]) pe_min[k] = en[k].esasa;
      if(en[k].esasa > pe_max[k]) pe_max[k] = en[k].esasa;
    }
    #endif
    #endif    

  }//end for i

i = 1;
#ifdef CHECKP
if(restart_equil != 2)          //Don't reinit production run if restarting production
#endif
{

  ioflush(1); 
  fprintf(stdout,"Equilibration completed!!!\n");

  /* ================================================================== */
  /* Apply pbc's, recalculate forces and properties of the system to    */
  /* prepare for production run.                                        */
  /* ================================================================== */
  for(int k=0; k<sim.NB; k++) {
    pbc_all(k);
    #ifndef GOLIK
    force_short(k); force_long(k); 
    #endif
    forces(k);   
    kinet(k);
    calcvalue(k);  
  }



/* **********************************************************************/
/*  PRODUCTION RUN STARTS                        */
/* **********************************************************************/
  #ifdef REM
  init_rem();
  #endif

  #ifdef TWHAM
  twham_init();
  #endif
}


#ifdef CHECKP
if(restart_prod)
{
        //Restart from production checkpoint file
        ReadCheckpoint(prodname[checkptr], &i, vscal, 1);
        for(int k=0; k<sim.NB; k++)
        {
                        #ifdef NLIST
                        #ifndef SASA    /* Neighborlist will be checkpointed if SASA is defined */
                                nblist(k);
			#endif
                        #endif
        }
}
#endif

  /* ================================================================== */
  /*                                                                    */
  /* Begin production iterations.                                       */
  /*                                                                    */
  /* ================================================================== */
  for( ; i<=sim.cyc_pr; i++) { // MC counter
	#ifdef CHECKP
	if(i%sim.blockc==0)
	{
		//Make sure to write all checkpoints sychronously (otherwise there are problems)
		#ifdef MPI
		MPI_Barrier(MPI_COMM_WORLD);
		#endif
		ioflush(1);
		WriteCheckpoint(prodname[checkptr], i, vscal, 1);
		checkptr = 1 - checkptr;
	}
	#endif

    for(int k=0; k<sim.NB; k++) {


    /* ----------------------------------------------- */
    /* Update the neighborlist with frequency          */
    /* sim.nlist.                                      */
    /* ----------------------------------------------- */
      #ifdef NLIST
      #ifndef DLIST
      if(i%sim.nlist==0) nblist(k);
      #else
      nl_check(0,box[k].boxns,k);
      if(nl_flag[k] == 1){
        nblist(k);
        //printf("nblist called at iteration %d\n",i);
      }
      #endif
      #endif  
      
    /* ----------------------------------------------- */
    /* Itegration of the Equations of Motion.          */
    /* ----------------------------------------------- */
    /* ----------------------------------------------- */
    /* Leap Frog Algorithm with Berendsen thermostat.  */
    /* ----------------------------------------------- */
      if (sim.ID2 ==1) integrate(k,vscal[k]);       
      
    /* ----------------------------------------------- */
    /* Velocity Verlet Algorithm with Berendsen        */
    /* thermostat.                                     */
    /* ----------------------------------------------- */
      else if (sim.ID2 ==2) integrate_vverlet(k, vscal[k]);  
    /* ----------------------------------------------- */
    /* Multiple Time Step (RESPA) Integration with     */
    /* Berendsen thermostat.                           */
    /* ----------------------------------------------- */
      else if (sim.ID2 ==3){                  
      /* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ */
        /* Because multiple time step algorithms use longer */
      /* time steps, adjust the number of iterations and  */
      /* update frequencies to maintain the desired run   */
      /* time.                                            */
      /* //////////////////////////////////////////////// */
        if(i==0 && k==0){
          sim.cyc_pr /=(long)sim.nsteps;
          sim.blockc /=(long)sim.nsteps;
          sim.blockd /=(long)sim.nsteps;
          sim.blockt /=(long)sim.nsteps;
          sim.drift  /=(long)sim.nsteps;
          sim.nlist  /=(int) sim.nsteps;
        }
        integrate_mts(k,vscal[k]);
      }
    /* ----------------------------------------------- */
    /* Nose Hoover Velocity Verlet Algorithm           */
    /* ----------------------------------------------- */
      else if (sim.ID2 ==4){
        integrate_nhc(sim.dt, k);
        integrate_vverlet(k, 1.0);
        integrate_nhc(sim.dt, k);
      }
    /* ----------------------------------------------- */
    /* Multiple Time Step Nose Hoover Integration      */
    /* (XO/XI-RESPA)                                   */
    /* ----------------------------------------------- */
      else if (sim.ID2 == 5){
      /* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ */
        /* Because multiple time step algorithms use longer */
      /* time steps, adjust the number of iterations and  */
      /* update frequencies to maintain the desired run   */
      /* time.                                            */
      /* //////////////////////////////////////////////// */
        if(i==0 && k==0){
          sim.cyc_pr /=(long)sim.nsteps;
          sim.blockc /=(long)sim.nsteps;
          sim.blockd /=(long)sim.nsteps;
          sim.blockt /=(long)sim.nsteps;
          sim.drift  /=(long)sim.nsteps;
          sim.nlist  /=(int) sim.nsteps;
        }
        integrate_mts_nhc(k,vscal[k]);
      }
      #ifdef WALL
      #ifndef SPHERE
      /*------------------------------------------*/
      /* check to make sure no sites have entered */
      /* the wall.                                */
      /*------------------------------------------*/
      for(int j=0; j<box[k].boxns; j++){
        if(atom[k][j].z < wall[k].z[0]){
          fprintf(stdout,"Particle %d (z coord = %f) is beyond the wall (z coord=%f) on iteration %lu\n",j, atom[k][j].z, wall[k].z[0], i);
	  exit(7923);
        }
      }
      #endif
      #endif

    /* ----------------------------------------------- */
    /* Calculate the properties of the system.         */
    /* ----------------------------------------------- */
      kinet(k);                          
      calcvalue(k);
      
    /* ----------------------------------------------- */
    /* If not an NVE ensemble, correct for drift,      */
    /* increment temperature, and determine the        */
    /* velocity scale.                                 */
    /* ----------------------------------------------- */
      if(sim.ID!=3 && sim.ID2 !=8){                    
        if((i%sim.drift) == 0) drift(k);        
        sim.T[k] += sim.dtemp[k];                
        if (sim.ID2 !=4 && sim.ID2!=5 && sim.ID2!=8) vscal[k] = velscale(k); 
      }

    /* ----------------------------------------------- */
    /* Accumulate the calculated values for            */
    /* block averages and error bars.                  */
    /* ----------------------------------------------- */
      svalues(k);
      #ifdef CONFIGT
      configtemp(k);
      #endif
    }// k boxes ends
    
    /* ----------------------------------------------- */
    /* Calculate the block averages accumulated in     */
    /* svalues() and reset the accumulators.           */
    /* ----------------------------------------------- */
    if((i%sim.blockd) == 0) {                
      vblock(i);
      for(int k=0; k<sim.NB; k++) init_block(k);
      outend(i);
      ioflush(0);

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
      for(int k=0; k<sim.NB; k++)  quatfit(k,rmsd_modes.file,rmsd_modes.atoms);
      #endif

      output(i,1);
      ioflush(0);
      #ifdef TWHAM
      for(int k=0; k<sim.NB; k++) twham_write(k,0);
      #endif
      #ifdef XWHAM
      for(int k=0; k<sim.NB; k++) xwham_write(k,0);
      #endif

    } //blockc ends  
    
    if((i%sim.blockt) == 0) {
      #ifdef TLATE
      for (int k =0; k<sim.NB; k++) translate(k);
      #endif

      #ifdef ROTN
      for (int k =0; k<sim.NB; k++) rotate_xyz(k,_AXIS);  
      #endif

      #ifdef TRR
      for (int k =0; k<sim.NB; k++) write_trr(k,i,1);
      #endif

      save_config(i,sim.cyc_pr);                
      save_veloc(i,sim.cyc_pr);                
      xyz_config(1,i);
      movie(1,i);
    }
#ifndef MPI    
    /* ----------------------------------------------- */
    /* Propose a swap move and check acceptance.     */
    /* ----------------------------------------------- */    
    if((i%sim_hyb.cyc_swap)==0){
      int n_success=0;
      struct{
        int box1;
        int box2;
      }swapped_boxes[NBOXES];

      if(sim_hyb.flag==0){
        int k=0;
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
        sim_hyb.flag=1;
      }
      else {
      int k=1;
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
        sim_hyb.flag=0;
      }

      if(n_success != 0){
        for(int j=0; j<n_success; j++) swap_ptr[swapped_boxes[j].box1] += sprintf(swap_ptr[swapped_boxes[j].box2],"%lu\t%d\n",i,swapped_boxes[j].box1);
        ioflush(0);
      }
    }// end swap counter here
#else
                /*------------------------------------------------------*/
                /* If MPI is defined then go to swap_mpi and later              */
                /* modify the value of flag from 0 to 1 or 1 to 0               */
                /* If flag = 0 the process sends the message if flag is */
                /* one the process receive a message from previous rank */
                /*------------------------------------------------------*/
	if((i%sim_hyb.cyc_swap)==0 )
	{
		int n_success = 0;
		int success=swap_box_mpi(mpi.flag);
		if (success ==1)
		{
			if (mpi.flag == 0)
				swap_ptr[0] += sprintf(swap_ptr[0],"%lu\t%d\n", i, mpi.my_rank-1); //only print if you are the "top" box (the one receiving the request
			n_success++;
			sim_hyb.swap_acc[0]++;
		}                      
		if (mpi.flag==0)
			mpi.flag = 1;
		else mpi.flag =0;
		if(n_success > 0)
			ioflush(0);
	}
#endif //MPI    

  }// end production run  MC counter i                            
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
