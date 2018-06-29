/* ======================================================================== */
/* nvt_md.cpp                                                               */
/*                                                                          */
/*		This subroutine is the main loop for simulaitons using an NVT       */
/* simulation.                                                              */
/* Revised by Thomas A. Knotts on 9 Sep 04 for I/O buffering and clean up.  */
/* ======================================================================== */

#include "defines.h"
  /* ================================================================== */
  /* Function prototypes                                                */
  /* ================================================================== */
void init_all     (void);
void init_block   (int);
void outend		  (unsigned long);
void xyz_config   (int,unsigned long);
void movie		  (int,unsigned long);
void vblock		  (int);
void pbc_all	  (int);
void forces		  (int);
void force_short  (int);
void force_long   (int);
void svalues	  (int);
void kinet		  (int);
void calcvalue	  (int);
void drift		  (int);
void nblist		  (int);

void save_config  (unsigned long,unsigned long);
void save_veloc   (unsigned long,unsigned long);
void output		  (unsigned long,int);
void ioflush	  (int);

void integrate			(int, double);
void integrate_vverlet	(int, double);
void integrate_mts		(int, double);
void integrate_mts_nhc	(int, double);
void integrate_nhc		(double, int);
		
double ran2		  (void);
double velscale   (int);
int	   ran_int	  (int, int);

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

#ifdef SMD
void isokin(int);
void integrate_smd(int,double);
#endif
#ifdef NSMD
void integrate_smd(int,double);
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

#ifdef TWHAM
void twham_init             (void);
void twham_write            (int,int);
#endif

#ifdef XWHAM
void xwham_write(int,int);
#endif

#ifdef ZHOU
void stressComputeInvariantsNow(int);
void stressSumInvariants(int);
void stressSumTensors(int);
void stressOutputInvariantsNow(int, unsigned int, int);
void stressOutputEigenvaluesNow(int, unsigned int, int);
void stressOutputTensorsNow(int, unsigned int, int);
void stressFinishBlock(int);
void stressResetBlock(int);
void stressResetTensors(int);
void stressResetAll(int);
void stressOutputStats(int);
        #ifdef TRR
        void stressInitMovie(int);
        void stressSumMovieStress(void);
        void stressOutputFrame(int, unsigned int, int);
        #endif
#endif

#ifdef CHECKP
int WriteCheckpoint(const char *filename, unsigned long i, double vscal[NBOXES]);
int ReadCheckpoint(const char *filename, unsigned long *i, double vscal[NBOXES]);
#endif

void nvt_md(){

	double vscal[NBOXES]; 

/* **********************************************************************/
/*	NVT EQUILIBRATION RUN STARTS										*/
/* **********************************************************************/

  #ifdef ZHOU
        for(int k=0;k<sim.NB;k++) stressResetAll(k);    //Initialize everything
  #ifdef TRR
        for(int k=0;k<sim.NB;k++) stressInitMovie(k);
  #endif
  #endif

	for(int k=0; k<sim.NB; k++) vscal[k]= 1.0;// initial scale for velocity
  /* ================================================================== */
#ifndef NSMD
#ifndef SMD
  /* ================================================================== */
  /*                                                                    */
  /* Begin equilibration iterations.                                    */
  /*                                                                    */
  /* ================================================================== */

	unsigned long i=1;

	//Determine if we should restart from a checkpoint.
	#ifdef CHECKP	
	if(restart_equil == 1)
	{
		ReadCheckpoint(equiname[checkptr], &i, vscal);		//Restart equilibration
	
		for(int k=0; k<sim.NB; k++)
		{	
			#ifdef NLIST
			#ifndef SASA	// Neighborlist will be checkpointed if SASA is defined
				nblist(k);
			#endif
			#endif	
		}	
	}
	if(restart_equil == 2)
		i = sim.cyc_eq + 1;	//Skip equilibration
	#endif

	for( ; i<=sim.cyc_eq; i++) {
		//printf("Eq Iteration %lu\n", i);
		n_iter=i;
		#ifdef CHECKP
		if(terminating)
		{
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
		for(int k=0; k<sim.NB; k++) {


//This is used to check the hessians
/*
hess();
forces(0);
double hessx =0.0;
double hessy =0.0;
double hessz =0.0;
for(int j=0; j<box[k].boxns; j++){
	hessx += hessian[3*j+0][3*j+0];
	hessy += hessian[3*j+1][3*j+1];
	hessz += hessian[3*j+2][3*j+2];
}
*/


  /* ----------------------------------------------- */
  /* Update the neighborlist with frequency          */
  /* sim.nlist.                                      */
  /* ----------------------------------------------- */
			#ifdef NLIST
			#ifndef DLIST
				if(i%sim.nlist==0) nblist(k);
			#endif
			#ifdef DLIST
				nl_check(0,box[k].boxns,k);
				if(nl_flag[k] == 1){
					nblist(k);
					//printf("nblist called at iteration %lu\n",i);
				}
			#endif
			#endif
	
                        #ifdef ZHOU
                            stressResetTensors(k);      //Zero-out tensors to prepare for accumulation
                        #endif
  /* ----------------------------------------------- */
  /* Integration of the Equations of Motion.          */
  /* ----------------------------------------------- */
			/* ----------------------------------------------- */
			/* Leap Frog Algorithm                             */
			/* ----------------------------------------------- */
			if (sim.ID2 ==1 ) integrate(k,vscal[k]); 

       /* ----------------------------------------------- */
       /* Velocity Verlet Algorithm                       */
	   /* Note: For the equilibration iterations, no      */
	   /* multiple time step algorthm is used, so the     */
	   /* same subroutine is used for both Berendsen and  */
	   /* mts Berendsen thermostats.                      */
       /* ----------------------------------------------- */
			else if (sim.ID2 ==2 || sim.ID2 ==3)	integrate_vverlet(k, vscal[k]);

       /* ----------------------------------------------- */
       /* Nose Hoover Velocity Verlet Algorithm           */
	   /* Note: For the equilibration iterations, no      */
	   /* multiple time step algorthm is used, so the     */
	   /* same subroutine is used for both Nose Hoover    */
	   /* and mts Nose Hoover thermostats.                */
	   /* Note2: For Nose Hoover algorithms, velocity     */
	   /* scaling is done in integrate_nhc() so           */
	   /* integrate_vverlet() is always called with       */
	   /* vscal[k] = 1.0.                                 */
       /* ----------------------------------------------- */
			else if (sim.ID2 ==4 || sim.ID2 ==5) 
			{
				integrate_nhc(sim.dt, k);
				integrate_vverlet(k, 1.0);
				integrate_nhc(sim.dt, k);
			}

			#ifdef SMD
			else if (sim.ID2 ==8) 
			{
				integrate_nhc(sim.dt, k);
				integrate_smd(k, 1.0);
				integrate_nhc(sim.dt, k);
			}
			#endif

			#ifdef NSMD
			else if (sim.ID2 ==8) 
			{
				integrate_nhc(sim.dt, k);
				integrate_smd(k, 1.0);
				integrate_nhc(sim.dt, k);
			}
			#endif

			#ifdef WALL
                        #ifndef SPHERE
		    /*------------------------------------------*/
			/* check to make sure no sites have entered */
			/* the wall.                                */
			/*------------------------------------------*/
			for(int j=0; j<box[k].boxns; j++){
				if(atnopbc[k][j].z < wall[k].z[0]){
					fprintf(stdout,"Particle %d (z coord=%f) is beyond the wall (z coord=%f) on iteration %lu\n",j,atnopbc[k][j].z,wall[k].z[0],i);
					exit(7923);
				}
			}
                        #endif
			#endif

			/* ----------------------------------------------- */
			/* Calculate the properties of the system.         */
			/* ----------------------------------------------- */
			kinet(k);    	            
			//svalues(k);//uncomment to accumulate properties in equil
			calcvalue(k);  
 			#ifdef CONFIGT
			configtemp(k);
			#endif
			
			/* ----------------------------------------------- */
			/* Correct for momentum drift with frequency       */
			/* sim.drift.                                      */
			/* ----------------------------------------------- */
			#ifndef SMD
			#ifndef NSMD
                        if((i%sim.drift) == 0) drift(k);

			#endif
			#endif

			/* ----------------------------------------------- */
			/* Increment tempererature by sim.dtemp.           */
			/* ----------------------------------------------- */
			sim.T[k] += sim.dtemp[k]; 
			
			/* ----------------------------------------------- */
			/* Determine the velocity scaling factor           */
			/* ----------------------------------------------- */
			if (sim.ID2!=4 && sim.ID2!=5 && sim.ID2!=8) vscal[k] = velscale(k); 
                        #ifdef ZHOU

                                #ifdef TRR
                                stressComputeInvariantsNow(k);  //Get I from S
                                stressSumMovieStress(); //Only for box 0, for now
                                #endif

                        #endif			
		}//end for k

		/* ----------------------------------------------- */
		/* Check to see if stop.txt is present. If it is   */
		/* set the counter to exit the next iteration.     */
		/* ----------------------------------------------- */
		if((i%sim.blockd) == 0) {							
	                       // vblock(i);
                        // for(int k=0; k<sim.NB; k++) init_block(k);
                        //outend(i);
                        //ioflush(0);  //uncomment to accumulate statistics in equil

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
                                #ifdef ZHOU
                                for(int k=0;k<sim.NB;k++)
                                {
                                        //stressComputeInvariantsNow(k);        //Get I from S
                                        stressOutputFrame(k,i,0);
                                }
                                #endif

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
			
	}//end for i
#endif//SMD
#endif //NSMD

i = 1;
#ifdef CHECKP
if(restart_equil != 2)		//Don't reinit production run if restarting production
#endif
{
	fprintf(stdout,"Equilibration completed!!!\n");

  /* ================================================================== */
  /* Apply pbc's, recalculate forces and properties of the system to    */
  /* prepare for production run.                                        */
  /* ================================================================== */
	for(int k=0; k<sim.NB; k++) {
		pbc_all(k);
#ifndef GOLIK
		force_short(k);force_long(k); 
#endif
		forces(k);    
		kinet(k);
		calcvalue(k);  
	}
/* **********************************************************************/
/*	PRODUCTION RUN STARTS												*/
/* **********************************************************************/
  #ifdef TWHAM
  twham_init();
  #endif

  #ifdef ZHOU
        for(int k=0;k<sim.NB;k++) stressResetAll(k);    //Initialize everything
  #endif
}

#ifdef CHECKP
if(restart_prod)
{
	//Restart from production checkpoint file
	ReadCheckpoint(prodname[checkptr], &i, vscal);
	for(int k=0; k<sim.NB; k++)
	{
			#ifdef NLIST
			#ifndef SASA /* Neighborlist will be checkpointed if SASA is defined */
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
	for( ; i<=sim.cyc_pr; i++) {
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
		for(int k=0; k<sim.NB; k++) {
		/* ----------------------------------------------- */
		/* Update the neighborlist with frequency          */
		/* sim.nlist.                                      */
		/* ----------------------------------------------- */
			#ifdef NLIST
			#ifndef DLIST
				if(i%sim.nlist==0) nblist(k);
			#endif
			#ifdef DLIST
				nl_check(0,box[k].boxns,k);
				if(nl_flag[k] == 1){
					nblist(k);
//					printf("nblist called at iteration %lu\n",i);
				}
			#endif
			#endif	

                        #ifdef ZHOU
                            stressResetTensors(k);      //Zero-out tensors to prepare for accumulation
                        #endif

			
		/* ----------------------------------------------- */
		/* Integration of the Equations of Motion.          */
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
#ifdef SMD
			else if (sim.ID2 ==8) 
			{
				integrate_nhc(sim.dt, k);
				integrate_smd(k, 1.0);
//				isokin(k);
				integrate_nhc(sim.dt, k);
			}
#endif			
#ifdef NSMD
			else if (sim.ID2 ==8) 
			{
				integrate_nhc(sim.dt, k);
				integrate_smd(k, 1.0);
				integrate_nhc(sim.dt, k);
			}
#endif
#ifdef WALL
#ifndef SPHERE
		    /*------------------------------------------*/
			/* check to make sure no sites have entered */
			/* the wall.                                */
			/*------------------------------------------*/
		for(int j=0; j<box[k].boxns; j++){
				if(atnopbc[k][j].z < wall[k].z[0]){
                                        fprintf(stdout,"Particle %d (z coord=%f) is beyond the wall (z coord=%f) on iteration %lu\n",j,atnopbc[k][j].z,wall[k].z[0],i);
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
                          //sim.T[k] += sim.dtemp[k];
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


                        #ifdef ZHOU

                                stressStepCount += 1;
                                stressTotalCount += 1;
                                stressComputeInvariantsNow(k);  //Get I from S
                                stressSumInvariants(k);         //Add S to S_sum and I to I_sum
                                stressSumTensors(k);
                                #ifdef TRR
                                stressSumMovieStress(); //Only for box 0, for now
                                #endif

                                if ((i%sim.blockc)==0)
                                {
                                        //Output instantaneous values (for error checking, mostly)
                                        stressOutputInvariantsNow(k, i, 1);
                                        stressOutputTensorsNow(k, i, 1);
                                        stressOutputEigenvaluesNow(k, i, 1);
                                        stressOutputCount += 1;

                                }

                                if ((i%sim.blockd)==0)
                                {
                                        //Do stuff for the statistics stuff
                                        stressFinishBlock(k);   //Get I_bar and I_bar^2 from I_sum, increment stressBlockCount
                                        stressResetBlock(k);    //Reset I_sum and stressStepCount
                                }

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
			for(int k=0; k<sim.NB; k++)	quatfit(k,rmsd_modes.file,rmsd_modes.atoms);
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
		
    }// end production run
  #ifdef ZHOU
       //Finish statistics and output.
       for(int k=0; k<sim.NB; k++)   stressOutputStats(k);
  #endif
	
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
