
#ifdef PRESSURE
/* ======================================================================== */
/* npt_md.cpp                                                               */
/*                                                                          */
/*		This subroutine is the main loop for simulaitons using an NPT       */
/* simulation.																*/
/* Revised by Thomas A. Knotts on 9 Sep 04 for I/O buffering and clean up.  */
/* ======================================================================== */

#include "defines.h"

void integrate_nhc(double, int);
void integrate_vverlet(int, double);
void integrate_mts_nhc(int, double);
void integrate_mts(int, double);
void boxscale(int);
double velscale (int);

void integrate_mts_npt_full(int);
void integrate_npt_full(int);
void pbc_npt_full(int);
void boxinv(int);
void force_short(int); 
void force_long(int); 
void forces(int);   
void kinet(int);
void calcvalue(int);
void nblist(int);
void output(unsigned long, int);
void ioflush (int); 
void outend(unsigned long);                 
void save_config(unsigned long,unsigned long);   
void save_veloc(unsigned long,unsigned long);   
void xyz_config(int,unsigned long);
void movie(int,unsigned long);
void quatfit(int, int, int);
void svalues(int);
void vblock(int);
void init_block(int);
void pbc_all(int);
#ifdef TLATE
void translate	(int);
#endif
#ifdef CONFIGT
void configtemp(int);
void hess (void);
#endif

#ifdef TRR
void write_trr(int,unsigned long,int);
#endif

#ifdef DLIST
void nl_check(int, int, int);
#endif
 
#ifdef CHECKP 
int WriteCheckpoint(const char *filename, unsigned long i, double vscal[NBOXES]);
int ReadCheckpoint(const char *filename, unsigned long *i, double vscal[NBOXES]);
#endif

void npt_md(void)
{
	double vscal[NBOXES];						// For Berendsen Thermostat
	for(int k=0; k<sim.NB; k++) vscal[k]= 1.0;	// initial scale for velocity

  /* ================================================================== */
  /* Calculate the the inverse box matrix.								*/
  /* ================================================================== */
	#ifdef PR_NPT
	for(int k=0; k<sim.NB; k++) boxinv(k);
	#endif

  /* ================================================================== */
  /*                                                                    */
  /* Begin equilibration iterations.                                    */
  /*                                                                    */
  /* ================================================================== */

	unsigned long i=1;

	#ifdef CHECKP
        //Determine if we should restart from a checkpoint.

        if(restart_equil == 1)
        {
                ReadCheckpoint(equiname[checkptr], &i, vscal);          //Restart equilibration

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

	for(; i<=sim.cyc_eq; i++) {

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
				if(nl_flag[k] == 1) nblist(k);
			#endif
			#endif	

			/* ----------------------------------------------- */
			/* Itegration of the Equations of Motion.          */
			/* ----------------------------------------------- */
			if(sim.ID2 == 6 || sim.ID2 == 7){
				integrate_nhc(sim.dt, k);
				integrate_vverlet(k, 1.0);
				integrate_nhc(sim.dt, k);

			}
			/* ----------------------------------------------- */
			/* If berendsen barostat is used; with or w/o MTS  */
			/* Remember no MTS in equilibration				   */
			/* ----------------------------------------------- */
			
			else integrate_vverlet(k, vscal[k]);
			
			
			/* ----------------------------------------------- */
			/* Calculate the properties of the system.         */
			/* ----------------------------------------------- */
			kinet(k);    	            
			calcvalue(k); 
			#ifndef PR_NPT
			boxscale(k);
			#endif

			/* ----------------------------------------------- */
			/* Increment tempererature by sim.dtemp.           */
			/* ----------------------------------------------- */
			sim.T[k] += sim.dtemp[k]; 

			/* ----------------------------------------------- */
			/* Determine the velocity scaling factor           */
			/* ----------------------------------------------- */
			if (sim.ID2!=6 && sim.ID2!=7) vscal[k] = velscale(k); 
						
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

	}//end for i
#ifdef CHECKP
if(restart_equil != 2)
#endif
{
	fprintf(stdout,"Equilibration completed!!!\n");

  /* ================================================================== */
  /* Apply pbc's, recalculate forces and properties of the system to    */
  /* prepare for production run.                                        */
  /* ================================================================== */
	for(int k=0; k<sim.NB; k++) {
		#ifdef PR_NPT
		pbc_npt_full(k);
		#else
		pbc_all(k);
		#endif
		force_short(k); force_long(k); 
		forces(k);   
		kinet(k);
		calcvalue(k);  
	}

	if(sim.ID2 == 7 || sim.ID2 == 10){					// Update frequencies for MTS
		sim.cyc_pr /=(long)sim.nsteps;
		sim.blockc /=(long)sim.nsteps;
		sim.blockd /=(long)sim.nsteps;
		sim.drift  /=(long)sim.nsteps;
		sim.nlist  /=(int) sim.nsteps;
	}
}

i = 1;

#ifdef CHECKP
if(restart_prod)
{
        //Restart from production checkpoint file
        ReadCheckpoint(prodname[checkptr], &i, vscal);
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


/* **********************************************************************/
/*	PRODUCTION RUN STARTS												*/
/* **********************************************************************/


  /* ================================================================== */
  /*                                                                    */
  /* Begin production iterations.                                       */
  /*                                                                    */
  /* ================================================================== */
	for( ; i<=sim.cyc_pr; i++) {
	
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
				if(nl_flag[k] == 1) nblist(k);
			#endif
			#endif	
			
			/* ----------------------------------------------- */
			/* Itegration of the Equations of Motion.          */
			/* ----------------------------------------------- */
		   /* *********************************************** */
		   /* Single Time Step Parinello Rahman NPT           */
		   /* *********************************************** */
			#ifdef PR_NPT
			if(sim.ID2 == 6) integrate_npt_full(k); 			
			
		   /* *********************************************** */
		   /* Mulitiple Time Step Parinello Rahman NPT        */
		   /* *********************************************** */
			else if (sim.ID2 == 7) integrate_mts_npt_full(k);

			#else

			/* ----------------------------------------------- */
			/* If berendsen barostat is used; with or w/o MTS  */
			/* ----------------------------------------------- */
			
			if (sim.ID2 == 1 || sim.ID2 == 2 || sim.ID2 == 3) integrate_vverlet(k, vscal[k]);
			else integrate_mts(k,vscal[k]);
			#endif		
			/* ----------------------------------------------- */
			/* Calculate the properties of the system.         */
			/* ----------------------------------------------- */
			kinet(k);    											
			calcvalue(k); 
			#ifndef PR_NPT
			boxscale(k);
			#endif			
			sim.T[k] += sim.dtemp[k]; 

			/* ----------------------------------------------- */
			/* Determine the velocity scaling factor           */
			/* ----------------------------------------------- */
			if (sim.ID2==9 || sim.ID2 ==10 || sim.ID2 == 1 || sim.ID2 ==2 || sim.ID2 ==3) vscal[k] = velscale(k); 

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
			for(int k=0; k<sim.NB; k++)	quatfit(k,rmsd_modes.file,rmsd_modes.atoms);
			#endif

			output(i,1);
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
			for (int k =0; k<sim.NB; k++) write_trr(k,i,1);
			#endif

			save_config(i,sim.cyc_pr);								
			save_veloc(i,sim.cyc_pr);								
			xyz_config(1,i);
			movie(1,i);
		}

      
    }//end production run
    /* ----------------------------------------------- */
    /* Call output on the last iteration everytime.    */
    /* ----------------------------------------------- */
	output(sim.cyc_pr,1);
	ioflush(0);
}

#endif //PRESSURE
