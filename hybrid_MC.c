/* ======================================================================== */
/* hybrid_MC.cpp	                                                        */
/* This suboutine performs hybrid monte carlo moves. The velocities are		*/
/* assigned to atoms in each box based on the box temperature. The regular  */
/* Molecular dynamics is executed and after certain steps the MD steps are	*/
/* accepted or rejected based on the metropolis criteria now applied to the */
/* total hamiltonian i.e. the sum of potential and kinetic energy			*/
/* ======================================================================== */
	
#include "defines.h"
void init_block   (int);
void vinit_hyb    (int);
void save_config  (unsigned long,unsigned long);
void save_veloc   (unsigned long,unsigned long);
void xyz_config   (int,unsigned long);
void movie		  (int,unsigned long);
void ofile		  (void);
void output		  (unsigned long,int);
void outend		  (unsigned long);
void ioflush	  (int);
void svalues	  (int);
void vblock		  (int);
void pbc_all	  (int);
void forces		  (int);
void integrate_vverlet	(int, double);
void integrate_nhc		(double, int);
void kinet		  (int);
void calcvalue	  (int);
void drift		  (int);
void nblist		  (int);
double ran2		  (void);
double velscale   (int);

#ifdef RMSD
int  quatfit (int, int, int);
#endif

#ifdef TLATE
void translate (int);
#endif

#ifdef ROTN
void rotate_xyz(int,int);
#endif

#ifdef DLIST
void nl_check(int, int, int);
#endif

#ifdef TRR
void write_trr(int, unsigned long, int);
#endif



void hybrid_MC()
{
	  double vscal[NBOXES];
	  double H_initial[NBOXES];
	  double H_final[NBOXES];
/* **********************************************************************/
/*	NVT EQUILIBRATION RUN STARTS										*/
/* **********************************************************************/
 if (sim.cyc_eq != 0) {

	for(int k=0; k<sim.NB; k++) vscal[k]= 1.0;// initial scale for velocity

  /* ================================================================== */
  /*                                                                    */
  /* Begin equilibration iterations.                                    */
  /*                                                                    */
  /* ================================================================== */
	for(unsigned long i=1; i<=sim.cyc_eq; i++) { // MC counter
		for(int k=0; k<sim.NB; k++) {
			vinit_hyb(k);
			kinet(k);
			calcvalue(k);
			for(int l =0; l< box[k].boxns; l++){
				atom_temp[k][l]		= atom[k][l];		/* Back up coordinates with pdb			*/
				atnopbc_temp[k][l]	= atnopbc[k][l];	/* Back up coordinates without pdb		*/
				ff_temp[k][l]		= ff[k][l];			/* Back up all the force  components	*/
			}
#ifdef PRESSURE
			pvir_temp[k]	= pvir[k];					/* Back up all the virial components	*/
#endif
			en_temp[k]		= en[k];					/* Back up all the energy components	*/
			H_initial[k]	= en[k].totals;
		}
		
		for (int j=0;j < sim_hyb.cyc_hybmd; j++){
			for(int k=0; k<sim.NB; k++) {
  /* ----------------------------------------------- */
  /* Update the neighborlist with frequency          */
  /* sim.nlist.                                      */
  /* ----------------------------------------------- */
				#ifdef NLIST
				#ifndef DLIST
				if(j%sim.nlist==0) nblist(k);
				#endif
				#ifdef DLIST
					nl_check(0, box[k].boxns,k);
					if(nl_flag[k] == 1) nblist(k);
				#endif
				#endif
			
  /* ----------------------------------------------- */
  /* Itegration of the Equations of Motion.          */
  /* Note: If an NVE ensemble is chosen, velscale()  */
  /* is never called and vscal[k] is always 1.0.     */
  /* ----------------------------------------------- */
			
  /* *********************************************** */
  /* Velocity Verlet Algorithm  in NVE mode			 */
  /*                                                 */
  /* *********************************************** */
				integrate_vverlet(k, vscal[k]);									
			}// 1st k boxes ends
		}// j loop ends here
  /* ----------------------------------------------- */
  /* Calculate the properties of the system.         */
  /* ----------------------------------------------- */
		for (int k=0; k<sim.NB; k++){
			kinet(k);    											
			calcvalue(k);
			H_final[k]=en[k].totals;
		}// 2nd k loop ends here


  /* ----------------------------------------------- */
  /* Propose a MC move and check for acceptance. If  */
  /* If the moves is accepted then accumulate the	 */
  /* calculated values for block averages and error  */
  /* bars else first return to original values and   */
  /* then accumulate statistics						 */
  /* ----------------------------------------------- */
		for (int k=0; k<sim.NB; k++){
			double H_delta = H_final[k] - H_initial[k];
//			fprintf(stdout,"exp (-betaU)= %lf,%lf!!!",H_delta,exp(-H_delta/sim.kT[k]));

			if (H_delta <0.0 || exp(-H_delta/sim.kT[k]) > ran2()) {
				sim_hyb.hyb_acc[k] ++;
			}// moves is accepted
			else {										/* if move not accepted					*/
				for(int l =0; l< box[k].boxns; l++){
					atom[k][l]		= atom_temp[k][l];	/* Back up coordinates with pdb			*/
					atnopbc[k][l]	= atnopbc_temp[k][l];/* Back up coordinates without pdb		*/
					ff[k][l]		= ff_temp[k][l];		/* Back up all the force  components	*/
				}
#ifdef PRESSURE
				pvir[k]		= pvir_temp[k];				/* Back to old  virial components		*/
#endif
				en[k]		= en_temp[k];				/* Back to old  energy components		*/

			}
			calcvalue(k);

		}// 3rd k loop ends here
    
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
 }//end  if (sim.cyc_eq != 0)
  
	fprintf(stdout,"Equilibration completed!!!\n");

  /* ================================================================== */
  /* Apply pbc's, recalculate forces and properties of the system to    */
  /* prepare for production run.                                        */
  /* ================================================================== */
	for(int k=0; k<sim.NB; k++) {
		pbc_all(k);
		forces(k);   
		kinet(k);
		calcvalue(k);
		sim_hyb.hyb_acc[k]=0;//reset counter
	}

/* **********************************************************************/
/*	PRODUCTION RUN STARTS												*/
/* **********************************************************************/
  sim.dt *=1.5;
  for(int k=0; k<sim.NB; k++) vscal[k] = 1.0;// initial scale for velocity

  /* ================================================================== */
  /*                                                                    */
  /* Begin production iterations.                                       */
  /*                                                                    */
  /* ================================================================== */
	for(unsigned long i=1; i<=sim.cyc_pr; i++) { // MC counter
		for(int k=0; k<sim.NB; k++) {
			vinit_hyb(k);
			kinet(k);
			calcvalue(k);
			for(int l =0; l< box[k].boxns; l++){
				atom_temp[k][l]	= atom[k][l];			/* Back up coordinates with pdb			*/
				atnopbc_temp[k][l]	= atnopbc[k][l];	/* Back up coordinates without pdb		*/
				ff_temp[k][l]		= ff[k][l];			/* Back up all the force  components	*/
			}
#ifdef PRESSURE
			pvir_temp[k]	= pvir[k];					/* Back up all the virial components	*/
#endif
			en_temp[k]		= en[k];					/* Back up all the energy components	*/
			H_initial[k]	= en[k].totals;
		}
		
		for (int j=0;j < sim_hyb.cyc_hybmd; j++){
			for(int k=0; k<sim.NB; k++) {
  /* ----------------------------------------------- */
  /* Update the neighborlist with frequency          */
  /* sim.nlist. This is done if move is accepted.    */
  /* ----------------------------------------------- */
			
  /* ----------------------------------------------- */
  /* Itegration of the Equations of Motion.          */
  /* Note: If an NVE ensemble is chosen, velscale()  */
  /* is never called and vscal[k] is always 1.0.     */
  /* ----------------------------------------------- */
			
  /* *********************************************** */
  /* Velocity Verlet Algorithm  in NVE mode			 */
  /*                                                 */
  /* *********************************************** */
				integrate_vverlet(k, vscal[k]);									
			}// 1st k boxes ends
		}// j loop ends here
  /* ----------------------------------------------- */
  /* Calculate the properties of the system.         */
  /* ----------------------------------------------- */
		for (int k=0; k<sim.NB; k++){
			kinet(k);    											
			calcvalue(k);
			H_final[k]=en[k].totals;
		}// 2nd k loop ends here
	
  /* ----------------------------------------------- */
  /* Propose a MC move and check for acceptance. If  */
  /* If the moves is accepted then accumulate the	 */
  /* calculated values for block averages and error  */
  /* bars else first return to original values and   */
  /* then accumulate statistics						 */
  /* ----------------------------------------------- */
		for (int k=0; k<sim.NB; k++){
			double H_delta = H_final[k] - H_initial[k];
			if (H_delta <0.0 || exp(-H_delta/sim.kT[k]) > ran2()) {
				sim_hyb.hyb_acc[k] ++;
				#ifdef NLIST
				#ifndef DLIST
				if(sim_hyb.hyb_acc[k] %sim.nlist==0) nblist(k);
				#endif
				#ifdef DLIST
					nl_check(0, box[k].boxns,k);
					if(nl_flag[k] == 1) nblist(k);
				#endif
				#endif

			}// moves is accepted
			else {											/* if move not accepted					*/
				for(int l =0; l< box[k].boxns; l++){
					atom[k][l]		= atom_temp[k][l];		/* Back up coordinates with pdb			*/
					atnopbc[k][l]	= atnopbc_temp[k][l];	/* Back up coordinates without pdb		*/
					ff[k][l]		= ff_temp[k][l];		/* Back up all the force  components	*/
				}
#ifdef PRESSURE
				pvir[k]		= pvir_temp[k];					/* Back to old  virial components		*/
#endif
				en[k]		= en_temp[k];					/* Back to old  energy components		*/
			}
			calcvalue(k);
			svalues(k);
		}// 3rd k loop ends here

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
		} //blockc ends  
		
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
      
    }// end production run	MC counter i														
  /* ----------------------------------------------- */
  /* Call output on the last iteration everytime.    */
  /* ----------------------------------------------- */
	output(sim.cyc_pr,1);
	ioflush(0);
  }
