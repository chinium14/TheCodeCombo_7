/*
TODO:

DONE  Make read_sim.c to read in Na_con and sim.epsRF[k]
DONE  Put DRSAM stuff in cgolik.c
DONE  Put DRSAM stuff in read_dnagolik.c
Find all other changes in the DNA code.... and put them in 06Aug2010
*/


/* ==================== Molecular Simulation Code V1.0 ==================== */
/*                                                                          */
/* Nitin Rathore & Tommy Knotts                                             */
/* Graduate Research Assistants                                             */
/* University of Wisconsin, Madison                                         */
/* Department of Chemical Engineering                                       */
/* Madison, WI 53706                                                        */
/* Completed May 31, 2002.                                                  */
/*                                                                          */
/* ======================================================================== */
/* Define options:                                                          */
/*   COULOMB    :       Coulombic interaction potential                     */
/*   NLIST      :       Use neighbor-list                                   */
/*   SASA       :       Implicit solvent model                              */
/*   NEUTRAL    :       Side chains of amino aicds neutralized              */
/*   RDIE       :       Distance dependent dielectric (with SASA)           */
/*   MOVIE      :       Sample pdb files to view trajectories               */
/*   PRESSURE   :       Enable pressure calculation                         */
/*                                                                          */
/* ==================== Molecular Simulation Code V2.0 ==================== */
/*                                                                          */ 
/*   RMSD       :       Enable superposition and rmsd computations          */
/*   LJ         :       Disable everything except Lennard Jones             */
/*                                                                          */
/* Hybrid MD MC and Replica exchange MD  added to Version 1.0               */
/* Multiple Box simulations enabled                                         */
/* Completed June 5, 2002.                                                  */
/*                                                                          */
/* ==================== Molecular Simulation Code V3.0 ==================== */
/*                                                                          */
/* This version of the code contains the improper torsion energy terms.     */
/* Also, some minor changes to the output files were made.                  */
/* Completed June 19, 2002.                                                 */
/*                                                                          */
/* ==================== Molecular Simulation Code V4.0 ==================== */
/*                                                                          */ 
/* This version of the code contains Density of States Algorithm. Also      */
/* vinit_hyb.cpp is added to be used by hybrid MC, replica exchange and DOS */
/* simulations. Also, the neighbor list update part need to be corercted    */
/* for hybrid MC, replica exchange and DOS simulations. To be done in V5.0  */
/* Define options:                                                          */
/*   DOS        :       Density of states simulation (sim.ID ==2)           */
/*                                                                          */
/* Completed July 10, 2002.                                                 */
/*                                                                          */
/* ==================== Molecular Simulation Code V5.0 ==================== */
/*                                                                          */ 
/* Some corrections to the LJ dynamics and properties were done. The        */
/* long-range correction to the virial was also corrected.  A pivot         */
/* move for DOS simualtions was also introduced.                            */
/* Define options:                                                          */
/*   FCC        :       Built an FCC lattice when LJ is defined.            */
/*   PIVOT      :       Attempt a pivot move of the protein.  Used          */
/*                      in MC with DOS.                                     */
/*                                                                          */
/* Completed July 30, 2002.                                                 */
/*                                                                          */
/* ==================== Molecular Simulation Code V6.0 ==================== */
/*                                                                          */ 
/* This version of the code added Ewald Sums for the calculation of         */
/* electrostatic interactions.                                              */
/* Define options:                                                          */
/*   EWALD      :       Enable calculation of electrostatic                 */
/*                      interactions by the use of Ewald Sums.              */
/*                                                                          */
/* Completed 29 Aug 02                                                      */
/*                                                                          */ 
/* ==================== Molecular Simulation Code V6.1 ==================== */
/*                                                                          */ 
/* This version of the code changes a few lines in Ewald subroutines to     */
/* reduce the computation time.                                             */
/* Completed 9 Sep 02                                                       */
/*                                                                          */ 
/* ==================== Molecular Simulation Code V6.2 ==================== */
/*                                                                          */ 
/* This version of the code now uses CHARMM22 parameters.                   */
/* Completed 10 Sep 02                                                      */
/*                                                                          */
/* ==================== Molecular Simulation Code V7.0 ==================== */
/*                                                                          */ 
/* This version of the code changed the read_ subroutines to use multiple   */
/* molecules. Also, some minor changes to the memory allocation for Ewald   */
/* sums was done to reduce memory issues.                                   */
/*                                                                          */
/* Completed 3 Oct 02                                                       */
/*                                                                          */
/* ==================== Molecular Simulation Code V9.0 ==================== */
/*                                                                          */ 
/* This version of the code has a ton of stuff added.                       */
/* Define option:                                                           */
/*   PR_NPT     :       Parinello-Rahman NPT                                */
/*                                                                          */
/* Completed 28 Feb 03                                                      */
/*                                                                          */ 
/* ======================================================================== */


#define MAIN
#include "defines.h"
  /* ================================================================== */
  /*                                                                    */
  /* Function prototypes                                                */
  /*                                                                    */
  /* ================================================================== */
void init_all     (void);
void init_param	  (void);
void init_block   (int);
void hybrid_MC	  (void);
void ljpot		  (void);
void ofile		  (void);
void outend		  (unsigned long);
void repexch	  (void);
void read_atom	  (void);
void read_solvent(void);
void read_config  (void);
void read_sim     (void);
void setlj		  (void);
void vinit        (void);
void weight       (void);
void nve_md		  (void);
void nvt_md		  (void);	
void nvt_mc     (void);

void read_veloc   (int*);
void xyz_config   (int,unsigned long);
void movie		  (int,unsigned long);
void vblock		  (int);
void pbc_all	  (int);
void forces		  (int);
void force_short  (int);
void force_long   (int);
void isokin       (int);
void kinet		  (int);
void calcvalue	  (int);
void drift		  (int);
void nblist		  (int);
void save_config  (unsigned long,unsigned long);
void save_veloc   (unsigned long,unsigned long);
void output		  (unsigned long,int);
void dos_awrite	  (unsigned long,int);
void dos_wwrite	  (unsigned long,int);
void ioflush	  (int);
void mem_free     (void);

void checks		  (void);
#ifdef WWALL 
void hydro_index_bgo(void);
#endif
#ifdef SPHERE
void hydro_index_bgo(void);
#endif
double ran2		  (void);
#ifdef NMA
double minim();
#endif

#ifdef PRESSURE
double virialcor (int);
#endif

#ifdef SASA
double csasa (int);
#ifdef SASAREX
double csasab (int);
#endif
#endif

#ifdef DOS
void dos(void);
#endif
#ifdef MMDOS
void mmdos(void);
#endif
#ifdef CTDOS
void ctdos(void);
#endif
#ifdef XEDOS
void xedos(void);
#endif //XEDOS
#ifdef FX_EDOS
void fx_edos(void);
#endif//FX_EDOS
#ifdef TDXEDOS
void tdxedos(void);
void tdxedos_init(void);
#endif
#ifdef EWALD
void ewald_setup(int);
#ifdef SPME
void pme_clean(void);
#endif //SPME
#endif
#ifdef PR_NPT
void boxinv(int);
void pbc_npt_full(int);

#endif
void npt_md(void);
#ifdef SMD	
void smd_setup (int);
#endif
#ifdef NSMD	
void smd_setup (int);
#endif
#ifdef CLIST
void init_cell(void);
#endif
#ifdef KONS
void read_cons(void);
#endif
#ifdef REST
void read_rest(void);
#endif
void write_param (void);
#ifdef GOLIK
#ifndef DNA_GOLIK
void read_golik(void);
#else
void read_dnagolik(void);
#endif
#endif
#ifdef WALL
void read_wall(void);
#endif
#ifdef PREEQUIL
void preequil(void);
#endif
#ifdef TRR
void write_trr(int, unsigned long, int);
#endif
#ifdef CHECKP
void InitCheckpoint();
#endif

/* ************************************************************************ */
/*                                                                          */
/* ====================== Begin main Subroutine =========================== */
/*                                                                          */
/* ************************************************************************ */
int main(int argc, char** argv)
{
#ifdef MPI
 MPI_Init(&argc,&argv);
 MPI_Comm_rank(MPI_COMM_WORLD,&mpi.my_rank);
 MPI_Comm_size(MPI_COMM_WORLD,&mpi.p);
#endif

#ifdef CHECKP
//Detect the presence of checkpoint files, register signal handler, etc.
InitCheckpoint();
#endif


  /* ================================================================== */
  /*                                                                    */
  /* Variables and file pointers needed in this subroutine.             */
  /*                                                                    */
  /* ================================================================== */
  int vtest; 
  FILE *sout;

  /* ================================================================== */
  /*                                                                    */
  /* Initialize the random number generator ran2 with negative seed     */
  /* value for idum.                                                    */
  /*                                                                    */
  /* ================================================================== */
#ifdef SEED
  time_t idum_clock;
  time(&idum_clock);
  idum = -1 * (int) idum_clock + getpid();
#else
  //idum = -6931;
  idum = -1192464077;
#endif 
  idum_seed = idum;
  ran2();

  /* ================================================================== */
  /*                                                                    */
  /* Initialization                                                     */
  /*                                                                    */
  /* ================================================================== */
  /* ----------------------------------------------- */
  /* Get the overall simulation parameters from the  */
  /* simul.input file.                               */
  /* ----------------------------------------------- */
  read_sim();
  /* ----------------------------------------------- */
  /* Initialize all the property variables to zero   */
  /* ----------------------------------------------- */
  init_all();
  /* ----------------------------------------------- */
  /* Read in the force field parameters from the     */
  /* simul19.input file.                             */
  /* ----------------------------------------------- */
  read_atom();
  /* ----------------------------------------------- */
  /* Read in the info on constrained sites.          */
  /* ----------------------------------------------- */
#ifdef KONS
  read_cons();
#endif
  /* ----------------------------------------------- */
  /* Read in the wall paramters.                     */
  /* ----------------------------------------------- */
#ifdef WALL
  read_wall();
#endif

  /* ----------------------------------------------- */
  /* Read in the x,y,z coordinates of the sites.     */
  /* ----------------------------------------------- */
  read_config();
  /* ----------------------------------------------- */
  /* Read in the info on constrained sites.          */
  /* ----------------------------------------------- */
#ifdef REST
  read_rest();
#endif
  /* ----------------------------------------------- */
  /* Check to make sure options don't conflict.      */
  /* ----------------------------------------------- */
  checks();
  /* ----------------------------------------------- */
  /* If a cell list is used, initialize the cells.   */
  /* ----------------------------------------------- */
#ifdef CLIST
  init_cell();
#endif

  /* ----------------------------------------------- */
  /* If EWALD is defined, set up the k-space         */
  /* vectors.                                        */
  /* ----------------------------------------------- */
#ifdef EWALD
  for(int k=0; k<sim.NB; k++) ewald_setup(k);
#endif

  /* ----------------------------------------------- */
  /* If PR_NPT is defined, calculate the inverse box */
  /* matrix.                                         */  
  /* ----------------------------------------------- */
#ifdef PR_NPT
  for(int k=0; k<sim.NB; k++){
	  boxinv(k);
  }
#endif

#ifdef SMD
  for(int k=0; k<sim.NB; k++) smd_setup(k);
  checks();
#endif
#ifdef NSMD
  for(int k=0; k<sim.NB; k++) smd_setup(k);
  checks();
#endif
  /* ----------------------------------------------- */
  /* If LJ is defined then over ride all the         */
  /* partameters like eps, sig, bonds, bends etc.    */  
  /* ----------------------------------------------- */
#ifdef LJ
  ljpot();
#endif
  /* ----------------------------------------------- */
  /* Calculate the LJ parameters for each site/site  */
  /* interaction.                                    */
  /* ----------------------------------------------- */
	setlj();

  /* ----------------------------------------------- */
  /* Initialize all the property variables to zero   */
  /* ----------------------------------------------- */
  init_param();
#ifndef XEDOS
  for(int k=0; k<sim.NB; k++) init_block(k);
#endif

  #ifdef TWHAM
  for(int k=0; k<sim.NB; k++){
    pe_min[k]=ELIMIT;
    pe_max[k]=-ELIMIT;
  }
  #endif

  /* ----------------------------------------------- */
  /* Read in Go-like model parameters.               */
  /* ----------------------------------------------- */
#ifdef GOLIK
#ifndef DNA_GOLIK
  read_golik();
#else
  read_dnagolik();
#endif
#endif

#ifdef BGO
#ifdef WWALL
hydro_index_bgo();
#endif
#ifdef SPHERE
hydro_index_bgo();
#endif
#endif

  /* ----------------------------------------------- */
  /* Read the velocities in from ./INPUT/simul.vel#  */
  /* file. If this file doesn't exist, generate      */
  /* random velocities.                              */
  /* ----------------------------------------------- */
  read_veloc(&vtest);
  /* ----------------------------------------------- */
  if(vtest != 0) vinit();

  /* ----------------------------------------------- */
  /* Calculate the total mass in each box.           */
  /* ----------------------------------------------- */
  weight();

  /* ----------------------------------------------- */
  /* Calculate the long range correction to the      */
  /* virial.                                         */
  /* ----------------------------------------------- */
#ifdef PRESSURE 
 for(int k=0; k<sim.NB; k++) pvir[k].lrc = virialcor(k);
#endif   
  /* ----------------------------------------------- */
  /* Variables to time the simulation.               */
  /* ----------------------------------------------- */
  time_t now0;
  time(&now0);
#ifndef MPI
  for(int k=0; k<sim.NB; k++){
#else
  for(int k=0; k<mpi.p; k++){
#endif
	  time(&simul_time[k]);
	  time(&ener_time[k]);
	  time(&ord_time[k]);
	  time(&swap_time[k]);
	  #if defined(SMD) || defined(NSMD)  
		time(&smd_time[k]);
	  #endif
    #ifdef TWHAM
    time(&wham_time[k]);
    #endif
  }
  /* ----------------------------------------------- */
  /* Create the output file containing the           */
  /* simulation data.                                */
  /* ----------------------------------------------- */

//Don't call ofile() if we are restarting!
#ifdef CHECKP
if(restart_equil == 0)
#endif
{
	  ofile();


	  /* ----------------------------------------------- */
	  /* Calculate the neighborlist, forces, energies,   */
	  /* and other properties.                           */
	  /* ----------------------------------------------- */
	  for(int k=0; k<sim.NB; k++) {
/*
        #ifndef MPI
          for(int k=0; k<sim.NB; k++){
        #else
          for(int k=0; k<mpi.p; k++){
        #endif
*/
	#ifdef NLIST
	    nblist(k);   
	#endif

	#ifndef GOLIK
	  if(sim.ID != 6){
			force_short(k); force_long(k); 
		}
	#endif

  	    forces(k);   
	    kinet(k);
	    calcvalue(k);  
	  }
	  xyz_config(0,0);
	  movie(0,0);
	  if(sim.ID ==2   || sim.ID ==7    ||sim.ID == 9 ||sim.ID == 10 || 
		  sim.ID == 11 || sim.ID == 12){
		  dos_awrite(0,0);
		  dos_wwrite(0,0);
	  }
	  else output(0,0);
	  write_param();
	  #ifdef TRR
	  write_trr(0,0,0);
	  #endif
	 ioflush(0);
}
  mem_free();
        

  fprintf(stdout,"Initialization completed!!!\n");

  /* ================================================================== */
  /* Preform pre-equilibration if it is defined.                        */
  /* ================================================================== */
  #ifdef PREEQUIL
  preequil();
  #endif

  /* ----------------------------------------------- */
  /*  The following modes can be assigned to sim.ID	 */
  /*												 */
  /*	Mode 1	:	NVT MD simulation				 */
  /*	Mode 2	:	DOS simulation					 */
  /*	Mode 3	:	NVE MD simulation				 */
  /*	Mode 4	:	NVT Hybrid MD MC simulation		 */
  /*	Mode 5	:	Replica Exchange simulation		 */
  /*	Mode 6	:	NPT Molecular Dynamics			 */
  /*	Mode 7	:	MMDOS simulation				 */
  /*	Mode 8  :	PE Minimization					 */
  /*	Mode 9	:	CTDOS simulation				 */
  /*	Mode 10	:	XEDOS simulation				 */
  /* ----------------------------------------------- */

#ifdef DOS
  if (sim.ID !=2) {
	  fprintf(stdout,"DOS simulation should be chosen if DOS is defined!!!\n");
	  exit(1);
  }
#endif

  if	  (sim.ID == 1)	nvt_md();
#ifdef DOS
  else if (sim.ID == 2) dos();
#endif
  else if (sim.ID == 3) nve_md();
  //else if (sim.ID == 4) hybrid_MC();
  else if (sim.ID == 4) nvt_mc();
  else if (sim.ID == 5) repexch();
#ifdef PRESSURE
  else if (sim.ID == 6)	npt_md();
#endif
#ifdef MMDOS
  else if (sim.ID == 7) mmdos();
#endif
#ifdef NMA
  else if (sim.ID == 8) double min_ener = minim();
#endif
#ifdef CTDOS
  else if (sim.ID == 9) ctdos();
#endif
#ifdef XEDOS
  else if (sim.ID == 10) xedos();
#endif
#ifdef FX_EDOS
  else if (sim.ID == 11) fx_edos();
#endif
#ifdef TDXEDOS
  else if (sim.ID ==12) tdxedos();
#endif
#ifdef SPME
  pme_clean();
#endif

  /* ================================================================== */
  /*                                                                    */
  /* Determine the simulation averages (average of the block values)    */
  /* and error bars and write them to file.  Then save the final        */
  /* velocities and configurations.                                     */
  /*                                                                    */
  /* ================================================================== */
  if(sim.ID !=2 && sim.ID !=7 && sim.ID != 9) outend(sim.cyc_pr);
  save_config(0,0);				/* (block, cycle)						*/
  save_veloc (0,0);				/* (block, cycle)						*/
  xyz_config (1,sim.cyc_pr);	/* flag==1 means simulfinal#.pdb saved	*/ 
  movie(1,sim.cyc_pr);
  ioflush(1);


  /* ================================================================== */
  /*                                                                    */
  /* Open the output file, calculate the simulation time, and write     */
  /* this data to the file and to the screen.							*/
  /*                                                                    */
  /* ================================================================== */
  for (int k=0; k<sim.NB; k++) {
#ifdef MPI
	sprintf(name, "./OUTPUT/BOX%d/simul%d.output",mpi.my_rank,mpi.my_rank); 
#endif
#ifndef MPI
	sprintf(name, "./OUTPUT/BOX%d/simul%d.output",k,k); 
#endif
	sout = fopen(name,"a");
	time_t now1;
	time(&now1);
	double total_time = (difftime(now1,now0))/60.0;
  
	fprintf(stdout,"\nElapsed run time [min.]:   %10.3f\n",total_time);
	fprintf(sout,  "\nElapsed run time [min.]:   %10.3f\n",total_time);
	fprintf(stdout,"Date:                      %.24s \n",ctime(&now1));
	fprintf(sout,  "Date:                      %.24s \n",ctime(&now1));
	char host[50];
	gethostname(host,50);
	fprintf(stdout,"Hostname:                  %s\n",host);
	fprintf(sout,  "Hostname:                  %s\n",host);
	fprintf(stdout,"End of simulation!!!\n\n");
	fprintf(sout,  "End of simulation!!!\n\n");
	fclose(sout);
  }
/* Finalize any pending messages for MPI*/

/*  	fprintf(stdout,"npr_1 is %lf \n",npr_1);
  	fprintf(stdout,"npr_2 is %lf \n",npr_2);
  	fprintf(stdout,"npr_3 is %lf \n",npr_3); 
*/
#ifdef MPI
  MPI_Finalize();
#endif

  exit(0);return 0;
}


