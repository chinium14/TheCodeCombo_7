/* ======================================================================== */
/* ofile.cpp                                                                */
/*                                                                          */
/*		This subroutine sets up the header for the output files, and        */
/* specifies what simulation modes were selected in simul.input.  It also   */
/* specifies what define options were used when the code was compiled.      */
/* The file is named simul#.output where # is 0,1,2,3,....  The number      */
/* If there is already a file in the OUTPUT folder named with a specific    */
/* number, the next number is used so as not to overwrite the other output  */
/* file.  The name of the file is stored in the global varible "name".      */
/*		This subroutine also specifies the name for the output file for     */
/* the energies.  This name is stored in the global variable "name1" and is */
/* called ener_box#.output, where # is the box number.                      */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

void output_headers(void);

void ofile (void)
{

  /* ================================================================== */
  /*                                                                    */
  /* Declare the varibles and file pointers needed for the output file. */
  /*                                                                    */
  /* ================================================================== */
  char host[100];  
  time_t now;
  FILE *sos;
  FILE *sos1;   

  /* ================================================================== */
  /*                                                                    */
  /* Create the output directories for each box.                        */
  /*                                                                    */
  /* ================================================================== */
#ifdef WIN

	_mkdir("./OUTPUT");

  for(int k=0; k<sim.NB; k++){
#ifdef MPI
	  sprintf(name,"./OUTPUT/BOX%d",mpi.my_rank);
#endif
#ifndef MPI
	  sprintf(name,"./OUTPUT/BOX%d",k);
#endif
	  
	  _mkdir(name);

#ifdef TWHAM
#ifdef MPI

	  sprintf(name,"./OUTPUT/BOX%d/TWHAM",mpi.my_rank);
#else
          sprintf(name,"./OUTPUT/BOX%d/TWHAM",k);
#endif
#endif
	  
	  _mkdir(name);

#ifdef XWHAM
	  sprintf(name,"./OUTPUT/BOX%d/XWHAM",k);
#endif
	  
	  _mkdir(name);

#ifdef ZHOU
#ifdef MPI
          sprintf(name,"./OUTPUT/BOX%d/STRESSES",mpi.my_rank);
#else
          sprintf(name,"./OUTPUT/BOX%d/STRESSES",k);
#endif
#endif
          _mkdir(name);

  }
#endif
#ifndef WIN

	mkdir("./OUTPUT", 0755);

  for(int k=0; k<sim.NB; k++){
#ifdef MPI
	  sprintf(name,"./OUTPUT/BOX%d",mpi.my_rank);
#endif
#ifndef MPI
	  sprintf(name,"./OUTPUT/BOX%d",k);
#endif
	  mkdir(name,0755);
#ifdef TWHAM
#ifdef MPI

          sprintf(name,"./OUTPUT/BOX%d/TWHAM",mpi.my_rank);
#else
          sprintf(name,"./OUTPUT/BOX%d/TWHAM",k);
#endif
#endif

    mkdir(name,0755);
#ifdef XWHAM
	  sprintf(name,"./OUTPUT/BOX%d/XWHAM",k);
#endif
    mkdir(name,0755);

#ifdef ZHOU
#ifdef MPI
          sprintf(name,"./OUTPUT/BOX%d/STRESSES",mpi.my_rank);
#else
          sprintf(name,"./OUTPUT/BOX%d/STRESSES",k);
#endif
#endif
    mkdir(name,0755);



  }
#endif

  /* ----------------------------------------------- */
  /* Make the directory   DOS if DOS is defined		 */
  /* ----------------------------------------------- */    
#ifdef DOS

#ifdef WIN
		sprintf(name,"./OUTPUT/DOS");
		_mkdir(name);
#endif
#ifndef WIN
		sprintf(name,"./OUTPUT/DOS");
		mkdir(name,0755);
#endif
#endif//DOS
#ifdef MMDOS

#ifdef WIN
		sprintf(name,"./OUTPUT/DOS");
		_mkdir(name);
#endif
#ifndef WIN
		sprintf(name,"./OUTPUT/DOS");
		mkdir(name,0755);
#endif
#endif//MMDOS
#ifdef CTDOS
#ifdef WIN
		sprintf(name,"./OUTPUT/DOS");
		_mkdir(name);
#endif
#ifndef WIN
		sprintf(name,"./OUTPUT/DOS");
		mkdir(name,0755);
#endif
#endif//CTDOS
#ifdef XEDOS
#ifdef WIN
		sprintf(name,"./OUTPUT/DOS");
		_mkdir(name);
#endif
#ifndef WIN
		sprintf(name,"./OUTPUT/DOS");
		mkdir(name,0755);
#endif
#endif//XEDOS
#ifdef FX_EDOS
#ifdef WIN
		sprintf(name,"./OUTPUT/DOS");
		_mkdir(name);
#endif
#ifndef WIN
		sprintf(name,"./OUTPUT/DOS");
		mkdir(name,0755);
#endif
#endif//FX_EDOS
#ifdef TDXEDOS
#ifdef WIN
		sprintf(name,"./OUTPUT/DOS");
		_mkdir(name);
#endif
#ifndef WIN
		sprintf(name,"./OUTPUT/DOS");
		mkdir(name,0755);
#endif
#endif//TDXEDOS
#ifdef WIN
	  sprintf(name,"./OUTPUT/GROMACS");
	  _mkdir(name);
#endif
#ifndef WIN
	  sprintf(name,"./OUTPUT/GROMACS");
	  mkdir(name,0755);
#endif



//#ifdef RMSD
  /* ================================================================== */
  /*                                                                    */
  /* Assign the name for the output file that contains the rmsd         */
  /* information for each box as rmsd#.output where # is the box        */
  /* number.                                                            */
  /*                                                                    */
  /* ================================================================== */
/*
	  char name2[100]; 
  for(int k=0; k<sim.NB; k++) {
	  FILE *statfile;               
#ifdef MPI
	  sprintf(name2,"./OUTPUT/BOX%d/rmsd%d.output",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	  sprintf(name2,"./OUTPUT/BOX%d/rmsd%d.output",k,k);
#endif
	  statfile= fopen(name2,"w");
	  fclose(statfile);
  }
#endif
*/
#ifdef EWALD
  /* ================================================================== */
  /*                                                                    */
  /* Assign the name for the output file that contains the different    */
  /* contributions to the coulombic energy as calculated by the EWALD   */
  /* method.  The file name is ewald#.output where # is the box number. */
  /*                                                                    */
  /* ================================================================== */
  char name3[100]; 
  for(int k=0; k<sim.NB; k++) {
	  FILE *ewaldfile;
#ifdef MPI
	  sprintf(name3,"./OUTPUT/BOX%d/ewald%d.output",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	  sprintf(name3,"./OUTPUT/BOX%d/ewald%d.output",k,k);
#endif
	  ewaldfile = fopen(name3,"w");
	  fclose(ewaldfile);
#ifdef EWALD_DEBUG
		FILE *ewaldfile2;
		char name4[100];
#ifdef MPI
	 sprintf(name4,"./OUTPUT/BOX%d/ewald_virial%d.output",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	  sprintf(name4,"./OUTPUT/BOX%d/ewald_virial%d.output",k,k);
#endif
		ewaldfile2= fopen(name4,"w");
		
		fclose(ewaldfile2);
#endif//EWALD_DEBUG
  }
#endif//EWALD

#ifdef PR_NPT
  char name5[100];
  FILE *nptfile;

  for(int k=0; k<sim.NB; k++){
#ifdef MPI
	sprintf(name5,"./OUTPUT/BOX%d/npt%d.output",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	sprintf(name5,"./OUTPUT/BOX%d/npt%d.output",k,k);
#endif
	nptfile= fopen(name5,"w");
	fclose(nptfile);
  }
#endif

/*#ifdef WALL_DEBUG
  char name7[100];
  FILE *wallfile;

  for(int k=0; k<sim.NB; k++){
	sprintf(name7,"./OUTPUT/BOX%d/angle%d.output",k,k);

	wallfile= fopen(name7,"w");
	fclose(wallfile);
  }
#endif
*/  
  /* ================================================================== */
  /*                                                                    */
  /* Assign the name for the output file that contains the energy       */
  /* information for each box as ener_box#.output where # is the box    */
  /* number.                                                            */
  /*                                                                    */
  /* ================================================================== */

char name6[50];
FILE *ordfile;

  for(int k=0; k<sim.NB; k++) {
#ifdef MPI
	 sprintf(name1,"./OUTPUT/BOX%d/ener_box%d.output",mpi.my_rank,mpi.my_rank);
 	 sprintf(name6,"./OUTPUT/BOX%d/order%d.output",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	 sprintf(name1,"./OUTPUT/BOX%d/ener_box%d.output",k,k);
	 sprintf(name6,"./OUTPUT/BOX%d/order%d.output",k,k);
#endif
    sos1 = fopen(name1,"w");
	ordfile = fopen(name6,"w");
	fclose(ordfile);
	fclose(sos1);
  }	
  
  for(int k=0; k<sim.NB; k++){
	  char swap_ptr[50];
	  FILE *swap1;
    #ifdef MPI
	  sprintf(swap_ptr,"./OUTPUT/swap%d.out",mpi.my_rank);
    #else
    sprintf(swap_ptr,"./OUTPUT/swap%d.out",k);
    #endif
    if((sim.ID == 2) || (sim.ID == 5) || (sim.ID == 7) || (sim.ID == 9) || (sim.ID == 10) 
      || (sim.ID == 11) || (sim.ID == 12)){
      swap1 = fopen(swap_ptr,"w");
      fclose(swap1);
	  }
  }
  
  time(&now);
  gethostname(host,50);

  /* ================================================================== */
  /*                                                                    */
  /* Assign the name for each simul#.output file, and write the         */
  /* information to this file.                                          */
  /*                                                                    */
  /* ================================================================== */
  for(int k = 0; k<sim.NB; k++){
	
	/* ----------------------------------------------- */
	/* Open the file.                                  */
	/* ----------------------------------------------- */    
#ifdef MPI
	  sprintf(name,"./OUTPUT/BOX%d/simul%d.output",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	 sprintf(name,"./OUTPUT/BOX%d/simul%d.output",k,k);
#endif
		sos = fopen(name,"w");

#ifdef SMD
		char nsmd[100];
		FILE *smd_ptr;
#ifdef MPI
	 sprintf(nsmd,"./OUTPUT/BOX%d/smd%d.output",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
 	 sprintf(nsmd,"./OUTPUT/BOX%d/smd%d.output",k,k);
#endif
	 smd_ptr = fopen(nsmd,"w");
	 fclose(smd_ptr);
#endif

#ifdef NSMD
		char nsmd[100];
		FILE *smd_ptr;
#ifdef MPI
	 sprintf(nsmd,"./OUTPUT/BOX%d/smd%d.output",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
 	 sprintf(nsmd,"./OUTPUT/BOX%d/smd%d.output",k,k);
#endif
	 smd_ptr = fopen(nsmd,"w");
	 fclose(smd_ptr);
#endif  
	/* ----------------------------------------------- */
	/* Output header.                                  */
	/* ----------------------------------------------- */    
	fprintf(sos,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n");
	fprintf(sos,"@@@@@@      Biological Molecule Simulation Code     @@@@@@ \n");
	fprintf(sos,"@@@@@@              by Nitin Rathore                @@@@@@ \n");
	fprintf(sos,"@@@@@@              & Thomas Knotts                 @@@@@@ \n");
	fprintf(sos,"@@@@@@       University of Wisconsin, Madison       @@@@@@ \n");
	fprintf(sos,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n");
	fprintf(sos,"\n");
	fprintf(sos,"Filename: %s \n",name);
	fprintf(sos,"Date and Time: %.24s \n",ctime(&now));
	fprintf(sos,"Hostname: %s \n",host);
	fprintf(sos,"\n");
	fprintf(sos,"%s\n",title);
  fprintf(sos,"Parameter File Name=%s\n",param_file_name);
	fprintf(sos,"\n");

	/* ----------------------------------------------- */
	/* Write the options selected in simul.input.      */
	/* ----------------------------------------------- */    
	if(sim.ID==1) fprintf(sos,"NVT MD simulation parameters: \n");
	if(sim.ID==2) fprintf(sos,"DOS    simulation parameters: \n");
	if(sim.ID==3) fprintf(sos,"NVE MD simulation parameters: \n");
	if(sim.ID==4) fprintf(sos,"Hybrid MD MC simulation parameters: \n");
	if(sim.ID==5) fprintf(sos,"Replica Exchange MD MC simulation parameters: \n");
	if(sim.ID==6) fprintf(sos,"NPT PR simulation parameters: \n");
	if(sim.ID==7) fprintf(sos,"MMDOS simulation parameters: \n");
	if(sim.ID==8) fprintf(sos,"Energy Minimization simulation parameters: \n");
	if(sim.ID==9) fprintf(sos,"CTDOS simulation parameters: \n");
	if(sim.ID==10)fprintf(sos,"XEDOS simulation parameters: \n\n");
	
	if(sim.ID2==1) fprintf(sos,"Leap Frog Integration/Berendson Thermostat: \n");
	if(sim.ID2==2) fprintf(sos,"Velocity Verlet Integration/Berendson Thermostat: \n");
	if(sim.ID2==3) fprintf(sos,"Multiple Time Step/Berendson Thermostat: \n");
	if(sim.ID2==4) fprintf(sos,"Nose-Hoover Chain Velocity Verlet: \n");
	if(sim.ID2==5) fprintf(sos,"MTS Nose-Hoover: \n");
	fprintf(sos,"Lengths of Box (x,y,z) in Angstroms			: %lf %lf %lf\n",box[k].lx, box[k].ly, box[k].lz);
  
	if(sim.ID2==1||sim.ID2==2||sim.ID2==3){
		fprintf(sos,"Temperature/coupling time/increment: %-6.2f K / %-4.2f ps / %-4.2f bar \n",sim.T[k],sim.Ttau[k],sim.dtemp[k]);
		if(sim.ID2 == 3) fprintf(sos, "Number of short time steps per long time step	: %d\n", sim.nsteps);
	}
	else if(sim.ID2 ==4 || sim.ID2 ==5 || sim.ID2==8)
	{	
		fprintf(sos,"Temperature							: %lf\n", sim.T[k]);
		fprintf(sos,"Number of NH Thermostats				: %d\n", nhc[k].M);
		fprintf(sos,"Mass of NH Thermostats(Q)				: %.2e\n", nhc[k].Q);
		fprintf(sos,"Number of NH Timesteps(Nc)				: %d\n", nhc[k].Nc);
		fprintf(sos,"Number of NH Integration Terms(Nys)			: %d\n", nhc[k].Nys);
	 

		if(sim.ID2 == 5){
			fprintf(sos, "Number of short time steps per long time step	: %d\n", sim.nsteps);
				if(sim.xRESPA[k]== 0) fprintf(sos, "XO-RESPA (MTS-NHC)\n\n");
				else fprintf(sos, "XI-RESPA (MTS-NHC)\n\n");
		}
	}
//	if(sim.ID==2 || sim.ID==4) fprintf(sos,"Pressure/coupling time/increment/compress.: %-5.2f bar / %-4.2f ps / %-4.2f bar / %-.3e \n",sim.P[k]/100.0,sim.Ptau[k],sim.dpress[k],sim.Pcomp[k]);
#ifdef DOS
	if(sim.ID==2){
		fprintf(sos,"Number of bins						: %d\n",sim_dos[k].e_bins);
		fprintf(sos,"Energy Range Begin					: %lf\n",sim_dos[k].e_begin);
		fprintf(sos,"Energy Range End						: %lf\n",sim_dos[k].e_end);
		fprintf(sos,"Temperature Range begin					: %lf\n",sim_dos[k].T_begin);
		fprintf(sos,"Temperature Range end					: %lf\n",sim_dos[k].T_end);
		fprintf(sos,"Bin Width							: %lf\n",sim_dos[k].e_width);
		fprintf(sos,"Number of MD steps per MC step				: %d\n", sim_hyb.cyc_hybmd);
		fprintf(sos,"Merge F							: %lf\n", MERGE_F);
		fprintf(sos,"Reset F							: %lf\n", RESET_F);
		fprintf(sos,"Stop F							: %lf\n", STOP_F);
	}
#endif
#ifdef MMDOS
	if(sim.ID==7){
		fprintf(sos,"Number of bins						: %d\n",sim_dos[k].e_bins);
		fprintf(sos,"Energy Range Begin					: %lf\n",sim_dos[k].e_begin);
		fprintf(sos,"Energy Range End						: %lf\n",sim_dos[k].e_end);
		fprintf(sos,"Temperature Range begin					: %lf\n",sim_dos[k].T_begin);
		fprintf(sos,"Temperature Range end				: %lf\n",sim_dos[k].T_end);
		fprintf(sos,"Bin Width						: %lf\n",sim_dos[k].e_width);
		fprintf(sos,"Number of MD steps per MC step				: %d\n", sim_hyb.cyc_hybmd);
	}
#endif
#ifdef CTDOS
	if(sim.ID==2){
		fprintf(sos,"Number of bins						: %d\n",sim_dos[k].e_bins);
		fprintf(sos,"Energy Range Begin					: %lf\n",sim_dos[k].e_begin);
		fprintf(sos,"Energy Range End						: %lf\n",sim_dos[k].e_end);
		fprintf(sos,"Temperature Range begin					: %lf\n",sim_dos[k].T_begin);
		fprintf(sos,"Temperature Range end					: %lf\n",sim_dos[k].T_end);
		fprintf(sos,"Bin Width							: %lf\n",sim_dos[k].e_width);
		fprintf(sos,"Number of MD steps per MC step				: %d\n", sim_hyb.cyc_hybmd);
		fprintf(sos,"Merge F							: %lf\n", MERGE_F);
		fprintf(sos,"Reset F							: %lf\n", RESET_F);
		fprintf(sos,"Stop F							: %lf\n", STOP_F);
	}
#endif
#ifdef XEDOS
	if(sim.ID==10){
		fprintf(sos,"Number of bins						: %d\n",sim_dos[k].l_bins);
		fprintf(sos,"Xi Range Begin                                 		: %lf\n",sim_dos[k].l_begin);
		fprintf(sos,"Xi Range End                                          	: %lf\n",sim_dos[k].l_end);
		fprintf(sos,"Temperature Range begin                                 : %lf\n",sim_dos[k].T_begin);
		fprintf(sos,"Temperature Range end                                   : %lf\n",sim_dos[k].T_end);
		fprintf(sos,"Bin Width                                               : %lf\n",sim_dos[k].l_width);
		fprintf(sos,"Pivot Move-Probability/Max Theta(rad)			: %lf/%lf\n",mc_pivot.PMS[k],mc_pivot.theta[k]);
		fprintf(sos,"Rigid Translation-Probability/Max Distance(Ang)      	: %lf/%lf\n",mc_trans.PMS[k],mc_trans.delta[k]);
                fprintf(sos,"Axial/Scale Move-Probability/Max strain(Ang)       	: %lf/%lf\n",mc_axial.PMS[k],mc_axial.strain[k]);
		fprintf(sos,"Random Atom Trans. Move-Probability/Max Distance(Ang)	: %lf/%lf\n",mc_rand.PMS[k],mc_rand.delta[k]);
		fprintf(sos,"Solvent Move-Probability/Max Distance(Ang)		: %lf/%lf\n",mc_solv.PMS[k],mc_solv.delta[k]);		
		fprintf(sos,"Number of MD steps per MC step                    	: %d\n", sim_hyb.cyc_hybmd);
		fprintf(sos, "Number of short time steps per long time step		: %d\n", sim.nsteps);
		fprintf(sos,"Merge F							: %lf\n", MERGE_F);
		fprintf(sos,"Reset F							: %lf\n", RESET_F);
		fprintf(sos,"Stop F							: %lf\n", STOP_F);
	}
#endif//XEDOS
#ifdef FX_EDOS
	if(sim.ID==11){
		fprintf(sos,"Number of bins						: %d\n",sim_dos[k].e_bins);
		fprintf(sos,"Energy Range Begin					: %lf\n",sim_dos[k].e_begin);
		fprintf(sos,"Energy Range End						: %lf\n",sim_dos[k].e_end);
		fprintf(sos,"Temperature Range begin					: %lf\n",sim_dos[k].T_begin);
		fprintf(sos,"Temperature Range end					: %lf\n",sim_dos[k].T_end);
		fprintf(sos,"Bin Width							: %lf\n",sim_dos[k].e_width);
		fprintf(sos,"Number of MD steps per MC step			: %d\n", sim_hyb.cyc_hybmd);
		fprintf(sos,"Target R value						: %lf\n",sim_dos[k].r_target);
		fprintf(sos,"Target round trips					: %lf\n",sim_dos[k].trips_target);
		fprintf(sos,"Merge F					: %lf\n", MERGE_F);
		fprintf(sos,"Reset F					: %lf\n", RESET_F);
		fprintf(sos,"Stop F					: %lf\n", STOP_F);
	}
#endif//FX_EDOS
	if(sim.ID==4 || sim.ID ==5){
		fprintf(sos,"Number of MD steps per MC step			: %ld\n", sim_hyb.cyc_hybmd);
	}
	fprintf(sos,"Number of equilibrium steps				: %lu\n",sim.cyc_eq);
	fprintf(sos,"Number of short production steps			: %lu\n",sim.cyc_pr);
	fprintf(sos,"Time step						: %-7.5f ps \n",sim.dt);
	fprintf(sos,"Cutoff for long range interactions			: %-4.2f A \n",sim.rc);
	fprintf(sos,"\n");
	fprintf(sos,"Number of components					: %d \n",sim.NC);
	fprintf(sos,"Number of simulation boxes				: %d \n",sim.NB);
	fprintf(sos,"Random Number Seed (idum)	 	  	      	: %ld\n",idum_seed);
	fprintf(sos,"\n");
	fprintf(sos,"Units for quantities: \n");
	fprintf(sos,"Temperature [=] K, Pressure [=] kPa, Density [=] kg/m3 \n");
	fprintf(sos,"Energy [=] kJ/mol, Time [=] ps \n");
	if(sim.ID2 ==4 || sim.ID2==5 || sim.ID2==8) fprintf(sos, "NHC Q [=] kg*ang^2\n");
	fprintf(sos,"\n");

	/* ----------------------------------------------- */
	/* Write the options defined when compiled.        */
	/* ----------------------------------------------- */    
	fprintf(sos, "Options Defined When Compiled.\n");

#ifdef COULOMB
	fprintf(sos,"COULOMB\n");
#endif
#ifdef WIN
	fprintf(sos,"WIN\n");
#endif
#ifdef NLIST
	fprintf(sos,"NLIST\n");
#endif
#ifdef SASA
	fprintf(sos,"SASA\n");
#endif
#ifdef NEUTRAL
	fprintf(sos,"NEUTRAL\n");
#endif
#ifdef RDIE
	 fprintf(sos,"RDIE\n");
#endif
#ifdef RFC
	fprintf(sos,"RFC\n");
#endif
#ifdef RMSD
	fprintf(sos, "RMSD\n");
#endif
#ifdef SHIFTF
	fprintf(sos,"SHIFTF\n");
#endif
#ifdef SHIFTV
	fprintf(sos,"SHIFTV\n");
#endif
#ifdef PRESSURE
	fprintf(sos,"PRESSURE\n");
#endif
#ifdef MOVIE
	fprintf(sos,"MOVIE\n");
#endif
#ifdef LJ
	fprintf(sos,"LJ\n");
#endif
#ifdef DOS
	fprintf(sos,"DOS\n");
#endif
#ifdef MMDOS
	fprintf(sos,"MMDOS\n");
#endif
#ifdef CTDOS
	fprintf(sos,"CTDOS\n");
#endif
#ifdef XEDOS
	fprintf(sos,"XEDOS\n");
#endif
#ifdef CONFIGT
	fprintf(sos,"CONFIGT\n");
#endif
#ifdef EWALD
	fprintf(sos,"EWALD\n");
#endif
#ifdef EWALD_DEBUG
	fprintf(sos,"EWALD_DEBUG\n");
#endif
#ifdef FCC
	fprintf(sos,"FCC\n");
#endif
#ifdef MPI
	fprintf(sos,"MPI\n");
#endif
#ifdef HESSIAN
	fprintf(sos,"HESSIAN\n");
#endif
#ifdef PR_NPT
	fprintf(sos,"PR_NPT\n");
#endif
#ifdef STYPE
	fprintf(sos,"STYPE\n");
#endif
#ifdef S14
	fprintf(sos,"S14\n");
#endif
#ifdef PHI
	fprintf(sos,"PHI\n");
#endif
#ifdef CRD
	fprintf(sos,"CRD\n");
#endif
#ifdef NMA
	fprintf(sos,"NMA\n");
#endif
#ifdef REST
	fprintf(sos,"REST\n");
#endif
#ifdef CLIST
	fprintf(sos,"CLIST\n");
#endif
#ifdef STATS
	fprintf(sos,"STATS\n");
#endif
#ifdef SMD
	fprintf(sos,"SMD\n");
#endif
#ifdef NSMD
	fprintf(sos,"NSMD\n");
#endif
#ifdef KONS
	fprintf(sos,"KONS\n");
#endif
#ifdef UMBP
	fprintf(sos,"UMBP\n");
#endif
#ifdef DLIST
	fprintf(sos,"DLIST\n");
#endif
#ifdef MC
	fprintf(sos,"MC\n");
#endif
#ifdef IONC
	fprintf(sos,"IONC\n");
#endif
#ifdef SLIST
	fprintf(sos,"SLIST\n");
#endif
#ifdef TLATE
	fprintf(sos,"TLATE\n");
#endif
#ifdef SPME
	fprintf(sos,"SPME\n");
#endif
#ifdef FLIM
	fprintf(sos,"FLIM\n");
#endif
	#ifdef SUG
	fprintf(sos,"SUG\n");
#endif
#ifdef TRR
	fprintf(sos,"TRR\n");
#endif
#ifdef FX_EDOS
	fprintf(sos,"FX_EDOS\n");
#endif
#ifdef GOLIK
	fprintf(sos, "GOLIK\n");
#endif
#ifdef GOBT
	fprintf(sos, "GOBT\n");
#endif
#ifdef DNA_GOLIK
	fprintf(sos, "DNA_GOLIK\n");
#endif
#ifdef WALL
	fprintf(sos, "WALL\n");
#endif
#ifdef AWALL
    fprintf(sos, "AWALL\n");
#endif
#ifdef NWALL
    fprintf(sos, "NWALL\n");
#endif
#ifdef RWALL
    fprintf(sos, "RWALL\n");
#endif
#ifdef HWALL
	fprintf(sos, "HWALL\n");
#endif
#ifdef CWALL
	fprintf(sos, "CWALL\n");
#endif
#ifdef WWALL
	fprintf(sos, "WWALL\n");
#endif
#ifdef SPHERE
	fprintf(sos, "SPHERE\n");
#endif
/*#ifdef WALL_DEBUG
	fprintf(sos, "WALL_DEBUG\n");
#endif
*/
#ifdef SEED
	fprintf(sos, "SEED\n");
#endif
#ifdef TDXEDOS
	fprintf(sos, "TDXEDOS\n");
#endif
#ifdef CAVITY
	fprintf(sos, "CAVITY\n");
#endif
#ifdef ACAVITY
	fprintf(sos, "ACAVITY\n");
#endif
#ifdef BCAVITY
	fprintf(sos, "BCAVITY\n");
#endif
#ifdef BETAPEP
	fprintf(sos, "BETAPEP\n");
#endif
#ifdef PREEQUIL
	fprintf(sos, "PREEQUIL\n");
#endif
#ifdef DEBHUCK
	fprintf(sos, "DEBHUCK\n");
#endif
#ifdef ZHOU
        fprintf(sos, "ZHOU\n");
#endif
#ifdef ZEROAM
        fprintf(sos, "ZEROAM\n");
#endif
#ifdef CHECKP
        fprintf(sos, "CHECKP\n");
#endif
#ifdef OHEADER
        fprintf(sos, "OHEADER\n");
#endif
#ifdef TWHAM
	fprintf(sos, "TWHAM\n");
#endif
#ifdef MM_ATTRACT
	fprintf(sos, "MM_ATTRACT\n");
#endif
#ifdef FSH
    fprintf(sos,"FSH\n");
#endif
#ifdef VSH
    fprintf(sos,"VSH\n");
#endif
#ifdef BGO
    fprintf(sos,"BGO\n");
#endif
	fprintf(sos, "\n");
#ifdef TRR
	FILE *trx;
	char nm[50];
	for(int i=0; i<100; i++) {
	#ifdef MPI
		  sprintf(nm,"./OUTPUT/GROMACS/gromacs%d_%d.trr",mpi.my_rank,i);
	#endif
	#ifndef MPI
		  sprintf(nm,"./OUTPUT/GROMACS/gromacs%d_%d.trr",k,i);
	#endif
		fileindex = i;
		if(i == 99) {
		  fprintf(stdout,"There are too many simul.ouput files (max allowed = 99) !!!\n");
		  exit(1);
		}
		else if( (trx=fopen(nm,"rb")) ) {
			fclose(trx);
			continue;
		}
		else {
		  trx = fopen(nm,"ab");
		  fclose(trx);
		  break;
		}
	}
#endif//TRR

#ifdef EWALD
#ifdef SPME
  fprintf(sos,"Long range electrostatic with Smooth Particle Mesh Ewald: Parameters in Box 0\n");
  fprintf(sos,"Kappa parameter: %lf / B-spline order: %d / Grid size: %d/%d/%d \n",
          sim.kappa[0],sim.order,sim.grid[0],sim.grid[1],sim.grid[2]);
#else
  fprintf(sos,"Long range electrostatic with standard Ewald summation : Parameters for Box 0\n");
  fprintf(sos,"Kappa parameter: %lf / Number of K-vectors: %d %d %d\n",sim.kappa[0],sim.kmax[0][0],sim.kmax[0][1],sim.kmax[0][2]);
#endif
#endif

	/* ----------------------------------------------- */
	/* Write the bond, bend, etc. info.                */
	/* ----------------------------------------------- */
	fprintf(sos,"Number of bonds			: %d	Number of bends			: %d\n",bondN[k], bendN[k]);
	fprintf(sos,"Number of torsions		: %d	Number of impropers		: %d\n",torsN[k], imprN[k]);
	fprintf(sos,"Number of H-bond acceptors	: %d	Number of H-bond donors		: %d\n", haccN[k],hdonN[k]);
	fprintf(sos,"Number of 1-4 Interactions	: %d	Number of NB exclusions		: %d\n",in14N[k], exNBN[k]);
        if(NUM_NBFIX == 0)
          fprintf(sos,"No NBFIX were detected.\n");
        else{
          fprintf(sos,"The Following NBFIX were used.\n");
          fprintf(sos,"Type 1, Type 2, epsilon, sigma\n");
          for(int i=0; i<NUM_NBFIX; i++)
            fprintf(sos,"%8d%8d%12.6lf   %.6lf\n",nbfix_prop[i].type1,nbfix_prop[i].type2,nbfix_prop[i].eps,nbfix_prop[i].sig);
        }
	fprintf(sos, "\n");

	
	fprintf(sos,"E14FAC = %f\n",E14FAC);

	/* ----------------------------------------------- */
	/* Write the spring coordinates if REST is         */
	/* defined.                                        */
	/* ----------------------------------------------- */
#ifdef REST
	fprintf(sos,"Restrained Parameters\n");
    for(int i=0; i<rest1[k].n; i++)
		fprintf(sos,"%d	%d	%lf	0.00000\n",rest1[k].site[i]+1,rest1[k].site[i]+1,rest1[k].k[i]/4.184);
	for(int i=0; i<rest2[k].n; i++)
		fprintf(sos,"%d	%d	%lf	%lf\n",rest2[k].site1[i]+1,rest2[k].site2[i]+1,rest2[k].k[i]/4.184,rest2[k].req[i]);
	fprintf(sos, "\n");
#endif

	/* ----------------------------------------------- */
	/* Write the native contacts.                      */
	/* ----------------------------------------------- */
	char con_filename[100];
	FILE *con_ptr;
	sprintf(con_filename, "./INPUT/contacts.inp");
	if( NULL != (con_ptr=fopen(con_filename,"r")) ){
		fprintf(sos,"Native Contacts\n");
		fprintf(sos,"%d\n",contacts[k].n);
		for(int i=0; i<contacts[k].n; i++){
				fprintf(sos,"%d	%d	%lf\n",contacts[k].a[i]+1, contacts[k].b[i]+1, contacts[k].d[i]);
		}
		fprintf(sos, "\n");
		fclose(con_ptr);
	}
	else fprintf(sos,"No native contacts were specified.\n\n");

#ifdef GOLIK
	/* ----------------------------------------------- */
	/* Write the Go-like native contacts.              */
	/* ----------------------------------------------- */
	fprintf(sos,"Go-like Model Parameters of Native Contacts\n");
	fprintf(sos,"Epsilon of beads =  %lf\n",golik[k].eps/4.184);
	fprintf(sos,"Equilibrium bond distance = %lf\n",golik[k].req);
	fprintf(sos,"Repulsive cut of = %lf\n",golik[k].d_repulsive);

	double factor = pow(2.0,-(1.0/6.0));
#ifndef DNA_GOLIK
  for(int a=0; a<golik[k].Psite; a++){
	  for(int b=a+1; b<golik[k].Psite; b++){
	    if(golik[k].contact[a][b] == 1){
	      double dr = ljset[k].pot[a][b].sig/factor;
              fprintf(sos,"%d	%d	%.13lf\n",a+1,b+1,dr);
	    } 
	  }
	}
#else
  for(int i=0; i<xintN[k]; i++){
   double sig = xint[k][i].sig/factor;
   if(sig>0.00001){
     fprintf(sos,"%d	%d	%.13lf\n",xint[k][i].a+1,xint[k][i].b+1,sig);
    }
  }
#endif
	fprintf(sos, "\n");
#endif

#ifdef WALL
	/* ----------------------------------------------- */
	/* Write the Wall parameters.                      */
	/* ----------------------------------------------- */
	fprintf(sos,"Wall Parameters\n");
	fprintf(sos,"Number of walls =  %d\n",wall[k].n);
	for(int i=0; i<wall[k].n; i++){
	  fprintf(sos,"Z location of Wall %d = %lf\n",i+1,wall[k].z[i]);
	}
	for(int i=0; i<wall[k].n; i++){
	  fprintf(sos,"Number of types of sites in Wall %d = %d\n",i+1,wall[k].nsites[i]);
	}
	for(int i=0; i<wall[k].n; i++){
	  for(int j=0; j<wall[k].nsites[i]; j++){
	    fprintf(sos,"Number density, sigma(ang), epsilon(kcal/mol), and hydropathy index of Wall %d = %lf, %lf, %.4e %lf\n",i+1,
			wall[k].num_density[i][j],wall[k].sig[i][j],wall[k].eps[i][j]/4.184,wall[k].hydro_index[i][j]);
	  }
	}

	fprintf(sos, "\n");
#endif


#ifdef XEDOS
	fprintf(sos,"Sites defining the reaction coordinate in XEDOS\n");
	fprintf(sos,"SITE1 = %d\nSITE2 = %d\nSITE3 = %d\n",SITE1+1,SITE2+1,SITE3+1);
	fprintf(sos,"\n");
#endif
	/* ----------------------------------------------- */
	/* Write the details of SMD simulation			   */
	/* ----------------------------------------------- */
#ifdef SMD
	fprintf(sos,"SMD Parameters\n");
	fprintf(sos,"Site 1 = Atom Number %d\n",steermd[k].site1);
	fprintf(sos,"Site 2 = Atom Number %d\n",steermd[k].site2);
	fprintf(sos,"Postion of Spring 1 (%lf,%lf,%lf)  : Spring constant 1(kJ/mol) :%lf\n",steermd[k].x1,
		steermd[k].y1,steermd[k].z1,steermd[k].k1);
	fprintf(sos,"Postion of Spring 2 (%lf,%lf,%lf)  : Spring constant 2(kJ/mol) :%lf\n",steermd[k].x2,
		steermd[k].y2,steermd[k].z2,steermd[k].k2);
	fprintf(sos,"Equilibrium end to end length of molecule		= %lf\n",steermd[k].r0);
	fprintf(sos,"Pulling Speed (A/ps)			= %lf\n",steermd[k].v);
	fprintf(sos,"Unit vector in direction of pull (%lf,%lf,%lf)\n",steermd[k].ex,steermd[k].ey,steermd[k].ez);
	
	fprintf(sos, "\n");
#endif//SMD
#ifdef NSMD
	fprintf(sos,"NSMD Parameters\n");
	fprintf(sos,"Site 1 = Atom Number %d\n",steermd[k].site1);
	fprintf(sos,"Site 2 = Atom Number %d\n",steermd[k].site2);
	fprintf(sos,"Postion of site 1 (%lf,%lf,%lf)  : \n",steermd[k].x1,steermd[k].y1,steermd[k].z1);
	fprintf(sos,"Postion of site 2 (%lf,%lf,%lf)  : \n",steermd[k].x2,steermd[k].y2,steermd[k].z2);
	fprintf(sos,"Equilibrium end to end length of molecule		= %lf\n",steermd[k].r0);
	fprintf(sos,"Pulling Speed (A/ps)			= %lf\n",steermd[k].v);
	fprintf(sos,"Unit vector in direction of pull (%lf,%lf,%lf)\n",steermd[k].ex,steermd[k].ey,steermd[k].ez);
	fprintf(sos, "\n");
#endif//NSMD

	/* ----------------------------------------------- */
	/* Write the constrained sites      			   */
	/* ----------------------------------------------- */
#ifdef KONS
	fprintf(sos,"Constrained Sites \n");
	for(int i=0; i<cons[k].n; i++) fprintf(sos,"%d\n",cons[k].site[i]+1);
	fprintf(sos, "\n");
#endif 
#ifdef CLIST
	fprintf(sos,"CELL_LENGTH_X	=	%lf\n",CELL_LENGTH_X);
	fprintf(sos,"CELL_LENGTH_Y	=	%lf\n",CELL_LENGTH_Y);
	fprintf(sos,"CELL_LENGTH_Z	=	%lf\n",CELL_LENGTH_Z);
#endif

	/* ----------------------------------------------- */
	/* Close the files.                                */
	/* ----------------------------------------------- */    
	fclose(sos);
  }//k

/********** Output the header files (ener_box, smd, order, prnpt_debug, ewald_debug) so we know which fields are what *********/
#ifdef OHEADER
	output_headers();
#endif


}
