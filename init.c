/* ======================================================================== */
/* init.cpp                                                                 */
/*                                                                          */
/*		This file contains two subroutines: init_all() and init_block.      */
/*                                                                          */
/* init_all(void):                                                          */
/*		This subroutine initializes the various variables that store the    */
/* the energy and other properties of the boxes.  It sets all the values to */
/* zero and is only called once at the beginning of main().                 */
/*                                                                          */
/* init_block(int):                                                        */
/*		This subroutine resets the block accumulators (res[k]) for the      */
/* properties calculated at each time step.  It is used to generate many    */
/* different "block averages" values for statistical purposes.              */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"
#ifdef TDXEDOS
void tdxedos_init(void);
#endif
#ifdef XWHAM
void xwham_init(void);
#endif
/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void init_all (void)
{

  /* ================================================================== */
  /*                                                                    */
  /* Initialize box variables.                                          */
  /*                                                                    */
  /* ================================================================== */
  for(int k=0; k<sim.NB; k++) {
    box[k].boxn		= 0;
    box[k].boxns	= 0;
	box[k].boxnres	= 0;
    box[k].vol		= 0.0;
    box[k].weight	= 0.0;
    box[k].dens		= 0.0;
    box[k].temp		= 0.0;
    box[k].press	= 0.0;
    box[k].pressmol = 0.0;
    box[k].nfree	= 0.0;
    for(int i=0; i<6; i++) {
      box[k].cpress[i]		= 0.0;
      box[k].cpressmol[i]	= 0.0;
    }

  /* ================================================================== */
  /*                                                                    */
  /* Initialize en variables.                                           */
  /*                                                                    */
  /* ================================================================== */
    en[k].total  = 0.0;
    en[k].totals = 0.0;
    en[k].poten  = 0.0;
    en[k].potens = 0.0;
    en[k].kinet  = 0.0;
    en[k].kinetmol  = 0.0;
    for(int i=0; i<6; i++) {
      en[k].ken[i] = 0.0;
      en[k].kenmol[i] = 0.0;
    }
    en[k].nbond  = 0.0;
    en[k].nbonds = 0.0;
    en[k].bond   = 0.0;
    en[k].bend   = 0.0;
#ifndef NEUTRAL
	en[k].urey	 = 0.0;
#endif
    en[k].tors   = 0.0;
	en[k].impr   = 0.0;
	en[k].coulomb= 0.0;
	en[k].H      = 0.0;
	sim_hyb.hyb_acc[k]	 = 0;
	sim_hyb.swap_acc[k]	 = 0;
	mc_pivot.pivot_acc[k]= 0;

	
  /* ================================================================== */
  /*                                                                    */
  /* Initialize oreder parameter variables								*/
  /*                                                                    */
  /* ================================================================== */
	ordparam[k].d_nc	 = 0.0;
	ordparam[k].d_nca	 = 0.0;
	ordparam[k].d_ncb	 = 0.0;
	ordparam[k].hel		 = 0.0;
	ordparam[k].hela	 = 0.0;
	ordparam[k].helb	 = 0.0;
	ordparam[k].con		 = 0.0;
	ordparam[k].cona	 = 0.0;
	ordparam[k].con_2	 = 0.0;
	ordparam[k].conb	 = 0.0;
	ordparam[k].force_1	 = 0.0;
	ordparam[k].force_2	 = 0.0;
	ordparam[k].gyr		 = 0.0;
	ordparam[k].gyra	 = 0.0;
	ordparam[k].gyrb	 = 0.0;
	ordparam[k].rmsd	 = 0.0;
	ordparam[k].rmsda	 = 0.0;
	ordparam[k].rmsdb	 = 0.0;
	ordparam[k].x1       = 0.0;
	ordparam[k].x1a      = 0.0;
	ordparam[k].x1b      = 0.0;
	ordparam[k].x2       = 0.0;
	ordparam[k].x2a      = 0.0;
	ordparam[k].x2b      = 0.0;


  /* ================================================================== */
  /*                                                                    */
  /* Initialize bp variables.                                           */
  /*                                                                    */
  /* ================================================================== */
    for(int i=0; i<sim.NC; i++) {
      bp[k][i].first = 0;
      bp[k][i].last  = 0;
      bp[k][i].nbox  = 0;
    }
 
 
  /* ================================================================== */
  /*                                                                    */
  /* Initialize pvir variables.                                         */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
    
	pvir[k].tnbond = 0.0;
    pvir[k].tnbondmol = 0.0;
    pvir[k].tbond  = 0.0;
    pvir[k].tbend  = 0.0;
#ifndef NEUTRAL
	pvir[k].turey  = 0.0;
#endif
	pvir[k].ttors  = 0.0;
	pvir[k].timpr  = 0.0;
    pvir[k].lrc    = 0.0;
#ifdef SASA
	pvir[k].tsasa  = 0.0;
#endif
    for(int i=0; i<6; i++) {
      pvir[k].nbond[i] = 0.0;
      pvir[k].bond[i]  = 0.0;
      pvir[k].bend[i]  = 0.0;
#ifndef NEUTRAL
	  pvir[k].urey[i]  = 0.0;
#endif
      pvir[k].tors[i]  = 0.0;
	  pvir[k].impr[i]  = 0.0;
      pvir[k].sum[i]   = 0.0;
      pvir[k].summol[i]= 0.0;
#ifdef SASA
	  pvir[k].sasa[i]  = 0.0;
#endif
#ifdef EWALD
	  pvir[k].ewald_real[i]=0.0;
	  pvir[k].ewald_kspace[i]=0.0;
	  pvir[k].ewald_intra[i]=0.0;
#endif//EWALD
    }

  /* ================================================================== */
  /*                                                                    */
  /* Initialize res variables.                                          */
  /*                                                                    */
  /* ================================================================== */
    res[k].pressa    = 0.0;	res[k].pressb    = 0.0;	res[k].pressc    = 0.0;

    for(int i=0; i<6; i++) {
      res[k].cpressa[i]  = 0.0;  
      res[k].cpressb[i]  = 0.0;  
      res[k].cpressc[i]  = 0.0;  
    }
#endif //endif for PRESSURE

	res[k].count     = 0;
    res[k].tempa     = 0.0;	res[k].tempb     = 0.0;	res[k].tempc     = 0.0;	
    res[k].densa     = 0.0; res[k].densb     = 0.0; res[k].densc     = 0.0;
    res[k].ebonda    = 0.0; res[k].ebondb    = 0.0; res[k].ebondc    = 0.0;
    res[k].ebenda    = 0.0; res[k].ebendb    = 0.0; res[k].ebendc    = 0.0;
#ifndef NEUTRAL
	res[k].eureya	 = 0.0; res[k].eureyb    = 0.0; res[k].eureyc    = 0.0;
#endif
    res[k].etorsa    = 0.0; res[k].etorsb    = 0.0; res[k].etorsc    = 0.0;
	res[k].eimpra    = 0.0; res[k].eimprb    = 0.0; res[k].eimprc	 = 0.0;
    res[k].enbonda   = 0.0; res[k].enbondb   = 0.0; res[k].enbondc   = 0.0;
    res[k].enbondsa  = 0.0; res[k].enbondsb  = 0.0; res[k].enbondsc  = 0.0;
    res[k].etpota	 = 0.0; res[k].etpotb    = 0.0; res[k].etpotc    = 0.0;
    res[k].etpotsa   = 0.0; res[k].etpotsb   = 0.0; res[k].etpotsc   = 0.0;
    res[k].etkina    = 0.0; res[k].etkinb    = 0.0; res[k].etkinc    = 0.0;
    res[k].etotala   = 0.0; res[k].etotalb   = 0.0; res[k].etotalc   = 0.0;
    res[k].etotalsa  = 0.0; res[k].etotalsb  = 0.0; res[k].etotalsc  = 0.0;
	res[k].Ha		 = 0.0;	res[k].Hb		 = 0.0; res[k].Hc        = 0.0;
#ifdef ELASTIC
	for(int i=0; i<21; i++){
		res[k].Cijklb[i]    = 0.0; res[k].Cijklc[i] = 0.0;
	}
#endif
	
  /* ================================================================== */
  /*                                                                    */
  /* Initialize resb variables.                                         */
  /*                                                                    */
  /* ================================================================== */
	resb[k].count    = 0;
	resb[k].tempa    = 0.0; resb[k].tempb    = 0.0; resb[k].tempc    = 0.0;
	resb[k].densa    = 0.0; resb[k].densb    = 0.0; resb[k].densc    = 0.0;
    resb[k].ebonda   = 0.0; resb[k].ebondb   = 0.0; resb[k].ebondc   = 0.0;
    resb[k].ebenda   = 0.0; resb[k].ebendb   = 0.0; resb[k].ebendc   = 0.0;
#ifndef NEUTRAL
	resb[k].eureya	 = 0.0; resb[k].eureyb   = 0.0; resb[k].eureyc  = 0.0;
#endif
    resb[k].etorsa   = 0.0; resb[k].etorsb   = 0.0; resb[k].etorsc   = 0.0;
    resb[k].eimpra   = 0.0; resb[k].eimprb   = 0.0; resb[k].eimprc   = 0.0;
    resb[k].enbonda  = 0.0; resb[k].enbondb  = 0.0; resb[k].enbondc  = 0.0;
    resb[k].enbondsa = 0.0; resb[k].enbondsb = 0.0; resb[k].enbondsc = 0.0;
    resb[k].etpota   = 0.0; resb[k].etpotb   = 0.0; resb[k].etpotc   = 0.0;
    resb[k].etpotsa  = 0.0; resb[k].etpotsb  = 0.0; resb[k].etpotsc  = 0.0;
    resb[k].etkina   = 0.0; resb[k].etkinb   = 0.0; resb[k].etkinc   = 0.0;
    resb[k].etotala  = 0.0; resb[k].etotalb  = 0.0; resb[k].etotalc  = 0.0;
    resb[k].etotalsa = 0.0; resb[k].etotalsb = 0.0; resb[k].etotalsc = 0.0;
	resb[k].Ha       = 0.0;	resb[k].Hb       = 0.0; resb[k].Hc       = 0.0;
#ifdef PRESSURE	
	resb[k].pressa   = 0.0; resb[k].pressb   = 0.0; resb[k].pressc   = 0.0;
    for(int i=0; i<6; i++) {
		resb[k].cpressa[i]  = 0.0;  
		resb[k].cpressb[i]  = 0.0;  
		resb[k].cpressc[i]  = 0.0;  
    }
#endif
#ifdef ELASTIC
	for(int i=0; i<21; i++){
		resb[k].Cijklb[i]   = 0.0; resb[k].Cijklc[i] = 0.0;
	}
#endif


	frames[k].count=0;
	}// loop k ends
	
  /* ================================================================== */
  /*                                                                    */
  /* Initialize NPT variables.                                          */
  /*                                                                    */
  /* ================================================================== */

#ifdef PR_NPT
	if(sim.ID == 6){
		if(sim.ID2==6){
			for(int k=0; k<sim.NB; k++){
				delta[0] = 1.0; 
				delta[1] = 0.0;
				delta[2] = 0.0;
				delta[3] = 0.0;
				delta[4] = 1.0;
				delta[5] = 0.0;
				delta[6] = 0.0;
				delta[7] = 0.0;
				delta[8] = 1.0;
				for(int i=0; i<9; i++){
					axesi[k][i]= 0.0;
					va[k][i]   = 0.0;
					Ga[k][i]   = 0.0;
					vvab[k][i] = 0.0;
					frab[k][i] = 0.0;
				}// for i
			}//for k
		}//if sim.ID2 == 6
	}
#endif
	#ifdef CONFIGT
	for(int k=0; k<sim.NB; k++){
		config[k].hesx =0.0;
		config[k].hesy =0.0;
		config[k].hesz =0.0;
		config[k].hesr =0.0;
		config[k].num = 0.0;
		config[k].den = 0.0;
		config[k].T	  = 0.0;	
	}
#endif


}// init_all ends


void init_param (void){
  for(int k=0; k<sim.NB; k++) {
	/* ----------------------------------------	*/
	/* Initialize phi and psi variable			*/
	/* ----------------------------------------	*/
	for(int i=0; i<torsN[k]; i++) {
		tors[k][i].phia  = 0.0;
		tors[k][i].phi	 = 0.0;
		tors[k][i].psia  = 0.0;
		tors[k][i].psi	 = 0.0;
		tors[k][i].thetaa= 0.0;
		tors[k][i].theta = 0.0;
		tors[k][i].count = 1;
	}
  }//end k loop

	/*--------------------------------------*/
	/* INITIALIZATION OF PIVOT MATRICES		*/
	/*--------------------------------------*/
	
	sim.dtlong	 = sim.dt*sim.nsteps;	//for MTS
#ifdef DOS
  /* ================================================================== */
  /*                                                                    */
  /* Initialize DOS variables.                                          */
  /*                                                                    */
  /* ================================================================== */
	for(int k=0; k<sim.NB; k++) {
		sim_dos[k].e_width		=  (int) (sim_dos[k].e_width*sim.kT[k] + 0.5);
		sim_dos[k].e_bins		=  (int) (0.5+ (sim_dos[k].e_end - sim_dos[k].e_begin)/sim_dos[k].e_width);
		if(sim_dos[k].e_bins >= NBINS){
			printf("Not enough memory allocated for number of bins. Increase NBINS to at least %i in defines.h\n", sim_dos[k].e_bins+1);
			exit(9);
		}
		sim_dos[k].dos_acc		= 0;
		mc_axial.axial_acc[k]   =0;
		mc_rand.rand_acc[k]     =0;
		mc_trans.trans_acc[k]	=0;
		char name6[100];
		FILE *dos3;

	#ifdef MPI
		sprintf(name6,"./INPUT/flathist1_%d.output",mpi.my_rank);
	#endif
	#ifndef MPI
		sprintf(name6,"./INPUT/flathist1_%d.output",k);
	#endif
		
		if( NULL == (dos3=fopen(name6,"r")) ) {
			sim_dos[k].mod_f		= 1.0;
			for ( int i=0; i<sim_dos[k].e_bins; i++){
				dos_hist[k][i].e_mid		= sim_dos[k].e_begin + sim_dos[k].e_width*(i+0.5);
				dos_hist[k][i].g_of_e		=	0.0;	
				dos_hist[k][i].h_of_e		=	0;
				average[k][i].count			=	1;
				average[k][i].ebond			=	0.0;
				average[k][i].ebend			=	0.0;
#ifndef NEUTRAL
				average[k][i].eurey			=	0.0;
#endif
				average[k][i].etors			=	0.0;
				average[k][i].eimpr			=	0.0;
				average[k][i].elj			=	0.0;
				average[k][i].ecoul			=	0.0;
				average[k][i].epotens		=	0.0;
				average[k][i].d_nc			=	0.0;
				average[k][i].hel			=	0.0;
				average[k][i].con			=	0.0;
				average[k][i].con_2			=	0.0;
				average[k][i].gyr			=	0.0;
				average[k][i].rmsd			=	0.0;
	#ifdef SASA
				average[k][i].esasa			=	0.0;
	#endif//SASA
			}
		}
		else {// if simulation is to be started from old flathist files.
			char name7[100];
			FILE *dos_sim;
			#ifdef MPI
				sprintf(name7,"./INPUT/simul%d.output",mpi.my_rank);
			#endif
			#ifndef MPI
				sprintf(name7,"./INPUT/simul%d.output",k);
			#endif
		
			if( NULL == (dos_sim=fopen(name7,"r")) ) {
				int dummy;
				for ( int i=0; i<sim_dos[k].e_bins; i++){
					fscanf(dos3, "%d	%lf	%lf	%lf	%d\n", &dummy, &sim_dos[k].mod_f, &dos_hist[k][i].e_mid, 
						&dos_hist[k][i].g_of_e, &dos_hist[k][i].h_of_e);
			//		dos_hist[k][i].h_of_e=0;
					dos_hist_temp[k][i]	=	dos_hist[k][i];
					if (dummy!= i) {
						fprintf(stdout,"incompatibility in reading flathist1_%d.output file\n",k);
						exit(10);
					}
					if(dos_hist[k][i].e_mid != (sim_dos[k].e_begin + sim_dos[k].e_width*(i+0.5))){
						fprintf(stdout,"incompatibility in assigning e_mid while reading flathist1.txt\n");
						exit(10);
					}
					average[k][i].count			=	1;
					average[k][i].ebond			=	0.0;
					average[k][i].ebend			=	0.0;
#ifndef NEUTRAL
					average[k][i].eurey			=	0.0;
#endif
					average[k][i].etors			=	0.0;
					average[k][i].eimpr			=	0.0;
					average[k][i].elj			=	0.0;
					average[k][i].ecoul			=	0.0;
					average[k][i].epotens		=	0.0;
					average[k][i].d_nc			=	0.0;
					average[k][i].hel			=	0.0;
					average[k][i].con			=	0.0;
					average[k][i].con_2			=	0.0;
					average[k][i].gyr			=	0.0;
					average[k][i].rmsd			=	0.0;
#ifdef SASA
					average[k][i].esasa			=	0.0;
#endif //SASA
				}
				sim_dos[k].mod_f =log(sim_dos[k].mod_f);
				fclose(dos3);
			}//NULL = dos_sim
			else{//old simul_output files to be read too

				int dummy;
				for ( int i=0; i<sim_dos[k].e_bins; i++){
					fscanf(dos3, "%d	%lf	%lf	%lf	%d\n", &dummy, &sim_dos[k].mod_f, &dos_hist[k][i].e_mid, 
						&dos_hist[k][i].g_of_e, &dos_hist[k][i].h_of_e);
					dos_hist_temp[k][i]	=	dos_hist[k][i];
					if (dummy!= i) {
						fprintf(stdout,"incompatibility in reading flathist1_%d.output file\n",k);
						exit(10);
					}
				#ifndef SASA
					#ifdef NEUTRAL
					fscanf(dos_sim, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",&average[k][i].epotensa, 
						&average[k][i].ebonda,&average[k][i].ebenda,&average[k][i].etorsa, &average[k][i].eimpra,
						&average[k][i].elja,&average[k][i].ecoula,&average[k][i].d_nca,&average[k][i].hela,&average[k][i].cona,
						&average[k][i].gyra,&average[k][i].rmsda,&average[k][i].count);
					
					#endif//NEUTRAL
					#ifndef NEUTRAL

					fscanf(dos_sim, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",&average[k][i].epotensa, 
						&average[k][i].ebonda,&average[k][i].ebenda,&average[k][i].eureya,&average[k][i].etorsa, &average[k][i].eimpra,
						&average[k][i].elja,&average[k][i].ecoula,&average[k][i].d_nca,&average[k][i].hela,&average[k][i].cona,
						&average[k][i].gyra,&average[k][i].rmsda,&average[k][i].count);
					
					#endif//NEUTRAL
				#endif//SASA
				#ifdef SASA
					#ifdef NEUTRAL
					fscanf(dos_sim, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",&average[k][i].epotensa,
						&average[k][i].ebonda,&average[k][i].ebenda,&average[k][i].etorsa,&average[k][i].eimpra,
						&average[k][i].elja,&average[k][i].ecoula,&average[k][i].esasaa,&average[k][i].d_nca,&average[k][i].hela,
						&average[k][i].cona,&average[k][i].gyra,
						&average[k][i].rmsda,&average[k][i].count);
					#endif
				#endif


					average[k][i].ebond			=	average[k][i].count*average[k][i].ebonda;
					average[k][i].ebend			=	average[k][i].count*average[k][i].ebenda;
#ifndef NEUTRAL
					average[k][i].eurey			=	average[k][i].count*average[k][i].eureya;
#endif
					average[k][i].etors			=	average[k][i].count*average[k][i].etorsa;
					average[k][i].eimpr			=	average[k][i].count*average[k][i].eimpra;
					average[k][i].elj			=	average[k][i].count*average[k][i].elja;
					average[k][i].ecoul			=	average[k][i].count*average[k][i].ecoula;
					average[k][i].epotens		=	average[k][i].count*average[k][i].epotensa;
					average[k][i].d_nc			=	average[k][i].count*average[k][i].d_nca;
					average[k][i].hel			=	average[k][i].count*average[k][i].hela;
					average[k][i].con			=	average[k][i].count*average[k][i].cona;
					average[k][i].gyr			=	average[k][i].count*average[k][i].gyra;
					average[k][i].rmsd			=	average[k][i].count*average[k][i].rmsda;
					#ifdef SASA
					average[k][i].esasa			=	average[k][i].count*average[k][i].esasaa;
					#endif //SASA
				}
				fclose(dos_sim);
				sim_dos[k].mod_f =log(sim_dos[k].mod_f);
				fclose(dos3);
			}//else if simul.output available

		}//else if flathist.output is available
	}//k
#endif//DOS

#ifdef MPI
	if((mpi.my_rank+1)%2==0) mpi.flag =0; else mpi.flag=1;

#endif

#ifdef MMDOS
  /* ================================================================== */
  /*                                                                    */
  /* Initialize MMDOS variables.                                        */
  /*                                                                    */
  /* ================================================================== */
	for(int k=0; k<sim.NB; k++) {
		sim_dos[k].e_width		=  (int) (sim_dos[k].e_width*sim.kT[k] + 0.5);
		sim_dos[k].e_bins		=  (int) (0.5+ (sim_dos[k].e_end - sim_dos[k].e_begin)/sim_dos[k].e_width);
		if(sim_dos[k].e_bins >= NBINS){
			printf("Not enough memory allocated for number of bins. Increase NBINS to at least %d in defines.h\n", sim_dos[k].e_bins+1);
			exit(9);
		}
		sim_dos[k].dos_acc		= 0;
		sim_dos[k].dos_acc		= 0;
		char name6[100];
		FILE *dos3;

	#ifdef MPI
		sprintf(name6,"./INPUT/flathist1_%d.output",mpi.my_rank);
	#endif
	#ifndef MPI
		sprintf(name6,"./INPUT/flathist1_%d.output",k);
	#endif
		sim.T[k]  = sim_dos[k].T_begin;
		sim.kT[k] = RG * sim.T[k]*.001;
		if( NULL == (dos3=fopen(name6,"r")) ) {
			sim_dos[k].mod_f		= 1.0;	//convergence criteria
			for ( int i=0; i<sim_dos[k].e_bins; i++){
				dos_hist[k][i].e_mid		= sim_dos[k].e_begin + sim_dos[k].e_width*(i+0.5);
				dos_hist[k][i].g_of_e		=	0.0;	
				dos_hist[k][i].h_of_e		=	1;	// as k_of_e not initialized to zero
				average[k][i].count			=	1;
				average[k][i].ebond			=	0.0;
				average[k][i].ebend			=	0.0;
#ifndef NEUTRAL
				average[k][i].eurey			=	0.0;
#endif
				average[k][i].etors			=	0.0;
				average[k][i].eimpr			=	0.0;
				average[k][i].elj			=	0.0;
				average[k][i].ecoul			=	0.0;
				average[k][i].epotens		=	0.0;
				average[k][i].d_nc			=	0.0;
				average[k][i].hel			=	0.0;
				average[k][i].con			=	0.0;
				average[k][i].con_2			=	0.0;
	#ifdef SASA
				average[k][i].esasa			=	0.0;
	#endif//SASA
			}
		}
		else {// if simulation is to be started from old flathist files.
			int dummy;
			for ( int i=0; i<sim_dos[k].e_bins; i++){
				fscanf(dos3, "%d	%lf	%lf	%lf	%d	%lf\n", &dummy, &sim_dos[k].mod_f, &dos_hist[k][i].e_mid, 
					&dos_hist[k][i].g_of_e, &dos_hist[k][i].h_of_e, &dos_hist[k][i].k_of_e);
				dos_hist[k][i].h_of_e=1;
				dos_hist_temp[k][i]	=	dos_hist[k][i];
				if (dummy!= i) {
					fprintf(stdout,"incompatibility in reading flathist1_%d.output file\n",k);
					exit(10);
				}
				if(dos_hist[k][i].e_mid != (sim_dos[k].e_begin + sim_dos[k].e_width*(i+0.5))){
					fprintf(stdout,"incompatibility in assigning e_mid while reading flathist1.txt\n");
					exit(10);
				}
				average[k][i].count			=	1;
				average[k][i].ebond			=	0.0;
				average[k][i].ebend			=	0.0;
#ifndef NEUTRAL
				average[k][i].eurey			=	0.0;
#endif
				average[k][i].etors			=	0.0;
				average[k][i].eimpr			=	0.0;
				average[k][i].elj			=	0.0;
				average[k][i].ecoul			=	0.0;
				average[k][i].epotens		=	0.0;
				average[k][i].d_nc			=	0.0;
				average[k][i].hel			=	0.0;
				average[k][i].con			=	0.0;
	#ifdef SASA
				average[k][i].esasa			=	0.0;
	#endif //SASA
			}
// claculate			sim_dos[k].mod_f =log(sim_dos[k].mod_f);
			fclose(dos3);	
		}
	}
#endif//MMDOS
#ifdef CTDOS
  /* ================================================================== */
  /*                                                                    */
  /* Initialize CTDOS variables.                                        */
  /*                                                                    */
  /* ================================================================== */
	for(int k=0; k<sim.NB; k++) {
		sim_dos[k].e_width		=  (int) (sim_dos[k].e_width*sim.kT[k] + 0.5);
		sim_dos[k].e_bins		=  (int) (0.5+ (sim_dos[k].e_end - sim_dos[k].e_begin)/sim_dos[k].e_width);
		if(sim_dos[k].e_bins >= NBINS){
			printf("Not enough memory allocated for number of bins. Increase NBINS to at least %i in defines.h\n", sim_dos[k].e_bins+1);
			exit(9);
		}
		sim_dos[k].dos_acc		= 0;
		sim_dos[k].dos_acc		= 0;
		mc_axial.axial_acc[k]   =0;
		mc_rand.rand_acc[k]     =0;
		mc_trans.trans_acc[k]	=0;
		char name6[100];
		FILE *dos3;

	#ifdef MPI
		sprintf(name6,"./INPUT/flathist1_%d.output",mpi.my_rank);
	#endif
	#ifndef MPI
		sprintf(name6,"./INPUT/flathist1_%d.output",k);
	#endif
		
		if( NULL == (dos3=fopen(name6,"r")) ) {
			sim_dos[k].mod_f		= 1.0;
			for ( int i=0; i<sim_dos[k].e_bins; i++){
				dos_hist[k][i].e_mid		= sim_dos[k].e_begin + sim_dos[k].e_width*(i+0.5);
				dos_hist[k][i].g_of_e		=	0.0;
				dos_hist[k][i].gT_of_e		=	0.0;
				dos_hist[k][i].h_of_e		=	0;
				dos_hist[k][i].ct_num  		=	0.0;
				dos_hist[k][i].ct_den  		=	0.0;
				average[k][i].count			=	1;
				average[k][i].ebond			=	0.0;
				average[k][i].ebend			=	0.0;
#ifndef NEUTRAL
				average[k][i].eurey			=	0.0;
#endif
				average[k][i].etors			=	0.0;
				average[k][i].eimpr			=	0.0;
				average[k][i].elj			=	0.0;
				average[k][i].ecoul			=	0.0;
				average[k][i].epotens		=	0.0;
				average[k][i].d_nc			=	0.0;
				average[k][i].hel			=	0.0;
				average[k][i].con			=	0.0;
				average[k][i].con_2			=	0.0;
				average[k][i].gyr			=	0.0;
				average[k][i].rmsd			=	0.0;
				#ifdef SASA
				average[k][i].esasa			=	0.0;
				#endif//SASA
			}
		}
		else {// if simulation is to be started from old flathist files.
			char name7[100];
			FILE *dos_sim;
			#ifdef MPI
				sprintf(name7,"./INPUT/simul%d.output",mpi.my_rank);
			#endif
			#ifndef MPI
				sprintf(name7,"./INPUT/simul%d.output",k);
			#endif
		
			if( NULL == (dos_sim=fopen(name7,"r")) ) {
				int dummy;
				for ( int i=0; i<sim_dos[k].e_bins; i++){
					fscanf(dos3, "%d	%lf	%lf	%lf	%d	%lf\n", &dummy, &sim_dos[k].mod_f, &dos_hist[k][i].e_mid, 
						&dos_hist[k][i].g_of_e, &dos_hist[k][i].h_of_e, &dos_hist[k][i].gT_of_e);
			//		dos_hist[k][i].h_of_e=0;
					dos_hist_temp[k][i]	=	dos_hist[k][i];
					if (dummy!= i) {
						fprintf(stdout,"incompatibility in reading flathist1_%d.output file\n",k);
						exit(10);
					}
					if(dos_hist[k][i].e_mid != (sim_dos[k].e_begin + sim_dos[k].e_width*(i+0.5))){
						fprintf(stdout,"incompatibility in assigning e_mid while reading flathist1.txt\n");
						exit(10);
					}
					average[k][i].count			=	1;
					average[k][i].ebond			=	0.0;
					average[k][i].ebend			=	0.0;
	#ifndef NEUTRAL
					average[k][i].eurey			=	0.0;
	#endif
					average[k][i].etors			=	0.0;
					average[k][i].eimpr			=	0.0;
					average[k][i].elj			=	0.0;
					average[k][i].ecoul			=	0.0;
					average[k][i].epotens		=	0.0;
					average[k][i].d_nc			=	0.0;
					average[k][i].hel			=	0.0;
					average[k][i].con			=	0.0;
					average[k][i].con_2			=	0.0;
					average[k][i].gyr			=	0.0;
					average[k][i].rmsd			=	0.0;
					#ifdef SASA
					average[k][i].esasa			=	0.0;
					#endif //SASA
				}
				sim_dos[k].mod_f =log(sim_dos[k].mod_f);
				fclose(dos3);
			}//if simul.output is not present
			else{//old simul_output files to be read too

				int dummy;
				for ( int i=0; i<sim_dos[k].e_bins; i++){
					fscanf(dos3, "%d	%lf	%lf	%lf	%d	%lf\n", &dummy, &sim_dos[k].mod_f, &dos_hist[k][i].e_mid, 
						&dos_hist[k][i].g_of_e, &dos_hist[k][i].h_of_e, &dos_hist[k][i].gT_of_e);
					dos_hist_temp[k][i]	=	dos_hist[k][i];
					if (dummy!= i) {
						fprintf(stdout,"incompatibility in reading flathist1_%d.output file\n",k);
						exit(10);
					}
				#ifndef SASA
					#ifdef NEUTRAL
					fscanf(dos_sim, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",&average[k][i].epotensa, 
						&average[k][i].ebonda,&average[k][i].ebenda,&average[k][i].etorsa, &average[k][i].eimpra,
						&average[k][i].elja,&average[k][i].ecoula,&average[k][i].d_nca,&average[k][i].hela,&dos_hist[k][i].config_T,&average[k][i].cona,
						&average[k][i].gyra,&average[k][i].rmsda,&average[k][i].count);
					
					#endif//NEUTRAL
					#ifndef NEUTRAL

					fscanf(dos_sim, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",&average[k][i].epotensa, 
						&average[k][i].ebonda,&average[k][i].ebenda,&average[k][i].eureya,&average[k][i].etorsa, &average[k][i].eimpra,
						&average[k][i].elja,&average[k][i].ecoula,&average[k][i].d_nca,&average[k][i].hela,&dos_hist[k][i].config_T,&average[k][i].cona,
						&average[k][i].gyra,&average[k][i].rmsda,&average[k][i].count);
					
					#endif//NEUTRAL
				#endif//SASA
				#ifdef SASA
					#ifdef NEUTRAL
					fscanf(dos_sim, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",&average[k][i].epotensa,
						&average[k][i].ebonda,&average[k][i].ebenda,&average[k][i].etorsa,&average[k][i].eimpra,
						&average[k][i].elja,&average[k][i].ecoula,&average[k][i].esasaa,&average[k][i].d_nca,&average[k][i].hela,&dos_hist[k][i].config_T,
						&average[k][i].cona,&average[k][i].gyra,
						&average[k][i].rmsda,&average[k][i].count);
					#endif
				#endif


					average[k][i].ebond			=	average[k][i].count*average[k][i].ebonda;
					average[k][i].ebend			=	average[k][i].count*average[k][i].ebenda;
#ifndef NEUTRAL
					average[k][i].eurey			=	average[k][i].count*average[k][i].eureya;
#endif
					average[k][i].etors			=	average[k][i].count*average[k][i].etorsa;
					average[k][i].eimpr			=	average[k][i].count*average[k][i].eimpra;
					average[k][i].elj			=	average[k][i].count*average[k][i].elja;
					average[k][i].ecoul			=	average[k][i].count*average[k][i].ecoula;
					average[k][i].epotens		=	average[k][i].count*average[k][i].epotensa;
					average[k][i].d_nc			=	average[k][i].count*average[k][i].d_nca;
					average[k][i].hel			=	average[k][i].count*average[k][i].hela;
					average[k][i].con			=	average[k][i].count*average[k][i].cona;
					average[k][i].gyr			=	average[k][i].count*average[k][i].gyra;
					average[k][i].rmsd			=	average[k][i].count*average[k][i].rmsda;
					#ifdef SASA
					average[k][i].esasa			=	average[k][i].count*average[k][i].esasaa;
					#endif //SASA
				}
				fclose(dos_sim);
				sim_dos[k].mod_f =log(sim_dos[k].mod_f);
				fclose(dos3);
			}//else if simul.output available

		}
	}
#endif//CTDOS
#ifdef XEDOS
  /* ================================================================== */
  /*                                                                    */
  /* Initialize XEDOS variables.                                        */
  /*                                                                    */
  /* ================================================================== */
	for(int k=0; k<sim.NB; k++) {
//		sim_dos[k].l_width		= (int) (sim_dos[k].l_width*sim_dos[k].l_begin + 0.5);
		sim_dos[k].l_bins		= (int) (0.5+ (sim_dos[k].l_end - sim_dos[k].l_begin)/sim_dos[k].l_width);
		if(sim_dos[k].l_bins >= NBINS){
			printf("Not enough memory allocated for number of bins. Increase NBINS to at least %i in defines.h\n", sim_dos[k].l_bins+1);
			exit(9);
		}
		sim_dos[k].dos_acc		= 0;
		sim_dos[k].dos_acc		= 0;
		mc_axial.axial_acc[k]   =0;
		mc_rand.rand_acc[k]     =0;
		mc_trans.trans_acc[k]	=0;
		char name6[100];
		FILE *dos3;

	#ifdef MPI
		sprintf(name6,"./INPUT/flathist1_%d.output",mpi.my_rank);
	#endif
	#ifndef MPI
		sprintf(name6,"./INPUT/flathist1_%d.output",k);
	#endif
		
		if( NULL == (dos3=fopen(name6,"r")) ) {
			sim_dos[k].mod_f		= 1.0;
			for ( int i=0; i<sim_dos[k].l_bins; i++){
				dos_hist[k][i].l_mid		= sim_dos[k].l_begin + sim_dos[k].l_width*(i+0.5);
				dos_hist[k][i].g_of_l		=	0.0;	
				dos_hist[k][i].h_of_l		=	0;
				average[k][i].count			=	1;
				average[k][i].ebond			=	0.0;
				average[k][i].ebend			=	0.0;
#ifndef NEUTRAL
				average[k][i].eurey			=	0.0;
#endif
				average[k][i].etors			=	0.0;
				average[k][i].eimpr			=	0.0;
				average[k][i].elj			=	0.0;
				average[k][i].ecoul			=	0.0;
				average[k][i].epotens		=	0.0;
				average[k][i].d_nc			=	0.0;
				average[k][i].hel			=	0.0;
				average[k][i].con			=	0.0;
				average[k][i].con_2			=	0.0;
				average[k][i].gyr			=	0.0;
				average[k][i].rmsd			=	0.0;
	#ifdef SASA
				average[k][i].esasa			=	0.0;
	#endif//SASA
			}
		}
		else {// if simulation is to be started from old flathist files.
			char name7[100];
			FILE *dos_sim;
			#ifdef MPI
				sprintf(name7,"./INPUT/simul%d.output",mpi.my_rank);
			#endif
			#ifndef MPI
				sprintf(name7,"./INPUT/simul%d.output",k);
			#endif
		
			if( NULL == (dos_sim=fopen(name7,"r")) ) {
				int dummy;
				for ( int i=0; i<sim_dos[k].l_bins; i++){
					fscanf(dos3, "%d	%lf	%lf	%lf	%d\n", &dummy, &sim_dos[k].mod_f, &dos_hist[k][i].l_mid, 
						&dos_hist[k][i].g_of_l, &dos_hist[k][i].h_of_l);
					// dos_hist[k][i].h_of_l=0;
					dos_hist_temp[k][i]	=	dos_hist[k][i];
					if (dummy!= i) {
						fprintf(stdout,"incompatibility in reading flathist1_%d.output file\n",k);
						exit(10);
					}
					/*
					if(dos_hist[k][i].l_mid != (sim_dos[k].l_begin + sim_dos[k].l_width*(i+0.5))){
						fprintf(stdout,"incompatibility in assigning l_mid while reading flathist1.txt\n");
						exit(10);
					}
					*/
					average[k][i].count			=	1;
					average[k][i].ebond			=	0.0;
					average[k][i].ebend			=	0.0;
#ifndef NEUTRAL
				average[k][i].eurey			=	0.0;
#endif
					average[k][i].etors			=	0.0;
					average[k][i].eimpr			=	0.0;
					average[k][i].elj			=	0.0;
					average[k][i].ecoul			=	0.0;
					average[k][i].epotens		=	0.0;
					average[k][i].d_nc			=	0.0;
					average[k][i].hel			=	0.0;
					average[k][i].con			=	0.0;
					average[k][i].con_2			=	0.0;
					average[k][i].gyr			=	0.0;
					average[k][i].rmsd			=	0.0;
					#ifdef SASA
					average[k][i].esasa			=	0.0;
					#endif //SASA
				}
				sim_dos[k].mod_f =log(sim_dos[k].mod_f);
				fclose(dos3);
			}
			else{//old simul_output files to be read too

				int dummy;
				for ( int i=0; i<sim_dos[k].l_bins; i++){
					fscanf(dos3, "%d	%lf	%lf	%lf	%d\n", &dummy, &sim_dos[k].mod_f, &dos_hist[k][i].l_mid, 
						&dos_hist[k][i].g_of_l, &dos_hist[k][i].h_of_l);
					dos_hist_temp[k][i]	=	dos_hist[k][i];
					if (dummy!= i) {
						fprintf(stdout,"incompatibility in reading flathist1_%d.output file\n",k);
						exit(10);
					}
				#ifndef SASA
					#ifdef NEUTRAL
					fscanf(dos_sim, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",&average[k][i].epotensa, 
						&average[k][i].ebonda,&average[k][i].ebenda,&average[k][i].etorsa, &average[k][i].eimpra,
						&average[k][i].elja,&average[k][i].ecoula,&average[k][i].d_nca,&average[k][i].hela,&average[k][i].cona,
						&average[k][i].force_1a,&average[k][i].force_2a,&average[k][i].gyra,&average[k][i].rmsda,&average[k][i].count);
					
					#endif//NEUTRAL
					#ifndef NEUTRAL

					fscanf(dos_sim, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",&average[k][i].epotensa, 
						&average[k][i].ebonda,&average[k][i].ebenda,&average[k][i].eureya,&average[k][i].etorsa, &average[k][i].eimpra,
						&average[k][i].elja,&average[k][i].ecoula,&average[k][i].d_nca,&average[k][i].hela,&average[k][i].cona,
						&average[k][i].force_1a,&average[k][i].force_2a,&average[k][i].gyra,&average[k][i].rmsda,&average[k][i].count);
					
					#endif//NEUTRAL
				#endif//SASA
				#ifdef SASA
					#ifdef NEUTRAL
					fscanf(dos_sim, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",&average[k][i].epotensa,
						&average[k][i].ebonda,&average[k][i].ebenda,&average[k][i].etorsa,&average[k][i].eimpra,
						&average[k][i].elja,&average[k][i].ecoula,&average[k][i].esasaa,&average[k][i].d_nca,&average[k][i].hela,
						&average[k][i].cona,&average[k][i].force_1a,&average[k][i].force_2a,&average[k][i].gyra,
						&average[k][i].rmsda,&average[k][i].count);
					#endif
				#endif

		//			average[k][i].count			=	dos_hist[k][i].h_of_l; // used when reading old format simul.output files
					average[k][i].ebond			=	average[k][i].count*average[k][i].ebonda;
					average[k][i].ebend			=	average[k][i].count*average[k][i].ebenda;
#ifndef NEUTRAL
					average[k][i].eurey			=	average[k][i].count*average[k][i].eureya;
#endif
					average[k][i].etors			=	average[k][i].count*average[k][i].etorsa;
					average[k][i].eimpr			=	average[k][i].count*average[k][i].eimpra;
					average[k][i].elj			=	average[k][i].count*average[k][i].elja;
					average[k][i].ecoul			=	average[k][i].count*average[k][i].ecoula;
					average[k][i].epotens		=	average[k][i].count*average[k][i].epotensa;
					average[k][i].d_nc			=	average[k][i].count*average[k][i].d_nca;
					average[k][i].hel			=	average[k][i].count*average[k][i].hela;
					average[k][i].con			=	average[k][i].count*average[k][i].cona;
					average[k][i].force_1		=	average[k][i].count*average[k][i].force_1a;
					average[k][i].force_2		=	average[k][i].count*average[k][i].force_2a;
					average[k][i].gyr			=	average[k][i].count*average[k][i].gyra;
					average[k][i].rmsd			=	average[k][i].count*average[k][i].rmsda;
					#ifdef SASA
					average[k][i].esasa			=	average[k][i].count*average[k][i].esasaa;
					#endif //SASA
				}
				fclose(dos_sim);
				sim_dos[k].mod_f =log(sim_dos[k].mod_f);
				fclose(dos3);
			}//else if simul.output available
		}//else if flathist available
	}
#endif//XEDOS
#ifdef FX_EDOS
  /* ================================================================== */
  /*                                                                    */
  /* Initialize FX_EDOS variables.                                      */
  /*                                                                    */
  /* ================================================================== */
	for(int k=0; k<sim.NB; k++) {
		sim_dos[k].e_width		=  (int) (sim_dos[k].e_width*sim.kT[k] + 0.5);
		sim_dos[k].e_bins		=  (int) (0.5+ (sim_dos[k].e_end - sim_dos[k].e_begin)/sim_dos[k].e_width);
		if(sim_dos[k].e_bins >= NBINS){
			printf("Not enough memory allocated for number of bins. Increase NBINS to at least %i in defines.h\n", sim_dos[k].e_bins+1);
			exit(9);
		}
		sim_dos[k].dos_acc		= 0;
		sim_dos[k].dos_acc		= 0;
		mc_axial.axial_acc[k]   =0;
		mc_rand.rand_acc[k]     =0;
		mc_trans.trans_acc[k]	=0;
		char name6[100];
		FILE *dos3;

	#ifdef MPI
		sprintf(name6,"./INPUT/hist_%d.output",mpi.my_rank);
	#endif
	#ifndef MPI
		sprintf(name6,"./INPUT/hist_%d.output",k);
	#endif
		
		if( NULL == (dos3=fopen(name6,"r")) ) {
			sim_dos[k].r_value		= 0.0;
			sim_dos[k].sign_walker	= 1;
			sim_dos[k].trips_num	= 0;
			sim_dos[k].trips_old	= 0;
			sim_dos[k].tol_w		= 1.0;

			for ( int i=0; i<sim_dos[k].e_bins; i++){
				dos_hist[k][i].e_mid		= sim_dos[k].e_begin + sim_dos[k].e_width*(i+0.5);
				dos_hist[k][i].g_of_e		=	0.0;	
				dos_hist[k][i].np_of_e		=	0;
				dos_hist[k][i].nm_of_e		=	0;
				dos_hist[k][i].nw_of_e		=	0;
				dos_hist[k][i].f_of_e		=	0.0;
				dos_hist[k][i].w_of_e		=	30.0*(sim_dos[k].e_end- dos_hist[k][i].e_mid)/(sim_dos[k].e_end- sim_dos[k].e_begin);
				dos_hist[k][i].m_of_e		=	0.0;

				average[k][i].count			=	1;
				average[k][i].ebond			=	0.0;
				average[k][i].ebend			=	0.0;
#ifndef NEUTRAL
				average[k][i].eurey			=	0.0;
#endif
				average[k][i].etors			=	0.0;
				average[k][i].eimpr			=	0.0;
				average[k][i].elj			=	0.0;
				average[k][i].ecoul			=	0.0;
				average[k][i].epotens		=	0.0;
				average[k][i].d_nc			=	0.0;
				average[k][i].hel			=	0.0;
				average[k][i].con			=	0.0;
				average[k][i].con_2			=	0.0;
				average[k][i].gyr			=	0.0;
				average[k][i].rmsd			=	0.0;
	#ifdef SASA
				average[k][i].esasa			=	0.0;
	#endif//SASA
				dos_hist_temp[k][i]= dos_hist[k][i];
			dos_hist_old[k][i] = dos_hist[k][i];
			}
			

		}
		else {// if simulation is to be started from old flathist files.

/*-----------------------------------------------------	*/
/* This need to be modified for reading in old files	*/
/* for assigning correct r-value, trips-num etc.		*/
/*-----------------------------------------------------	*/

//			int dummy;
			sim_dos[k].r_value		= 0.0;
			sim_dos[k].sign_walker	= 1;
			sim_dos[k].trips_num	= 0;
			sim_dos[k].trips_old	= 0;
			sim_dos[k].tol_w		= 1.0;
			for ( int i=0; i<sim_dos[k].e_bins; i++){
/*				fscanf(dos3, "%d	%lf	%lf	%d	%d	%d	%lf	%lf	%lf\n", &dummy, &dos_hist[k][i].e_mid,&dos_hist[k][i].g_of_e, 
					&dos_hist[k][i].np_of_e,&dos_hist[k][i].nm_of_e,&dos_hist[k][i].nw_of_e,&dos_hist[k][i].f_of_e,
					&dos_hist[k][i].w_of_e,&dos_hist[k][i].m_of_e);
				dos_hist_temp[k][i]	=	dos_hist[k][i];
				if (dummy!= i) {
					fprintf(stdout,"incompatibility in reading hist_%d.output file\n",k);
					exit(10);
				}
				if(dos_hist[k][i].e_mid != (sim_dos[k].e_begin + sim_dos[k].e_width*(i+0.5))){
					fprintf(stdout,"incompatibility in assigning e_mid while reading hist.txt\n");
					exit(10);
				}
*/
				dos_hist[k][i].e_mid		= sim_dos[k].e_begin + sim_dos[k].e_width*(i+0.5);
				dos_hist[k][i].g_of_e		=	0.0;	
				dos_hist[k][i].np_of_e		=	0;
				dos_hist[k][i].nm_of_e		=	0;
				dos_hist[k][i].nw_of_e		=	0;
				dos_hist[k][i].f_of_e		=	0.0;
				fscanf(dos3, "%lf\n",&dos_hist[k][i].w_of_e);
				dos_hist[k][i].m_of_e		=	0.0;	
			average[k][i].count			=	1;
				average[k][i].ebond			=	0.0;
				average[k][i].ebend			=	0.0;
#ifndef NEUTRAL
				average[k][i].eurey			=	0.0;
#endif
				average[k][i].etors			=	0.0;
				average[k][i].eimpr			=	0.0;
				average[k][i].elj			=	0.0;
				average[k][i].ecoul			=	0.0;
				average[k][i].epotens		=	0.0;
				average[k][i].d_nc			=	0.0;
				average[k][i].hel			=	0.0;
				average[k][i].con			=	0.0;
				average[k][i].con_2			=	0.0;
				average[k][i].gyr			=	0.0;
				average[k][i].rmsd			=	0.0;
	#ifdef SASA
				average[k][i].esasa			=	0.0;
	#endif //SASA
				dos_hist_temp[k][i]= dos_hist[k][i];
				dos_hist_old[k][i] = dos_hist[k][i];
			}
			fclose(dos3);	
		}
	}
#endif//FX_EDOS
#ifdef TDXEDOS
  tdxedos_init();
#endif
#ifdef XWHAM
  xwham_init();
#endif

}

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

void init_block (int ibox)
{
  int k = ibox;
//  for(int k=0; k<sim.NB; k++) {
  /* ================================================================== */
  /*                                                                    */
  /* Reset res accumulators.                                            */
  /*                                                                    */
  /* ================================================================== */
	res[k].count    = 0;
	res[k].tempb    = 0.0;	res[k].tempc    = 0.0;
	res[k].densb    = 0.0;  res[k].densc    = 0.0;
    res[k].ebondb   = 0.0;	res[k].ebondc   = 0.0;
    res[k].ebendb   = 0.0;	res[k].ebendc   = 0.0;
#ifndef NEUTRAL
	res[k].eureyb	= 0.0;	res[k].eureyc   = 0.0;
#endif
    res[k].etorsb   = 0.0;	res[k].etorsc   = 0.0;
	res[k].eimprb   = 0.0;  res[k].eimprc   = 0.0;
    res[k].enbondb  = 0.0;	res[k].enbondc  = 0.0;
   	res[k].enbondsb = 0.0;	res[k].enbondsc = 0.0;
    res[k].etpotb   = 0.0;	res[k].etpotc   = 0.0;
    res[k].etpotsb  = 0.0;	res[k].etpotsc  = 0.0;
    res[k].etkinb   = 0.0;	res[k].etkinc   = 0.0;
    res[k].etotalb  = 0.0;	res[k].etotalc  = 0.0;
    res[k].etotalsb = 0.0;	res[k].etotalsc = 0.0;

  /* ----------------------------------------- */
  /* If the Nose-Hoover thermostat is being    */
  /* used, resets its values.                  */
  /* ----------------------------------------- */
	if(sim.ID == 6 || ((sim.ID2 == 4 || sim.ID2 == 5 || sim.ID2 == 8) && sim.ID ==1))
	{
		res[k].Hb   = 0.0;	res[k].Hc       = 0.0;
	} 

  /* ----------------------------------------- */
  /* If PRESSURE is defined, reset the         */
  /* pressure values.                          */
  /* ----------------------------------------- */
#ifdef PRESSURE    
	res[k].pressb   = 0.0; res[k].pressc    = 0.0;
	for(int i=0; i<6; i++) {
		res[k].cpressa[i] = 0.0;
		res[k].cpressb[i] = 0.0;
		res[k].cpressc[i] = 0.0;
    }
#endif
	ordparam[k].d_ncb	 = 0.0;
	ordparam[k].helb	 = 0.0;
	ordparam[k].conb     = 0.0;
	ordparam[k].gyrb	 = 0.0;
	ordparam[k].rmsdb	 = 0.0;
	for(int i=0; i<torsN[k]; i++) {
		tors[k][i].phia  = tors[k][i].phia/tors[k][i].count;
		tors[k][i].psia  = tors[k][i].psia/tors[k][i].count;
		tors[k][i].thetaa = tors[k][i].thetaa/tors[k][i].count;
		tors[k][i].count = 1;
	}

#ifdef DOS
	  for(int i=0; i<sim_dos[k].e_bins; i++){
		average[k][i].count		= 1;
		average[k][i].ebond		= 0.0;
		average[k][i].ebend		= 0.0;
#ifndef NEUTRAL
		average[k][i].eurey		= 0.0;
#endif
		average[k][i].etors		= 0.0;
		average[k][i].eimpr		= 0.0;
		average[k][i].elj		= 0.0;
		average[k][i].ecoul		= 0.0;
		average[k][i].epotens	= 0.0;
		average[k][i].d_nc		= 0.0;
		average[k][i].hel		= 0.0;
		average[k][i].con		= 0.0;
		average[k][i].con_2			=	0.0;
		average[k][i].gyr		= 0.0;
		average[k][i].rmsd		= 0.0;
		#ifdef SASA
		average[k][i].esasa		= 0.0;
		#endif
	  }
#endif //DOS
#ifdef MMDOS
	  for(int i=0; i<sim_dos[k].e_bins; i++){
		average[k][i].count		= 1;
		average[k][i].ebond		= 0.0;
		average[k][i].ebend		= 0.0;
#ifndef NEUTRAL
		average[k][i].eurey		= 0.0;
#endif
		average[k][i].etors		= 0.0;
		average[k][i].eimpr		= 0.0;
		average[k][i].elj		= 0.0;
		average[k][i].ecoul		= 0.0;
		average[k][i].epotens	= 0.0;
		average[k][i].d_nc		= 0.0;
		average[k][i].hel		= 0.0;
		average[k][i].con		= 0.0;
		average[k][i].con_2			=	0.0;
		#ifdef SASA
		average[k][i].esasa		= 0.0;
		#endif
	  }
#endif //MMDOS
#ifdef CTDOS
	  for(int i=0; i<sim_dos[k].e_bins; i++){
		average[k][i].count		= 1;
		average[k][i].ebond		= 0.0;
		average[k][i].ebend		= 0.0;
#ifndef NEUTRAL
		average[k][i].eurey		= 0.0;
#endif
		average[k][i].etors		= 0.0;
		average[k][i].eimpr		= 0.0;
		average[k][i].elj		= 0.0;
		average[k][i].ecoul		= 0.0;
		average[k][i].epotens	= 0.0;
		average[k][i].d_nc		= 0.0;
		average[k][i].hel		= 0.0;
		average[k][i].con		= 0.0;
		average[k][i].con_2		= 0.0;
		average[k][i].gyr		= 0.0;
		average[k][i].rmsd		= 0.0;
		#ifdef SASA
		average[k][i].esasa		= 0.0;
		#endif
	  }
#endif //CTDOS  
#ifdef XEDOS
	  for(int i=0; i<sim_dos[k].l_bins; i++){
		average[k][i].count		= 1;
		average[k][i].ebond		= average[k][i].ebonda;
		average[k][i].ebend		= average[k][i].ebenda;
#ifndef NEUTRAL
		average[k][i].eurey		= average[k][i].eureya;
#endif
		average[k][i].etors		= average[k][i].etorsa;
		average[k][i].eimpr		= average[k][i].eimpra;
		average[k][i].elj		= average[k][i].elja;
		average[k][i].ecoul		= average[k][i].ecoula;
		average[k][i].epotens	= average[k][i].epotensa;
		average[k][i].d_nc		= average[k][i].d_nca;
		average[k][i].hel		= average[k][i].hela;
		average[k][i].con		= average[k][i].cona;
		average[k][i].con_2		= average[k][i].con_2a;
		average[k][i].force_1	= average[k][i].force_1a;
		average[k][i].force_2	= average[k][i].force_2a;
		average[k][i].gyr		= average[k][i].gyra;
		average[k][i].rmsd		= average[k][i].rmsda;
		#ifdef SASA
		average[k][i].esasa		= average[k][i].esasaa;
		#endif
	  }
#endif //XEDOS
#ifdef FX_EDOS
	  for(int i=0; i<sim_dos[k].e_bins; i++){
		average[k][i].count		= 1;
		average[k][i].ebond		= average[k][i].ebonda;
		average[k][i].ebend		= average[k][i].ebenda;
#ifndef NEUTRAL
		average[k][i].eurey		= average[k][i].eureya;
#endif
		average[k][i].etors		= average[k][i].etorsa;
		average[k][i].eimpr		= average[k][i].eimpra;
		average[k][i].elj		= average[k][i].elja;
		average[k][i].ecoul		= average[k][i].ecoula;
		average[k][i].epotens	= average[k][i].epotensa;
		average[k][i].d_nc		= average[k][i].d_nca;
		average[k][i].hel		= average[k][i].hela;
		average[k][i].con		= average[k][i].cona;
		average[k][i].con_2		= average[k][i].con_2a;
		average[k][i].gyr		= average[k][i].gyra;
		average[k][i].rmsd		= average[k][i].rmsd;
		#ifdef SASA
		average[k][i].esasa		= average[k][i].esasaa;
		#endif
	  }
#endif //FX_EDOS

//  }//end k loop
}
