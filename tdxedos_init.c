#ifdef TDXEDOS
/* ======================================================================== */
/* tdxedos_init.cpp		                                                    */
/* This suboutine initializes the TDXEDOS variables.  It is analogous to    */
/* the operations done in init_all and is called from that subroutine.      */
/* It is in a separate subroutine to simply help in bookkeeping.            */
/*                                                                          */
/* Written by Thomas Knotts 8 June 04                                       */	
/* ======================================================================== */
#include "defines.h"
void read_swap(void);

void tdxedos_init(void){

    char name[100], name2[100];
	FILE *flathist_ptr, *simul_ptr;

  /* ================================================================== */
  /* Read the swap#.input file.                                         */
  /* ================================================================== */
  read_swap();

  /* ================================================================== */
  /* Initialize XEDOS variables.                                        */
  /* ================================================================== */
  for(int k=0; k<sim.NB; k++) {
	int flag = 0;
    sim_dos[k].x1_bins		= (int)(0.5+(sim_dos[k].x1_end - sim_dos[k].x1_begin)/sim_dos[k].x1_width);
	sim_dos[k].x2_bins		= (int)(0.5+(sim_dos[k].x2_end - sim_dos[k].x2_begin)/sim_dos[k].x2_width);
	if(sim_dos[k].x1_bins >= NBINS1){
		printf("Not enough memory allocated for number of bins.  Increase NBINS1 to at least %i in defines.h\n", sim_dos[k].x1_bins+1);
		flag = 1;
	}else if(sim_dos[k].x2_bins >= NBINS2){
        printf("Not enough memory allocated for number of bins.  Increase NBINS2 to at least %i in defines.h\n", sim_dos[k].x2_bins+1);
        flag = 1;
	}
	if(flag == 1) exit(9);
    sim_dos[k].dos_acc		= 0;
    mc_axial.axial_acc[k]   = 0;
    mc_rand.rand_acc[k]     = 0;
    mc_trans.trans_acc[k]	= 0;

    /* --------------------------------------------------- */
	/* If flathist#.output files are present, read these   */
	/* to initialize, otherwise, start a fresh simulation. */
	/* --------------------------------------------------- */

	#ifdef MPI
		sprintf(name,"./INPUT/flathist1_%d.output",mpi.my_rank);
	#endif
	#ifndef MPI
		sprintf(name,"./INPUT/flathist1_%d.output",k);
	#endif
		
    if(NULL == (flathist_ptr=fopen(name,"r")) ){
      sim_dos[k].mod_f = 1.0;
      for(int i=0; i<sim_dos[k].x1_bins; i++){
		for(int j=0; j<sim_dos[k].x2_bins; j++){
          dos_hist[k][i][j].x1_mid  = sim_dos[k].x1_begin + sim_dos[k].x1_width*(i+0.5);
		  dos_hist[k][i][j].x2_mid  = sim_dos[k].x2_begin + sim_dos[k].x2_width*(j+0.5);
          dos_hist[k][i][j].g_of_l  = 0.0;	
          dos_hist[k][i][j].h_of_l  = 0;
          average[k][i][j].count    = 1;
          average[k][i][j].ebond    = 0.0;
          average[k][i][j].ebend    = 0.0;
          #ifndef NEUTRAL
          average[k][i][j].eurey    = 0.0;
          #endif
          average[k][i][j].etors    = 0.0;
          average[k][i][j].eimpr    = 0.0;
          average[k][i][j].elj      = 0.0;
          average[k][i][j].ecoul    = 0.0;
          average[k][i][j].epotens  = 0.0;
          average[k][i][j].d_nc     = 0.0;
          average[k][i][j].hel      = 0.0;
          average[k][i][j].con      = 0.0;
          average[k][i][j].gyr      = 0.0;
          average[k][i][j].rmsd     = 0.0;
          #ifdef SASA
          average[k][i][j].esasa    = 0.0;
          #endif
		  average[k][i][j].x1       = 0.0;
		  average[k][i][j].x2       = 0.0;
		}
	  }
	}
    else{// if simulation is to be started from old flathist files.
      int dummy1,dummy2;
	  for(int i=0; i<sim_dos[k].x1_bins; i++){
        for(int j=0; j<sim_dos[k].x2_bins; j++){
          fscanf(flathist_ptr, "%d	%lf	%lf	%lf	%lf	%lf	%d\n", &dummy1, &dummy2, &sim_dos[k].mod_f, &dos_hist[k][i][j].x1_mid, 
            &dos_hist[k][i][j].x2_mid,&dos_hist[k][i][j].g_of_l, &dos_hist[k][i][j].h_of_l);
          dos_hist_temp[k][i][j]	=	dos_hist[k][i][j];
          if (dummy1!= i && dummy2 != j){
            fprintf(stdout,"incompatibility in reading flathist1_%d.output file\n",k);
            exit(10);
		  }
		}
	  }        
	  sim_dos[k].mod_f =log(sim_dos[k].mod_f);
      fclose(flathist_ptr);

    /* --------------------------------------------------- */
	/* If simul#.output files are present, read these      */
	/* to initialize, otherwise, start a fresh simulation. */
	/* --------------------------------------------------- */
      #ifdef MPI
      sprintf(name2,"./INPUT/simul%d.output",mpi.my_rank);
      #endif
      #ifndef MPI
      sprintf(name2,"./INPUT/simul%d.output",k);
      #endif
		
      if(NULL == (simul_ptr=fopen(name2,"r")) ) {
        for(int i=0; i<sim_dos[k].x1_bins; i++){
		  for(int j=0; j<sim_dos[k].x2_bins; j++){
            average[k][i][j].count    = 1;
            average[k][i][j].ebond    = 0.0;
            average[k][i][j].ebend    = 0.0;
            #ifndef NEUTRAL
            average[k][i][j].eurey    = 0.0;
            #endif
            average[k][i][j].etors    = 0.0;
            average[k][i][j].eimpr    = 0.0;
            average[k][i][j].elj      = 0.0;
            average[k][i][j].ecoul    = 0.0;
            average[k][i][j].epotens  = 0.0;
            average[k][i][j].d_nc     = 0.0;
            average[k][i][j].hel      = 0.0;
            average[k][i][j].con      = 0.0;
            average[k][i][j].gyr      = 0.0;
            average[k][i][j].rmsd     = 0.0;
            #ifdef SASA
            average[k][i][j].esasa    = 0.0;
            #endif
		    average[k][i][j].x1       = 0.0;
		    average[k][i][j].x2       = 0.0;
		  }
		}		
	  }
      else{//old simul_output files to be read too
        for(int i=0; i<sim_dos[k].x1_bins; i++){
		  for(int j=0; j<sim_dos[k].x2_bins; j++){
            #ifndef SASA
            #ifdef NEUTRAL
            fscanf(simul_ptr, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",&average[k][i][j].epotensa, 
              &average[k][i][j].ebonda,&average[k][i][j].ebenda,&average[k][i][j].etorsa, &average[k][i][j].eimpra,
              &average[k][i][j].elja,&average[k][i][j].ecoula,&average[k][i][j].d_nca,&average[k][i][j].hela,&average[k][i][j].cona,
              &average[k][i][j].force_1a,&average[k][i][j].force_2a,&average[k][i][j].gyra,&average[k][i][j].rmsda,%average[k][i][j].x1a,
			  &average[k][i][j].x2a,&average[k][i][j].count);
            #else
            fscanf(simul_ptr, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",&average[k][i][j].epotensa, 
              &average[k][i][j].ebonda,&average[k][i][j].ebenda,&average[k][i][j].eureya,&average[k][i][j].etorsa, &average[k][i][j].eimpra,
              &average[k][i][j].elja,&average[k][i][j].ecoula,&average[k][i][j].d_nca,&average[k][i][j].hela,&average[k][i][j].cona,
              &average[k][i][j].force_1a,&average[k][i][j].force_2a,&average[k][i][j].gyra,&average[k][i][j].rmsda,&average[k][i][j].x1a,
              &average[k][i][j].x2a,&average[k][i][j].count);	
            #endif//NEUTRAL
            #endif//SASA

            #ifdef SASA
            #ifdef NEUTRAL
            fscanf(dos_sim, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",&average[k][i][j].epotensa,
              &average[k][i][j].ebonda,&average[k][i][j].ebenda,&average[k][i][j].etorsa,&average[k][i][j].eimpra,
              &average[k][i][j].elja,&average[k][i][j].ecoula,&average[k][i][j].esasaa,&average[k][i][j].d_nca,&average[k][i][j].hela,
              &average[k][i][j].cona,&average[k][i][j].force_1a,&average[k][i][j].force_2a,&average[k][i][j].gyra,
              &average[k][i][j].rmsda,&average[k][i][j].x1a,&average[k][i][j].x2a,&average[k][i][j].count);
            #endif
            #endif
            average[k][i][j].ebond    = average[k][i][j].count*average[k][i][j].ebonda;
            average[k][i][j].ebend    = average[k][i][j].count*average[k][i][j].ebenda;
            #ifndef NEUTRAL
            average[k][i][j].eurey    = average[k][i][j].count*average[k][i][j].eureya;
            #endif
            average[k][i][j].etors    = average[k][i][j].count*average[k][i][j].etorsa;
            average[k][i][j].eimpr    = average[k][i][j].count*average[k][i][j].eimpra;
            average[k][i][j].elj      = average[k][i][j].count*average[k][i][j].elja;
            average[k][i][j].ecoul    = average[k][i][j].count*average[k][i][j].ecoula;
            average[k][i][j].epotens  = average[k][i][j].count*average[k][i][j].epotensa;
            average[k][i][j].d_nc     = average[k][i][j].count*average[k][i][j].d_nca;
            average[k][i][j].hel      = average[k][i][j].count*average[k][i][j].hela;
            average[k][i][j].con      = average[k][i][j].count*average[k][i][j].cona;
            average[k][i][j].force_1  = average[k][i][j].count*average[k][i][j].force_1a;
            average[k][i][j].force_2  = average[k][i][j].count*average[k][i][j].force_2a;
            average[k][i][j].gyr      = average[k][i][j].count*average[k][i][j].gyra;
            average[k][i][j].rmsd     = average[k][i][j].count*average[k][i][j].rmsda;
            #ifdef SASA
            average[k][i][j].esasa    = average[k][i][j].count*average[k][i][j].esasaa;
           #endif //SASA
			average[k][i][j].x1       = average[k][i][j].count*average[k][i][j].x1a;
 			average[k][i][j].x2       = average[k][i][j].count*average[k][i][j].x2a;
 		  }//j
		}//i 
	  }//else if simul#.output
      fclose(simul_ptr);
	}//if flathist#.output
  }//k


}
#endif
