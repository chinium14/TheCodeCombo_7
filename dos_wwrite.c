/* ======================================================================== */
/* dos_write.cpp                                                           */
/*                                                                          */
/*		This subroutine writes the instantaneous property values to the     */
/* output file periodically during the simulation.  The frequency is        */
/* is controled by sim.blockc which is specified in simul.input.  The       */
/* subroutine is called in main.  Two files are written to in this          */
/* subroutine ener_box#.output, and simul#.output.                          */
/*                                                                          */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						icycle:		The current iteration number            */
/*                      flag:       Flag used to denote if the subroutine   */
/*                                  is being called during equilibration(0) */
/*                                  or production(1).                       */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void dos_wwrite (unsigned long icycle, int flag)
{

  /* ================================================================== */
  /*                                                                    */
  /* Declare the file pointers needed for the output file.  "sout" is   */
  /* the file pointer for simul#.output, and sout1 is the file pointer  */
  /* for ener_box#.output.                                              */
  /*                                                                    */
  /* ================================================================== */
  FILE *sout;

  /* ================================================================== */
  /*                                                                    */
  /* Calculate the current simulation time.                             */
  /*                                                                    */
  /* ================================================================== */
  double tt;
  if (flag==0) tt = sim.dt * 1.0*icycle;
  else if (flag==1 && (sim.ID2==3 || sim.ID2==5)) tt= sim.dt*sim.cyc_eq+ sim.dtlong*1.0*icycle; // for multiple time step
  else tt = sim.dt * 1.0*(sim.cyc_eq+icycle);

  /* ================================================================== */
  /*                                                                    */
  /* Write the current values to simul#.output.                         */
  /*                                                                    */
  /* ================================================================== */
  for(int k=0; k<sim.NB; k++) {
#ifdef MPI
	 sprintf(name, "./OUTPUT/DOS/simul%d.output",mpi.my_rank);
#endif
#ifndef MPI
	 sprintf(name, "./OUTPUT/DOS/simul%d.output",k);
#endif
	sout = fopen(name,"w");
#ifdef NEUTRAL
#ifdef DOS
	for(int i=0; i<sim_dos[k].e_bins; i++){
		#ifndef SASA
		fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",average[k][i].epotensa, average[k][i].ebonda,average[k][i].ebenda,
			average[k][i].etorsa, average[k][i].eimpra,average[k][i].elja,average[k][i].ecoula,average[k][i].d_nca,
			average[k][i].hela,average[k][i].cona,average[k][i].gyra,average[k][i].rmsda,average[k][i].con_2a,average[k][i].count);
		#endif
		#ifdef SASA
		fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",average[k][i].epotensa, average[k][i].ebonda,average[k][i].ebenda,
			 average[k][i].etorsa, average[k][i].eimpra,average[k][i].elja,average[k][i].ecoula,average[k][i].esasaa,
			average[k][i].d_nca,average[k][i].hela,average[k][i].cona,average[k][i].gyra,average[k][i].rmsda,average[k][i].con_2a,average[k][i].count);
		#endif
	}
#endif // DOS
#ifdef MMDOS
	for(int i=0; i<sim_dos[k].e_bins; i++){
		#ifndef SASA
		fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",average[k][i].epotensa, average[k][i].ebonda,average[k][i].ebenda,
			average[k][i].etorsa, average[k][i].eimpra,average[k][i].elja,average[k][i].ecoula,average[k][i].d_nca,
			average[k][i].hela,average[k][i].cona,average[k][i].rmsda);
		#endif
		#ifdef SASA
		fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",average[k][i].epotensa, average[k][i].ebonda,average[k][i].ebenda,
			average[k][i].etorsa, average[k][i].eimpra,average[k][i].elja,average[k][i].ecoula,average[k][i].esasaa,
			average[k][i].d_nca,average[k][i].hela,average[k][i].cona,average[k][i].rmsda);
		#endif
	}
#endif // MMDOS
#ifdef CTDOS
	for(int i=0; i<sim_dos[k].e_bins; i++){
		#ifndef SASA
		fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",average[k][i].epotensa, average[k][i].ebonda,average[k][i].ebenda,
			average[k][i].etorsa, average[k][i].eimpra,average[k][i].elja,average[k][i].ecoula,average[k][i].d_nca,
			average[k][i].hela,dos_hist[k][i].config_T,average[k][i].cona,average[k][i].gyra,average[k][i].rmsda,average[k][i].con_2a,average[k][i].count);
		#endif
		#ifdef SASA
		fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",average[k][i].epotensa, average[k][i].ebonda,average[k][i].ebenda,
			average[k][i].etorsa, average[k][i].eimpra,average[k][i].elja,average[k][i].ecoula,average[k][i].esasaa,
			average[k][i].d_nca,average[k][i].hela,dos_hist[k][i].config_T,average[k][i].cona,average[k][i].gyra,average[k][i].rmsda,
			average[k][i].con_2a,average[k][i].count);
		#endif
	}
#endif // CTDOS
#ifdef XEDOS
	for(int i=0; i<sim_dos[k].l_bins; i++){
		#ifndef SASA
		fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",average[k][i].epotensa, average[k][i].ebonda,
			average[k][i].ebenda,average[k][i].etorsa, average[k][i].eimpra,average[k][i].elja,average[k][i].ecoula,
			average[k][i].d_nca,average[k][i].hela,average[k][i].cona,average[k][i].force_1a,average[k][i].force_2a,average[k][i].gyra
			,average[k][i].rmsda,average[k][i].con_2a,average[k][i].count);
		#endif
		#ifdef SASA
		fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",average[k][i].epotensa, average[k][i].ebonda,
			average[k][i].ebenda,average[k][i].etorsa, average[k][i].eimpra,average[k][i].elja,average[k][i].ecoula,
			average[k][i].esasaa,average[k][i].d_nca,average[k][i].hela,average[k][i].cona,average[k][i].force_1a,average[k][i].force_2a,
			average[k][i].gyra,average[k][i].rmsda,average[k][i].con_2a,average[k][i].count);
		#endif
	}
#endif // XEDOS
#ifdef FX_EDOS
	for(int i=0; i<sim_dos[k].e_bins; i++){
		#ifndef SASA
		fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",average[k][i].epotensa, average[k][i].ebonda,average[k][i].ebenda,
			average[k][i].etorsa, average[k][i].eimpra,average[k][i].elja,average[k][i].ecoula,average[k][i].d_nca,
			average[k][i].hela,average[k][i].cona,average[k][i].gyra,average[k][i].rmsda,average[k][i].count);
		#endif
		#ifdef SASA
		fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",average[k][i].epotensa, average[k][i].ebonda,average[k][i].ebenda,
			 average[k][i].etorsa, average[k][i].eimpra,average[k][i].elja,average[k][i].ecoula,average[k][i].esasaa,
			average[k][i].d_nca,average[k][i].hela,average[k][i].cona,average[k][i].gyra,average[k][i].rmsda,average[k][i].count);
		#endif
	}
#endif // FX_EDOS
#ifdef TDXEDOS
	for(int i=0; i<sim_dos[k].x1_bins; i++){
		for(int j=0; j<sim_dos[k].x2_bins; j++){
			#ifndef SASA
			fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		%lf %lf	%d\n",average[k][i][j].epotensa, average[k][i][j].ebonda,
			average[k][i][j].ebenda,average[k][i][j].etorsa, average[k][i][j].eimpra,average[k][i][j].elja,average[k][i][j].ecoula,average[k][i][j].d_nca,
			average[k][i][j].hela,average[k][i][j].cona,average[k][i][j].gyra,average[k][i][j].rmsda,average[k][i][j].x1a,
			average[k][i][j].x2a,average[k][i][j].count);
			#endif
			#ifdef SASA
			fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",average[k][i][j].epotensa, average[k][i][j].ebonda,
			average[k][i][j].ebenda,average[k][i][j].etorsa, average[k][i][j].eimpra,average[k][i][j].elja,average[k][i][j].ecoula,average[k][i][j].esasaa,
			average[k][i][j].d_nca,average[k][i][j].hela,average[k][i][j].cona,average[k][i][j].gyra,average[k][i][j].rmsda,
			average[k][i][j].x1a,average[k][i][j].x2a,average[k][i][j].count);
			#endif
		}
	}
#endif //TDXEDOS
#endif//  NEUTRAL
#ifndef NEUTRAL
#ifdef DOS
	for(int i=0; i<sim_dos[k].e_bins; i++){
		fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",average[k][i].epotensa, average[k][i].ebonda,average[k][i].ebenda,
			average[k][i].eureya, average[k][i].etorsa, average[k][i].eimpra,average[k][i].elja,average[k][i].ecoula,average[k][i].d_nca,
			average[k][i].hela,average[k][i].cona,average[k][i].gyra,average[k][i].rmsda,average[k][i].con_2a,average[k][i].count);
	}
#endif // DOS
#ifdef MMDOS
	for(int i=0; i<sim_dos[k].e_bins; i++){
		fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",average[k][i].epotensa, average[k][i].ebonda,average[k][i].ebenda,
			average[k][i].eureya,average[k][i].etorsa, average[k][i].eimpra,average[k][i].elja,average[k][i].ecoula,average[k][i].d_nca,
			average[k][i].hela,average[k][i].cona,average[k][i].rmsda);
	}
#endif // MMDOS
#ifdef CTDOS
	for(int i=0; i<sim_dos[k].e_bins; i++){
		fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",average[k][i].epotensa, average[k][i].ebonda,average[k][i].ebenda,
			average[k][i].eureya,average[k][i].etorsa, average[k][i].eimpra,average[k][i].elja,average[k][i].ecoula,average[k][i].d_nca,
			average[k][i].hela,dos_hist[k][i].config_T,average[k][i].cona,average[k][i].gyra,average[k][i].rmsda,average[k][i].con_2a,average[k][i].count);
	}
#endif // CTDOS
#ifdef XEDOS
	for(int i=0; i<sim_dos[k].l_bins; i++){
		fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",average[k][i].epotensa, average[k][i].ebonda,
			average[k][i].ebenda,average[k][i].eureya,average[k][i].etorsa, average[k][i].eimpra,average[k][i].elja,average[k][i].ecoula,
			average[k][i].d_nca,average[k][i].hela,average[k][i].cona,average[k][i].force_1a,average[k][i].force_2a,average[k][i].gyra,
			average[k][i].rmsda,average[k][i].con_2a,average[k][i].count);
	}
#endif // XEDOS
#ifdef FX_EDOS
	for(int i=0; i<sim_dos[k].e_bins; i++){
		fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",average[k][i].epotensa, average[k][i].ebonda,average[k][i].ebenda,
			average[k][i].eureya, average[k][i].etorsa, average[k][i].eimpra,average[k][i].elja,average[k][i].ecoula,average[k][i].d_nca,
			average[k][i].hela,average[k][i].cona,average[k][i].gyra,average[k][i].rmsda,average[k][i].count);	}
#endif // FX_EDOS
#ifdef TDXEDOS
	for(int i=0; i<sim_dos[k].x1_bins; i++){
		for(int j=0; j<sim_dos[k].x2_bins; j++){
			fprintf(sout, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%d\n",average[k][i][j].epotensa, average[k][i][j].ebonda,
				average[k][i][j].ebenda,average[k][i][j].eureya,average[k][i][j].etorsa, average[k][i][j].eimpra,average[k][i][j].elja,average[k][i][j].ecoula,
				average[k][i][j].d_nca,average[k][i][j].hela,average[k][i][j].cona,average[k][i][j].force_1a,average[k][i][j].force_2a,average[k][i][j].gyra,
				average[k][i][j].rmsda,average[k][i][j].x1a,average[k][i][j].x2a,average[k][i][j].count);
		}
	}
#endif // TDXEDOS
#endif// ndef NEUTRAL
	fclose(sout);
  }// k loop ends

  /* ------------------------------------------------ */
  /* Density of States MD MC without thermostat		  */
  /* Print out the current histogram				  */	
  /* ------------------------------------------------ */
#ifdef DOS
  char name3[100];
  FILE * fptr_hist1;
  for(int k=0; k<sim.NB; k++) {	
	#ifdef MPI
	sprintf(name3, "./OUTPUT/DOS/flathist1_%d.output",mpi.my_rank); 
	#endif
	#ifndef MPI
    sprintf(name3, "./OUTPUT/DOS/flathist1_%d.output",k);
	#endif
	fptr_hist1 = fopen(name3,"w");
	for (int i =0; i<sim_dos[k].e_bins; i++){
			fprintf(fptr_hist1, "%d	%6.12lf	%lf	%lf	%d\n", i, exp(sim_dos[k].mod_f), dos_hist[k][i].e_mid, dos_hist[k][i].g_of_e, dos_hist[k][i].h_of_e);
	}
	fclose(fptr_hist1);
  }
#endif //DOS

#ifdef MMDOS
  char name3[100];
  FILE * fptr_hist1;
  for(int k=0; k<sim.NB; k++) {
	#ifdef MPI
	sprintf(name3, "./OUTPUT/DOS/flathist1_%d.output",mpi.my_rank); 
	#endif
	#ifndef MPI
	 sprintf(name3, "./OUTPUT/DOS/flathist1_%d.output",k);
	#endif
	fptr_hist1 = fopen(name3,"w");
	for (int i =0; i<sim_dos[k].e_bins; i++){
			fprintf(fptr_hist1, "%d	%6.12lf	%lf	%lf	%d	%lf\n", i, sim_dos[k].mod_f, dos_hist[k][i].e_mid, dos_hist[k][i].g_of_e, dos_hist[k][i].h_of_e, (dos_hist[k][i].k_of_e/dos_hist[k][i].h_of_e));
	}
	fclose(fptr_hist1);
  }
#endif //MMDOS
#ifdef CTDOS
  char name3[100];
  FILE * fptr_hist1;
  for(int k=0; k<sim.NB; k++) {	
	#ifdef MPI
	sprintf(name3, "./OUTPUT/DOS/flathist1_%d.output",mpi.my_rank); 
	#endif
	#ifndef MPI
    sprintf(name3, "./OUTPUT/DOS/flathist1_%d.output",k);
	#endif
	fptr_hist1 = fopen(name3,"w");
	for (int i =0; i<sim_dos[k].e_bins; i++){
			fprintf(fptr_hist1, "%d	%6.12lf	%lf	%lf	%d	%lf\n", i, exp(sim_dos[k].mod_f), dos_hist[k][i].e_mid, dos_hist[k][i].g_of_e, dos_hist[k][i].h_of_e,dos_hist[k][i].gT_of_e);
	}
	fclose(fptr_hist1);
  }
#endif//CTDOS
#ifdef XEDOS
  char name3[100];
  FILE * fptr_hist1;
  for(int k=0; k<sim.NB; k++) {	
	#ifdef MPI
	sprintf(name3, "./OUTPUT/DOS/flathist1_%d.output",mpi.my_rank); 
	#endif
	#ifndef MPI
    sprintf(name3, "./OUTPUT/DOS/flathist1_%d.output",k);
	#endif
	fptr_hist1 = fopen(name3,"w");
	for (int i =0; i<sim_dos[k].l_bins; i++){
			fprintf(fptr_hist1, "%d	%6.12lf	%lf	%lf	%d\n", i, exp(sim_dos[k].mod_f), dos_hist[k][i].l_mid, dos_hist[k][i].g_of_l, dos_hist[k][i].h_of_l);
	}
	fclose(fptr_hist1);
  }
#endif //XEDOS
#ifdef FX_EDOS
  char name3[100];
  FILE * fptr_hist1;
  for(int k=0; k<sim.NB; k++) {	
	#ifdef MPI
	sprintf(name3, "./OUTPUT/DOS/hist_%d.output",mpi.my_rank); 
	#endif
	#ifndef MPI
    sprintf(name3, "./OUTPUT/DOS/hist_%d.output",k);
	#endif
	fptr_hist1 = fopen(name3,"w");
	for (int i =0; i<sim_dos[k].e_bins; i++){
			fprintf(fptr_hist1, "%d	%lf	%d	%d	%d	%f	%f	%f	%f	%d	%d	%d	%f\n", i, dos_hist[k][i].e_mid, 
			dos_hist[k][i].np_of_e,	dos_hist[k][i].nm_of_e,dos_hist[k][i].nw_of_e,dos_hist[k][i].f_of_e,
			dos_hist[k][i].m_of_e,dos_hist[k][i].w_of_e,dos_hist[k][i].g_of_e,
			sim_dos[k].trips_num-sim_dos[k].trips_old,sim_dos[k].trips_num,n_iter, sim_dos[k].r_value);
	}
	fclose(fptr_hist1);
  }
#endif //FX_EDOS
#ifdef TDXEDOS
  char name3[100];
  FILE * fptr_hist1;
  for(int k=0; k<sim.NB; k++) {	
	#ifdef MPI
	sprintf(name3, "./OUTPUT/DOS/flathist1_%d.output",mpi.my_rank); 
	#endif
	#ifndef MPI
    sprintf(name3, "./OUTPUT/DOS/flathist1_%d.output",k);
	#endif
	fptr_hist1 = fopen(name3,"w");
	for(int i=0; i<sim_dos[k].x1_bins; i++){
	  for(int j=0; j<sim_dos[k].x2_bins; j++)
        fprintf(fptr_hist1, "%d	%d	%6.12lf	%lf	%lf	%lf	%d\n", i,j, exp(sim_dos[k].mod_f), dos_hist[k][i][j].x1_mid, 
		  dos_hist[k][i][j].x2_mid*180.0/PI, dos_hist[k][i][j].g_of_l, dos_hist[k][i][j].h_of_l);
	}
	fclose(fptr_hist1);
  }
#endif //TDXEDOS

}
