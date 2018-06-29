/* ======================================================================== */
/* output.cpp                                                               */
/*                                                                          */
/*		This subroutine writes the instantaneous property values to the     */
/* output file periodically during the simulation.  The frequency is        */
/* is controled by sim.blockc which is specified in simul.input.  The       */
/* subroutine is called in main.  Two files are written to in this          */
/* subroutine ener_box#.output, and simul#.output.                          */
/*                                                                          */
/*                                                                          */
/* Passed Parameters:                                                       */
/*			icycle:		The current iteration number            */
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
#ifdef SMD
void init_smd(int);
#endif
#ifdef NSMD
void init_smd(int);
#endif
double tors_av(int,int);
void output (unsigned long icycle, int flag)
{
  /* ================================================================== */
  /* Calculate the current simulation time.                             */
  /* ================================================================== */
  double tt;
  if (flag==0) tt = sim.dt * 1.0*icycle;
  else if (flag==1 && (sim.ID2==3 || sim.ID2==5)) tt= sim.dt*sim.cyc_eq+ sim.dtlong*1.0*icycle; // for multiple time step
  else tt = sim.dt * 1.0*(sim.cyc_eq+icycle);
 
  /* ================================================================== */
  /* Write the current values to simul#.output.                         */
  /* ================================================================== */
	for(int k=0; k<sim.NB; k++) {
		double theta_av	=	tors_av(k,0);
		double psi_av	=	tors_av(k,1);
		double phi_av	=	tors_av(k,2);
  

		/* ================================================================== */
		/* Write the energy values to order#.output.                          */
		/* ================================================================== */
		#ifdef CONFIGT
		if(flag!=0){
			ord_ptr[k] += sprintf(ord_ptr[k],"%lu	%lf	%lf	%lf	",icycle, ordparam[k].hel,
				ordparam[k].d_nc,ordparam[k].con);
			for(int i=0; i<contacts[k].n; i++)	ord_ptr[k] += sprintf(ord_ptr[k],"%d",contacts[k].code[i]);
			ord_ptr[k] += sprintf(ord_ptr[k],"	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",config[k].T, 
				box[k].temp,ordparam[k].force_1,ordparam[k].force_2,ordparam[k].gyr,ordparam[k].rmsd,
				theta_av*180.0/PI,psi_av*180.0/PI,phi_av*180.0/PI,ordparam[k].x1,ordparam[k].x2,ordparam[k].con_2);
		 }
		#endif
		#ifndef CONFIGT
		if(flag!=0){
			ord_ptr[k] += sprintf(ord_ptr[k],"%lu	%lf	%lf	%lf	",icycle, ordparam[k].hel,
				ordparam[k].d_nc,ordparam[k].con);
			for(int i=0; i<contacts[k].n; i++)	ord_ptr[k] += sprintf(ord_ptr[k],"%d",contacts[k].code[i]);
			ord_ptr[k] += sprintf(ord_ptr[k],"	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",box[k].temp,ordparam[k].force_1,
				ordparam[k].force_2,ordparam[k].gyr,ordparam[k].rmsd,theta_av*180.0/PI,psi_av*180.0/PI,phi_av*180.0/PI,
				ordparam[k].x1,ordparam[k].x2,ordparam[k].con_2);
		}

		/* ================================================================== */
		/* Write the energy values to smd#.output.                            */
		/* ================================================================== */
		#ifdef SMD
		init_smd(k);
		smd_ptr[k] += sprintf(smd_ptr[k],"%lu	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",icycle,tt,steermd[k].effr1_av,
			ordparam[k].d_nc,steermd[k].effr2_av,steermd[k].r2,steermd[k].dw1,steermd[k].dw2,en[k].H);
		#endif
		#ifdef NSMD
		init_smd(k);
		smd_ptr[k] += sprintf(smd_ptr[k],"%lu	%lf	%lf	%lf	%lf	\n",icycle,tt,steermd[k].effr1_av,steermd[k].f2_av,en[k].H);
		#endif
		#endif

		/* ================================================================== */
		/* Write the energy values to ener_box#.output.                       */
		/* ================================================================== */
		
		/* ----------------------------------------------- */
		/* If a Nose-Hoover thermostat is being used, the  */
		/* extended Hamiltonian is also written.NVE/NVT MD */
		/* ----------------------------------------------- */    
		if(sim.ID == 6 || ((sim.ID2 == 4 || sim.ID2 == 5 || sim.ID2==8) && sim.ID ==1))
		{
			#ifdef NEUTRAL
			#ifndef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf %lf		 \n", icycle,tt,en[k].bond,en[k].bend,en[k].tors,
				en[k].impr,en[k].nbonds-en[k].coulomb,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals, en[k].H, box[k].press);
			#endif
			#ifdef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf %lf		 \n", icycle,tt,en[k].bond,en[k].bend,en[k].tors,
				en[k].impr,en[k].nbonds-en[k].coulomb,en[k].coulomb,en[k].potens,en[k].kinet, en[k].esasa,en[k].totals, en[k].H, box[k].press);
			#endif
			#endif//NEUTRAL

			#ifndef NEUTRAL
                        #ifndef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf %lf		 \n", icycle,tt,en[k].bond,en[k].bend,en[k].urey,
				en[k].tors,en[k].impr,en[k].nbonds-en[k].coulomb,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals, en[k].H, box[k].press);
                        #endif
			#ifdef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf %lf %lf			 \n", icycle,tt,en[k].bond,en[k].bend,en[k].urey,
                                en[k].tors,en[k].impr,en[k].nbonds-en[k].coulomb,en[k].coulomb,en[k].potens,en[k].kinet,en[k].esasa,en[k].totals,en[k].H,box[k].press);
			#endif
			#endif//NEUTRAL
		}

		/* ------------------------------------------------ */
		/* Either NVE or NVT MD without Nose Hoover			*/
		/* ------------------------------------------------ */
		else if ((sim.ID ==1) || (sim.ID ==3) || (sim.ID==8)) 
		{
			#ifdef NEUTRAL
			#ifndef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf %lf		 \n", icycle,tt,en[k].bond,en[k].bend,en[k].tors,
			en[k].impr,en[k].nbonds-en[k].coulomb,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals,box[k].press);
			#endif
			#ifdef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf %lf %lf			 \n", icycle,tt,en[k].bond,en[k].bend,en[k].tors,
			en[k].impr,en[k].nbonds-en[k].coulomb,en[k].coulomb,en[k].potens,en[k].kinet,en[k].esasa,en[k].totals,box[k].press);
			#endif
			#endif//NEUTRAL

			#ifndef NEUTRAL
                        #ifndef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf %lf		 \n", icycle,tt,en[k].bond,en[k].bend,en[k].urey,
			en[k].tors,en[k].impr,en[k].nbonds-en[k].coulomb,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals,box[k].press);
                        #endif
			#ifdef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf %lf %lf			 \n", icycle,tt,en[k].bond,en[k].bend,en[k].urey,
                   	en[k].tors,en[k].impr,en[k].nbonds-en[k].coulomb,en[k].coulomb,en[k].potens,en[k].kinet,en[k].esasa,en[k].totals,box[k].press);
			#endif
			#endif//NEUTRAL
		}
		/* ------------------------------------------------ */
		/* Hybrid MD MC without thermostat					*/
		/* ------------------------------------------------ */
		else if (sim.ID ==4) 
		{
			#ifdef NEUTRAL
			#ifndef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%ld	%ld	%ld	%ld	%ld	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf %lf		 \n", icycle,mc_pivot.pivot_acc[k],mc_trans.trans_acc[k],mc_rand.rand_acc[k],mc_hmc_acc[k],sim_hyb.swap_acc[k],en[k].bond,en[k].bend,en[k].tors,
			en[k].impr,en[k].nbonds-en[k].coulomb,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals,box[k].press);
			#endif
			#ifdef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%ld	%ld	%ld	%ld	%ld	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf %lf		 \n", icycle,mc_pivot.pivot_acc[k],mc_trans.trans_acc[k],mc_rand.rand_acc[k],mc_hmc_acc[k],sim_hyb.swap_acc[k],en[k].bond,en[k].bend,en[k].tors,
			en[k].impr,en[k].nbonds-en[k].coulomb,en[k].coulomb,en[k].potens,en[k].kinet,en[k].esasa,en[k].totals,box[k].press);
			#endif
			#endif//NEUTRAL

			#ifndef NEUTRAL
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%ld	%ld	%ld	%ld	%ld	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf %lf		 \n", icycle,mc_pivot.pivot_acc[k],mc_trans.trans_acc[k],mc_rand.rand_acc[k],mc_hmc_acc[k],sim_hyb.swap_acc[k],en[k].bond,en[k].bend,en[k].urey,
			en[k].tors,en[k].impr,en[k].nbonds-en[k].coulomb,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals,box[k].press);
			#endif//NEUTRAL
		}
		/* ------------------------------------------------ */
		/* Replica Exchange MD without thermostat			*/
		/* ------------------------------------------------ */
		else if (sim.ID ==5) 
		{
			#ifdef NEUTRAL
			#ifndef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf %lf		 \n", icycle,sim_hyb.swap_acc[k],
				en[k].bond,en[k].bend,en[k].tors,en[k].impr,en[k].nbonds-en[k].coulomb,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals,box[k].press);
			#endif
			#ifdef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%ld	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf %lf	 \n", icycle,sim_hyb.swap_acc[k],
				en[k].bond,en[k].bend,en[k].tors,en[k].impr,en[k].nbonds-en[k].coulomb,en[k].coulomb,en[k].potens,en[k].kinet,en[k].esasa,en[k].totals,box[k].press);
			#endif
			#endif//NEUTRAL

			#ifndef NEUTRAL
			#ifndef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%ld	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf %lf		\n", icycle,sim_hyb.swap_acc[k],
				en[k].bond,en[k].bend,en[k].urey,en[k].tors,en[k].impr,en[k].nbonds-en[k].coulomb,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals,box[k].press);
			#endif
			#ifdef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf %lf %lf			 \n", icycle,sim_hyb.swap_acc[k],
                                en[k].bond,en[k].bend,en[k].urey,en[k].tors,en[k].impr,en[k].nbonds-en[k].coulomb,en[k].coulomb,en[k].potens,en[k].kinet,en[k].esasa,en[k].totals,box[k].press);
			#endif
			#endif//NEUTRAL
		}

		/* ================================================================== */
		/* Write the contributions to the coumbic energy values calculated by */
		/* the Ewald method to ewald#.output.                                 */
		/* ================================================================== */
		#ifdef EWALD_DEBUG
		FILE *ewald_fp;
		char ewald_name[100];
		#ifdef MPI
		sprintf(ewald_name,"./OUTPUT/BOX%d/ewald%d.output",mpi.my_rank,mpi.my_rank);
		#endif
		#ifndef MPI
		sprintf(ewald_name,"./OUTPUT/BOX%d/ewald%d.output",k,k);
		#endif
	
		ewald_fp= fopen(ewald_name,"a");
		fprintf(ewald_fp, "%lu	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n", icycle,tt,en[k].rewald,en[k].rewalds,en[k].kewald,
			en[k].sewald,en[k].iewald,en[k].tewald);
		fclose(ewald_fp);
		#endif

		#ifdef PR_NPT_DEBUG
		FILE npt_fp;
		char npt_name[100];

		sprintf(npt_name,"./OUTPUT/BOX%d/npt%d.output",k,k);
		npt_fp= fopen(npt_name,"a");
		fprintf(npt_fp, "%lf	%lf	%lf	%lf	%lf %lf	%lf	%lf	%lf\n", tt, en[k].poten,en[k].kinet, ke_thermostat[k], 
		sum_ktzeta_i[k], cnf9ktzeta1[k],pr_pv[k],cka[k],en[k].H);
		fclose(npt_fp);
		#endif

  /* ================================================================== */
  /* Write the angle the molecule makes with the wall in angle#.output. */
  /* ================================================================== */
/*#ifdef WALL_DEBUG
  FILE *wall_ptr;
  char wall_name[100];
  double x = atom[k][wall[k].angle_site].x;
  double y = atom[k][wall[k].angle_site].y;
  double z = atom[k][wall[k].angle_site].z;
  double angle = atan(z/sqrt(x*x + y*y))*180.0/PI;
  sprintf(wall_name,"./OUTPUT/BOX%d/angle%d.output",k,k);
  wall_ptr= fopen(wall_name,"a");
  fprintf(wall_ptr, "%lu	%d	%lf	%lf	%lf	%lf\n",icycle,wall[k].angle_site,x,y,z,angle);
  fclose(wall_ptr);
#endif
*/

  }// k loop ends

}
