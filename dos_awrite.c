/* ======================================================================== */
/* dos_output.cpp                                                           */
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
double tors_av(int,int);
void dos_awrite (unsigned long icycle, int flag)
{
 
	for(int k=0; k<sim.NB; k++) {

		double theta_av	=	tors_av(k,0);
		double psi_av	=	tors_av(k,1);
		double phi_av	=	tors_av(k,2);
  

		/* ================================================================== */
		/* Write the energy values to ener_box#.output.                       */
		/* ================================================================== */

		#ifdef CONFIGT
		ord_ptr[k] += sprintf(ord_ptr[k],"%lu	%lf	%lf	%lf	",icycle, ordparam[k].hel,ordparam[k].d_nc,
			ordparam[k].con);
		for(int i=0; i<contacts[k].n; i++)	ord_ptr[k] += sprintf(ord_ptr[k],"%d",contacts[k].code[i]);
		ord_ptr[k] += sprintf(ord_ptr[k],"	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",config[k].T, box[k].temp,ordparam[k].force_1,
			ordparam[k].force_2,ordparam[k].gyr,ordparam[k].rmsd,theta_av*180.0/PI,psi_av*180.0/PI,phi_av*180.0/PI,ordparam[k].x1, 
			ordparam[k].x2,ordparam[k].con_2);
		#endif
		#ifndef CONFIGT
		ord_ptr[k] += sprintf(ord_ptr[k],"%lu	%lf	%lf	%lf	",icycle, ordparam[k].hel,ordparam[k].d_nc,	ordparam[k].con);
		for(int i=0; i<contacts[k].n; i++)	ord_ptr[k] += sprintf(ord_ptr[k],"%d",contacts[k].code[i]);
		ord_ptr[k] += sprintf(ord_ptr[k],"	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",ordparam[k].force_1,ordparam[k].force_2,ordparam[k].gyr,
			ordparam[k].rmsd,theta_av*180.0/PI,psi_av*180.0/PI,phi_av*180.0/PI,ordparam[k].x1, ordparam[k].x2,ordparam[k].con_2);
		#endif


		/* ================================================================== */
		/* Write the energy values to ener_box#.output.                       */
		/* ================================================================== */
	
		/* ------------------------------------------------ */
		/* Density of States MD MC without thermostat		*/
		/* ------------------------------------------------ */
		#ifdef NEUTRAL

		#ifdef DOS 
		if (sim.ID ==2) 
		{
			#ifndef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%ld	%ld	%ld	%ld	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,sim_dos[k].dos_acc,mc_pivot.pivot_acc[k],
				mc_trans.trans_acc[k],mc_rand.rand_acc[k],,en[k].bond,en[k].bend,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals);
			#else
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%d	%d	%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,sim_dos[k].dos_acc,mc_pivot.pivot_acc[k],
				mc_trans.trans_acc[k],en[k].bond,en[k].bend,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].esasa,en[k].totals);
			#endif
		}
		#endif //DOS
	
		#ifdef MMDOS 
		else if (sim.ID ==7) 
		{
			#ifndef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%d	%d	%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,sim_dos[k].dos_acc,mc_pivot.pivot_acc[k],
				sim_dos[k].dke_acc,en[k].bond,en[k].bend,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals);
			#else
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%d	%d	%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,sim_dos[k].dos_acc,mc_pivot.pivot_acc[k],
				sim_dos[k].dke_acc,en[k].bond,en[k].bend,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].esasa,en[k].totals);
			#endif
		}
		#endif //MMDOS

		#ifdef CTDOS 
		if (sim.ID ==9) 
		{
			#ifndef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%d		%d	%d	%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,sim_dos[k].dos_acc,mc_pivot.pivot_acc[k],
				mc_trans.trans_acc[k],mc_rand.rand_acc[k],en[k].bond,en[k].bend,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals);
			#else
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%d	%d	%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,sim_dos[k].dos_acc,mc_pivot.pivot_acc[k],
				mc_trans.trans_acc[k],en[k].bond,en[k].bend,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].esasa,en[k].totals);
			#endif
		}
		#endif //CTDOS

		#ifdef XEDOS 
		if (sim.ID ==10) 
		{
			#ifndef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%d	%d	%d	%d	%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,mc_pivot.pivot_acc[k],
				mc_trans.trans_acc[k],mc_axial.axial_acc[k],mc_rand.rand_acc[k],sim_dos[k].dos_acc,en[k].bond,en[k].bend,en[k].tors,
				en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals);
			#else
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%d	%d	%d	%d %d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,mc_pivot.pivot_acc[k],
				mc_trans.trans_acc[k],mc_axial.axial_acc[k],mc_rand.rand_acc[k],sim_dos[k].dos_acc,en[k].bond,en[k].bend,en[k].tors,en[k].impr,
				en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].esasa,en[k].totals);
			#endif
		}
		#endif // XEDOS

		#ifdef FX_EDOS 
		if (sim.ID ==11) 
		{
			#ifndef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%d	%d	%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,sim_dos[k].dos_acc,mc_pivot.pivot_acc[k],
				mc_trans.trans_acc[k],en[k].bond,en[k].bend,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals);
			#else
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%d	%d	%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,sim_dos[k].dos_acc,mc_pivot.pivot_acc[k],
				mc_trans.trans_acc[k],en[k].bond,en[k].bend,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].esasa,en[k].totals);
			#endif
		}
		#endif //FX_EDOS

		#ifdef TDXEDOS 
		if (sim.ID ==12) 
		{
			#ifndef SASA
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%d	%d	%d	%d	%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,mc_pivot.pivot_acc[k],
				mc_trans.trans_acc[k],mc_axial.axial_acc[k],mc_rand.rand_acc[k],sim_dos[k].dos_acc,en[k].bond,en[k].bend,en[k].tors,
				en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals);
			#else
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%d	%d	%d	%d %d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,mc_pivot.pivot_acc[k],
				mc_trans.trans_acc[k],mc_axial.axial_acc[k],mc_rand.rand_acc[k],sim_dos[k].dos_acc,en[k].bond,en[k].bend,en[k].tors,en[k].impr,
				en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].esasa,en[k].totals);
			#endif
		}
		#endif // TDXEDOS

		#else //NEUTRAL

		#ifdef DOS 
		if (sim.ID ==2) 
		{
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%ld	%ld	%ld	%ld	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,sim_dos[k].dos_acc,mc_pivot.pivot_acc[k],
				mc_trans.trans_acc[k],mc_rand.rand_acc[k],en[k].bond,en[k].bend,en[k].urey,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals);
		}
		#endif //DOS
			
		#ifdef MMDOS 
		if (sim.ID ==7) 
		{
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%d	%d	%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,sim_dos[k].dos_acc,mc_pivot.pivot_acc[k],
				sim_dos[k].dke_acc,en[k].bond,en[k].bend,en[k].urey,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals);
		}
		#endif //MMDOS

		#ifdef CTDOS 
		if (sim.ID ==9) 
		{
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%ld	%ld	%ld	%ld	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,sim_dos[k].dos_acc,mc_pivot.pivot_acc[k],
				mc_trans.trans_acc[k],mc_rand.rand_acc[k],en[k].bond,en[k].bend,en[k].urey,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals);
		}
		#endif //CTDOS

		#ifdef XEDOS 
		if (sim.ID ==10) 
		{
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%d	%d	%d	%d	%d	%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,mc_pivot.pivot_acc[k],
				mc_trans.trans_acc[k],mc_axial.axial_acc[k],mc_rand.rand_acc[k],mc_solv.solv_acc[k],sim_dos[k].dos_acc,en[k].bond,en[k].bend,
				en[k].urey,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals);
		}
		#endif //XEDOS

		#ifdef FX_EDOS 
		if (sim.ID ==11) 
		{
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%d	%d	%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,sim_dos[k].dos_acc,mc_pivot.pivot_acc[k],
				mc_trans.trans_acc[k],en[k].bond,en[k].bend,en[k].urey,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals);
		}
		#endif //FX_EDOS

		#ifdef TDXEDOS 
		if (sim.ID ==12) 
		{
			ener_ptr[k] += sprintf(ener_ptr[k], "%lu	%d	%d	%d	%d	%d	%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf		 \n", icycle,mc_pivot.pivot_acc[k],
				mc_trans.trans_acc[k],mc_axial.axial_acc[k],mc_rand.rand_acc[k],mc_solv.solv_acc[k],sim_dos[k].dos_acc,en[k].bond,en[k].bend,
				en[k].urey,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals);
		}
		#endif //TDXEDOS

		#endif //NEUTRAL


	}// k loop ends

}
