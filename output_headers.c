/* ======================================================================== */
/* output_headers.c                                                         */
/* ======================================================================== */

#ifdef OHEADER

#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

void output_headers ()
{
  /* ================================================================== */
  /* Write the current values to simul#.header.                         */
  /* ================================================================== */


	for(int k=0; k<sim.NB; k++) {

		/* ================================================================== */
		/* Write the energy values to order#.header.                          */
		/* ================================================================== */

		char ord_ptr_hdr[32768];

		#ifdef CONFIGT
			sprintf(ord_ptr_hdr,"%s	%s	%s	%s	","icycle", "ordparam[k].hel","ordparam[k].d_nc","ordparam[k].con");
			for(int i=0; i<contacts[k].n; i++)	sprintf(ord_ptr_hdr,"%s","contacts[k].code[i]");
			sprintf(ord_ptr_hdr,"	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s\n","config[k].T", 
				"box[k].temp","ordparam[k].force_1","ordparam[k].force_2","ordparam[k].gyr","ordparam[k].rmsd",
				"theta_av*180.0/PI","psi_av*180.0/PI","phi_av*180.0/PI","ordparam[k].x1","ordparam[k].x2","ordparam[k].con_2");
		#endif
		#ifndef CONFIGT
			sprintf(ord_ptr_hdr,"%s	%s	%s	%s	","icycle", "ordparam[k].hel",
				"ordparam[k].d_nc","ordparam[k].con");
			for(int i=0; i<contacts[k].n; i++)	sprintf(ord_ptr_hdr,"%s","contacts[k].code[i]");
			sprintf(ord_ptr_hdr,"	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s\n","box[k].temp","ordparam[k].force_1",
				"ordparam[k].force_2","ordparam[k].gyr","ordparam[k].rmsd","theta_av*180.0/PI","psi_av*180.0/PI","phi_av*180.0/PI",
				"ordparam[k].x1","ordparam[k].x2","ordparam[k].con_2");

			/* ================================================================== */
			/* Write the energy values to smd#.header.                            */
			/* ================================================================== */
			char smd_ptr_hdr[32768];
			#ifdef SMD
				init_smd(k);
				sprintf(smd_ptr_hdr,"%s	%s	%s	%s	%s	%s	%s	%s	%s\n","icycle","tt","steermd[k].effr1_av",
					"ordparam[k].d_nc","steermd[k].effr2_av","steermd[k].r2","steermd[k].dw1","steermd[k].dw2","en[k].H");
			#endif
			#ifdef NSMD
				init_smd(k);
				sprintf(smd_ptr_hdr,"%s	%s	%s	%s	%s	\n","icycle","tt","steermd[k].effr1_av","steermd[k].f2_av","en[k].H");
			#endif

			char smd_name[100];
			FILE *fp_smd;
			#ifdef MPI
			sprintf(smd_name,"./OUTPUT/BOX%d/smd%d.header",mpi.my_rank,mpi.my_rank);
			#endif
			#ifndef MPI
			sprintf(smd_name,"./OUTPUT/BOX%d/smd%d.header",k,k);
			#endif
			if ( !(fp_smd = fopen(smd_name, "w")) )
			{
				fprintf(stderr, "Could not open %s for input!\n", smd_name);
				exit(1);
			}		
			fprintf(fp_smd, "%s", smd_ptr_hdr);
			fclose(fp_smd);

		#endif

		char ord_name[100];
		FILE *fp_ord;
		#ifdef MPI
		sprintf(ord_name,"./OUTPUT/BOX%d/order%d.header",mpi.my_rank,mpi.my_rank);
		#endif
		#ifndef MPI
		sprintf(ord_name,"./OUTPUT/BOX%d/order%d.header",k,k);
		#endif
		if ( !(fp_ord = fopen(ord_name, "w")) )
		{
			fprintf(stderr, "Could not open %s for input!\n", ord_name);
			exit(1);
		}		
		fprintf(fp_ord, "%s", ord_ptr_hdr);
		fclose(fp_ord);
		

		/* ================================================================== */
		/* Write the energy values to ener_box#.header.                       */
		/* ================================================================== */
		char ener_ptr_hdr[32768];	
		/* ----------------------------------------------- */
		/* If a Nose-Hoover thermostat is being used, the  */
		/* extended Hamiltonian is also written.NVE/NVT MD */
		/* ----------------------------------------------- */    
		if(sim.ID == 6 || ((sim.ID2 == 4 || sim.ID2 == 5 || sim.ID2==8) && sim.ID ==1))
		{
			#ifdef NEUTRAL
			#ifndef SASA
				sprintf(ener_ptr_hdr, "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s %s		 \n", "icycle","tt","en[k].bond","en[k].bend","en[k].tors",
					"en[k].impr","en[k].nbonds-en[k].coulomb","en[k].coulomb","en[k].potens","en[k].kinet","en[k].totals"," en[k].H"," box[k].press");
			#endif
			#ifdef SASA
			sprintf(ener_ptr_hdr, "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s %s		 \n","icycle","tt","en[k].bond","en[k].bend","en[k].tors",
				"en[k].impr","en[k].nbonds-en[k].coulomb","en[k].coulomb","en[k].potens","en[k].kinet"," en[k].esasa","en[k].totals"," en[k].H"," box[k].press");
			#endif
			#endif//NEUTRAL

			#ifndef NEUTRAL
			sprintf(ener_ptr_hdr, "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s %s		 \n", "icycle","tt","en[k].bond","en[k].bend","en[k].urey",
				"en[k].tors","en[k].impr","en[k].nbonds-en[k].coulomb","en[k].coulomb","en[k].potens","en[k].kinet","en[k].totals"," en[k].H"," box[k].press");
			#endif//NEUTRAL
		}

		/* ------------------------------------------------ */
		/* Either NVE or NVT MD without Nose Hoover			*/
		/* ------------------------------------------------ */
		else if ((sim.ID ==1) || (sim.ID ==3) || (sim.ID==8)) 
		{
			#ifdef NEUTRAL
				#ifndef SASA
				sprintf(ener_ptr_hdr, "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s %s		 \n", "icycle","tt","en[k].bond","en[k].bend","en[k].tors",
				"en[k].impr","en[k].nbonds-en[k].coulomb","en[k].coulomb","en[k].potens","en[k].kinet","en[k].totals","box[k].press");
				#endif
				#ifdef SASA
				sprintf(ener_ptr_hdr, "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s %s %s			 \n", "icycle","tt","en[k].bond","en[k].bend","en[k].tors",
				"en[k].impr","en[k].nbonds-en[k].coulomb","en[k].coulomb","en[k].potens","en[k].kinet","en[k].esasa","en[k].totals","box[k].press");
				#endif
			#endif//NEUTRAL

			#ifndef NEUTRAL
			sprintf(ener_ptr_hdr, "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s %s		 \n", "icycle","tt","en[k].bond","en[k].bend","en[k].urey",
			"en[k].tors","en[k].impr","en[k].nbonds-en[k].coulomb","en[k].coulomb","en[k].potens","en[k].kinet","en[k].totals","box[k].press");
			#endif//NEUTRAL
		}
		/* ------------------------------------------------ */
		/* Hybrid MD MC without thermostat					*/
		/* ------------------------------------------------ */
		else if (sim.ID ==4) 
		{
			#ifdef NEUTRAL
				#ifndef SASA
				sprintf(ener_ptr_hdr, "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s %s		 \n", "icycle","mc_pivot.pivot_acc[k]","mc_trans.trans_acc[k]","mc_rand.rand_acc[k]","mc_hmc_acc[k]","sim_hyb.swap_acc[k]","en[k].bond","en[k].bend","en[k].tors",
				"en[k].impr","en[k].nbonds-en[k].coulomb","en[k].coulomb","en[k].potens","en[k].kinet","en[k].totals","box[k].press");
				#endif
				#ifdef SASA
				sprintf(ener_ptr_hdr, "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s %s		 \n", "icycle","mc_pivot.pivot_acc[k]","mc_trans.trans_acc[k]","mc_rand.rand_acc[k]","mc_hmc_acc[k]","sim_hyb.swap_acc[k]","en[k].bond","en[k].bend","en[k].tors",
				"en[k].impr","en[k].nbonds-en[k].coulomb","en[k].coulomb","en[k].potens","en[k].kinet","en[k].esasa","en[k].totals","box[k].press");
				#endif
			#endif//NEUTRAL

			#ifndef NEUTRAL
			sprintf(ener_ptr_hdr, "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s %s		 \n", "icycle","mc_pivot.pivot_acc[k]","mc_trans.trans_acc[k]","mc_rand.rand_acc[k]","mc_hmc_acc[k]","sim_hyb.swap_acc[k]","en[k].bond","en[k].bend","en[k].urey",
			"en[k].tors","en[k].impr","en[k].nbonds-en[k].coulomb","en[k].coulomb","en[k].potens","en[k].kinet","en[k].totals","box[k].press");
			#endif//NEUTRAL
		}
		/* ------------------------------------------------ */
		/* Replica Exchange MD without thermostat			*/
		/* ------------------------------------------------ */
		else if (sim.ID ==5) 
		{
			#ifdef NEUTRAL
				#ifndef SASA
				sprintf(ener_ptr_hdr, "%s	%d	%s	%s	%s	%s	%s	%s	%s	%s	%s %s		 \n", "icycle","sim_hyb.swap_acc[k]",
					"en[k].bond","en[k].bend","en[k].tors","en[k].impr","en[k].nbonds-en[k].coulomb","en[k].coulomb","en[k].potens","en[k].kinet","en[k].totals","box[k].press");
				#endif
				#ifdef SASA
				sprintf(ener_ptr_hdr, "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s %s	 \n", "icycle","sim_hyb.swap_acc[k]",
					"en[k].bond","en[k].bend","en[k].tors","en[k].impr","en[k].nbonds-en[k].coulomb","en[k].coulomb","en[k].potens","en[k].kinet","en[k].esasa","en[k].totals","box[k].press");
				#endif
			#endif//NEUTRAL

			#ifndef NEUTRAL
			sprintf(ener_ptr_hdr, "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s %s		 \n", "icycle","sim_hyb.swap_acc[k]",
				"en[k].bond","en[k].bend","en[k].urey","en[k].tors","en[k].impr","en[k].nbonds-en[k].coulomb","en[k].coulomb","en[k].potens","en[k].kinet","en[k].totals","box[k].press");
			#endif//NEUTRAL
		}

		char ener_name[100];
		FILE *fp_ener;
		#ifdef MPI
		sprintf(ener_name,"./OUTPUT/BOX%d/ener_box%d.header",mpi.my_rank,mpi.my_rank);
		#endif
		#ifndef MPI
		sprintf(ener_name,"./OUTPUT/BOX%d/ener_box%d.header",k,k);
		#endif
		if ( !(fp_ener = fopen(ener_name, "w")) )
		{
			fprintf(stderr, "Could not open %s for input!\n", ener_name);
			exit(1);
		}		
		fprintf(fp_ener, "%s", ener_ptr_hdr);
		fclose(fp_ener);

		/* ================================================================== */
		/* Write the contributions to the coumbic energy values calculated by */
		/* the Ewald method to ewald#.header.                                 */
		/* ================================================================== */
		#ifdef EWALD_DEBUG
			FILE *ewald_fp;
			char ewald_name[100];
			#ifdef MPI
			sprintf(ewald_name,"./OUTPUT/BOX%d/ewald%d.header",mpi.my_rank,mpi.my_rank);
			#endif
			#ifndef MPI
			sprintf(ewald_name,"./OUTPUT/BOX%d/ewald%d.header",k,k);
			#endif
		
			ewald_fp= fopen(ewald_name,"a");
			fprintf(ewald_fp, "%s	%s	%s	%s	%s	%s	%s	%s\n", "icycle","tt","en[k].rewald","en[k].rewalds","en[k].kewald",
				"en[k].sewald","en[k].iewald","en[k].tewald");
			fclose(ewald_fp);
		#endif

		#ifdef PR_NPT_DEBUG
			FILE npt_fp;
			char npt_name[100];

			sprintf(npt_name,"./OUTPUT/BOX%d/npt%d.header",k,k);
			npt_fp= fopen(npt_name,"a");
			fprintf(npt_fp, "%s	%s	%s	%s	%s %s	%s	%s	%s\n", "tt"," en[k].poten","en[k].kinet"," ke_thermostat[k]", 
			"sum_ktzeta_i[k]"," cnf9ktzeta1[k]","pr_pv[k]","cka[k]","en[k].H");
			fclose(npt_fp);
		#endif

  /* ================================================================== */
  /* Write the angle the molecule makes with the wall in angle#.header. */
  /* ================================================================== */
/*#ifdef WALL_DEBUG
  FILE *wall_ptr;
  char wall_name[100];
  double x = atom[k][wall[k].angle_site].x;
  double y = atom[k][wall[k].angle_site].y;
  double z = atom[k][wall[k].angle_site].z;
  double angle = atan(z/sqrt(x*x + y*y))*180.0/PI;
  sprintf(wall_name,"./OUTPUT/BOX%d/angle%d.header",k,k);
  wall_ptr= fopen(wall_name,"a");
  fprintf(wall_ptr, "%s	%d	%s	%s	%s	%s\n",icycle,wall[k].angle_site,x,y,z,angle);
  fclose(wall_ptr);
#endif
*/

  }// k loop ends

}

#endif
