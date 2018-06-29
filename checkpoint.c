/********************************/
/* checkpoint.c                 */
/* Writes and reads checkpoint  */
/* files in binary, allowing    */
/* restart of simulation        */
/*				*/
/* Eric Mansfield		*/
/********************************/

#ifdef CHECKP

#include <stdlib.h>
#include <stdio.h>

#include "defines.h"

/*
#if defined(TWHAM) || defined(XWHAM)
#include "tak_histogram.h"
#endif
*/
#if defined(TWHAM) || defined(XWHAM)
void WriteHistogram(tak_histogram *hist, FILE *fout);
void ReadHistogram(tak_histogram *hist, FILE *fout);
#endif

void Terminate()
{
	printf("SIGTERM received.\n");
	printf("Exiting--please wait....\n");
	terminating = 1;
}

void Interrupt()
{
	printf("SIGINT received.\n");
	printf("Exiting--please wait....\n");
	terminating = 1;
}

int file_exists(const char *filename)
{
	FILE *file;
    if ( (file = fopen(filename, "r")) )
    {
        fclose(file);
        return 1;
    }
    return 0;
}

void InitCheckpoint()
{

	//Register the signal handler
	//Don't do signal checkpointing for parallel jobs--the behavior is not well defined
	#ifndef MPI
		signal(SIGINT, Interrupt);
		signal(SIGTERM, Terminate);
	#endif

	restart_equil = 0;
	restart_prod = 0;

	mkdir ("CHECKP", 0755);

	#ifdef MPI
	sprintf(equiname[0], "CHECKP/restart_equil%d.chkA", mpi.my_rank);
	sprintf(prodname[0], "CHECKP/restart_prod%d.chkA", mpi.my_rank);
	sprintf(equiname[1], "CHECKP/restart_equil%d.chkB", mpi.my_rank);
	sprintf(prodname[1], "CHECKP/restart_prod%d.chkB", mpi.my_rank);
	#else
	sprintf(equiname[0], "CHECKP/restart_equil.chkA");
	sprintf(prodname[0], "CHECKP/restart_prod.chkA");
	sprintf(equiname[1], "CHECKP/restart_equil.chkB");
	sprintf(prodname[1], "CHECKP/restart_prod.chkB");
	#endif

	//Detect and read the first value from each of the checkpoint files (should be iteration number)
	//And set checkptr to point to the file that is more recent (based on iteration number)
	
	unsigned long equi[2] = {0, 0};
	unsigned long prod[2] = {0, 0};

	for(int alt = 0; alt <= 1; alt++)
	{
		if(file_exists(equiname[alt]))
		{
			restart_equil = 1;
			//printf("Detected equilibration checkpoint file:%s\n", equiname[alt]);
			FILE *test;
			test = fopen(equiname[alt], "rb");
			fread(&equi[alt], sizeof(int), 1, test);
			fclose(test);
		}
	}
	for(int alt = 0; alt <= 1; alt++)
	{
		if(file_exists(prodname[alt]))
		{
			restart_equil = 2;
			restart_prod = 1;
			//printf("Detected production checkpoint file: %s\n", prodname[alt]);
			FILE *test;
			test = fopen(prodname[alt], "rb");
			fread(&prod[alt], sizeof(int), 1, test);
			fclose(test);
		}
	}
	
	if(restart_equil == 0)
	{
		checkptr = 0;	//Just start with 'A'
	}
	else if(restart_equil == 1)
	{
		//Figure out which of the checkpoint files is the most recent
		int mostrecent = 0;
		if(equi[1] > equi[0])
			mostrecent = 1;
		checkptr = mostrecent;
		#ifdef MPI
		printf("Restarting equilibration on job %i,\t iteration %lu\n", mpi.my_rank, equi[mostrecent]);
		#else
		printf("Restarting equilibration, iteration %lu\n", equi[mostrecent]);
		#endif
	}
	else if(restart_equil == 2)
	{
		int mostrecent = 0;
		if(prod[1] > prod[0])
			mostrecent = 1;
		checkptr = mostrecent;
		#ifdef MPI
		printf("Restarting production on job %i,\t iteration %lu\n", mpi.my_rank, prod[mostrecent]);
		#else
		printf("Restarting production, iteration %lu\n", prod[mostrecent]);
		#endif
	}
	//Now checkptr is set to the checkpoint file we will read next.
	//It is also the next file to which we will write.	
	
}


int WriteCheckpoint(const char *filename, unsigned long i, double vscal[NBOXES], int state)
{
	//Writes a binary checkpoint file. This routine should only
	//be called at the beginning of an iteration loop.
	//It may be called from production or equilibration.

	//Write all local variables (iteration number, etc)
	//Write all global variables that cannot be initialized in the
	//usual way


	FILE *fout;
	if(!(fout = fopen(filename, "w")))
	{
		fprintf(stderr, "ERROR: Cannot open %s for output!\n", filename);
		return 1;
	}
	else
	{
//		printf("Writing checkpoint file: %s\n", filename);
	}
	// ALL LOCAL VARIABLESa
	
	//Places where this routine may be called
/*
	I will hold off on hybrid_MC.c for now...
	hybrid_MC.c
		double vscal[NBOXES]
//		double H_initial[NBOXES]
//		double H_final[NBOXES]
		ulong i (n_iter is NOT set)
	npt_md.c
		double vscal[NBOXES]
		ulong i
	nve_md.c
		double vscal[NBOXES] (but not important)
		ulong i
	nvt_md.c
		double vscal[NBOXES]
		ulong i
	repexch.c
		double vscal[NBOXES]
		ulong i (n_iter NOT set)
*/	

	//Generate random numbers to use as hash codes to ensure a non-corrupt file

	

	/******** Write local variables ********/
	fwrite(&i,	sizeof(int),	1,	fout);	
	fwrite(vscal,	sizeof(double),	NBOXES,	fout);


	/******** idum_seed should be unique for each job. Write it here ********/

	fwrite(&idum_seed,	sizeof(idum_seed), 1, fout);
	
	/******** Write global variables *******/

	// struct atoms
	// Each struct is allocated with a separate calloc, so 
	// it is probably wise to write one struct at a time

	for(int k=0; k<sim.NB; k++)
	{
		int na = box[k].boxns + MEXTRA;
		for(int n=0; n<na; n++)
		{
			fwrite(&atom[k][n],		sizeof(struct atoms), 1, fout);
			fwrite(&atnopbc[k][n],		sizeof(struct atoms), 1, fout);
			fwrite(&atom_temp[k][n],		sizeof(struct atoms), 1, fout);
			fwrite(&atnopbc_temp[k][n],		sizeof(struct atoms), 1, fout);
			#if defined(SASA) && defined(DLIST)  // See below
			fwrite(&rx0[k][n], sizeof(double), 1, fout);
			fwrite(&ry0[k][n], sizeof(double), 1, fout);
			fwrite(&rz0[k][n], sizeof(double), 1, fout);
			#endif
		}
	}

	#ifdef SASA	// Because of the strong dependence of neighborlist on the SASA potential,
			// it is necessary to write all of the neighbor data when checkpointing.

	for(int k=0; k<sim.NB; k++)
	{
		int size = box[k].boxns
		#ifdef WALL
		+ 2
		#endif
		;
		int max_num_neighbors = 0;
    		if(size < 2000) max_num_neighbors = size;
		    else max_num_neighbors = 1500;

		#ifdef DNA_GOLIK
		    max_num_neighbors = 1000;
		#endif

		#ifndef SPME
			if(size > 2000) 	max_num_neighbors = box[k].boxns/2;
		#endif
		
		for(int n=0; n<size; n++)
		{
			fwrite(&nlist[k].count_sasa[n],	sizeof(int), 1, fout);
			fwrite(&nlist[k].count[n],	sizeof(int), 1, fout);

			for(int m=0; m<max_num_neighbors; m++)
			{
				fwrite(&nlist[k].list_sasa[n][m], sizeof(int), 1, fout);
				fwrite(&nlist[k].xflag_sasa[n][m], sizeof(int), 1, fout);
				fwrite(&nlist[k].list[n][m], sizeof(int), 1, fout);
				fwrite(&nlist[k].xflag[n][m], sizeof(int), 1, fout);
			}
		}
	}

	#endif
	
//	struct veloc **vv, **uu, **ff, **vcm, **ff_short, **ff_long, **ff_temp, **vv_temp
	for(int k=0; k<sim.NB; k++)
	{
		int na = box[k].boxns + MEXTRA;
		for(int n=0; n<na; n++)
		{
			fwrite(&vv[k][n],	sizeof(struct veloc), 1, fout);
			fwrite(&uu[k][n],	sizeof(struct veloc), 1, fout);
			fwrite(&ff[k][n],	sizeof(struct veloc), 1, fout);
/*			fwrite(&ff_short[k][n],	sizeof(struct veloc), 1, fout);
			fwrite(&ff_long[k][n],	sizeof(struct veloc), 1, fout);
			fwrite(&ff_temp[k][n],	sizeof(struct veloc), 1, fout);
			fwrite(&vv_temp[k][n],	sizeof(struct veloc), 1, fout);
*/		}
	}

	fwrite(en, sizeof(struct energy), NBOXES, fout);
	fwrite(en_temp, sizeof(struct energy), NBOXES, fout);

	#ifdef PRESSURE
	fwrite(pvir, sizeof(struct virial), NBOXES, fout);
	fwrite(pvir_temp, sizeof(struct virial), NBOXES, fout);
	#ifdef MC
	fwrite(pviro, sizeof(struct virial), NBOXES, fout);
	#endif
	#endif

//	struct struct_quant ordparam[NBOXES]
	fwrite(ordparam,	sizeof(struct struct_quant), NBOXES, fout);		
	
//	fwrite(&sim,	sizeof(struct SimData), 1, fout);

	fwrite(box,	sizeof(struct BoxData), NBOXES,	fout);

//	struct results res[NBOXES], resb[NBOXES]
	fwrite(res,	sizeof(struct results), NBOXES,	fout);
	fwrite(resb,	sizeof(struct results), NBOXES, fout);
	
	fwrite(&idum,		sizeof(idum), 1, fout);	
	fwrite(&idum2,		sizeof(idum2), 1, fout);	
	fwrite(&iy,		sizeof(iy), 1, fout);	
	fwrite(iv,		sizeof(iv[0]), 32, fout);	
	fwrite(&idum_seed,	sizeof(idum_seed), 1, fout);

	//Nose Hoover Chain Stuff

	for(int k=0; k<sim.NB; k++)
	{
		for(int nw=0; nw<nhc[k].Nys; nw++)
		{
			fwrite(&nhc[k].w[nw], sizeof(double), 1, fout);
		}
		for(int nzeta=0; nzeta<nhc[k].M; nzeta++)
		{
			fwrite(&nhc[k].zeta[nzeta], sizeof(struct NH_therm), 1, fout);
		}
	}

#ifdef PR_NPT
	fwrite(axes,	sizeof(double), NBOXES * 9, fout);
	fwrite(axesi,	sizeof(double), NBOXES * 9, fout);
	fwrite(va,	sizeof(double), NBOXES * 9, fout);
	fwrite(Ga,	sizeof(double), NBOXES * 9, fout);
	fwrite(vvab,	sizeof(double), NBOXES * 9, fout);
	fwrite(frab,	sizeof(double), NBOXES * 9, fout);
	fwrite(delta,	sizeof(double), 9, fout);
	fwrite(keaxes,	sizeof(double), NBOXES, fout);
	fwrite(sum_ktzeta_i,	sizeof(double), NBOXES, fout);
	fwrite(ke_thermostat,	sizeof(double), NBOXES, fout);
	fwrite(cnf9ktzeta1,	sizeof(double), NBOXES, fout);
	fwrite(cka,	sizeof(double), NBOXES, fout);
	fwrite(pr_pv,	sizeof(double), NBOXES, fout);
#endif //PR_NPT

//	sasa
#ifdef SASA
	for(int k=0; k<sim.NB; k++)
	{
		int na = box[k].boxns + MEXTRA;
		for(int n=0; n<na; n++)
		{
			fwrite(&sasa[k][n], sizeof(struct imp_sasa), 1, fout);
		}
	}
#endif

//	struct hybrid sim_hyb;
	fwrite(&sim_hyb, sizeof(struct hybrid), 1, fout);
	
//	double PIVOT[3][2][3][3];
	fwrite(PIVOT,	sizeof(double), 3*2*3*3, fout);
//	struct pivot mc_pivot;
	fwrite(&mc_pivot, sizeof(struct pivot), 1, fout);
//	struct trans mc_trans;
	fwrite(&mc_trans, sizeof(struct trans), 1, fout);
//	struct axial mc_axial;
	fwrite(&mc_axial, sizeof(struct axial), 1, fout);
//	struct arandoml mc_rand;
	fwrite(&mc_rand, sizeof(struct arandoml), 1, fout);
//	long mc_hmc_acc[NBOXES];
	fwrite(mc_hmc_acc, sizeof(long), NBOXES, fout);
//	struct solv_rand mc_solv;
	fwrite(&mc_solv, sizeof(struct solv_rand), 1, fout);

#ifdef MPI
  //only need mpi.flag -- DO NOT checkpoint things like rank; it is assigned
  //by the scheduler
  // ----------------------------------------- /
  // This is for running the code on parallel  /
  // processors using MPI					   /
  // ----------------------------------------- /
//	struct mpi_def mpi;

	fwrite(&mpi.flag, sizeof(int), 1, fout);

#endif


//	struct counter frames[NBOXES];
	fwrite(frames, sizeof(struct counter), NBOXES, fout);

#ifdef CONFIGT
//	struct hessn config[NBOXES],config_temp[NBOXES];
	fwrite(config,		sizeof(struct hessn), NBOXES, fout);
#endif


#ifdef ELASTIC
//	struct elasticity_tensor elast[NBOXES];
	fwrite(elast, sizeof(struct elasticity_tensor), NBOXES, fout);
#endif //ELASTIC

#ifdef SMD
//	struct smd_struct *steermd;
	for(int n=0; n<sim.NB; n++)
	{
		fwrite(steermd[n], sizeof(struct smd_struct), 1, fout);
	}
#endif //SMD

#ifdef NSMD
//	struct smd_struct *steermd;
	for(int n=0; n<sim.NB; n++)
	{
		fwrite(steermd[n], sizeof(struct smd_struct), 1, fout);
	}
#endif //NSMD

//	int fileindex;
	fwrite(&fileindex, sizeof(int), 1, fout);

//	unsigned long n_iter;				// This is a global variable to keep track of number of iterations
	fwrite(&n_iter, sizeof(unsigned long), 1, fout);

#ifdef REM
//	struct replica rep_exc[NBOXES];
	fwrite(rep_exc, sizeof(struct replica), NBOXES, fout);
//	struct junction rep_jnc[NBOXES];
	fwrite(rep_jnc, sizeof(struct replica), NBOXES, fout);
//	struct replica_opt rep_opt[NBOXES];
	fwrite(rep_opt, sizeof(struct replica), NBOXES, fout);
#endif//REM


// We DO need these if we want to be safe (in case the program crashes and the
// buffers are lost
/*
char **simul_buff;
char **simul_ptr;
time_t *simul_time;
char **ener_buff;
char **ener_ptr;
time_t *ener_time;
char **ord_buff;
char **ord_ptr;
time_t *ord_time;
char **swap_buff;
char **swap_ptr;
time_t *swap_time;


#if defined(SMD) || defined(NSMD)
	char **smd_buff;
	char **smd_ptr;
	time_t *smd_time;
#endif

#if defined(TWHAM) || defined(XWHAM)
  time_t *wham_time;
#endif
*/

//Write the buffers
for(int k=0; k<sim.NB; k++)
{
	fwrite(simul_buff[k],	1, MAXBUFF, fout);
	fwrite(simul_ptr[k],	sizeof(simul_ptr[k]), 1, fout);
	fwrite(&simul_time[k],	sizeof(time_t), 1, fout);
	fwrite(ener_buff[k],	1, MAXBUFF, fout);
	fwrite(ener_ptr[k],	sizeof(simul_ptr[k]), 1, fout);
	fwrite(&ener_time[k],	sizeof(time_t), 1, fout);
	fwrite(ord_buff[k],	1, MAXBUFF, fout);
	fwrite(ord_ptr[k],	sizeof(simul_ptr[k]), 1, fout);
	fwrite(&ord_time[k],	sizeof(time_t), 1, fout);
	fwrite(swap_buff[k],	1, MAXBUFF, fout);
	fwrite(swap_ptr[k],	sizeof(simul_ptr[k]), 1, fout);
	fwrite(&swap_time[k],	sizeof(time_t), 1, fout);
	#if defined(SMD) || defined(NSMD)
		fwrite(smd_buff[k],	1, MAXBUFF, fout);
		fwrite(smd_ptr[k],	sizeof(smd_ptr[k]), MAXBUFF, fout);
		fwrite(&smd_time[k],	sizeof(time_t), 1, fout);
	#endif
	#if defined(TWHAM) || defined(XWHAM)		
		fwrite(&wham_time[k],	sizeof(time_t), 1, fout);
	#endif

}	
	


// ========================================== //
// variable needed for wham.                  //
// ========================================== //
#if defined(TWHAM) || defined(XWHAM)
//double *pe_max;
//double *pe_min;
#ifdef TWHAM

	for(int k=0; k<sim.NB; k++)
	{
		fwrite(&pe_max[k], sizeof(double), 1, fout);
		fwrite(&pe_min[k], sizeof(double), 1, fout);
	}
#endif //TWHAM

	if(state)
	{
		for(int k=0; k<sim.NB; k++)
		{
			WriteHistogram(h_ebond[k],	fout);
			WriteHistogram(h_ebend[k],	fout);
			WriteHistogram(h_eurey[k],	fout);
			WriteHistogram(h_etors[k],	fout);
			WriteHistogram(h_eimpr[k],	fout);
			WriteHistogram(h_elj[k],	fout);
			WriteHistogram(h_eqq[k],	fout);
			WriteHistogram(h_pe[k],	fout);
			WriteHistogram(h_ke[k],	fout);
			WriteHistogram(h_etot[k],	fout);
			WriteHistogram(h_hel[k],	fout);
			WriteHistogram(h_dnc[k],	fout);
			WriteHistogram(h_con[k],	fout);
			WriteHistogram(h_con_2[k],	fout);
			WriteHistogram(h_x1[k],	fout);
			WriteHistogram(h_x2[k],	fout);
			WriteHistogram(h_rmsd[k],	fout);
			WriteHistogram(h_gyr[k],	fout);
		}
	}
#endif //TWHAM || XWHAM

	/***** Write idum_seed one final time *****/	
	fwrite(&idum_seed,	sizeof(idum_seed), 1, fout);

	fclose(fout);
	return 0;

}

int ReadCheckpoint(const char *filename, int *i, double vscal[NBOXES], int state)
{
	//TODO: Test the local variables!!!!!!!
	
	//Reads a binary checkpoint file. This routine should only
	//be called at the beginning of an iteration loop.
	//It may be called from production or equilibration.

	//Read all local variables (iteration number, etc)
	//Read all global variables that cannot be initialized in the
	//usual way

	FILE *fout;
	if(!(fout = fopen(filename, "rb")))
	{
		fprintf(stderr, "Error reading %s: Cannot open for input.\n", filename);
		exit(1);
	}
	else
	{
		//printf("Reading checkpoint file: %s\n", filename);
	}

	// ALL LOCAL VARIABLES
	
	//Places where this routine may be called
/*
	I will hold off on hybrid_MC.c for now...
	hybrid_MC.c
		double vscal[NBOXES]
//		double H_initial[NBOXES]
//		double H_final[NBOXES]
		ulong i (n_iter is NOT set)
	npt_md.c
		double vscal[NBOXES]
		ulong i
	nve_md.c
		double vscal[NBOXES] (but not important)
		ulong i
	nvt_md.c
		double vscal[NBOXES]
		ulong i
	repexch.c
		double vscal[NBOXES]
		ulong i (n_iter NOT set)
*/	

	/******** Write local variables ********/
	fread(i,	sizeof(int),	1,	fout);	
	fread(vscal,	sizeof(double),	NBOXES,	fout);

	/******** Read idum_seed to ensure the file is not corrupt ******/

	long idum_seed_check1 = 1, idum_seed_check2 = 2;	
	fread(&idum_seed_check1,	sizeof(long), 1, fout);


	/******** Write global variables *******/

	// struct atoms
	// Each struct is allocated with a separate calloc, so 
	// it is probably wise to write one struct at a time
	for(int k=0; k<sim.NB; k++)
	{
		int na = box[k].boxns + MEXTRA;
		for(int n=0; n<na; n++)
		{
			fread(&atom[k][n],		sizeof(struct atoms), 1, fout);
			fread(&atnopbc[k][n],		sizeof(struct atoms), 1, fout);
			fread(&atom_temp[k][n],		sizeof(struct atoms), 1, fout);
			fread(&atnopbc_temp[k][n],		sizeof(struct atoms), 1, fout);
			#if defined(SASA) && defined(DLIST)
			fread(&rx0[k][n], sizeof(double), 1, fout);
			fread(&ry0[k][n], sizeof(double), 1, fout);
			fread(&rz0[k][n], sizeof(double), 1, fout);
			#endif
		}
	}

	#ifdef SASA	// Because of the strong dependence of neighborlist on the SASA potential,
			// it is necessary to write all of the neighbor data when checkpointing.

	for(int k=0; k<sim.NB; k++)
	{
		int size = box[k].boxns
		#ifdef WALL
		+ 2
		#endif
		;
		int max_num_neighbors = 0;
    		if(size < 2000) max_num_neighbors = size;
		    else max_num_neighbors = 1500;

		#ifdef DNA_GOLIK
		    max_num_neighbors = 1000;
		#endif

		#ifndef SPME
			if(size > 2000) 	max_num_neighbors = box[k].boxns/2;
		#endif
		
		for(int n=0; n<size; n++)
		{
			fread(&nlist[k].count_sasa[n],	sizeof(int), 1, fout);
			fread(&nlist[k].count[n],	sizeof(int), 1, fout);

			for(int m=0; m<max_num_neighbors; m++)
			{
				fread(&nlist[k].list_sasa[n][m], sizeof(int), 1, fout);
				fread(&nlist[k].xflag_sasa[n][m], sizeof(int), 1, fout);
				fread(&nlist[k].list[n][m], sizeof(int), 1, fout);
				fread(&nlist[k].xflag[n][m], sizeof(int), 1, fout);
			}
		}
	}

	#endif

//	struct veloc **vv, **uu, **ff, **vcm, **ff_short, **ff_long, **ff_temp, **vv_temp


	for(int k=0; k<sim.NB; k++)
	{
		int na = box[k].boxns + MEXTRA;
		for(int n=0; n<na; n++)
		{
			fread(&vv[k][n],	sizeof(struct veloc), 1, fout);
			fread(&uu[k][n],	sizeof(struct veloc), 1, fout);
			fread(&ff[k][n],	sizeof(struct veloc), 1, fout);
/*			fread(&ff_short[k][n],	sizeof(struct veloc), 1, fout);
			fread(&ff_long[k][n],	sizeof(struct veloc), 1, fout);
			fread(&ff_temp[k][n],	sizeof(struct veloc), 1, fout);
			fread(&vv_temp[k][n],	sizeof(struct veloc), 1, fout);
*/		}
	}

	fread(en, sizeof(struct energy), NBOXES, fout);
	fread(en_temp, sizeof(struct energy), NBOXES, fout);

	#ifdef PRESSURE
	fread(pvir, sizeof(struct virial), NBOXES, fout);
	fread(pvir_temp, sizeof(struct virial), NBOXES, fout);
	#ifdef MC
	fread(pviro, sizeof(struct virial), NBOXES, fout);
	#endif
	#endif

//	struct struct_quant ordparam[NBOXES]
	fread(ordparam,	sizeof(struct struct_quant), NBOXES, fout);		

//	fread(&sim,	sizeof(struct SimData), 1, fout);
//	struct BoxData box[NBOXES]
	fread(box,	sizeof(struct BoxData), NBOXES,	fout);

//	struct results res[NBOXES], resb[NBOXES]
	fread(res,	sizeof(struct results), NBOXES,	fout);
	fread(resb,	sizeof(struct results), NBOXES, fout);
	
	fread(&idum,		sizeof(idum), 1, fout);	
	fread(&idum2,		sizeof(idum2), 1, fout);	
	fread(&iy,		sizeof(iy), 1, fout);	
	fread(&iv,		sizeof(iv[0]), 32, fout);	
	fread(&idum_seed,	sizeof(idum_seed), 1, fout);


	//Nose Hoover Chain Stuff

	for(int k=0; k<sim.NB; k++)
	{
		for(int nw=0; nw<nhc[k].Nys; nw++)
		{
			fread(&nhc[k].w[nw], sizeof(double), 1, fout);
		}
		for(int nzeta=0; nzeta<nhc[k].M; nzeta++)
		{
			fread(&nhc[k].zeta[nzeta], sizeof(struct NH_therm), 1, fout);
		}
	}

#ifdef PR_NPT
	fread(axes,	sizeof(double), NBOXES * 9, fout);
	fread(axesi,	sizeof(double), NBOXES * 9, fout);
	fread(va,	sizeof(double), NBOXES * 9, fout);
	fread(Ga,	sizeof(double), NBOXES * 9, fout);
	fread(vvab,	sizeof(double), NBOXES * 9, fout);
	fread(frab,	sizeof(double), NBOXES * 9, fout);
	fread(delta,	sizeof(double), 9, fout);
	fread(keaxes,	sizeof(double), NBOXES, fout);
	fread(sum_ktzeta_i,	sizeof(double), NBOXES, fout);
	fread(ke_thermostat,	sizeof(double), NBOXES, fout);
	fread(cnf9ktzeta1,	sizeof(double), NBOXES, fout);
	fread(cka,	sizeof(double), NBOXES, fout);
	fread(pr_pv,	sizeof(double), NBOXES, fout);
#endif //PR_NPT

//	sasa
#ifdef SASA
	for(int k=0; k<sim.NB; k++)
	{
		int na = box[k].boxns + MEXTRA;
		for(int n=0; n<na; n++)
		{
			fread(&sasa[k][n], sizeof(struct imp_sasa), 1, fout);
		}
	}
#endif

//	struct hybrid sim_hyb;
	fread(&sim_hyb, sizeof(struct hybrid), 1, fout);
	
//	double PIVOT[3][2][3][3];
	fread(PIVOT,	sizeof(double), 3*2*3*3, fout);
//	struct pivot mc_pivot;
	fread(&mc_pivot, sizeof(struct pivot), 1, fout);
//	struct trans mc_trans;
	fread(&mc_trans, sizeof(struct trans), 1, fout);
//	struct axial mc_axial;
	fread(&mc_axial, sizeof(struct axial), 1, fout);
//	struct arandoml mc_rand;
	fread(&mc_rand, sizeof(struct arandoml), 1, fout);
//	long mc_hmc_acc[NBOXES];
	fread(mc_hmc_acc, sizeof(long), NBOXES, fout);
//	struct solv_rand mc_solv;
	fread(&mc_solv, sizeof(struct solv_rand), 1, fout);


#ifdef MPI
  //only need mpi.flag -- DO NOT checkpoint things like rank; it is assigned
  //by the scheduler
  // ----------------------------------------- /
  // This is for running the code on parallel  /
  // processors using MPI					   /
  // ----------------------------------------- /
//	struct mpi_def mpi;

	fread(&mpi.flag, sizeof(int), 1, fout);

#endif

//	struct counter frames[NBOXES];
	fread(frames, sizeof(struct counter), NBOXES, fout);

#ifdef CONFIGT
//	struct hessn config[NBOXES],config_temp[NBOXES];
	fread(config,		sizeof(struct hessn), NBOXES, fout);
#endif

#ifdef ELASTIC
//	struct elasticity_tensor elast[NBOXES];
	fread(elast, sizeof(struct elasticity_tensor), NBOXES, fout);
#endif //ELASTIC

#ifdef SMD
//	struct smd_struct *steermd;
	for(int n=0; n<sim.NB; n++)
	{
		fread(steermd[n], sizeof(struct smd_struct), 1, fout);
	}
#endif //SMD

#ifdef NSMD
//	struct smd_struct *steermd;
	for(int n=0; n<sim.NB; n++)
	{
		fread(steermd[n], sizeof(struct smd_struct), 1, fout);
	}
#endif //NSMD

//	int fileindex;
	fread(&fileindex, sizeof(int), 1, fout);

//	unsigned long n_iter;				// This is a global variable to keep track of number of iterations
	fread(&n_iter, sizeof(unsigned long), 1, fout);

#ifdef REM
//	struct replica rep_exc[NBOXES];
	fread(rep_exc, sizeof(struct replica), NBOXES, fout);
//	struct junction rep_jnc[NBOXES];
	fread(rep_jnc, sizeof(struct replica), NBOXES, fout);
//	struct replica_opt rep_opt[NBOXES];
	fread(rep_opt, sizeof(struct replica), NBOXES, fout);
#endif//REM


//We don't need these (I don't think) as long as we call ioflush(1) before exiting
/*
char **simul_buff;
char **simul_ptr;
time_t *simul_time;
char **ener_buff;
char **ener_ptr;
time_t *ener_time;
char **ord_buff;
char **ord_ptr;
time_t *ord_time;
char **swap_buff;
char **swap_ptr;
time_t *swap_time;
*/

//Read the buffers
for(int k=0; k<sim.NB; k++)
{
	fread(simul_buff[k],	1, MAXBUFF, fout);
	fread(simul_ptr[k],	sizeof(simul_ptr[k]), 1, fout);
	fread(&simul_time[k],	sizeof(time_t), 1, fout);
	fread(ener_buff[k],	1, MAXBUFF, fout);
	fread(ener_ptr[k],	sizeof(simul_ptr[k]), 1, fout);
	fread(&ener_time[k],	sizeof(time_t), 1, fout);
	fread(ord_buff[k],	1, MAXBUFF, fout);
	fread(ord_ptr[k],	sizeof(simul_ptr[k]), 1, fout);
	fread(&ord_time[k],	sizeof(time_t), 1, fout);
	fread(swap_buff[k],	1, MAXBUFF, fout);
	fread(swap_ptr[k],	sizeof(simul_ptr[k]), 1, fout);
	fread(&swap_time[k],	sizeof(time_t), 1, fout);
	#if defined(SMD) || defined(NSMD)
		fread(smd_buff[k],	1, MAXBUFF, fout);
		fread(smd_ptr[k],	sizeof(smd_ptr[k]), MAXBUFF, fout);
		fread(&smd_time[k],	sizeof(time_t), 1, fout);
	#endif
	#if defined(TWHAM) || defined(XWHAM)		
		fread(&wham_time[k],	sizeof(time_t), 1, fout);
	#endif

}	
// ========================================== //
// variable needed for wham.                  //
// ========================================== //
#if defined(TWHAM) || defined(XWHAM)
//double *pe_max;
//double *pe_min;
#ifdef TWHAM
	for(int k=0; k<sim.NB; k++)
	{
		fread(&pe_max[k], sizeof(double), 1, fout);
		fread(&pe_min[k], sizeof(double), 1, fout);
	}
#endif //TWHAM

	if(state)
	{
		//Init the histograms given pe_max, pe_min
		for(int k=0; k<sim.NB; k++) twham_init(k);
		
		for(int k=0; k<sim.NB; k++)
		{
			ReadHistogram(h_ebond[k],	fout);
			ReadHistogram(h_ebend[k],	fout);
			ReadHistogram(h_eurey[k],	fout);
			ReadHistogram(h_etors[k],	fout);
			ReadHistogram(h_eimpr[k],	fout);
			ReadHistogram(h_elj[k],	fout);
			ReadHistogram(h_eqq[k],	fout);
			ReadHistogram(h_pe[k],	fout);
			ReadHistogram(h_ke[k],	fout);
			ReadHistogram(h_etot[k],	fout);
			ReadHistogram(h_hel[k],	fout);
			ReadHistogram(h_dnc[k],	fout);
			ReadHistogram(h_con[k],	fout);
			ReadHistogram(h_con_2[k],	fout);
			ReadHistogram(h_x1[k],	fout);
			ReadHistogram(h_x2[k],	fout);
			ReadHistogram(h_rmsd[k],	fout);
			ReadHistogram(h_gyr[k],	fout);
		}
	}
#endif //TWHAM || XWHAM

	fread(&idum_seed_check2,	sizeof(long), 1, fout);

	fclose(fout);

	if(idum_seed != idum_seed_check1 || idum_seed != idum_seed_check2)
	{
		printf("idum_seed = %li, idum_seed_check1 = %li, idum_seed_check2 = %li\n", idum_seed, idum_seed_check1, idum_seed_check2);
		fprintf(stderr, "Error reading %s: Corrupt checkpoint file. Please delete!\n", filename);
		exit(1);
	}
	
	return 0;
}

#if defined(TWHAM) || defined(XWHAM)
void WriteHistogram(tak_histogram *hist, FILE *fout)
{
	//In theory, we should only need to write the bin data.
//	printf("Entering WriteHistogram\n");
//	printf("hist = %i", hist);
	for(int i=0; i<hist->n; i++)
	{
//		printf("Going to output hist->bin[i] (i = %i)\n", i);
//		printf("%f\n", hist->bin[i]);
		fwrite(&(hist->bin[i]), sizeof(double), 1, fout);
	}
}

void ReadHistogram(tak_histogram *hist, FILE *fout)
{
	//In theory, we should only need to write the bin data.
	for(int i=0; i<hist->n; i++)
	{
		fread(&(hist->bin[i]), sizeof(double), 1, fout);
	}
}
#endif


#endif /* CHECKP */
