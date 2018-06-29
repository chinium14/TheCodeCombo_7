/* ======================================================================== */
/* read_sim.cpp                                                             */
/*                                                                          */
/*		This subroutine reads in the simulation parameters. These           */
/* parameters include the ensemble choice(NVT,NVE,NPT), the integration     */
/* mode, the number of steps, the frequency of the neighborlist update,     */
/* etc.  The input file name is simul.input.                                */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

void read_sim (void)
{

/* ======================================================================= */
/*                                                                         */
/* BEGIN READING simul.input AND repex#.input                              */
/*                                                                         */
/* ======================================================================= */


	/* ================================================================== */
	/* Variables and pointers needed in the subroutine.                   */
	/* ================================================================== */
	char tt[150];
	int boxnum;
	FILE *das;

	/* ================================================================== */
	/* Set the file pointers to the correct input file.                   */
	/* ================================================================== */
	if(NULL == (das=fopen("./INPUT/simul.input","r")) )
	{
		fprintf(stdout,"input file simul.input does not exist!!!\n");
		exit(1);
	}
	//das = fopen("./INPUT/simul.input","r");

	/* ================================================================== */
	/* Start reading the file information into the sim sturcture.         */
	/* ================================================================== */

	#ifdef SCRIPT

	/* =================================================== */
	/*     SIM VARIABLES                                   */
	/* =================================================== */

	char buff[256];

	// Default variables (some will cast warnings if not set)
	strcpy(title, "Untitled Simulation");
	sim.cyc_eq = 0;
	sim.cyc_pr = 0;
	sim.blockd = ~0;	//Highest unsigned long possible
	sim.blockc = ~0;
	sim.blockt = ~0;
	sim_hyb.cyc_hybmd = 0;
	sim_hyb.cyc_swap = 0;
	sim.drift = ~0;
	sim.nsteps = 1;
	sim.rc = 0;
	sim.rclist = 0;
	sim.rc_ew = 0;
	sim.rclist_ew = 0;
	sim.nlist = 0;

	double epsRF_temp = 0;

	// User must supply these variables.
	strcpy(param_file_name, "");
	sim.NB = 0;	// (Will be assumed to be 1, with a warning)
	sim.dt = -1;
	sim.ID = 0;
	sim.ID2 = 0;
	

	while ( !feof(das) )
	{
		// Read a single line from the file stream
		fgets(buff, 256, das);
		
		// Discard any whitespace
		int c=0;
		while(isspace(buff[c])) c++;
		memcpy(buff, buff+c, 256-c);		

		// Discard text appearing after the first # character
		char *pound = strchr(buff, '#');	//Return pointer to first '#'
		if (pound) *pound = '\0';		//Terminate string here
	
		if (strlen(buff) > 2)			//Sometimes small junk will come through--skip that
		{
			// Read in the name of this variable (until we hit whitespace) and convert to lowercase
			char varname_full[64];
			sscanf(buff, "%s", varname_full);

			// Is this a sim variable?
			if( strstr(varname_full, "sim_") == varname_full )	// The string begins with "sim"
			{
				
				//Here's a little hack to extract the rest of the string (the variable value)
				char *variable = buff + strlen(varname_full);
				char *varname = varname_full + 4;   //Get rid of the "sim_"



				for(size_t cc=0; cc < strlen(varname); cc++)
					varname[cc] = tolower(varname[cc]);

				if ( !strcmp(varname, "title") )
				{
					char title_first[80];
					sscanf (variable, "%s", title_first);	//Identify the beginning of the title (discard preceeding whitespace)
					strcpy(title, strstr(variable, title_first));   //Copy from that onward into title
					title[strlen(title)-2] = '\0';	//Remove junk at end of string
				}
				else if ( !strcmp(varname, "paramfile") )
					sscanf (variable, "%s", param_file_name);
				else if ( !strcmp(varname, "mode") )
					sscanf (variable, "%i", &sim.ID);
				else if ( !strcmp(varname, "intmode") )
					sscanf (variable, "%i", &sim.ID2);
				else if ( !strcmp(varname, "eq_steps") )
					sscanf (variable, "%lu", &sim.cyc_eq);
				else if ( !strcmp(varname, "pr_steps") )
					sscanf (variable, "%lu", &sim.cyc_pr);
				else if ( !strcmp(varname, "mdmc_steps") )
					sscanf (variable, "%li", &sim_hyb.cyc_hybmd);
				else if ( !strcmp(varname, "swap_freq") )
					sscanf (variable, "%li", &sim_hyb.cyc_swap);
				else if ( !strcmp(varname, "blockd") )
					sscanf (variable, "%lu", &sim.blockd);
				else if ( !strcmp(varname, "blockc") )
					sscanf (variable, "%lu", &sim.blockc);
				else if ( !strcmp(varname, "blockt") )
					sscanf (variable, "%lu", &sim.blockt);
				else if ( !strcmp(varname, "drift_freq") )
					sscanf (variable, "%lu", &sim.drift);
				else if ( !strcmp(varname, "timestep") )
					sscanf (variable, "%lf", &sim.dt);
				else if ( !strcmp(varname, "mts_steps") )
					sscanf (variable, "%i", &sim.nsteps);
				else if ( !strcmp(varname, "cutoff") )
					sscanf (variable, "%lf", &sim.rc);
				else if ( !strcmp(varname, "cutoff_nblist") )
					sscanf (variable, "%lf", &sim.rclist);
				else if ( !strcmp(varname, "cutoff_ewald") )
					sscanf (variable, "%lf", &sim.rc_ew);
				else if ( !strcmp(varname, "cutoff_ewald_nblist") )
					sscanf (variable, "%lf", &sim.rclist_ew);
				else if ( !strcmp(varname, "nblist_freq") )
					sscanf (variable, "%lu", &sim.nlist);
				else if ( !strcmp(varname, "nboxes") )
					sscanf (variable, "%i", &sim.NB);
				else if ( !strcmp(varname, "dielectric") )
					sscanf (variable, "%lf", &epsRF_temp);
				else
				{
					fprintf(stderr, "Unrecognized sim_variable name: %s\n", varname);
				}
			}
			else if( strstr(varname_full, "box_") == varname_full )	// The string begins with "sim"
			{
				//Ok
			}
			else 
			{
				fprintf(stderr, "Unrecognized command: %s\n", varname_full);
			}
		}
	}
	fclose(das);

	if(sim.ID == 0)
	{
		fprintf(stderr, "ERROR: You must supply sim_mode in simul.input.\n");
		exit(2);
	}
	if(sim.ID2 == 0)
	{
		fprintf(stderr, "ERROR: You must supply sim_intmode in simul.input.\n");
		exit(2);
	}
	if( !strcmp(param_file_name, "") )
	{
		fprintf(stderr, "ERROR: You must supply a parameter file name.\n");
		exit(2);
	}
	if(sim.dt == -1)
	{
		fprintf(stderr, "ERROR: You must supply a timestep.\n");
		exit(2);
	}

	if(sim.cyc_eq == 0)
		printf("WARNING: No equilibration steps will be performed.\n");
	if(sim.cyc_pr == 0)
		printf("WARNING: No production steps will be performed.\n");
	if(sim.blockd == ~0)
		printf("WARNING: sim_blockd not specified.\n");
	if(sim.blockc == ~0)
		printf("WARNING: sim_blockc not specified.\n");
	if(sim.blockt == ~0)
		printf("WARNING: sim_blockt not specified.\n");
	if(sim.drift == ~0)
		printf("WARNING: sim_drift not specified.\n");
		
	if(sim.NB == 0)
	{
		fprintf(stdout, "INFO: sim_nboxes was not supplied in simul.input; assuming = 1\n");
		sim.NB = 1;
	}

	// Default variables (some will cast warnings if not set)
/*
	sim_hyb.cyc_hybmd = 0;
	sim_hyb.cyc_swap = 0;
	sim.nsteps = 1;
	sim.rc = 0;
	sim.rclist = 0;
	sim.rc_ew = 0;
	sim.rclist_ew = 0;
	sim.nlist = 0;
	epsRF = 0;
*/

	/* =================================================== */
	/*     BOX VARIABLES                                   */
	/* =================================================== */

	for(int k=0; k<NBOXES; k++) sim.epsRF[k] = epsRF_temp;	

	int k;

	#ifdef MPI
		// Open the repex#.input file.
		char repex_name[100];
		sprintf(repex_name,"./INPUT/repex%d.input",mpi.my_rank);
		if(NULL == (das=fopen(repex_name,"r")) )
		{
			fprintf(stdout,"input file repex%d.input does not exist!!!\n",mpi.my_rank);
			exit(1);
		}
	#else
		// Just reopen the same file to read the box information
		
		k=0;
		if(NULL == (das=fopen("./INPUT/simul.input","r")) )
		{
			fprintf(stdout,"input file simul.input does not exist!!!\n");
			exit(1);
		}
	#endif

	// Set defaults here
	for(k=0; k<sim.NB; k++)		// Action if user does not specify:
	{
		box[k].lx = 0;		// Will causer ERROR and abort
		box[k].ly = 0;
		box[k].lz = 0;

		nhc[k].Nc = 1;		// Cannot be changed
		nhc[k].Nys = 1;		// Cannot be changed
		sim.T[k] = -1;		// Will cause ERROR and abort
		sim.dtemp[k] = 0;	
		sim.Ttau[k] = 0;	// Will cause warning and set to default: 0.2
		nhc[k].M = 0;		// Will cause warning and set to default: 8
		nhc[k].Q = 0;		// Will cause warning and set to default: 1.0e-18
		for(int d = 0; d<3; d++)
		{
			sim.P[k][d] = -1;	// Will cause warning and set to default: 0
			sim.dpress[k][d] = 0;	
			sim.Ptau[k][d] = 0;	// Will cause warning and set to default: 0.5
			sim.Pcomp[k][d] = 0;	// Will cause warning and set to default: 1.0e-6
			sim.kmax[k][d] = 0;	// Will cause warning and set to default: 6
			sim.grid[d] = 0;	// Will cause warning and set to default: 35
		}

		Wg[k] = 0;		// Will cause warning and set to default: 7.813381e-22
		sim.xRESPA[k] = 0;	
		sim.kappa[k] = 0;	
		sim.order = 0;	// Will cause warning and set to default: 6

		sim.Na_con[k] = -1;	// Will be set to 0 with a warning if DRSAM is defined
		
	}

	

	k = -1;	//Force user to supply a box number
	while ( !feof(das) )
	{
		// Read a single line from the file stream
		fgets(buff, 256, das);

		// Discard any whitespace
		int c=0;
		while(isspace(buff[c])) c++;
		memcpy(buff, buff+c, 256-c);		

		// Discard text appearing after the first # character

		char *pound = strchr(buff, '#');	//Return pointer to first '#'
		if (pound) *pound = '\0';		//Terminate string here
	
		if (strlen(buff) > 2)	//Sometimes junk stays behind--skip this
		{
			// Read in the name of this variable (until we hit whitespace) and convert to lowercase
			char varname_full[64];
			sscanf(buff, "%s", varname_full);
			for(size_t c=0; c < strlen(varname_full); c++)
				if ('A' <= varname_full[c] && varname_full[c] <= 'Z') varname_full[c] += 'a' - 'A';

			// Is this a box variable?
			if( strstr(varname_full, "box_") == varname_full )	// The string begins with "box_"
			{
				//Here's a little hack to extract the rest of the string (the variable value)
				char *variable = buff + strlen(varname_full);
				char *varname = varname_full + 4;


				if      ( !strcmp(varname, "number") )
				{
					sscanf (variable, "%i", &k);
					if(k >= sim.NB)
					{
						fprintf(stderr, "ERROR: Invalid box number: %i. Box numbers begin at 0, with maximum box number = nboxes-1 = %i.\n", k, sim.NB-1);
						exit(2);
					}
				}
				else
				{
					if(k == -1)
					{
						fprintf(stderr, "ERROR: You must supply a box number before any other box variables!\n");
						exit(2);
					}
					else if ( !strcmp(varname, "size_xyz") )
						sscanf (variable, "%lf%lf%lf", &box[k].lx, &box[k].ly, &box[k].lz);
					else if ( !strcmp(varname, "temp") )
					{
						sscanf (variable, "%lf", &sim.T[k]);
						sim.kT[k] = RG * sim.T[k]*.001; // This is modified to get the units in KJ/mol
					}
					else if ( !strcmp(varname, "dtemp") )
						sscanf (variable, "%lf", &sim.dtemp[k]);
					else if ( !strcmp(varname, "ber_tstat_coupling") )
						sscanf (variable, "%lf", &sim.Ttau[k]);
					else if ( !strcmp(varname, "nh_tstat_n") )
						sscanf (variable, "%i", &nhc[k].M);
					else if ( !strcmp(varname, "nh_tstat_mass") )
						sscanf (variable, "%lf", &nhc[k].Q);
					else if ( !strcmp(varname, "press_xyz") )
						sscanf (variable, "%lf%lf%lf", &sim.P[k][0],&sim.P[k][1],&sim.P[k][2]);
					else if ( !strcmp(varname, "dpress_xyz") )
						sscanf (variable, "%lf%lf%lf", &sim.dpress[k][0],&sim.dpress[k][1],&sim.dpress[k][2]);
					else if ( !strcmp(varname, "ber_bstat_coupling_xyz") )
						sscanf (variable, "%lf%lf%lf", &sim.Ptau[k][0],&sim.Ptau[k][1],&sim.Ptau[k][2]);
					else if ( !strcmp(varname, "isothermal_comp_xyz") )
						sscanf (variable, "%lf%lf%lf", &sim.Pcomp[k][0],&sim.Pcomp[k][1],&sim.Pcomp[k][2]);
					else if ( !strcmp(varname, "pr_bstat_mass") )
						sscanf (variable, "%lf", &Wg[k]);
					else if ( !strcmp(varname, "respa") )
						sscanf (variable, "%i", &sim.xRESPA[k]);
					else if ( !strcmp(varname, "kappa") )
					{
						#ifdef DRSAM
							printf("INFO: DRSAM is defined, therefore the explicit value specified for kappa in simul.input or repex#.input in box %i will be ignored.\n", k);
						#else
							sscanf (variable, "%lf", &sim.kappa[k]);
						#endif
					}
					else if ( !strcmp(varname, "kspace_vectors_xyz") )
						sscanf (variable, "%i%i%i",&sim.kmax[k][0],&sim.kmax[k][1],&sim.kmax[k][2]);
					else if ( !strcmp(varname, "spme_spline_order") )
						sscanf (variable, "%i",&sim.order);
					else if ( !strcmp(varname, "grid_pts_xyz") )
						sscanf (variable, "%i%i%i",&sim.grid[0],&sim.grid[1],&sim.grid[2]);
					else if ( !strcmp(varname, "na_conc") )
					{
						#ifdef DRSAM
						if(sim.T[k] == -1)
						{
							printf("ERROR: DRSAM is defined, therefore in simul.input or repex#.input in box %i, you must specify temp before Na_conc!\n", k);
							exit(2);
						}
						sscanf (variable, "%lf", &sim.Na_con[k]);
						sim.epsRF[k] =(249.4-0.788*sim.T[k]+7.2E-4*sim.T[k]*sim.T[k]) * (1.0-0.2551*sim.Na_con[k]+0.05151*sim.Na_con[k]*sim.Na_con[k]-6.889E-3*sim.Na_con[k]*sim.Na_con[k]*sim.Na_con[k]);
						if(sim.Na_con[k]==0)
							sim.kappa[k] = 78;
						else
							sim.kappa[k] = sqrt((sim.epsRF[k] * sim.kT[k]) / (8 * PI * ELEC * NA * NA * ELEQ * ELEQ * sim.Na_con[k]) * 1.0E30);

						#else
						printf("INFO: DRSAM is NOT defined, therefore the value specified for Na_conc in simul.input or repex#.input for box %i will be ignored.\n", k);
						#endif
					}
					else if ( !strcmp(varname, "dielectric") )
					{
						#ifdef DRSAM
							printf("INFO: DRSAM is defined, therefore the explicit value for dielectric constant specified in simul.input or repex#.input for box %i will be ignored.\n", k);
						#else
							sscanf (variable, "%lf",&sim.epsRF[k]);
						#endif
					}
					else
					{
						fprintf(stderr, "Unrecognized box_variable name in box %i: %s\n", k, varname);
					}
				}
			}
			else if( strstr(varname_full, "sim_") == varname_full )	// The string begins with "sim"
			{
				//Ok
			}
			else 
			{
				fprintf(stderr, "Unrecognized command: %s\n", varname_full);
			}
		}
	}
	fclose(das);

	// Check inputs here
	for(k=0; k<sim.NB; k++)		// Action if user does not specify:
	{
		if(box[k].lx == 0)   { fprintf(stderr, "ERROR: You must specify box_size in box %i!\n", k); exit(2); }
		if(box[k].ly == 0)   { fprintf(stderr, "ERROR: You must specify box_size in box %i!\n", k); exit(2); }
		if(box[k].lz == 0)   { fprintf(stderr, "ERROR: You must specify box_size in box %i!\n", k); exit(2); }
		if(sim.T[k] == -1)   { fprintf(stderr, "ERROR: You must specify a temperature in box %i!\n", k); exit(2); }
		if(sim.Ttau[k] == 0) { fprintf(stdout, "WARNING: ber_Tstat_coupling not specified in box %i; default = 0.2\n", k); sim.Ttau[k] = 0.2; }
		if(nhc[k].M == 0)    { fprintf(stdout, "WARNING: nh_Tstat_N not specified in box %i; default = 8\n", k); nhc[k].M = 8; }
		if(nhc[k].Q == 0)    { fprintf(stdout, "WARNING: nh_Tstat_mass not specified in box %i, default = 1.0e-18\n", k); nhc[k].Q = 1.0e-18; }
		for(int d = 0; d<3; d++)
		{
			if(sim.P[k][d] == -1)    {printf("WARNING: press[%i] not specified in box %i; default = 0\n", d, k); sim.P[k][d] = 0; }
			if(sim.Ptau[k][d] == 0)  {printf("WARNING: ber_Bstat_coupling_time[%i] not specified in box %i; default = 0.5\n", d, k); sim.Ptau[k][d] = 0.5; }
			if(sim.Pcomp[k][d] == 0) {printf("WARNING: isothermal_comp[%i] not specified in box %i; default = 1.0e-6\n", d, k); sim.Pcomp[k][d] = 1.0e-6; }
			if(sim.kmax[k][d] == 0)  {printf("WARNING: kspace_vectors[%i] not specified in box %i; default = 6\n", d, k); sim.kmax[k][d] = 6; }
			if(sim.grid[d] == 0)  {printf("WARNING: grid_pts[%i] not specified in box %i; default = 35\n", d, k); sim.grid[d] = 35; }
		}

		if(Wg[k] == 0) { fprintf(stdout, "WARNING: pr_Bstat_mass not specified in box %i; default = 7.813381e-22\n", k); Wg[k] = 7.813381e-22; }
		if(sim.order == 0)    { fprintf(stdout, "WARNING: spme_spline_order not specified in box %i, default = 6\n", k); sim.order = 6; }
		#ifdef DRSAM
		if(sim.Na_con[k] == -1)
		{
			fprintf(stdout, "WARNING: DRSAM is defined but Na_conc was not specified in box %i; default = 0\n", k);
			sim.Na_con[k] = 0;
			sim.epsRF[k] =(249.4-0.788*sim.T[k]+7.2E-4*sim.T[k]*sim.T[k]) * (1.0-0.2551*sim.Na_con[k]+0.05151*sim.Na_con[k]*sim.Na_con[k]-6.889E-3*sim.Na_con[k]*sim.Na_con[k]*sim.Na_con[k]);
			if(sim.Na_con[k]==0)
				sim.kappa[k] = 78;
			else
				sim.kappa[k] = sqrt((sim.epsRF[k] * sim.kT[k]) / (8 * PI * ELEC * NA * NA * ELEQ * ELEQ * sim.Na_con[k]) * 1.0E30);

		}
		#endif
	}
	for(int k=0; k<sim.NB; k++)
	{
		sim.max_k_square[k] = sim.kmax[k][0]*sim.kmax[k][0];
	}
	#else  //SCRIPT

	fgets(tt,150,das);
	fgets(tt,150,das);
	fgets(tt,150,das);
	fgets(tt,150,das);
	fgets(tt,150,das);
	fgets(tt,150,das);
	fgets(tt,150,das);
	fgets(tt,150,das);
	fgets(tt,150,das);
	fgets(tt,150,das);
	fgets(tt,150,das);
	fgets(tt,150,das);
	fgets(title,150,das);
	int len = strlen(title); title[len-2] = '\0';
	fscanf(das,"Parameter File Name=%s",param_file_name);
	fgets(tt,150,das);
	fscanf(das,"%d", &sim.ID);				fgets(tt,150,das);
	fscanf(das,"%d", &sim.ID2);				fgets(tt,150,das);
	fscanf(das,"%lu",&sim.cyc_eq);			fgets(tt,150,das);
	fscanf(das,"%lu",&sim.cyc_pr);			fgets(tt,150,das);
	fscanf(das,"%ld",&sim_hyb.cyc_hybmd);	fgets(tt,150,das);
	fscanf(das,"%ld",&sim_hyb.cyc_swap);	fgets(tt,150,das);
	fscanf(das,"%ld",&sim.blockd);			fgets(tt,150,das);
	fscanf(das,"%ld	%ld",&sim.blockc,&sim.blockt);	fgets(tt,150,das);
	fscanf(das,"%ld",&sim.drift);			fgets(tt,150,das);
	fscanf(das,"%lf",&sim.dt);				fgets(tt,150,das);
	fscanf(das,"%d", &sim.nsteps);			fgets(tt,150,das);
	fscanf(das,"%lf",&sim.rc);			
	fscanf(das,"%lf",&sim.rclist);			fgets(tt,150,das);
	fscanf(das,"%lf",&sim.rc_ew);			
	fscanf(das,"%lf",&sim.rclist_ew);		fgets(tt,150,das);
	fscanf(das,"%lu", &sim.nlist);			fgets(tt,150,das);
	fscanf(das,"%d", &sim.NB);				fgets(tt,150,das);
	#ifndef DRSAM
		double epsRF_temp;
		fscanf(das,"%lf",&epsRF_temp);				fgets(tt,150,das);
		for(int k=0; k<sim.NB; k++) sim.epsRF[k] = epsRF_temp;
	#endif

	fgets(tt,150,das);


	/* ================================================================== */
	/* Loop around number of boxes to get data on box properties.         */
	/* ================================================================== */
	#ifndef MPI
	for (int k=0; k<sim.NB; k++)
	{
		#else
		fclose(das);
		char repex_name[100];
		sprintf(repex_name,"./INPUT/repex%d.input",mpi.my_rank);
		if(NULL == (das=fopen(repex_name,"r")) )
		{
			fprintf(stdout,"input file repex%d.input does not exist!!!\n",mpi.my_rank);
			exit(1);
		}
		int k=0;
		#endif

		fscanf(das, "%d", &boxnum);			fgets(tt,150,das);
		if (k != boxnum)
		{
			fprintf(stdout,"ERROR: Box Number (%i) does not match in simul.input!!! (should be %i) \n", boxnum, k);
			exit(1);
		}
		fscanf(das,"%lf %lf %lf",&box[k].lx,&box[k].ly,&box[k].lz); fgets(tt,80,das);
		fscanf(das,"%lf",&sim.T[k]);		fgets(tt,150,das);
		fscanf(das,"%lf",&sim.dtemp[k]);	fgets(tt,150,das);
		fscanf(das,"%lf",&sim.Ttau[k]);		fgets(tt,150,das);
		fscanf(das,"%d",&nhc[k].M);			fgets(tt,150,das);
		nhc[k].Nc=1;
		nhc[k].Nys=1;
		fscanf(das,"%lf",&nhc[k].Q);		fgets(tt,150,das);
		fscanf(das,"%lf %lf %lf",&sim.P[k][0],&sim.P[k][1],&sim.P[k][2]);					fgets(tt,150,das);
		fscanf(das,"%lf %lf %lf",&sim.dpress[k][0],&sim.dpress[k][1],&sim.dpress[k][2]);	fgets(tt,150,das);
		fscanf(das,"%lf %lf %lf",&sim.Ptau[k][0],&sim.Ptau[k][1],&sim.Ptau[k][2]);			fgets(tt,150,das);
		fscanf(das,"%lf %lf %lf",&sim.Pcomp[k][0],&sim.Pcomp[k][1],&sim.Pcomp[k][2]);		fgets(tt,150,das);
		fscanf(das,"%lf",&Wg[k]);			fgets(tt,150,das);
		fscanf(das,"%d",&sim.xRESPA[k]);	fgets(tt,150,das);
		fscanf(das,"%lf",&sim.kappa[k]);	fgets(tt,150,das);
		fscanf(das,"%d	%d	%d",&sim.kmax[k][0],&sim.kmax[k][1],&sim.kmax[k][2]); 			fgets(tt,150,das);
		sim.max_k_square[k] = sim.kmax[k][0]*sim.kmax[k][0];
		fscanf(das,"%d",&sim.order);		fgets(tt,150,das);
		fscanf(das,"%d",&sim.grid[0]);
		fscanf(das,"%d",&sim.grid[1]);
		fscanf(das,"%d",&sim.grid[2]);		fgets(tt,150,das);
		fscanf(das,"%lf",&sim.epsRF[k]);
		fgets(tt,150,das);

		sim.kT[k] = RG * sim.T[k]*.001; // This is modified to get the units in KJ/mol

                #ifdef DRSAM
			sim.Na_con[k] = sim.kappa[k];
			sim.epsRF[k] =(249.4-0.788*sim.T[k]+7.2E-4*sim.T[k]*sim.T[k]) * (1.0-0.2551*sim.Na_con[k]+0.05151*sim.Na_con[k]*sim.Na_con[k]-6.889E-3*sim.Na_con[k]*sim.Na_con[k]*sim.Na_con[k]);

			if(sim.Na_con[k]==0)
			{
				sim.kappa[k] = 78;
			}       
			else
			{   
				sim.kappa[k] = sqrt((sim.epsRF[k] * sim.kT[k]) / (8 * PI * ELEC * NA * NA * ELEQ * ELEQ * sim.Na_con[k]) * 1.0E30);
			}       
                #endif


	#ifndef MPI
	}//for k
	#endif

	fclose(das);

	#endif // #ifdef SCRIPT


/* ======================================================================= */
/*                                                                         */
/* END READING simul.input AND repex#.input                                */
/*                                                                         */
/* ======================================================================= */



	#ifdef RMSD
	FILE *rmsd_ptr;
	int nat_r;		     /* number of all atoms in reference molecule */
	char temp[5];
	char temp1[5]="END";
	char temp2[7]="REMARK";

	if( NULL == (rmsd_ptr=fopen("./INPUT/ref.pdb","r")) )
	{		// open reference file
		fprintf(stdout,"input file ref.pdb does not exist!!!\n");
		exit(10);
	} 

	int i =0;	

	while (!feof(rmsd_ptr))
	{
		fscanf(rmsd_ptr, "%s",   temp);
		if (strcmp(temp,temp1)==0 || strcmp(temp,temp2)==0) {
		fgets(tt,150,rmsd_ptr);
			continue;
		}
		fscanf(rmsd_ptr, "%d",   &i); nat_r =i;
		if(nat_r >= MAXPOINTS) {
			fprintf(stderr,"Error: Molecule too big. Recompile program with larger MAXPOINTS (at least %i)\n", nat_r+1);
			exit(12);
		}
		fscanf(rmsd_ptr,"%4s",  ref[i].name);
		fscanf(rmsd_ptr,"%3s",  ref[i].res);
		fscanf(rmsd_ptr,"%d",   &ref[i].nres);
		fscanf(rmsd_ptr,"%lf",  &ref[i].x);
		fscanf(rmsd_ptr,"%lf",  &ref[i].y);
		fscanf(rmsd_ptr,"%lf",  &ref[i].z);
		fgets(tt,150,rmsd_ptr);
	}
	refN=nat_r;
	fclose(rmsd_ptr);

	/* Read in the excluded residues */
	count_excl=0;

	FILE *exc1;
	if( NULL == (exc1=fopen("./INPUT/rmsd.inp","r")) )
	{
		fprintf(stdout,"input file rmsd.inp does not exist!!!\n");
		exit(10);
	} 
	fgets(tt,150,exc1);    
	fgets(tt,150,exc1);    
	fgets(tt,150,exc1);    
	fgets(tt,150,exc1);    
	fgets(tt,150,exc1);    
	fgets(tt,150,exc1);    
	fgets(tt,150,exc1);    
	fgets(tt,150,exc1);    
	fgets(tt,150,exc1);     
	fscanf(exc1,"%d %d",&rmsd_modes.file,&rmsd_modes.atoms); fgets(tt,150,exc1);
	while (!feof(exc1))
	{
		count_excl++;
		fscanf(exc1, "%d\n",&res_excl[count_excl]);
		res_excl[count_excl]=res_excl[count_excl]-1;
	}

	fclose(exc1);
	#endif
  
  
	/*----------------------------------*/
	/* If DOS is defined read in the	*/
	/* DOS parameters.					*/
	/*----------------------------------*/

#ifdef DOS
	char name3[100];
	FILE *dos1;
	for(int k=0; k<sim.NB; k++)
	{
		#ifdef MPI
			sprintf(name3, "./INPUT/dos%d.input",mpi.my_rank);
		#endif
		#ifndef MPI
			sprintf(name3, "./INPUT/dos%d.input",k);
		#endif
		dos1 = fopen(name3,"r");
		if( NULL == dos1)
		{
			fprintf(stdout,"input file dos.input (%s) does not exist!!!\n", name3);
			exit(1);
		}
		fgets(tt,150,dos1);
		fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].e_begin);		    fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].e_end);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].e_width);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].T_begin);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].T_end);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].flat);			fgets(tt,150,dos1);
		fscanf(dos1,"%f",&MERGE_F);				fgets(tt,150,dos1);
		fscanf(dos1,"%d",&MERGE_N);				fgets(tt,150,dos1);
		fscanf(dos1,"%f",&RESET_F);				fgets(tt,150,dos1);
		fscanf(dos1,"%f",&STOP_F);				fgets(tt,150,dos1);
		fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_pivot.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_pivot.theta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_trans.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_trans.delta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_axial.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_axial.strain[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_rand.PMS[k]);				fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_rand.delta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_solv.PMS[k]);				fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_solv.delta[k]);			fgets(tt,150,dos1);
		fclose (dos1);
	}// kloop ends here
	#endif//DOS

	/*----------------------------------*/
	/* If MMDOS is defined read in the	*/
	/* MMDOS parameters.				*/
	/*----------------------------------*/

	#ifdef MMDOS
	char name3[100];
	FILE *dos1;
	for(int k=0; k<sim.NB; k++)
	{
		#ifdef MPI
			sprintf(name3, "./INPUT/dos%d.input",mpi.my_rank);
		#endif
		#ifndef MPI
			sprintf(name3, "./INPUT/dos%d.input",k);
		#endif
		dos1 = fopen(name3,"r");
		if( NULL == dos1)
		{
			fprintf(stdout,"input file dos.input (%s) does not exist!!!\n", name3);
			exit(1);
		}
		fgets(tt,150,dos1);
		fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].e_begin);		    fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].e_end);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].e_width);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].T_begin);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].T_end);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_pivot.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_pivot.theta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_trans.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_trans.delta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].p_dke);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].d_max);			fgets(tt,150,dos1);
		fclose (dos1);
	}// kloop ends here
	#endif//MMDOS

	#ifdef CTDOS
	char name3[100];
	FILE *dos1;
	for(int k=0; k<sim.NB; k++)
	{
		#ifdef MPI
		sprintf(name3, "./INPUT/dos%d.input",mpi.my_rank);
		#endif
		#ifndef MPI
		sprintf(name3, "./INPUT/dos%d.input",k);
		#endif
		dos1 = fopen(name3,"r");
		if( NULL == dos1)
		{
			fprintf(stdout,"input file dos.input (%s) does not exist!!!\n", name3);
			exit(1);
		}
		fgets(tt,150,dos1);
		fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].e_begin);		    fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].e_end);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].e_width);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].T_begin);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].T_end);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].flat);			fgets(tt,150,dos1);
		fscanf(dos1,"%f",&MERGE_F);				fgets(tt,150,dos1);
		fscanf(dos1,"%d",&MERGE_N);				fgets(tt,150,dos1);
		fscanf(dos1,"%f",&RESET_F);				fgets(tt,150,dos1);
		fscanf(dos1,"%f",&STOP_F);				fgets(tt,150,dos1);
		fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_pivot.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_pivot.theta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_trans.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_trans.delta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_axial.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_axial.strain[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_rand.PMS[k]);				fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_rand.delta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_solv.PMS[k]);				fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_solv.delta[k]);			fgets(tt,150,dos1);
		fclose (dos1);
	}// kloop ends here
	#endif//CTDOS

	#ifdef XEDOS
	char name3[100];
	FILE *dos1;
	for(int k=0; k<sim.NB; k++)
	{
		#ifdef MPI
		sprintf(name3, "./INPUT/dos%d.input",mpi.my_rank);
		#endif
		#ifndef MPI
		sprintf(name3, "./INPUT/dos%d.input",k);
		#endif
		dos1 = fopen(name3,"r");
		if( NULL == dos1)
		{
			fprintf(stdout,"input file dos.input (%s) does not exist!!!\n", name3);
			exit(1);
		}
		fgets(tt,150,dos1);
		fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].l_begin);		    fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].l_end);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].l_width);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].T_begin);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].T_end);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].flat);			fgets(tt,150,dos1);
		fscanf(dos1,"%f",&MERGE_F);				fgets(tt,150,dos1);
		fscanf(dos1,"%d",&MERGE_N);				fgets(tt,150,dos1);
		fscanf(dos1,"%f",&RESET_F);				fgets(tt,150,dos1);
		fscanf(dos1,"%f",&STOP_F);				fgets(tt,150,dos1);
		fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_pivot.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_pivot.theta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_trans.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_trans.delta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_axial.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_axial.strain[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_rand.PMS[k]);				fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_rand.delta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_solv.PMS[k]);				fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_solv.delta[k]);			fgets(tt,150,dos1);
		fclose (dos1);
	}// kloop ends here
	#endif//XEDOS

	/*----------------------------------*/
	/* If FX_EDOS is defined read in 	*/
	/* FX_EDOS parameters.				*/
	/*----------------------------------*/

	#ifdef FX_EDOS
	char name3[100];
	FILE *dos1;
	for(int k=0; k<sim.NB; k++)
	{
		#ifdef MPI
		sprintf(name3, "./INPUT/fxdos%d.input",mpi.my_rank);
		#endif
		#ifndef MPI
		sprintf(name3, "./INPUT/fxdos%d.input",k);
		#endif
		dos1 = fopen(name3,"r");
		if( NULL == dos1)
		{
			fprintf(stdout,"input file fxdos.input (%s) does not exist!!!\n", name3);
			exit(1);
		}
		fgets(tt,150,dos1);
		fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].e_begin);		    fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].e_end);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].e_width);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].T_begin);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].T_end);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].r_target);		fgets(tt,150,dos1);
		fscanf(dos1,"%d", &sim_dos[k].trips_target);	fgets(tt,150,dos1);
		fscanf(dos1,"%f", &MERGE_F);					fgets(tt,150,dos1);
		fscanf(dos1,"%d", &MERGE_N);					fgets(tt,150,dos1);
		fscanf(dos1,"%f", &RESET_F);					fgets(tt,150,dos1);
		fscanf(dos1,"%f", &STOP_F);						fgets(tt,150,dos1);
														fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_pivot.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_pivot.theta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_trans.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_trans.delta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_axial.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_axial.strain[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_rand.PMS[k]);				fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_rand.delta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_solv.PMS[k]);				fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_solv.delta[k]);			fgets(tt,150,dos1);
		fclose (dos1);
	}// kloop ends here
	#endif//FX_EDOS

	
	#ifdef TDXEDOS
	char name3[100];
	FILE *dos1;
	for(int k=0; k<sim.NB; k++)
	{
		#ifdef MPI
		sprintf(name3, "./INPUT/tdxedos%d.input",mpi.my_rank);
		#endif
		#ifndef MPI
		sprintf(name3, "./INPUT/tdxedos%d.input",k);
		#endif
		dos1 = fopen(name3,"r");
		if( NULL == dos1)
		{
			fprintf(stdout,"input file tdxedos.input (%s) does not exist!!!\n", name3);
			exit(1);
		}
		fgets(tt,150,dos1);
		fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].x1_begin);		fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].x1_end);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].x1_width);		fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].x2_begin);		fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].x2_end);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].x2_width);		fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].T_begin);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].T_end);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&sim_dos[k].flat);			fgets(tt,150,dos1);
		fscanf(dos1,"%f",&MERGE_F);				fgets(tt,150,dos1);
		fscanf(dos1,"%d",&MERGE_N);				fgets(tt,150,dos1);
		fscanf(dos1,"%f",&RESET_F);				fgets(tt,150,dos1);
		fscanf(dos1,"%f",&STOP_F);				fgets(tt,150,dos1);
		fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_pivot.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_pivot.theta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_trans.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_trans.delta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_axial.PMS[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_axial.strain[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_rand.PMS[k]);				fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_rand.delta[k]);			fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_solv.PMS[k]);				fgets(tt,150,dos1);
		fscanf(dos1,"%lf",&mc_solv.delta[k]);			fgets(tt,150,dos1);
		fclose (dos1);
		sim_dos[k].x2_begin *= PI/180.0;
		sim_dos[k].x2_end   *= PI/180.0;
		sim_dos[k].x2_width *= PI/180.0;

	}// kloop ends here
	#endif//TDXEDOS
	
	/* ================================================================== */
	/* Read in information about the native contacts.                     */
	/* ================================================================== */
	char con_filename[100];
	FILE *con_ptr;
	sprintf(con_filename, "./INPUT/contacts.inp");
 	/* ---------------------------------------------------------- */
	/* Allocate memory for restraint sturcture and read in the    */
	/* parameters of each box from restrain.input.                */
	/* ---------------------------------------------------------- */
	contacts = (struct struct_contacts*) calloc(sim.NB,sizeof(struct struct_contacts));
	if(contacts == NULL)
	{
		fprintf(stdout, "ERROR: cannot allocate memory for contacts\n");
		exit(11);
	}
 
	if( NULL != (con_ptr=fopen(con_filename,"r")) )
	{

		fgets(tt,150,con_ptr);
		fgets(tt,150,con_ptr);
		fgets(tt,150,con_ptr);
		
		for(int k=0; k<sim.NB; k++)
		{
			fscanf(con_ptr, "%d", &boxnum);		fgets(tt,150,con_ptr);
			if (k != boxnum)
			{
				fprintf(stdout,"ERROR: Box Number (%i) does not match in contacts.input!!! (should be %i)\n", boxnum, k);
				exit(1);
			}
			fscanf(con_ptr, "%d", &contacts[k].n);		fgets(tt,150,con_ptr);
			contacts[k].a = (int*) calloc(contacts[k].n,sizeof(int));
				if(contacts[k].a == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for contacts[k].a\n"); exit(11);}
			contacts[k].b = (int*) calloc(contacts[k].n,sizeof(int));
				if(contacts[k].b == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for contacts[k].b\n"); exit(11);}
			contacts[k].d = (double*) calloc(contacts[k].n,sizeof(double));
				if(contacts[k].d == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for contacts[k].d\n"); exit(11);}
			contacts[k].code = (int*) calloc(contacts[k].n,sizeof(int));
				if(contacts[k].code == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for contacts[k].code\n"); exit(11);}

			for(int i=0; i<contacts[k].n; i++)
			{
				fscanf(con_ptr, "%d",  &contacts[k].a[i]);
				fscanf(con_ptr, "%d",  &contacts[k].b[i]);
				fscanf(con_ptr, "%lf", &contacts[k].d[i]); fgets(tt,150,con_ptr);
			
				//subtract off one from the site number so they begin at 0
				contacts[k].a[i] -= 1;
				contacts[k].b[i] -= 1;
			
			}
			fgets(tt,150,con_ptr);
		}//for k
	  
		fclose(con_ptr); 
	}//if( NULL !=...)
	else
	{
		for(int k=0; k<sim.NB; k++)
			contacts[k].n=0;
	}

	/* ================================================================== */
	/* Read in information about the native contacts. 2                   */
	/* ================================================================== */
	
	FILE *con_ptr_2;
	sprintf(con_filename, "./INPUT/contacts_2.inp");
 	/* ---------------------------------------------------------- */
	/* Allocate memory for contacts 2 sturcture and read in the   */
	/* parameters of each box from contacts_2.input.              */
	/* ---------------------------------------------------------- */
	contacts_2 = (struct struct_contacts*) calloc(sim.NB,sizeof(struct struct_contacts));
		if(contacts_2 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for contacts_2\n"); exit(11);}
 
	if( NULL != (con_ptr_2=fopen(con_filename,"r")) )
	{

		fgets(tt,150,con_ptr_2);
		fgets(tt,150,con_ptr_2);
		fgets(tt,150,con_ptr_2);
		
		for(int k=0; k<sim.NB; k++){
			fscanf(con_ptr_2, "%d", &boxnum);		fgets(tt,150,con_ptr_2);
			if (k != boxnum) {
				fprintf(stdout,"ERROR: Box Number (%i) does not match in contacts_2.input!!! (should be %i)\n", boxnum, k);
				exit(1);
			}
			fscanf(con_ptr_2, "%d", &contacts_2[k].n);		fgets(tt,150,con_ptr_2);
			contacts_2[k].a = (int*) calloc(contacts_2[k].n,sizeof(int));
				if(contacts_2[k].a == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for contacts_2[k].a\n"); exit(11);}
			contacts_2[k].b = (int*) calloc(contacts_2[k].n,sizeof(int));
				if(contacts_2[k].b == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for contacts_2[k].b\n"); exit(11);}
			contacts_2[k].d = (double*) calloc(contacts_2[k].n,sizeof(double));
				if(contacts_2[k].d == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for contacts_2[k].d\n"); exit(11);}
			contacts_2[k].code = (int*) calloc(contacts_2[k].n,sizeof(int));
				if(contacts_2[k].code == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for contacts_2[k].code\n"); exit(11);}

			for(int i=0; i<contacts_2[k].n; i++){
				fscanf(con_ptr_2, "%d",  &contacts_2[k].a[i]);
				fscanf(con_ptr_2, "%d",  &contacts_2[k].b[i]);
				fscanf(con_ptr_2, "%lf", &contacts_2[k].d[i]); fgets(tt,150,con_ptr_2);
			
				//subtract off one from the site number so they begin at 0
				contacts_2[k].a[i] -= 1;
				contacts_2[k].b[i] -= 1;
			
			}
			fgets(tt,150,con_ptr_2);
		}//for k
	  
		fclose(con_ptr_2); 
	}//if( NULL !=...)
	else
	{
		for(int k=0; k<sim.NB; k++)
			contacts_2[k].n=0;
	}

	/* ================================================================== */
	/* Read in information about the radius of gyrations                  */
	/* ================================================================== */
	char gyr_filename[100];
	FILE *gyr_ptr;
	sprintf(gyr_filename, "./INPUT/gyration.inp");
 	/* ---------------------------------------------------------- */
	/* Allocate memory											  */
	/* ---------------------------------------------------------- */

	gyration = (struct struct_gyr*) calloc(sim.NB,sizeof(struct struct_gyr));
	if(gyration == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for gyration\n"); exit(11);}
 
	if( NULL != (gyr_ptr=fopen(gyr_filename,"r")) )
	{

		fgets(tt,150,gyr_ptr);
		fgets(tt,150,gyr_ptr);
		fgets(tt,150,gyr_ptr);
		fgets(tt,150,gyr_ptr);
		
		for(int k=0; k<sim.NB; k++)
		{
			fscanf(gyr_ptr, "%d", &boxnum);		fgets(tt,150,gyr_ptr);
			if (k != boxnum) {
				fprintf(stdout,"ERROR: Box Number (%i) does not match in gyration.input!!! (should be %i)\n", boxnum, k);
				exit(1);
			}
			fscanf(gyr_ptr, "%d", &gyration[k].n);		fgets(tt,150,gyr_ptr);
			gyration[k].sites = (int*) calloc(gyration[k].n,sizeof(int));

			if(gyration[k].sites == NULL){
				fprintf(stdout, "ERROR: cannot allocate memory for gyration[k].sites\n"); exit(11);
			}
			int count=0;
			for(int i=0; i<gyration[k].n; i++){
				fscanf(gyr_ptr, "%d",  &gyration[k].sites[i]);
				count++;
				if (count==10) {
					fscanf(gyr_ptr,"\n");
					count=0;
				}
	      else if (i==gyration[k].n-1) fscanf(gyr_ptr,"\n");
				gyration[k].sites[i] -= 1; //subtract off one from the site number so they begin at 0
			}
			//fgets(tt,150,gyr_ptr); //this line was wrong.  If you had a power of 10 number of gyration
				     //sites, the fgets actually read in the next box number rather than
				     //a line.  To fix it, I added the else if two lines up analagous
				     //to read_topology with bonds, bends, etc.
		}//for k
	  
	   fclose(gyr_ptr); 
	}//if( NULL !=...)
	else
	{
		for(int k=0; k<sim.NB; k++) gyration[k].n=0;
	}

	#ifdef SMD	
	char name7[100];
	FILE *das7;
	sprintf(name7, "./INPUT/smd.input");
	das7 = fopen(name7,"r");

	/* ================================================================== */
	/* Allocate memory for steermd sturcture and read in the              */
	/* parameters of each box from smd.input.		                        */
	/* ================================================================== */
	steermd = (struct smd_struct*) calloc(sim.NB,sizeof(struct smd_struct));
	if(steermd == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for steermd\n"); exit(11);}

	fgets(tt,150,das7);
	fgets(tt,150,das7);
	fgets(tt,150,das7);
	
	for(int k=0; k<sim.NB; k++)
	{
		fscanf(das7, "%d", &boxnum);		fgets(tt,150,das7);
		if (k != boxnum) {
			fprintf(stdout,"ERROR: Box Number (%i) does not match in smd.input!!! (should be %i) \n", boxnum, k);
			exit(1);
		}
		fscanf(das7,"%d", &steermd[k].site1);		fgets(tt,150,das7);
		fscanf(das7,"%d", &steermd[k].site2);		fgets(tt,150,das7);
		steermd[k].site1	-=	1;
		steermd[k].site2	-=	1;
		fscanf(das7,"%lf",&steermd[k].k1);			fgets(tt,150,das7);
		fscanf(das7,"%lf",&steermd[k].k2);			fgets(tt,150,das7);
		fscanf(das7,"%lf",&steermd[k].v);			fgets(tt,150,das7);
		fgets(tt,150,das7);

		/* Change the units on k1 and k2 to kJ/mol */
		steermd[k].k1 *= 4.184;
		steermd[k].k2 *= 4.184;

	}
	fclose(das7);
#endif //SMD
#ifdef NSMD	
	char name7[100];
	FILE *das7;
	sprintf(name7, "./INPUT/smd.input");
	das7 = fopen(name7,"r");

	/* ================================================================== */
	/* Allocate memory for steermd sturcture and read in the              */
	/* parameters of each box from smd.input.		                        */
	/* ================================================================== */
	steermd = (struct smd_struct*) calloc(sim.NB,sizeof(struct smd_struct));
	if(steermd == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for steermd\n"); exit(11);}

	fgets(tt,150,das7);
	fgets(tt,150,das7);
	fgets(tt,150,das7);
	
	for(int k=0; k<sim.NB; k++)
	{
		fscanf(das7, "%d", &boxnum);		fgets(tt,150,das7);
		if (k != boxnum) {
			fprintf(stdout,"ERROR: Box Number (%i) does not match in smd.input!!! (should be %i)\n", boxnum, k);
			exit(1);
		}
		fscanf(das7,"%d", &steermd[k].site1);		fgets(tt,150,das7);
		fscanf(das7,"%d", &steermd[k].site2);		fgets(tt,150,das7);
		steermd[k].site1	-=	1;
		steermd[k].site2	-=	1;
		fscanf(das7,"%lf",&steermd[k].k1);			fgets(tt,150,das7);
		fscanf(das7,"%lf",&steermd[k].k2);			fgets(tt,150,das7);
		fscanf(das7,"%lf",&steermd[k].v);			fgets(tt,150,das7);
		fgets(tt,150,das7);

		/* Change the units on k1 and k2 to kJ/mol */
		steermd[k].k1 *= 4.184;
		steermd[k].k2 *= 4.184;

	}
	fclose(das7);
#endif //NSMD
  /* ==================================================== */
  /* Read the mc move info if a Monte Carlo simulation    */
  /* is specified.                                        */
  /* ==================================================== */
	if(sim.ID == 4)
	{
		for(int k=0; k<sim.NB; k++)
		{
			char mc_fn[80];
			FILE* mc_fp;
			sprintf(mc_fn, "./INPUT/mc%d.input",k);
			if(NULL!=(mc_fp=fopen(mc_fn,"r")))
			{
				fgets(tt,150,mc_fp);
				fgets(tt,150,mc_fp);
				fscanf(mc_fp,"%lf",&mc_pivot.PMS[k]);     fgets(tt,150,mc_fp);
				fscanf(mc_fp,"%lf",&mc_pivot.theta[k]);   fgets(tt,150,mc_fp);
				fscanf(mc_fp,"%lf",&mc_trans.PMS[k]);     fgets(tt,150,mc_fp);
				fscanf(mc_fp,"%lf",&mc_trans.delta[k]);   fgets(tt,150,mc_fp);
				fscanf(mc_fp,"%lf",&mc_rand.PMS[k]);      fgets(tt,150,mc_fp);
				fscanf(mc_fp,"%lf",&mc_rand.delta[k]);    fgets(tt,150,mc_fp);
			}
			else
			{
				fprintf(stdout,"The file mc%d.input does not exist!\n",k);
				exit(1);
			}
			fclose(mc_fp);
		}//for k
	}//if(sim.ID==4)

  /* =================================================== */
  /* Check to make sure that blockd is not greater than  */
  /* production steps or errors won't work in outend().  */
  /* [vblock() will never be called.]                    */
  /* =================================================== */
  if((unsigned int)sim.blockd > sim.cyc_pr && sim.cyc_pr != 0) sim.blockd = sim.cyc_pr;
  


}  


