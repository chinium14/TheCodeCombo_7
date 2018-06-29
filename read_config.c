/* ======================================================================== */
/* read_config.cpp                                                          */
/*                                                              		    */
/*		This subroutine reads in the x,y,z coordinates of the sites from    */
/* the simul#.crd file.  It also does some consistancy checks to make sure  */
/* all the files contain the same number of molecules and sites.  After     */
/* the box information is read in, variables like the box volume, Nose      */
/* Hoover parameters, and SASA parameters are calculated.                   */
/*                                                                          */
/* Note: There is a separate simul#.crd file for each box.                  */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

  /* ================================================================== */
  /*                                                                    */
  /* Function Prototypes                                                */
  /*                                                                    */
  /* ================================================================== */
void allocm (void);
void pbc_init (int);
void read_topology(int);
#ifdef BGO
#ifdef SASA
void vdw_radii_bgo(void);
void read_solvent(void);
#endif
#endif
/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void read_config (void) 
{
  /* ================================================================== */
  /* Variables and pointers needed in the subroutine.                   */
  /* ================================================================== */
 int qw;
  char tt[80], namec[50];
  FILE *das;
  FILE *molinput_ptr;

  /* ================================================================== */
  /* Open and read mol.input.                                           */
  /* ================================================================== */
  if(NULL == (molinput_ptr=fopen("./INPUT/mol.input","r"))){
	  fprintf(stdout,"ERROR: Input file mol.input does not exist!!!\n");
	  exit(1);
  }

    fgets(tt,80,molinput_ptr);
    fgets(tt,80,molinput_ptr);

	fscanf(molinput_ptr,"%d", &sim.NC);		fgets(tt,80,molinput_ptr);
	fscanf(molinput_ptr,"%d", &Nmolall);	fgets(tt,80,molinput_ptr);

	fgets(tt,80,molinput_ptr);
  
  /* ================================================================== */
  /* Reset the counter for the total number of molecules i in all the   */
  /* boxes.  (This is used below as a consistancy check.)               */
  /* ================================================================== */
 // for (int i=0; i<sim.NC; i++) mol[i].Ni = 0;

  /* ----------------------------------------------- */
  /* Read in the number of each type of molecule in  */
  /* the box.                                        */
  /* ----------------------------------------------- */    
    for (int i=0; i<sim.NC; i++) {
      fscanf(molinput_ptr,"%d",&bp[0][i].nbox); fgets(tt,80,molinput_ptr);
    }
	for(int k=1; k<sim.NB; k++){
		for(int i=0; i<sim.NC; i++) bp[k][i].nbox = bp[0][i].nbox;
	}
    fgets(tt,80,molinput_ptr);


  /* ================================================================== */
  /* Allocate memory for mol sturcture and read in molecular weights    */
  /* of each component.                                                 */
  /* ================================================================== */
  mol = (struct mol_info*) calloc(sim.NC,sizeof(struct mol_info));
	if(mol == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for mol\n"); exit(11);}

  /* ================================================================== */
  /* Reset the counter for the total number of molecules i in all the   */
  /* boxes.  (This is used below as a consistancy check.)               */
  /* ================================================================== */
  for (int i=0; i<sim.NC; i++) mol[i].Ni = 0;


  /* ----------------------------------------------- */
  /* Read in the number of sites on one moluecule    */
  /* for each type  in of moleculethe box.           */
  /* ----------------------------------------------- */    
    for (int i=0; i<sim.NC; i++) {
      fscanf(molinput_ptr,"%d",&mol[i].Nsite);  fgets(tt,80,molinput_ptr);
    }
    fgets(tt,80,molinput_ptr);


  /* ----------------------------------------------- */
  /* Calculate the total number of molecules of each */
  /* component in all the boxes.					 */
  /* ----------------------------------------------- */    
    for (int k=0; k<sim.NB; k++){
		for (int i=0; i<sim.NC; i++)  mol[i].Ni += bp[k][i].nbox;
	}
	
  /* ----------------------------------------------- */
  /* Read in the number of residues in one moluecule */
  /* for each type in of molecule the box.           */
  /* ----------------------------------------------- */    
	for (int i=0; i<sim.NC; i++) {
      fscanf(molinput_ptr,"%d",&mol[i].Nres);  fgets(tt,80,molinput_ptr);
    }
    fgets(tt,80,molinput_ptr);

  /* ----------------------------------------------- */
  /* Read in the number of molecules to perform      */
  /* fancy moves on.                                 */
  /* ----------------------------------------------- */    
	fscanf(molinput_ptr, "%d",&molec.Nsolute); fgets(tt,80,molinput_ptr);
	fscanf(molinput_ptr, "%d",&SITE1); fgets(tt,80,molinput_ptr);
	fscanf(molinput_ptr, "%d",&SITE2); fgets(tt,80,molinput_ptr);
	fscanf(molinput_ptr, "%d",&SITE3); fgets(tt,80,molinput_ptr);
	SITE1 = SITE1 - 1;
    	SITE2 = SITE2 - 1;
	SITE3 = SITE3 - 1;

	fgets(tt,80,molinput_ptr);

	fclose(molinput_ptr);    
  /* ================================================================== */
  /*                                                                    */
  /* Calculate the total number of molecules, in each box by summing up */
  /* the number of each type of molecule in each box.  Also, calculate  */
  /* the total number of sites and residues in each box by multiplying  */
  /* the number of each type of molecule in each box by the number of   */
  /* sites in each molecule and summing them up.                        */
  /*                                                                    */
  /* ================================================================== */
    for (int k=0; k<sim.NB; k++){
		box[k].boxn  = 0;
		box[k].boxns = 0;
		for(int i=0; i<sim.NC; i++) {
		box[k].boxn  += bp[k][i].nbox;// number of molecules in each box
		box[k].boxns += bp[k][i].nbox * mol[i].Nsite;// number of sites in each box
		box[k].boxnres+= bp[k][i].nbox * mol[i].Nres;// number of residues in each box
		}
	}
    
  /* ----------------------------------------------- */
  /* Check to make sure the total number of          */
  /* molecules is correct. (This check compares      */
  /* numbers from the mol.input file with the        */
  /* number given in simul.input.                    */
  /* ----------------------------------------------- */    
  int qt = 0;
  for (int i=0; i<sim.NC; i++)  
	  qt += mol[i].Ni;
  if (qt != Nmolall) {
    fprintf(stdout,"total number of molecules do not match in mol.input!!! (Total number of molecules = %i, sum of individual molecules = %i) \n", Nmolall, qt); exit(1);}

  
  /* ================================================================== */
  /* Allocate the molec variable and calculate the beginning site of    */
  /* each molecule.                                                     */
  /* ================================================================== */
  molec.fsite = (int*) calloc(box[0].boxn,sizeof(int));
	if(molec.fsite == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for molec.fsite\n"); exit(11);}
  molec.lsite = (int*) calloc(box[0].boxn,sizeof(int));
	if(molec.lsite == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for molec.lsite\n"); exit(11);}
  molec.fres  = (int*) calloc(box[0].boxn,sizeof(int));
    if(molec.fres == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for molec.fres\n"); exit(11);}
  molec.lres  = (int*) calloc(box[0].boxn,sizeof(int));
    if(molec.lres == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for molec.lres\n"); exit(11);}
/* //This is not right
  int nmolec = 0;
  int prev_moliNsite = mol[0].Nsite;
  int prev_moliNres  = mol[0].Nres;
  for(int i=0; i<sim.NC; i++){
	  for(int j=0; j<bp[0][i].nbox; j++){
		  if(i==0 && j==0){
			  molec.fsite[0] = 0;
			  molec.lsite[0] = mol[0].Nsite;
			  molec.fres[0]  = 0;
			  molec.lres[0]  = mol[0].Nres;
		  }
		  else{
			  molec.fsite[nmolec] = molec.fsite[nmolec-1]+prev_moliNsite;
			  molec.lsite[nmolec] = molec.lsite[nmolec-1]+prev_moliNsite;
			  molec.fres[nmolec] = molec.fres[nmolec-1]+prev_moliNres;
			  molec.lres[nmolec] = molec.lres[nmolec-1]+prev_moliNres;
		  }
		  nmolec++;
		  prev_moliNsite = mol[i].Nsite;
		  prev_moliNres  = mol[i].Nres;
	  }

  }
*/

  int nmolec = 0;
  for(int i=0; i<sim.NC; i++){
	  for(int j=0; j<bp[0][i].nbox; j++){
		  if(i==0 && j==0){
			  molec.fsite[0] = 0;
			  molec.lsite[0] = mol[0].Nsite;
			  molec.fres[0]  = 0;
			  molec.lres[0]  = mol[0].Nres;
		  }
		  else{
			  molec.fsite[nmolec] = molec.lsite[nmolec-1];
			  molec.lsite[nmolec] = molec.fsite[nmolec] +mol[i].Nsite;
			  molec.fres[nmolec]  = molec.lres[nmolec-1];
			  molec.lres[nmolec]  = molec.lres[nmolec-1]+mol[i].Nres;
		  }
		  nmolec++;
	  }
  }

  /* ================================================================== */
  /* Now that the box numbers are known, allocate the needed memory.    */
  /* ================================================================== */
    allocm();

  /* ================================================================== */
  /* Read in the connectivity.                                          */
  /* ================================================================== */
  for(int k=0; k<sim.NB; k++) read_topology(k);
	

  /* ================================================================== */
  /* Now, read in the x,y,z coordinates from the simul#.crd file.       */
  /* ================================================================== */
  for(int k=0; k<sim.NB; k++) {
  /* ----------------------------------------------- */
  /* Set the file pointers to the correct input      */
  /* file.                                           */
  /* ----------------------------------------------- */    
#ifdef MPI
	sprintf(namec,"./INPUT/simul%d.crd",mpi.my_rank);
#endif
#ifndef MPI
	sprintf(namec,"./INPUT/simul%d.crd",k);
#endif
    das = fopen(namec,"r");

  /* ----------------------------------------------- */
  /* Variables needed to read the files.             */
  /* ----------------------------------------------- */    
	int rescount =0; char temp[5]; char tempc[5];

  /* ----------------------------------------------- */
  /* Check to make sure the number of sites and      */
  /* number of coordinate sets match.                */
  /* ----------------------------------------------- */    
		fscanf(das, "%d",   &qw); 
		if (qw != box[k].boxns) {
        fprintf(stdout,"ERROR: number of sites in box (%i) do not match with number of coordinate sets (%i) in %s!!!\n", box[k].boxns, qw, namec);
        exit(1);
		}
  /* ----------------------------------------------- */
  /* Read in the site name and residue number of     */
  /* each site.  Compare them to those read in       */
  /* read_topology().  If they are different, exit.  */
  /* ----------------------------------------------- */    
	for (int i=0; i<box[k].boxns; i++){
		  fscanf(das, "%d",   &qw);	
		  fscanf(das, "%d",   &rescount);
		  fscanf(das, "%s",  temp);
		  fscanf(das, "%4s",  tempc);
		  if ( (strcmp(atnopbc[k][i].name,tempc)) != 0) {

		 	fprintf(stdout,"ERROR: Name of site %d (%s) in crd does not match name in psf %d (%s).\n",qw,temp,i,atnopbc[k][i].name);
			exit(1);
		  }
			  
  /* ----------------------------------------------- */
  /* Read in the site x,y,z coordinates.             */
  /* ----------------------------------------------- */    
	  fscanf(das,"%lf",  &atnopbc[k][i].x);
	  fscanf(das,"%lf",  &atnopbc[k][i].y);
	  fscanf(das,"%lf",  &atnopbc[k][i].z);fgets(tt,80,das);
	}// end i
    fclose(das);    
  }// end k


  /* ================================================================== */
  /*                                                                    */
  /* With the information above, initialize some variables.             */
  /*                                                                    */
  /* ================================================================== */
  /* ----------------------------------------------- */
  /* Assign the x,y,z c coordinates for use with     */
  /* periodic boundary conditions.                   */
  /* ----------------------------------------------- */    
  for(int k=0; k<sim.NB; k++) 
    for(int i=0; i<box[k].boxns; i++) atom[k][i] = atnopbc[k][i];


  /* ----------------------------------------------- */
  /* Apply periodic boundary conditions.             */
  /* ----------------------------------------------- */    
  #ifndef WALL
  for(int k=0; k<sim.NB; k++) pbc_init(k);
  #endif
  
  /* ----------------------------------------------- */
  /* Calculate the box height and volume.            */
  /* ----------------------------------------------- */    
  for(int k=0; k<sim.NB; k++){
  	box[k].hx = box[k].lx/2.0;
	box[k].hy = box[k].ly/2.0;
	box[k].hz = box[k].lz/2.0;
	box[k].vol  = box[k].lx*box[k].ly*box[k].lz;

	double least_boxl = box[k].lx;
	if(least_boxl > box[k].ly) least_boxl = box[k].ly;
	if(least_boxl > box[k].lz) least_boxl = box[k].lz;
	double least_boxh = least_boxl/2.0;
	
	box[k].least_boxl = least_boxl;
	box[k].least_boxh = least_boxh;
    if(sim.rc > 0.0) box[k].rc2 = sim.rc * sim.rc;
    else box[k].rc2 = least_boxl * least_boxl;
	box[k].rc4 = box[k].rc2 * box[k].rc2;	
  }


  /* ----------------------------------------------- */
  /* Calculate the number of degrees of freedom in   */
  /* each box.                                       */
  /* ----------------------------------------------- */    
  for(int k=0; k<sim.NB; k++) {
	  box[k].nfree = 3.0 * box[k].boxns;
#ifdef SMD
	  box[k].nfree = box[k].nfree -3.0;
#endif
#ifdef NSMD
	  box[k].nfree = box[k].nfree -6.0;
#endif
#ifdef KONS
	  box[k].nfree = box[k].nfree -3.0*(double)cons[k].n;
#endif
  
  }

  /* ----------------------------------------------- */
  /* Calculate the variables needed for higher order */
  /* integration of Nose-Hoover thermostats.         */
  /* ----------------------------------------------- */    
  if(sim.ID2 == 4 || sim.ID2 == 5 || sim.ID==6 || sim.ID2==8){
	for(int k=0; k<sim.NB;k++)
	{

			double kT = KB*sim.T[k]*1E-4; //kg*ang^2/ps^2
			double Q = nhc[k].Q;

			for (int i=0; i<nhc[k].M; i++)
			{
				nhc[k].zeta[i].r = 0.0;
				nhc[k].zeta[i].v = 0.0;
			}
			for (int i=0; i<nhc[k].M-1; i++)		
				nhc[k].zeta[i+1].G = (Q*nhc[k].zeta[i].v*nhc[k].zeta[i].v-kT)/Q;

			if(nhc[k].Nys%2 == 0)
			{
				fprintf(stdout, "The number of higher order terms for Nose-Hoover (%i) must be odd!\n", nhc[k].Nys);
				exit(1);
			}
			else
			{
				double num = (double)nhc[k].Nys-1.0;
				int cent_num = nhc[k].Nys/2;
		
				if(nhc[k].Nys==1)
					nhc[k].w[0]=1;
				else
				{	
					for(int i=0; i<nhc[k].Nys; i++)
					{
						if(i!=cent_num)
						nhc[k].w[i] = 1.0/(num-pow(num,1.0/3.0));
						else
							nhc[k].w[i] = 1.0-num*1.0/(num-pow(num,1.0/3.0));
					}
				}
			}

	}
  }
  /* ================================================================== */
  /*                                                                    */
  /* If SASA is defined, read in the information from sasa.inp.         */
  /*                                                                    */
  /* ================================================================== */
#ifdef SASA
#ifdef BGO
  /* ----------------------------------------------- */
  /* Read in the info on solvent for BGO solvent     */
  /* energy.                                         */
  /* ----------------------------------------------- */
  read_solvent();
  // If BGO is defined, read in SASA parameters from code not by sasa.inp
  vdw_radii_bgo();
#endif 
#ifndef BGO
  /* ----------------------------------------------- */
  /* Variables and pointers needed to read the       */
  /* files.                                          */
  /* ----------------------------------------------- */    
  int sa; double sa1, sa2, sa3, sa4;
  FILE *das1;

  /* ----------------------------------------------- */
  /* Open the file.                                  */
  /* ----------------------------------------------- */    
  if( NULL == (das1=fopen("./INPUT/sasa.inp","r")) ) {
    fprintf(stdout,"input file sasa.inp does not exist!!!\n");
    exit(1);
  }  
 
  /* ----------------------------------------------- */
  /* Zero out the parameters.                        */
  /* ----------------------------------------------- */    
  for(int k=0; k<sim.NB; k++){
		  for(int i =0; i<box[k].boxns; i++){//move this?
			  sasa[k][i].A		= 0.0;
			  sasa[k][i].p		= 0.0;
			  sasa[k][i].r		= 0.0;
			  sasa[k][i].S		= 0.0;
			  sasa[k][i].sigma	= 0.0;
		  }
	  }

  /* ----------------------------------------------- */
  /* Read in the parameters.                         */
  /* ----------------------------------------------- */    
  while (!feof(das1)){


	  fscanf(das1, "%d",       &sa);
	  fscanf(das1, "%lf",     &sa1);
	  fscanf(das1, "%lf",     &sa2);
	  fscanf(das1, "%lf",     &sa3);
	  fscanf(das1, "%lf\n",   &sa4);
    
		for(int k=0; k<sim.NB; k++){
		  for(int i =0; i<box[k].boxns; i++){
			  if (atom[k][i].atomid == sa){
				  sasa[k][i].A		= 0.0;
				  sasa[k][i].p		= sa3;
				  sasa[k][i].r		= sa2;
				  sasa[k][i].S		= 4.0*PI*(RPROBE+sasa[k][i].r)*(RPROBE+sasa[k][i].r);
				  sasa[k][i].sigma	= sa4*4.184;
			  }
		  }
	  }
  }
  fclose(das1);
  /* ----------------------------------------------- */
  /* Verify that all atom types have valid sasa      */
  /* parameters.                                     */
  /* ----------------------------------------------- */    
	for(int k=0; k<sim.NB; k++){	
		  for(int i =0; i<box[k].boxns; i++){
			  if (sasa[k][i].r == 0.0 || sasa[k][i].S ==0.0){
				 fprintf(stdout,"no sasa parameter for atom number %i (atom type = %i) !!!\n", i+1, atom[k][i].atomid);
				 exit(1);
			  }
		  }
	  }
#endif
#endif

  /* ----------------------------------------------- */
  /* If running PR NPT, set the initial box matrix.  */
  /* ----------------------------------------------- */    
#ifdef PR_NPT
	if(sim.ID==6){
		if(sim.ID2 == 6){
			for(int k=0; k<sim.NB; k++){
				axes[k][0] = box[k].lx; 
				axes[k][1] = 0.0;
				axes[k][2] = 0.0;
				axes[k][3] = 0.0;
				axes[k][4] = box[k].ly;
				axes[k][5] = 0.0;
				axes[k][6] = 0.0;
				axes[k][7] = 0.0;
				axes[k][8] = box[k].lz;
			}
		}
		//put other NPT routine here
	}
  
#endif//PR_NPT

  
}
