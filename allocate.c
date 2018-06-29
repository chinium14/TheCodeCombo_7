/* ======================================================================== */
/* allocate.cpp                                                             */
/*                                                                          */
/*       This subroutine allocates memory by using the calloc command on    */
/* structures found in defines.h.                                           */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"


/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

void allocm (void) 
{

  // NOTE: extra memory is allocated in case the number of molecules
  //       in each box changes

  /* ================================================================== */
  /* Allocate memory for all sites.                                     */
  /*                                                                    */
  /*	Here, the sturcture called "atom" has two indices.  The first   */
  /* is the box number and the second is the site number.               */
  /*                                                                    */
  /* ================================================================== */
  atom    = (struct atoms**) calloc(sim.NB,sizeof(struct atoms*));
		if(atom == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for atom\n"); exit(11);}
  atnopbc = (struct atoms**) calloc(sim.NB,sizeof(struct atoms*));
		if(atnopbc == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for atnopbc\n"); exit(11);}
  atom_temp = (struct atoms**) calloc(sim.NB,sizeof(struct atoms*));
		if(atom_temp == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for atom_temp\n"); exit(11);}
  atnopbc_temp = (struct atoms**) calloc(sim.NB,sizeof(struct atoms*));
		if(atnopbc_temp == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for atnopbc_temp\n"); exit(11);}
  nlist = (struct nlist_def*) calloc(sim.NB,sizeof(struct nlist_def));
		if(nlist == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist\n"); exit(11);}
#ifdef SLIST
  slist = (struct slist_def*) calloc(sim.NB,sizeof(struct slist_def));
		if(slist == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for slist\n"); exit(11);}
#endif

#ifdef EWALD
  nlist_ew = (struct nlist_def*) calloc(sim.NB,sizeof(struct nlist_def));
		if(nlist_ew == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist_ew\n"); exit(11);}
#endif
		
#ifdef DLIST
  rx0	=	(double**) calloc(sim.NB,sizeof(double*));
		if(rx0 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rx0\n"); exit(11);}
  ry0	=	(double**) calloc(sim.NB,sizeof(double*));
		if(ry0 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ry0\n"); exit(11);}
  rz0	=	(double**) calloc(sim.NB,sizeof(double*));
		if(rz0 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rz0\n"); exit(11);}
  nl_flag=	(int*)calloc(sim.NB,sizeof(int));
		if(nl_flag == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nl_flag\n"); exit(11);}
#endif
  for(int i=0; i<sim.NB; i++) {
    int na = box[i].boxns + MEXTRA;
    atom[i]    = (struct atoms*) calloc(na,sizeof(struct atoms));
		if(atom[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for atom[i]\n"); exit(11);}
	atom_temp[i]    = (struct atoms*) calloc(na,sizeof(struct atoms));
		if(atom_temp[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for atom[i]_temp\n"); exit(11);}
    atnopbc[i] = (struct atoms*) calloc(na,sizeof(struct atoms));
		if(atnopbc[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for atnopbc[i]\n"); exit(11);}
	atnopbc_temp[i] = (struct atoms*) calloc(na,sizeof(struct atoms));
		if(atnopbc_temp[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for atnopbc_temp[i]\n"); exit(11);}

#ifdef DLIST
    rx0[i]	=	(double*) calloc(na,sizeof(double));
		if(rx0[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rx0[i]\n"); exit(11);}
    ry0[i]	=	(double*) calloc(na,sizeof(double));
		if(ry0[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ry0[i]\n"); exit(11);}
    rz0[i]	=	(double*) calloc(na,sizeof(double));
		if(rz0[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for rz0[i]\n"); exit(11);}
#endif
  }

  /* ================================================================== */
  /* Allocate memory for the SASA implicit solvent model.               */
  /*                                                                    */
  /*	Here, the sturcture called "sasa" has two indices.  The first   */
  /* is the box number and the second is the site number.               */
  /*                                                                    */
  /* ================================================================== */  
#ifdef SASA
#ifndef BGO
  sasa    = (struct imp_sasa**) calloc(sim.NB,sizeof(struct imp_sasa*));
		if(sasa == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for sasa\n"); exit(11);}
  for(int i=0; i<sim.NB; i++) {
    int na = box[i].boxns + MEXTRA;
    sasa[i]    = (struct imp_sasa*) calloc(na,sizeof(struct imp_sasa));
		if(sasa[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for sasa[i]\n"); exit(11);}
  }
#endif
#ifdef BGO
  sasa    = (struct bgo_sasa**) calloc(sim.NB,sizeof(struct bgo_sasa*));
		if(sasa == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for sasa\n"); exit(11);}
  for(int i=0; i<sim.NB; i++) {
    int na = box[i].boxns + MEXTRA;
    sasa[i]    = (struct bgo_sasa*) calloc(na,sizeof(struct bgo_sasa));
		if(sasa[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for sasa[i]\n"); exit(11);}
  }
  tfe    = (struct tfe_sasa**) calloc(sim.NB,sizeof(struct tfe_sasa*));
		if(tfe == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for sasa\n"); exit(11);}
  for(int i=0; i<sim.NB; i++) {
    int na = box[i].boxns + MEXTRA;
    tfe[i]    = (struct tfe_sasa*) calloc(na,sizeof(struct tfe_sasa));
		if(tfe[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for sasa[i]\n"); exit(11);}
  }
  beta    = (double **) calloc(sim.NB,sizeof(double*));
		if(beta == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for sasa\n"); exit(11);}
  for(int i=0; i<sim.NB; i++) {
    int na = box[i].boxns + MEXTRA;
    beta[i]    = (double*) calloc(na,sizeof(double));
		if(beta[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for sasa[i]\n"); exit(11);}
  }
#ifdef MEMBRANE
  #ifndef MPI
    concent_local = (double *) calloc(box[0].boxns,sizeof(double));
  #else
    concent_local = (double *) calloc(box[0].boxns,sizeof(double));
  #endif
#endif

#ifndef MPI
//  solvent_type    = (int *) calloc(sim.NB,sizeof(int));
  concentration    = (double *) calloc(sim.NB,sizeof(double));
#else
//  solvent_type    = (int *) calloc(mpi.p,sizeof(int));
  concentration    = (double *) calloc(mpi.p,sizeof(double));
#endif

#endif
#endif

  /* ================================================================== */
  /* Allocate memory for the residue structure.                         */
  /*                                                                    */
  /*	Here, the sturcture called "residue" has two indices.  The      */
  /* first is the box number and the second is the residue number.      */
  /*                                                                    */
  /* ================================================================== */
  residue = (struct res**) calloc(sim.NB,sizeof(struct res*));
		if(residue == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for residue"); exit(11);}
  for(int i=0; i<sim.NB; i++) {
    int na = box[i].boxnres + MEXTRA;						// boxnres = total number of amino acid residue in box i
    residue[i] = (struct res*) calloc(na,sizeof(struct res));
		if(residue[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for residue[i]\n"); exit(11);}
  }

  /* ================================================================== */
  /* Allocate memory for the velocity structures.                       */
  /*                                                                    */
  /*	Here, the sturctures called "uu" & "vv" have two indices.  The  */
  /* first is the box number and the second is the site/atom number.    */
  /*                                                                    */
  /* ================================================================== */
  vv  = (struct veloc**) calloc(sim.NB,sizeof(struct veloc*));
		if(vv == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for vv\n"); exit(11);}
  uu  = (struct veloc**) calloc(sim.NB,sizeof(struct veloc*));
		if(uu == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for uu\n"); exit(11);}
  vv_temp=(struct veloc**) calloc(sim.NB,sizeof(struct veloc*));
		if(vv_temp == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for vv_temp\n"); exit(11);}
  for(int i=0; i<sim.NB; i++) {
    int vn = box[i].boxns + MEXTRA;
    vv[i]  = (struct veloc*) calloc(vn,sizeof(struct veloc));
		if(vv[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for vv[i]\n"); exit(11);}
    uu[i]  = (struct veloc*) calloc(vn,sizeof(struct veloc));
		if(uu[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for uu[i]\n"); exit(11);}
	vv_temp[i]  = (struct veloc*) calloc(vn,sizeof(struct veloc));
		if(vv_temp[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for vv_temp[i]\n"); exit(11);}
  }

  /* ================================================================== */
  /* Allocate memory for the force structures.                          */
  /*                                                                    */
  /*	Here, the sturctures called "ff", "ff_short", & ff_long" have   */
  /* two indices.  The first is the box number and the second is the    */
  /* site/atom number.                                                  */
  /* ================================================================== */
  ff = (struct veloc**) calloc(sim.NB,sizeof(struct veloc*));
		if(ff == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ff\n"); exit(11);}
  for(int i=0; i<sim.NB; i++) {
    int fn = box[i].boxns + MEXTRA;
    ff[i] = (struct veloc*) calloc(fn,sizeof(struct veloc));
		if(ff[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ff[i]\n"); exit(11);}
  }
  ff_short = (struct veloc**) calloc(sim.NB,sizeof(struct veloc*));
		if(ff_short == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ff_short\n"); exit(11);}
  for(int i=0; i<sim.NB; i++) {
    int fn = box[i].boxns + MEXTRA;
    ff_short[i] = (struct veloc*) calloc(fn,sizeof(struct veloc));
		if(ff_short[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory ff_short[i]\n"); exit(11);}
  }
  ff_long = (struct veloc**) calloc(sim.NB,sizeof(struct veloc*));
		if(ff_long == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ff_long\n"); exit(11);}
  for(int i=0; i<sim.NB; i++) {
    int fn = box[i].boxns + MEXTRA;
    ff_long[i] = (struct veloc*) calloc(fn,sizeof(struct veloc));
		if(ff_long[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ff_long[i]\n"); exit(11);}
  }
  ff_temp=(struct veloc**) calloc(sim.NB,sizeof(struct veloc*));
		if(ff_temp == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ff_temp\n"); exit(11);}
  for(int i=0; i<sim.NB; i++) {
    int fn = box[i].boxns + MEXTRA;
    ff_temp[i] = (struct veloc*) calloc(fn,sizeof(struct veloc));
		if(ff_temp[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ff_temp[i]\n"); exit(11);}
  }
#ifdef MC
  ffox = (double**) calloc(sim.NB,sizeof(double*));
		if(ffox == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ffox\n"); exit(11);}
  for(int i=0; i<sim.NB; i++) {
    int fn = box[i].boxns + MEXTRA;
    ffox[i] = (double*) calloc(fn,sizeof(double));
		if(ffox[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ffox[i]\n"); exit(11);}
  }
    ffoy = (double**) calloc(sim.NB,sizeof(double*));
		if(ffoy == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ffoy\n"); exit(11);}
  for(int i=0; i<sim.NB; i++) {
    int fn = box[i].boxns + MEXTRA;
    ffoy[i] = (double*) calloc(fn,sizeof(double));
		if(ffoy[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ffoy[i]\n"); exit(11);}
  }
    ffoz = (double**) calloc(sim.NB,sizeof(double*));
		if(ffoz == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ffoz\n"); exit(11);}
  for(int i=0; i<sim.NB; i++) {
    int fn = box[i].boxns + MEXTRA;
    ffoz[i] = (double*) calloc(fn,sizeof(double));
		if(ffoz[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ffoz[i]\n"); exit(11);}
  }

/*#ifdef CONFIGT
  hesox = (double*) calloc(sim.NB,sizeof(double));
		if(hesox == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for hesox\n"); exit(11);}
  hesoy = (double*) calloc(sim.NB,sizeof(double));
		if(hesoy == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for hesoy\n"); exit(11);}
  hesoz = (double*) calloc(sim.NB,sizeof(double));
		if(hesoz == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for hesoz\n"); exit(11);}
  hesor = (double*) calloc(sim.NB,sizeof(double));
		if(hesor == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for hesor\n"); exit(11);}
#endif
*/
#ifdef PRESSURE

#endif

#endif//MC
  /* ================================================================== */
  /*                                                                    */
  /* Allocate memory for the force field parameter.                     */
  /*                                                                    */
  /* ================================================================== */
  pott  = (struct inter**) calloc(sim.NB,sizeof(struct inter*));
		if(pott == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for pott\n"); exit(11);}
#ifdef S14	
  pott14  = (struct inter**) calloc(sim.NB,sizeof(struct inter*));
		if(pott14 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for pott14\n"); exit(11);}
#endif
#ifndef STYPE
  ljset = (struct interb*) calloc(sim.NB,sizeof(struct interb));
		if(ljset == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ljset\n"); exit(11);}
#ifdef S14
  ljset14 = (struct interb*) calloc(sim.NB,sizeof(struct interb));
		if(ljset14 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ljset14\n"); exit(11);}
#endif
#endif
  bond  = (struct bonds**) calloc(sim.NB,sizeof(struct bonds*));
		if(bond == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for bond\n"); exit(11);}
  bend  = (struct bends**) calloc(sim.NB,sizeof(struct bends*));
		if(bend == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for bend\n"); exit(11);}
#ifndef NEUTRAL
  urey  = (struct ureys**) calloc(sim.NB,sizeof(struct ureys*));
		if(urey == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for urey\n"); exit(11);}
#endif

  tors  = (struct torsions**) calloc(sim.NB,sizeof(struct torsions*));
		if(tors == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for tors\n"); exit(11);}
  impr  = (struct impropers**) calloc(sim.NB, sizeof(struct impropers*));
		if(impr == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for impr\n"); exit(11);}
#ifdef STYPE
  xint  = (struct excinter**) calloc(sim.NB,sizeof(struct excinter*));
 		if(xint == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for xint\n"); exit(11);} 
#endif
  in14  = (struct interactions1_4**) calloc(sim.NB, sizeof(struct interactions1_4*));
		if(in14 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for in14\n"); exit(11);}
  hdon  = (struct donors**) calloc(sim.NB, sizeof(struct donors*));
		if(hdon == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for hdon\n"); exit(11);}
  hacc  = (struct acceptors**) calloc(sim.NB, sizeof(struct acceptors*));
		if(hacc == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for hacc\n"); exit(11);}
  exNB  = (struct NBexclusions**) calloc(sim.NB, sizeof(struct NBexclusions*));
		if(exNB == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for exNB\n"); exit(11);}
#ifdef STYPE
  table = (struct hashlist**) calloc(sim.NB,sizeof(struct hashlist*));
		if(table == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for table\n"); exit(11);}
#endif

		
		
		//  table = (struct hashlist**) calloc(sim.NB,sizeof(struct hashlist*));

  /* ================================================================== */
  /* Allocate memory for the neighbor list structure.                   */
  /*                                                                    */
  /*	Here, the sturcture called "nlist" has one indice for the       */
  /* box number.  The "count" member of the structure has one indice    */
  /* for the site/atom number and the member "list" has two indices.    */
  /* The first is the site/atom number and the second is the site       */
  /* number of its neighbor.                                            */
  /*                                                                    */
  /* ================================================================== */
  for(int k=0; k<sim.NB; k++) {

    int size = box[k].boxns;
#ifdef WALL
    size += 2;
#endif   
    if(size < 2000) max_num_neighbors = size;
    else max_num_neighbors = 1500;

#ifdef DNA_GOLIK
    max_num_neighbors = 1000;
#endif

    nlist[k].count = (int*) calloc(size,sizeof(int));
	if(nlist[k].count == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist[k].count\n"); exit(11);}
    nlist[k].list  = (int**) calloc(size,sizeof(int*));
	if(nlist[k].list == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist[k].list\n"); exit(11);}
#if defined(STYPE) || defined(DNA_GOLIK)
    nlist[k].xflag = (int**) calloc(size,sizeof(int*));
	if(nlist[k].xflag == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist[k].xflag\n"); exit(11);}
#endif
#ifdef DNA_GOLIK
    nlist[k].sig = (double**) calloc(size,sizeof(double*));
	if(nlist[k].sig == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist[k].sig\n"); exit(11);}
#endif
    for(int i=0; i<size; i++){ 
      	nlist[k].list[i]   = (int*) calloc(max_num_neighbors,sizeof(int));
		if(nlist[k].list[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist[k].list[i]\n"); exit(11);}
#if defined(STYPE) || defined(DNA_GOLIK)
	nlist[k].xflag[i]  = (int*) calloc(max_num_neighbors,sizeof(int));
		if(nlist[k].xflag[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist[k].xflag[i]\n"); exit(11);}
#endif
#ifdef DNA_GOLIK
	nlist[k].sig[i]  = (double*) calloc(max_num_neighbors,sizeof(double));
		if(nlist[k].sig[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist[k].sig[i]\n"); exit(11);}
#endif
    }
#ifdef EWALD    
#ifndef SPME
	if(size > 2000){
		max_num_neighbors = box[k].boxns/2;
	}
#endif
    nlist_ew[k].count = (int*) calloc(size,sizeof(int));
	if(nlist_ew[k].count == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist_ew[k].count\n"); exit(11);}
    nlist_ew[k].list  = (int**) calloc(size,sizeof(int*));
	if(nlist_ew[k].list == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist_ew[k].list\n"); exit(11);}
#ifdef STYPE
    nlist_ew[k].xflag = (int**) calloc(size,sizeof(int*));
	if(nlist_ew[k].xflag == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist_ew[k].xflag\n"); exit(11);}
#endif
    for(int i=0; i<size; i++){ 
      	nlist_ew[k].list[i]   = (int*) calloc(max_num_neighbors,sizeof(int));
		if(nlist_ew[k].list[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist_ew[k].list[i]\n"); exit(11);}
#ifdef STYPE
	nlist_ew[k].xflag[i]  = (int*) calloc(max_num_neighbors,sizeof(int));
		if(nlist_ew[k].xflag[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist_ew[k].xflag[i]\n"); exit(11);}
#endif
    }  
#endif //EWALD
    
#ifdef SLIST
	int size_slist	= molec.lsite[molec.Nsolute-1];
	int max_neibor	= max_num_neighbors/3;	
	slist[k].count = (int*) calloc(size_slist,sizeof(int));
		if(slist[k].count == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for slist[k].count\n"); exit(11);}
    slist[k].list  = (int**) calloc(size_slist,sizeof(int*));
		if(slist[k].list == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for slist[k].list\n"); exit(11);}
	for(int i=0; i<size_slist; i++){ 
		slist[k].list[i]   = (int*) calloc(max_neibor,sizeof(int));
		  if(slist[k].list[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for slist[k].list[i]\n"); exit(11);}
	}
#endif//SLIST



  /* ================================================================== */
  /* Allocate memory for the neighbor list of the implicit solvent      */
  /* SASA model.                                                        */
  /*                                                                    */
  /*	This is the same as the structure just above except it is used  */
  /* when employing the implicit solvent model.                         */
  /*                                                                    */
  /* ================================================================== */
#ifdef STYPE
#ifdef SASA
	nlist[k].count_sasa = (int*) calloc(size,sizeof(int));
    nlist[k].list_sasa  = (int**) calloc(size,sizeof(int*));
    nlist[k].xflag_sasa      = (int**) calloc(size,sizeof(int*));
	for(int i=0; i<size; i++){ 
      nlist[k].list_sasa[i] = (int*) calloc(max_num_neighbors,sizeof(int));
	  nlist[k].xflag_sasa[i]     = (int*) calloc(max_num_neighbors,sizeof(int));
	}
#endif
#endif //STYPE

#ifndef STYPE
#ifdef SASA
	nlist[k].count_sasa = (int*) calloc(size,sizeof(int));
		if(nlist[k].count_sasa == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist[k].count_sasa\n"); exit(11);}
    nlist[k].list_sasa  = (int**) calloc(size,sizeof(int*));
		if(nlist[k].list_sasa == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist[k].list_sasa\n"); exit(11);}
    for(int i=0; i<size; i++){ 
      nlist[k].list_sasa[i] = (int*) calloc(size,sizeof(int));
		if(nlist[k].list_sasa[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nlist[k].list_sasa\n"); exit(11);}}
#endif
#endif// STYPE
  }

  /* ================================================================== */
  /* Allocate memory for the Nose Hoover Chain thermostat method.       */
  /*                                                                    */
  /*	The nhc sturcture has one indice for the box number.  The       */
  /* "zeta" member of this structure has one indice for the thermostat  */
  /* number, and the "w" member has one indice for the higher order     */
  /* integration terms.                                                 */
  /*                                                                    */
  /* ================================================================== */
  for(int k=0; k<sim.NB; k++){
	int size_1 = nhc[k].M;
	int size_2 = nhc[k].Nys;
	nhc[k].zeta = (struct NH_therm*) calloc(size_1, sizeof(struct NH_therm));
		if(nhc[k].zeta == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nhc[k].zeta\n"); exit(11);}
	nhc[k].w = (double*) calloc(size_2, sizeof(double));
		if(nhc[k].w == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nhc[k].w\n"); exit(11);}
  }
  /* ================================================================== */
  /*                                                                    */
  /* Allocate memory for Ewald Sums k vectors.                          */
  /*                                                                    */
  /* ================================================================== */
#ifdef EWALD
#ifdef SPME
  int ORDER = sim.order;

  int bx = (int)(box[0].lx);
  int by = (int)(box[0].ly);
  int bz = (int)(box[0].lz);

  if(oddeven(bx)) bx += 1;
  if(oddeven(by)) by += 1;
  if(oddeven(bz)) bz += 1;

  //sim.grid[0] = bx;
  //sim.grid[1] = by;
  //sim.grid[2] = bz;

  int GRIDX = sim.grid[0];
  int GRIDY = sim.grid[1];
  int GRIDZ = sim.grid[2];

  int k = 0;
  int ns = box[k].boxns+1;
  ss = (struct pmedat*) calloc(ns,sizeof(struct pmedat));
  if(ss == NULL) {fprintf(stderr,"ERROR: cannot allocate memory for SPME::ss\n"); exit(1);}

  mm = (struct pmedat**) calloc(ns,sizeof(struct pmedat*));
  if(mm == NULL) {fprintf(stderr,"ERROR: cannot allocate memory for SPME::mm\n"); exit(1);}
  dm = (struct pmedat**) calloc(ns,sizeof(struct pmedat*));
  if(dm == NULL) {fprintf(stderr,"ERROR: cannot allocate memory for SPME::dm\n"); exit(1);}
  for(int i=0; i<ns; i++) {
    mm[i] = (struct pmedat*) calloc(ORDER+1,sizeof(struct pmedat));
    if(mm[i] == NULL) {fprintf(stderr,"ERROR: cannot allocate memory for SPME::mm*\n"); exit(1);}
    dm[i] = (struct pmedat*) calloc(ORDER+1,sizeof(struct pmedat));
    if(dm[i] == NULL) {fprintf(stderr,"ERROR: cannot allocate memory for SPME::dm*\n"); exit(1);}
  }

  mn = (double*) calloc(ORDER+2,sizeof(double));
  if(mn == NULL) {fprintf(stderr,"ERROR: cannot allocate memory for SPME::mn\n"); exit(1);}
  dmn = (double*) calloc(ORDER+2,sizeof(double));
  if(dmn == NULL) {fprintf(stderr,"ERROR: cannot allocate memory for SPME::dmn\n"); exit(1);}

  bsp_modx = (double*) calloc(GRIDX+5,sizeof(double));
  if(bsp_modx == NULL) {fprintf(stderr,"ERROR: cannot allocate memory for SPME::bsp_modx\n"); exit(1);}
  bsp_mody = (double*) calloc(GRIDY+5,sizeof(double));
  if(bsp_mody == NULL) {fprintf(stderr,"ERROR: cannot allocate memory for SPME::bsp_mody\n"); exit(1);}
  bsp_modz = (double*) calloc(GRIDZ+5,sizeof(double));
  if(bsp_modz == NULL) {fprintf(stderr,"ERROR: cannot allocate memory for SPME::bsp_modz\n"); exit(1);}

  gridmx = greater(GRIDX,greater(GRIDY,GRIDZ))+5;

  bsp_arr = (double*) calloc(gridmx,sizeof(double));
  if(bsp_arr == NULL) {fprintf(stderr,"ERROR: cannot allocate memory for SPME::bsp_arr\n"); exit(1);}

  int ngrid = (GRIDX+5)*(GRIDY+5)*(GRIDZ+5);
  qgrida = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ngrid);
  if(qgrida == NULL) {fprintf(stderr,"ERROR: cannot allocate memory for SPME::qgrida\n"); exit(1);}

  for(int i=0; i<ngrid; i++) {
    qgrida[i][0] = 0.0;   qgrida[i][1] = 0.0;
  }

#else
  kvec    = (double**) calloc(sim.NB,sizeof(double*));
		if(kvec == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for kvec\n"); exit(11);}
  for(int k=0; k<sim.NB; k++){
	  kvec[k] = (double*) calloc((sim.kmax[k][0]+1)*(2*sim.kmax[k][1]+1)*(2*sim.kmax[k][2]+1),sizeof(double));
		if(kvec[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for kvec[k]\n"); exit(11);}
  }
  total_k_vectors = (int*) calloc(sim.NB,sizeof(int));
		if(total_k_vectors == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for total_k_vectors\n"); exit(11);}
#endif
#endif

  /* ================================================================== */
  /* Allocate memory for Hessian matrix                                 */
  /* ================================================================== */
#ifdef HESSIAN
//hessian = (double ***) calloc(sim.NB,sizeof(double**));
//if(hessian == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for hessian\n"); exit(11);}
//for(int k=0; k<sim.NB; k++){
		hessian = (double**) calloc(3*box[0].boxns,sizeof(double*));
		if(hessian == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for hessian\n"); exit(11);}
		for(int i=0; i<(3*box[0].boxns); i++){
			hessian[i] = (double*) calloc(3*box[0].boxns,sizeof(double));
				if(hessian[i] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for hessian[i]\n"); exit(11);}
		}
//	}
#endif

  /* ================================================================== */
  /* Allocate memory for I/O buffers.                                   */
  /* ================================================================== */
  simul_buff = (char**) calloc(sim.NB,sizeof(char*));
	if(simul_buff == NULL){fprintf(stderr,"ERROR: cannot allocate memory for simul_buff\n"); exit(1);}
  simul_ptr = (char**) calloc(sim.NB,sizeof(char*));
	if(simul_ptr == NULL){fprintf(stderr,"ERROR: cannot allocate memory for simul_ptr\n"); exit(1);}

  ener_buff = (char**) calloc(sim.NB,sizeof(char*));
	if(ener_buff == NULL){fprintf(stderr,"ERROR: cannot allocate memory for ener_buff\n"); exit(1);} 
  ener_ptr = (char**) calloc(sim.NB,sizeof(char*));
	if(ener_ptr == NULL){fprintf(stderr,"ERROR: cannot allocate memory for ener_ptr\n"); exit(1);} 

  ord_buff = (char**) calloc(sim.NB,sizeof(char*));
	if(ord_buff == NULL){fprintf(stderr,"ERROR: cannot allocate memory for ord_buff\n"); exit(1);}
  ord_ptr = (char**) calloc(sim.NB,sizeof(char*));
	if(ord_ptr == NULL){fprintf(stderr,"ERROR: cannot allocate memory for ord_ptr\n"); exit(1);}
  
  swap_buff = (char**) calloc(sim.NB,sizeof(char*));
	if(swap_buff == NULL){fprintf(stderr,"ERROR: cannot allocate memory for swap_buff\n"); exit(1);}
  swap_ptr = (char**) calloc(sim.NB,sizeof(char*));
	if(swap_ptr == NULL){fprintf(stderr,"ERROR: cannot allocate memory for swap_ptr\n"); exit(1);}

  #if defined(SMD) || defined(NSMD)  
	smd_buff = (char**) calloc(sim.NB,sizeof(char*));
	  if(smd_buff == NULL){fprintf(stderr,"ERROR: cannot allocate memory for smd_buff\n"); exit(1);}
	smd_ptr = (char**) calloc(sim.NB,sizeof(char*));
	  if(smd_ptr == NULL){fprintf(stderr,"ERROR: cannot allocate memory for smd_ptr\n"); exit(1);}
  #endif

  simul_time = (time_t*) calloc(sim.NB,sizeof(time_t));
	if(simul_time == NULL){fprintf(stderr,"ERROR: cannot allocate memory for simul_time\n"); exit(1);}
  ener_time = (time_t*) calloc(sim.NB,sizeof(time_t));
	if(ener_time == NULL){fprintf(stderr,"ERROR: cannot allocate memory for ener_time\n"); exit(1);}
  ord_time = (time_t*) calloc(sim.NB,sizeof(time_t));
	if(ord_time == NULL){fprintf(stderr,"ERROR: cannot allocate memory for ord_time\n"); exit(1);}
  swap_time = (time_t*) calloc(sim.NB,sizeof(time_t));
	if(swap_time == NULL){fprintf(stderr,"ERROR: cannot allocate memory for swap_time\n"); exit(1);}

  #if defined(SMD) || defined(NSMD)  
    smd_time = (time_t*) calloc(sim.NB,sizeof(time_t));
	  if(smd_time == NULL){fprintf(stderr,"ERROR: cannot allocate memory for smd_time\n"); exit(1);}
  #endif

  #if defined(TWHAM) || defined(XWHAM)
   wham_time = (time_t*) calloc(sim.NB,sizeof(time_t));
	  if(wham_time == NULL){fprintf(stderr,"ERROR: cannot allocate memory for wham_time\n"); exit(1);}
  #endif

  for(int k=0; k<sim.NB; k++){
	simul_buff[k] = (char*) calloc(MAXBUFF,sizeof(char));
	  if(simul_buff[k] == NULL){fprintf(stderr,"ERROR: cannot allocate memory for simul_buff[k]\n"); exit(1);}
    simul_ptr[k] = simul_buff[k];

	ener_buff[k] = (char*) calloc(MAXBUFF,sizeof(char));
	  if(ener_buff[k] == NULL){fprintf(stderr,"ERROR: cannot allocate memory for ener_buff[k]\n"); exit(1);}
    ener_ptr[k] = ener_buff[k];

	ord_buff[k] = (char*) calloc(MAXBUFF,sizeof(char));
	  if(ord_buff[k] == NULL){fprintf(stderr,"ERROR: cannot allocate memory for ord_buff[k]\n"); exit(1);}
    ord_ptr[k] = ord_buff[k];

	swap_buff[k] = (char*) calloc(MAXBUFF,sizeof(char));
		if(swap_buff[k] == NULL){fprintf(stderr,"ERROR: cannot allocate memory for swap_buff[k]\n"); exit(1);}
	swap_ptr[k] = swap_buff[k];

    #if defined(SMD) || defined(NSMD)  
	  smd_buff[k] = (char*) calloc(MAXBUFF,sizeof(char));
	    if(smd_buff[k] == NULL){fprintf(stderr,"ERROR: cannot allocate memory for smd_buff[k]\n"); exit(1);}
  	  smd_ptr[k] = smd_buff[k];
    #endif
  }

  #if defined(TWHAM) || defined(XWHAM)
  #ifdef TWHAM
  pe_max = (double*)calloc(sim.NB,sizeof(double));
    if(pe_max == NULL){fprintf(stderr,"ERROR: cannot allocate memory for pe_max.\n"); exit(1);}

  pe_min = (double*)calloc(sim.NB,sizeof(double));
    if(pe_min == NULL){fprintf(stderr,"ERROR: cannot allocate memory for pe_min.\n"); exit(1);}
  #endif
  h_ebond = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_ebond == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_ebond.\n"); exit(1);}

  h_ebend = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_ebend == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_ebend.\n"); exit(1);}

  h_eurey = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_eurey == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_eurey.\n"); exit(1);}

  h_etors = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_etors == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_etors.\n"); exit(1);}

  h_eimpr = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_eimpr == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_eimpr.\n"); exit(1);}

  h_elj = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_elj == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_elj.\n"); exit(1);}

  h_eqq = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_eqq == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_eqq.\n"); exit(1);}

  h_pe = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_pe == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_pe.\n"); exit(1);}

  h_ke = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_ke == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_ke.\n"); exit(1);}

  h_etot = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_etot == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_etot.\n"); exit(1);}

  h_hel = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_hel == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_hel.\n"); exit(1);}

  h_dnc = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_dnc == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_dnc.\n"); exit(1);}

  h_con = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_con == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_con.\n"); exit(1);}

  h_con_2 = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_con_2 == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_con_2.\n"); exit(1);}

  h_x1 = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_x1 == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_x1.\n"); exit(1);}

  h_x2 = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_x2 == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_x2.\n"); exit(1);}

  h_rmsd = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_rmsd == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_rmsd.\n"); exit(1);}

  h_gyr = (tak_histogram**)calloc(sim.NB,sizeof(tak_histogram*));
    if(h_gyr == NULL){fprintf(stderr,"ERROR: cannot allocate memory for h_gyr.\n"); exit(1);}
  #endif //WHAM

  #ifdef ZHOU
  /*
  stresses = (double***)calloc(box[0].boxns,sizeof(double**));
                if(stresses == NULL){fprintf(stdout, "ERROR: cannot allocate memory for stresses\n"); exit(11);}
        for (int kk=0;kk<box[0].boxns;kk++)
                {stresses[kk] = (double**)calloc(3,sizeof(double*));
                        if(stresses[kk] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for stresses\n"); exit(11);}
                }
        for (int kk=0;kk<box[0].boxns;kk++)
                {for (int jj=0;jj<3;jj++)
                 {stresses[kk][jj] = (double*)calloc(3,sizeof(double));
                        if(stresses[kk][jj] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for stresses\n"); exit(11);}}
                }
  */

  //Allocate memory for stress invariants (output instaneous stress every block, for inspection only, really, to make sure everything looks OK)
  stress_i = (struct stress_type**)calloc(sim.NB,sizeof(struct stress_type*));
  if(stress_i == NULL) {fprintf(stdout, "ERROR: cannot allocate memory for stress_i\n"); exit(11);}

  for (int nn=0;nn<sim.NB;nn++)
        {stress_i[nn] = (struct stress_type*)calloc(box[nn].boxns,sizeof(struct stress_type));
                if(stress_i[nn] == NULL){fprintf(stdout, "ERROR: cannot allocate memory for stress_i\n"); exit(11);}
        }               

  //Allocate memory for cumulative stress invariants (summed, then zeroed after each block)
  stress_isum = (struct stress_type**)calloc(sim.NB,sizeof(struct stress_type*));
  if(stress_isum == NULL) {fprintf(stdout, "ERROR: cannot allocate memory for stress_isum\n"); exit(11);}

  for (int nn=0;nn<sim.NB;nn++)
        {stress_isum[nn] = (struct stress_type*)calloc(box[nn].boxns,sizeof(struct stress_type));
                if(stress_isum[nn] == NULL){fprintf(stdout, "ERROR: cannot allocate memory for stress_isum\n"); exit(11);}
        }

  //Allocate memory for stress invariants: SUM of (AVERAGES over each block)
  stress_ibarsum = (struct stress_type**)calloc(sim.NB,sizeof(struct stress_type*));
  if(stress_ibarsum == NULL) {fprintf(stdout, "ERROR: cannot allocate memory for stress_ibarsum\n"); exit(11);}

  for (int nn=0;nn<sim.NB;nn++)
        {stress_ibarsum[nn] = (struct stress_type*)calloc(box[nn].boxns,sizeof(struct stress_type));
                if(stress_ibarsum[nn] == NULL){fprintf(stdout, "ERROR: cannot allocate memory for stress_ibarsum\n"); exit(11);}
        }

  //Allocate memory for stress invariants: SUM of (AVERAGES over each block)^2
  stress_ibarsumsq = (struct stress_type**)calloc(sim.NB,sizeof(struct stress_type*));
  if(stress_ibarsumsq == NULL) {fprintf(stdout, "ERROR: cannot allocate memory for stress_ibarsumsq\n"); exit(11);}

  for (int nn=0;nn<sim.NB;nn++)
        {stress_ibarsumsq[nn] = (struct stress_type*)calloc(box[nn].boxns,sizeof(struct stress_type));
                if(stress_ibarsumsq[nn] == NULL){fprintf(stdout, "ERROR: cannot allocate memory for stress_ibarsumsq\n"); exit(11);}
        }

  //Allocate memory for stresses tensor [0..boxes-1][0..sites-1][0..2][0..2]
  stresses = (double****)calloc(sim.NB,sizeof(double***));
  if(stresses == NULL) {fprintf(stdout, "ERROR: cannot allocate memory for stresses\n"); exit(11);}

  for (int nn=0;nn<sim.NB;nn++)
        {stresses[nn] = (double***)calloc(box[nn].boxns,sizeof(double**));
                if(stresses[nn] == NULL){fprintf(stdout, "ERROR: cannot allocate memory for stresses\n"); exit(11);}
        }

  for (int nn=0;nn<sim.NB;nn++)
        {for (int kk=0;kk<box[nn].boxns;kk++)
                {stresses[nn][kk] = (double**)calloc(3,sizeof(double*));
                        if(stresses[nn][kk] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for stresses\n"); exit(11);}
                }
        }

  for (int nn=0;nn<sim.NB;nn++)
        {for (int kk=0;kk<box[nn].boxns;kk++)
                {for (int jj=0;jj<3;jj++)
                        {stresses[nn][kk][jj] = (double*)calloc(3,sizeof(double));
                                if(stresses[nn][kk][jj] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for stresses\n"); exit(11);}
                        }
                }
        }

  //Allocate memory for cumulative stresses tensor [0..boxes-1][0..sites-1][0..2][0..2]
  stresses_sum = (double****)calloc(sim.NB,sizeof(double***));
  if(stresses_sum == NULL) {fprintf(stdout, "ERROR: cannot allocate memory for stresses_sum\n"); exit(11);}

  for (int nn=0;nn<sim.NB;nn++)
        {stresses_sum[nn] = (double***)calloc(box[nn].boxns,sizeof(double**));
                if(stresses_sum[nn] == NULL){fprintf(stdout, "ERROR: cannot allocate memory for stresses_sum\n"); exit(11);}
        }

  for (int nn=0;nn<sim.NB;nn++)
        {for (int kk=0;kk<box[nn].boxns;kk++)
                {stresses_sum[nn][kk] = (double**)calloc(3,sizeof(double*));
                        if(stresses_sum[nn][kk] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for stresses_sum\n"); exit(11);}
                }
        }

  for (int nn=0;nn<sim.NB;nn++)
        {for (int kk=0;kk<box[nn].boxns;kk++)
                {for (int jj=0;jj<3;jj++)
                        {stresses_sum[nn][kk][jj] = (double*)calloc(3,sizeof(double));
                                if(stresses_sum[nn][kk][jj] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for stresses_sum\n"); exit(11);}
                        }
                }
        }

        #ifdef TRR

        //Allocate memory for accumulating mean stress for use in movies (to be averaged over each movie frame)
        //Currently only allocated for box 0 only!
        stress_movie_psum = (double*)calloc(box[0].boxns,sizeof(double));
        #endif /* TRR */


#endif /* ZHOU */


}
