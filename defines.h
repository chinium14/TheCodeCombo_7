/* ======================================================================== */
/* defines.h                                                                */
/*                                                                          */
/*		This header file contains the constants, structures, and libraries  */
/* used in the program.  The structures defined here are global variables   */
/* accessible to any subroutine where #include defines.h is put at the top  */
/* of the file.  							    */	
/*                                                                          */
/* Note: To make a variable global in this manner, #ifdef MAIN is used to   */
/* include the key word "extern" before the structure is defined.  MAIN is  */
/* only defined in the subroutine main(). Thereafter, extern tells the      */
/* compiler that this variable has been declared elsewhere, but to make it  */
/* accessible to the current subroutine.                                    */
/* ======================================================================== */



  /* ================================================================== */
  /*                                                                    */
  /* Libraries and Constants                                            */
  /*                                                                    */
  /* ================================================================== */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef MKL
#include <mkl.h>
#endif
#include <math.h>
#include <time.h>
//#include "malloc.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <signal.h>

#ifdef MPI
#include <mpi.h>
#endif

#ifdef SPME
#include <fftw3.h>
#endif

#ifdef SCRIPT
#include <ctype.h>
#endif

#ifndef WIN
#include <unistd.h>
//#include <sys/stat.h>
#endif
#ifdef WIN
#include <Winsock2.h>
#include <direct.h>
#define for if(0);else for
#endif

#if defined(TWHAM) || defined(XWHAM)
#include "wham.h"
#endif

#define PI      3.141592653589793238462643383			/* pi			                   */
#define KB      1.380662E-23					/* Boltzmann constant [J/K]	           */
#define NA      6.022045E23					/* Avogadro  constant [/mol]	           */
#define RG      (KB*NA)
#define AU      1.6605655E-27					/* atomic massunit    [kg]           	   */
#define PLANCK  6.626176E-34					/* Planck constant    [Js]                 */
#define ELEC    8.987551789e19					/* [J.A/C^2] (4.Pi.e)^-1                   */
#define ELEQ    1.602176462e-19					/* [C]                                     */
#define DIPOLE  3.33564e-20					/* [C.A] 1 Debye = 3.33564e-20 C.A         */
#define QFACTOR ELEC*ELEQ*ELEQ*NA*1.0E-3                        /* conversion factor for charges to kJ/mol */
#ifndef BGO
#define E14FAC  0.4
#else
#define E14FAC 1.0
#endif	
			
#define DIM     3								/* dimension                               */
#ifndef MPI
#define NBOXES  31								/* max number of boxes                     */
#endif
#ifdef MPI
#define NBOXES	1
#endif
#define NTYPES  8								/* max type of molecules                   */
#define MEXTRA  100								/* size for extra memory allocation        */
#define MAXT    80								/* max number of characters                */
#define NPRIME  1009							/* number of entries in the hash table     */
#ifdef SASA
#define RPROBE  1.40							/* probe radius for implicit solvent SASA  */ 
#endif
#define ELIMIT  1.0E10							/* arbitrary large number - energy test    */
#define MAXBUFF 40000							/* Max size of buffer to hold output before writing to disk */
#define OUTBUFF 200							/* Size at which the buffer will be written to disk */

#ifdef SPME
#define index3D(I,J,K)  ((K) + (J)*GRIDZ + (I)*GRIDY*GRIDZ)
#endif

#define greater(A,B) ((A)>(B) ? (A):(B))
#define lesser(A,B) ((A)<(B) ? (A):(B))
#define oddeven(A)  ((A%2) ? 1 : 0)              /* even = 0, odd = 1*/



/*
#ifndef NEUTRAL
#define NUM_ATOM	30
#define NUM_BOND	70
#define	NUM_BEND	191
#define NUM_TORS	41
#define NUM_IMPR    43  
#endif
*/

#define NUM_LJSET	300//200							/* This number controls the size of ljset  */

/*#if defined(NEUTRAL)
  #define NUM_ATOM	31							
  #define NUM_BOND	87
  #define NUM_BEND	212
  #define NUM_TORS	47
  #define NUM_IMPR  50
#elif defined(SUG)
  #define NUM_ATOM	161						
  #define NUM_BOND	297
  #define NUM_BEND	738
  #define NUM_TORS	1147
  #define NUM_IMPR  76
  #define NUM_UREY	199 
#elif defined(DNA_GOLIK)
  #define NUM_ATOM	7					
  #define NUM_BOND	5
  #define NUM_BEND	6
  #define NUM_TORS	6
  #define NUM_IMPR  1
  #define NUM_UREY	1
#else
  #define NUM_ATOM	157		
  #define NUM_BOND	295
  #define NUM_BEND	740
  #define NUM_TORS	1156
  #define NUM_IMPR  76
  #define NUM_UREY	202  
#endif
*/

//#ifndef NEUTRAL
//#ifndef SUG									/* If NEUTRAL is not defined, the CHARMM22 */
//#define NUM_ATOM	157//155//131							/* parameters are used.                    */
//#define NUM_BOND	295//291//250
//#define	NUM_BEND	741//728//725//622
//#define NUM_TORS	1152//1145//1143//1049
//#define NUM_IMPR    76//73
//#define NUM_UREY	202//199//145  
//#endif
//#ifdef SUG									/* If NEUTRAL is not defined, the CHARMM22 */
//#define NUM_ATOM	161//131							/* parameters are used.                    */
//#define NUM_BOND	297//250
//#define	NUM_BEND	738//622
//#define NUM_TORS	1147//1049
//#define NUM_IMPR    76//73
//#define NUM_UREY	199//145  
//#endif
//#endif



//#ifdef NEUTRAL									/* new atom type CR is added and hence the  */
//#define NUM_ATOM	31							/* correspomding bonds, bends and torsions	*/
//#define NUM_BOND	87
//#define	NUM_BEND	212
//#define NUM_TORS	47
//#define NUM_IMPR    50
//#endif

#ifdef DOS
#define NBINS		1000
#endif
#ifdef MMDOS
#define NBINS		1000
#endif
#ifdef CTDOS
#define NBINS		1000
#endif
#ifdef XEDOS
#define NBINS		1000
#endif

#ifdef FX_EDOS
#define NBINS		1000
#endif

#ifdef TDXEDOS
#define NBINS1		200
#define NBINS2		50
#endif

#ifndef MAIN
extern
#endif
float STOP_F;

#ifndef MAIN
extern
#endif
float RESET_F;

#ifndef MAIN
extern
#endif
float MERGE_F;

#ifndef MAIN
extern
#endif
int MERGE_N;




#ifdef ROTN
#define _AXIS 0                                 /* If rotatation is defined, align to x (0),y (1), or z (2) axis. */
#endif
#ifdef CLIST
#define CELL_LENGTH_X 21.0
#define CELL_LENGTH_Y 21.0
#define CELL_LENGTH_Z 21.0
#endif
#ifdef FLIM
#define FORCE_LIMIT 3000   /*kJ/mol/ang*/
#endif
#define MINIMG(x,box_xh,box_x) (((x) > (box_xh)) ? ((x) - (box_x)) : (((x) < (-(box_xh))) ? ((x) + (box_x)) : x))
#define CROSS(A,B,C)  (C[0] = (((A[1]) * (B[2])) - ((B[1]) * (A[2])))); \
                      (C[1] = (((B[0]) * (A[2])) - ((B[2]) * (A[0])))); \
                      (C[2] = (((A[0]) * (B[1])) - ((B[0]) * (A[1]))))
#define DET(A) (((A[0][0]) * (A[1][1]) * (A[2][2])) - ((A[0][0]) * (A[1][2]) * (A[2][1])) - ((A[1][0]) * (A[0][1]) * (A[2][2])) + ((A[1][0]) * (A[0][2]) * (A[2][1])) + ((A[2][0]) * (A[0][1]) * (A[1][2])) - ((A[2][0]) * (A[0][2]) * (A[1][1])))

#define INV(I,II)   (II[0][0]  = (((I[1][1]) * (I[2][2]) - (I[1][2]) * (I[2][1])) / DET(I)));  \
                    (II[0][1] = (-((I[0][1]) * (I[2][2]) - (I[0][2]) * (I[2][1])) / DET(I)));  \
					(II[0][2] = (-((-I[0][1]) * (I[1][2]) + (I[0][2]) * (I[1][1])) / DET(I))); \
					(II[1][0] = (((-I[1][0]) * (I[2][2]) + (I[1][2]) * (I[2][0])) / DET(I)));  \
					(II[1][1] = (((I[0][0]) * (I[2][2]) - (I[0][2]) * (I[2][0])) / DET(I)));   \
					(II[1][2] = (-((I[0][0]) * (I[1][2]) - (I[0][2]) * (I[1][0])) / DET(I)));  \
					(II[2][0] = (-((-I[1][0]) * (I[2][1]) + (I[1][1]) * (I[2][0])) / DET(I))); \
					(II[2][1] = (((-I[0][0]) * (I[2][1]) + (I[0][1]) * (I[2][0])) / DET(I)));  \
				    (II[2][2] = (((I[0][0]) * (I[1][1]) - (I[0][1]) * (I[1][0])) / DET(I)))
  /* ================================================================== */
  /*                                                                    */
  /* Strutures of Molecular Variables                                   */
  /*                                                                    */
  /* ================================================================== */
  /* ----------------------------------------- */
  /* Globals for number of atoms types, bonds  */
  /* bends, ureys, torsions, impropers that    */
  /* are read in read_atom from the param file */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
char param_file_name[50];

#ifndef MAIN
extern
#endif
int NUM_ATOM;

#ifndef MAIN
extern
#endif
int NUM_BOND;

#ifndef MAIN
extern
#endif
int NUM_BEND;

#ifndef MAIN
extern
#endif
int NUM_TORS;

#ifndef MAIN
extern
#endif
int NUM_IMPR;

#ifndef MAIN
extern
#endif
int NUM_UREY;

#ifndef MAIN
extern
#endif
int NUM_NBFIX;

  /* ----------------------------------------- */
  /* atoms                                     */
  /* Contains information about each site.     */
  /* The site symbol, name, position, charge   */
  /* molecule number, and id number are here.  */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct atoms {
  char   t[5];          /* symbol for site           */
  char   name[5];
  double x;             /* x-position of site        */
  double y;             /* y-position of site        */
  double z;             /* z-position of site        */
  double q;             /* site charge               */
  int n;                /* index for residue number  */
  int molec;            /* index for site number     */
  int atomid;						
} **atom, **atnopbc, **atom_temp, **atnopbc_temp; 

  /* ----------------------------------------- */
  /* res                                       */
  /* Contains information about residue: the   */
  /* residue type, number of sites in the      */
  /* residue, and the number of chains.        */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct res {
  int type;
  int Nsite;
  int chain;  
} **residue;

  /* ----------------------------------------- */
  /* mol_info                                  */
  /* Contains information about each type of   */
  /* molecule: total number of molecule i,     */
  /* the number of sites in molecule i, the    */
  /* number of residues in molecule i, and the */
  /* molecular weight of molecule i.           */
  /* ----------------------------------------- */
#ifndef MAIN
extern 
#endif
struct mol_info {
  int    Ni;                           /* total number of molecules i      */
  int    Nsite;                        /* number of sites of molecule i    */
  int	 Nres;                         /* number of residues in molecule i */
  double MW;                           /* molecular weight of molecule i   */
} *mol;

  /* ----------------------------------------- */
  /* mol_id                                    */
  /* Contains information about each specific  */
  /* molecule in the system.  It gives the     */
  /* beginning site number of each molecule    */
  /* and the number of molecules of interest   */
  /* to perform special moves on, e.g. the     */
  /* proteins in the system rather than the    */
  /* water might want to have a pivot move.    */
  /* ----------------------------------------- */
#ifndef MAIN
extern 
#endif
struct mol_id {
  int   *fsite;                          /* first site number of molecule i                 */
  int   *lsite;							 /* last site number of molecule i                  */
  int   *fres;							 /* beginning residue number of each molecule i     */
  int	*lres;							 /* last residue number in each molecule            */
  int	Nsolute;						 /* Number of molecules to perform fancy MC move on */

} molec;

  /* ----------------------------------------- */
  /* Nmolall = Total number of molecules in    */
  /*	       all the boxes.                  */
  /* Nsitall = Total number of sites in all    */
  /*		   the boxes.                      */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
int Nmolall, Nsitall;


  /* ----------------------------------------- */
  /* nlist                                     */
  /* The neighborlist structure. The index on  */
  /* count runs over all sites.  The indicies  */
  /* on list run over all sites and all        */
  /* neighbors.                                */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct nlist_def{
  int *count;
  int **list;
#if defined(STYPE) || defined(DNA_GOLIK)
  int **xflag;
#endif
#ifdef DNA_GOLIK
  double **sig;
#endif
#ifdef SASA
  int *count_sasa;
  int **list_sasa;
#ifdef STYPE
  int **xflag_sasa;
#endif //STYPE
#endif//SASA
} *nlist
#ifdef EWALD
, *nlist_ew
#endif
;



#ifndef MAIN

extern

#endif

  int max_num_neighbors;

#ifdef SLIST
#ifndef MAIN
extern
#endif
struct slist_def{
  int *count;
  int **list;
} *slist;

#endif//SLIST
#ifdef DLIST
  /* ----------------------------------------- */
  /* rx0 ry0 rz0                               */
  /* coordiantes when neighborlist was last    */
  /* updated.                                  */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
	double **rx0;
#ifndef MAIN
extern
#endif
	double **ry0;
#ifndef MAIN
extern
#endif
	double **rz0;
#ifndef MAIN
extern
#endif
	int *nl_flag;

#endif

  /* ----------------------------------------- */
  /* clist                                     */
  /* The cell list structures                  */
  /* ----------------------------------------- */
#ifdef CLIST
#ifndef MAIN
extern
#endif
int **head;

#ifndef MAIN
extern
#endif
int **clist;

#ifndef MAIN
extern
#endif
int **cell_no;
 
#ifndef MAIN
extern
#endif
struct clistn_struct{
  int Mx;
  int My;
  int Mz;
  int ncell;
  int *count;
  int **list;
}*clistn;


#endif //CLIST

  /* ----------------------------------------- */
  /* veloc                                     */
  /* This structure contains the x,y,z,        */
  /* components of the velocities and forces.  */
  /* ff_short and ff_long are used for mts     */
  /* integration.                              */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif 
struct veloc {
  double x;
  double y;
  double z;
} **vv, **uu, **ff, **vcm, **ff_short, **ff_long, **ff_temp, **vv_temp;             

  /* ----------------------------------------- */
  /* energy                                    */
  /* This structure contains the information   */
  /* on all the energy of the system.  Both    */
  /* kinetic and potential energy are          */
  /* included as well as the contributions to  */
  /* the potential energy and the extended     */
  /* hamiltonian for Nose Hoover dynamics.     */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct energy {
  double total;
  double totals;
  double poten;
  double potens;
  double kinet;
  double kinetmol;
  double ken[6];
  double kenmol[6];
  double nbond;
  double nbonds;
  double bond;
  double bend;
  double tors;
  double impr;
  double coulomb;
  double H;
#ifdef SASA
  double esasa;
  double esasab;
#endif
#ifdef EWALD
  double rewald; //real part of Ewald electrostatics
  double rewalds; //shifted real part of Ewald electrostatics
  double kewald; //recipricol part of Ewald electrostatics
  double sewald; //self interaction part of Ewald electrostatics
  double iewald; //intra-molecular self energy
  double tewald; //total Ewald electrostatics
#endif
#ifndef NEUTRAL
  double urey;
#endif
} en[NBOXES], en_temp[NBOXES];

  /* ----------------------------------------- */
  /* virial                                    */
  /* This structure contains the information   */
  /* on all the contributions to the atomic    */
  /* virial as well as the total viral.  The   */
  /* virial tensor is also here.               */
  /* ----------------------------------------- */
#ifdef PRESSURE
#ifndef MAIN
extern
#endif
struct virial {
  double nbond[6];
  double nbondmol[6];
  double bond[6];
  double bend[6];
#ifndef NEUTRAL
  double urey[6];
#endif  
  double tors[6];
  double impr[6];
  double sum[6];
  double summol[6];
  double tnbond;
  double tnbondmol;
  double tbond;
  double tbend;
#ifndef NEUTRAL
  double turey;
#endif
  double ttors;
  double timpr;
#ifdef SASA
  double tsasa;
#endif
#ifdef EWALD
  double tewald;
#endif
  double total;
  double totalmol;
  double lrc;
#ifdef SASA
  double sasa[6];
#endif
  #ifdef EWALD
  double ewald_real[6];
  double ewald_kspace[6];
  double ewald_intra[6];
#endif

} pvir[NBOXES], pvir_temp[NBOXES]
#ifdef MC
,pviro[NBOXES]
#endif
;
#endif//Pressure


  /* ----------------------------------------- */
  /* inter                                     */
  /* This structure contains the information   */
  /* on the LJ parameters, mass, and charge    */
  /* of each site.  The indices on the arrays  */
  /* run over the box number and the site      */
  /* number.                                   */
  /* ----------------------------------------- */

#ifndef MAIN
extern
#endif
struct inter {
  double mas;                       /* mass of each atom/site              */
  double sig;                       /* sigma parameter for LJ potential    */
  double eps;                       /* epsilon parameter for LJ potential  */
  double qq;                        /* site charge for Coulombic potential */  
  int   id;                     /* symbol for atom/site                */ 
#ifdef SASA
  double pab;  
#endif
#ifndef	STYPE
#ifdef EWALD
  int bond_flag;
  int bend_flag;
  int tors_flag;
  int exNB_flag;
#endif
#endif//STYPE
#ifdef DNA_GOLIK
  double rc2_non;
#endif

}
#ifndef STYPE
#ifdef S14
**pott14,
#endif
**pota, **pott;
#endif// STYPE

#ifdef STYPE
#ifdef S14
**pott14, **ljset14,
#endif
**pota, **pott, **ljset;
#endif// STYPE


  /* ================================================================== */
  /*                                                                    */
  /* The following 9 structures contain the information used to         */
  /* calculate the energies and forces on the system. They are assigned */
  /* by comparing the atoms/sites involved in system with the different */
  /* read in read_atom from the simul22.param file.  The sturctures     */
  /* that contain the infor from the simul22.param file are the 5       */
  /* just below the 9 here.                                             */
  /*                                                                    */
  /* ================================================================== */
  /* ----------------------------------------- */
  /* interb                                    */
  /* This structure contains the specific site */
  /* site interactions. The indice on ljset    */
  /* runs over the box number, while the       */
  /* indicies on its pot member run over sites */
  /* a and b. This structure holds the sigma   */
  /* and epsilon values calculated from the    */
  /* mixing rules for the values stored in     */
  /* pott above.  Values are set in setlj().   */
  /* ----------------------------------------- */
#ifndef STYPE
#ifndef MAIN   
extern
#endif
struct interb {
 struct  inter **pot;
} *ljset
#ifdef S14
 ,*ljset14
#endif
;
#endif


  /* ----------------------------------------- */
  /* bonds                                     */
  /* This structure contains the information   */
  /* on each bond in the system.  The site     */
  /* number of the atoms involved in the bond  */
  /* and the force field parameters are here.  */
  /* the indicies run over the number of       */
  /* boxes and the number of bonds.            */
  /* Values are set in read_topology.          */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct bonds {
  int a;
  int b;
  double req;
  double krbond;
} **bonda, **bond;

  /* ----------------------------------------- */
  /* bondN[k] = Total number of bonds in box k */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif 
int bondNa[NTYPES], bondN[NBOXES];


  /* ----------------------------------------- */
  /* urey-bradleys                             */
  /* This structure contains the information   */
  /* on each Urey-Bradley 1-3 distance.        */
  /* ----------------------------------------- */
#ifndef NEUTRAL
#ifndef MAIN
extern
#endif
struct ureys {
  int a;
  int b;
  double Seq;
  double k;
} **ureya, **urey;

  /* ----------------------------------------- */
  /* ureyN[k] = Total number of ureys in box k */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif 
int ureyNa[NTYPES], ureyN[NBOXES];
#endif//NEUTRAL

  /* ----------------------------------------- */
  /* bends                                     */
  /* This structure contains the information   */
  /* on each bend in the system.  The site     */
  /* number of the atoms involved in the bend  */
  /* and the force field parameters are here.  */
  /* the indicies run over the number of       */
  /* boxes and the number of bends.            */
  /* Values are set in read_topology.          */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct bends {
  int a;
  int b;
  int c;
  double angeq;
  double krbend;
} **benda, **bend;

  /* ----------------------------------------- */
  /* bendN[k] = Total number of bends in box k */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
int bendNa[NTYPES], bendN[NBOXES];


  /* ----------------------------------------- */
  /* torsions                                  */
  /* This structure contains the information   */
  /* on each torsion in the system.  The site  */
  /* number of the atoms involved in the       */
  /* dihedral and the force field parameters   */
  /* are here. The indicies run over the       */
  /* number of boxes and the number of         */
  /* torsions.                                 */
  /* Values are set in read_topology.          */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct torsions {
  int a;
  int b;
  int c;
  int d;
  int count;
  double phi;
  double psi;
  double theta;
  double phia;
  double psia;
  double thetaa;
  char phitag;
  char psitag;
  char thetag;
  double kphi[6];
  double nphi[6];
  double delphi[6];
} **torsa, **tors;

  /* ----------------------------------------- */
  /* torsN[k] = Total number of torsions in    */
  /* box k.                                    */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
int torsNa[NTYPES], torsN[NBOXES];

  /* ----------------------------------------- */
  /* improper torsions                         */
  /* This structure contains the information   */
  /* on each improper torsion in the system.   */
  /* The site number of the atoms involved in  */
  /* the improper torsion and the force field  */
  /* parameters are here. The indicies run     */
  /* over the     number of boxes and the      */
  /* number of improper torsions.              */
  /* Values are set in read_topology.          */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct impropers {
  int a;
  int b;
  int c;
  int d;
  double kimpr;
  double angeq;
} **impra, **impr;

  /* ----------------------------------------- */
  /* imprN[k] = Total number of improper       */
  /* torsions in  box k.                       */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
int imprNa[NTYPES], imprN[NBOXES];

  /* ----------------------------------------- */
  /* excinter                                  */
  /* This structure contains the information   */
  /* on each modified interaction in the       */
  /* system.  Note: These interactions are not */
  /* used in the current version.              */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct excinter {
  int a;
  int b;
#ifdef STYPE
  int xflag;
#endif
  double eps;
  double sig;
  double qa;
  double qb;
} **xinta, **xint;

  /* ----------------------------------------- */
  /* xintN[k] = Total number of modified       */
  /* interactions in box k.                    */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
int xintNa[NTYPES], xintN[NBOXES];

  /* ----------------------------------------- */
  /* 1-4 Interactions                          */
  /* This structure contains the information   */
  /* on each 1-4 interaction in the system.    */
  /* The CHARMM force field does not have      */
  /* parameters for every torsion so this is   */
  /* then the complete list of 1-4             */
  /* interactions.                             */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct interactions1_4 {
  int a;
  int b;
  int c;
  int d;
} **in14a, **in14;

  /* ----------------------------------------- */
  /* in14N[k] = Total number of 1-4            */
  /* interactions in box k.                    */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
int in14Na[NTYPES], in14N[NBOXES];

  /* ----------------------------------------- */
  /* Hydrogen bond donors                      */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct donors {
  int a;
  int b;
} **hdona, **hdon;

  /* ----------------------------------------- */
  /* ndonN[k] = Total number donors.           */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
int hdonNa[NTYPES], hdonN[NBOXES];

  /* ----------------------------------------- */
  /* Hydrogen bond acceptors                   */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct acceptors {
  int a;
  int b;
} **hacca, **hacc;

  /* ----------------------------------------- */
  /* haccN[k] = Total number acceptors.        */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
int haccNa[NTYPES], haccN[NBOXES];

  /* ----------------------------------------- */
  /* Fixed non-bonded exclusions(ring members) */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct NBexclusions {
  int a;
  int b;
} **exNBa, **exNB;

  /* ----------------------------------------- */
  /* Fixed non-bonded exclusions(ring members) */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct nbfixes{
  int eps;
  int sig;
}**nbfix;

  /* ----------------------------------------- */
  /* exNBN[k] = Total number fixed exclusions. */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
int exNBNa[NTYPES], exNBN[NBOXES];
  /* ----------------------------------------- */
  /* bondprop                                  */
  /* This structure contains the information   */
  /* on each type of bond in the system as     */
  /* read in read_atom() from the              */
  /* simul22.param file.                       */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct bondprop {
  int type1;
  int type2;
  double kbond;
  double eqbond;
} *bond_prop;


  /* ----------------------------------------- */
  /* ureyprop                                  */
  /* This structure contains the information   */
  /* on each type of bond in the system as     */
  /* read in read_atom() from the              */
  /* simul22.param file.                       */
  /* ----------------------------------------- */
#ifndef NEUTRAL
#ifndef MAIN
extern
#endif
struct ureyprop {
  int type1;
  int type2;
  int type3;
  double k;
  double Seq;
} *urey_prop;
#endif //NEUTRAL

  /* ----------------------------------------- */
  /* bendprop                                  */
  /* This structure contains the information   */
  /* on each type of bend in the system as     */
  /* read in read_atom() from the              */
  /* simul22.param file.                       */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct bendprop {
  int type1;
  int type2;
  int type3;
  double kbend;
  double eqbend;
}*bend_prop;


  /* ----------------------------------------- */
  /* dihedralprop                              */
  /* This structure contains the information   */
  /* on each type of torsion in the system as  */
  /* read in read_atom() from the              */
  /* simul22.param file.                       */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct dihedralprop {
  int type1;
  int type2;
  int type3;
  int type4;
  double kphi;
  double nphi;
  double delphi;
} *dihedral_prop;

  /* ----------------------------------------- */
  /* improperprop                              */
  /* This structure contains the information   */
  /* on each type of improper torsion in the   */
  /* system as read in read_atom() from the    */
  /* simul22.param file.                       */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct improperprop {
  int type1;
  int type2;
  int type3;
  int type4;
  double kimpr;
  double angeq;
} *improper_prop;

#ifndef MAIN
extern
#endif
struct nbfixprop{
  int type1;
  int type2;
  double sig;                       /* sigma parameter for LJ potential    */
  double eps;                       /* epsilon parameter for LJ potential  */
}*nbfix_prop;



  /* ----------------------------------------- */
  /* structure                                 */
  /* This structure contains the information   */
  /* about the end to end distance, helicity,  */
  /* and hydrogen bonds.                       */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct struct_quant{
  double d_nc;
  double d_nca;
  double d_ncb;
  double hel;
  double hela;
  double helb;
  int hbonds;
  double con; //native contacts
  double cona;
  double conb;
  double con_2;
  double force_1;
  double force_2;
  double gyr;
  double gyra;
  double gyrb;
  double rmsd;
  double rmsda;
  double rmsdb;
  double x1;
  double x1a;
  double x1b;
  double x2;
  double x2a;
  double x2b;
} ordparam[NBOXES];

  /* ----------------------------------------- */
  /* structure                                 */
  /* This structure contains the information   */
  /* about the native contacts.                */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct struct_contacts{
	int n;
	int *a;		//site 1
	int *b;		//site 2
	double *d;	//distance
	int *code; // native contact code
} *contacts, *contacts_2;
  /* ----------------------------------------- */
  /* structure                                 */
  /* This structure contains the information   */
  /* about the radius of gyration.             */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct struct_gyr{
	int n;
	int *sites;		//sites to be included in radius of gyration calculation
} *gyration;

  /* ----------------------------------------- */
  /* epstmp                                    */
  /* This structure contains the epsilon value */
  /* of the LJ potential for each type of site */
  /* as read in reac_atom() from the           */
  /* simul22.param file.                       */
  /* The first index:                          */
  /*					0 = Regular NBOND      */
  /*					1 = Special 14 NBOND   */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
double epstmp[2][2000];

  /* ----------------------------------------- */
  /* sigtmp                                    */
  /* This structure contains the sigma value   */
  /* of the LJ potential for each type of site */
  /* as read in reac_atom() from the           */
  /* simul22.param file.                       */
  /* The first index:                          */
  /*					0 = Regular NBOND      */
  /*					1 = Special 14 NBOND   */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
double sigtmp[2][2000];

  /* ----------------------------------------- */
  /* atomtype                                  */
  /* This structure contains each type of      */
  /* atom/site read in read_atom.              */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct atomtype {
	char name[5];
	int atomid;
	double mass;
} *atom_type;


  /* ----------------------------------------- */
  /* residuetype                               */
  /* This structure contains each type of      */
  /* residue read in read_atom.                */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct residuetype {
	char name[4];
	char aminoid[2];
	int residueid;
} residue_type[20];


  /* ----------------------------------------- */
  /* SimData                                   */
  /* This structure contains information on    */
  /* the simulation as read in read sim.       */
  /* ----------------------------------------- */
#ifndef MAIN
extern 
#endif
struct SimData {
  int    ID;                   /* ensemble id mode                     */
  int	 ID2;                  /* integration id mode                  */
  int    NB;                   /* number of simulation boxes           */
  int    NC;                   /* number of components                 */
  double T[NBOXES];            /* temperature [K]                      */
  double Ttau[NBOXES];         /* thermostat coupling time [ps]        */
  double dtemp[NBOXES];        /* temperature increment per time step  */
  double P[NBOXES][3];         /* pressure  [kPa]                      */
  double Ptau[NBOXES][3];      /* barostat coupling time [ps]          */
  double Pcomp[NBOXES][3];     /* isothermal compressibility [1/bar]   */
  double dpress[NBOXES][3];    /* pressure increment per time step     */
  double kT[NBOXES];           /* kT = KB * T                          */
  unsigned long   cyc_eq;               /* equilibration cycles                 */
  unsigned long   cyc_pr;               /* production cycles                    */
  unsigned long   blockc;               /* blocking of configurations           */
  unsigned long   blockd;               /* blocking of data                     */
  unsigned long	 blockt;			   /* saving trajectory					   */
  unsigned long   drift;                /* check for drift of momentum          */ 
  double rc;                   /* cutoff radius [A]                    */
  double rclist;               /* cutoff radius for neighbor-list[A]   */

  double rc_ew;				   /* cutoff radius for ewald [A]		   */

  double rclist_ew;			   /* cutoff radius for ewald nlist [A]    */
  unsigned long    nlist;                /* frequency for updating neighbor list */
  double dt;                   /* time step [fs]                       */
  double dtlong;               /* time of long time step [fs]          */
  int    nsteps;               /* number of short time steps/long step */
  int	 xRESPA[NBOXES];       /* flag for extended RESPA              */
  double kappa[NBOXES];        /* Ewald parameter for Gaussian width   */
  double Na_con[NBOXES];       /* Sodium concentration */
  double epsRF[NBOXES];        /* Dielectric constant */
  int	 kmax[NBOXES][3];      /* Ewald parameter for reciprical space */
  int	 max_k_square[NBOXES]; /* Ewald parameter for reciprical space */
  int    order;          /* order for B-spline in PME calculation      */
  int    grid[3];        /* XYZ dimensions of grid for PME calculation */
} sim;


  /* ----------------------------------------- */
  /* BoxData                                   */
  /* This structure contains information on    */
  /* each box.                                 */
  /* ----------------------------------------- */
#ifndef MAIN
extern 
#endif
struct BoxData {
  int    boxn;					/* total number of molecules in box i     */
  int    boxns;					/* total number of sites in box i         */
  int    boxnres;				/* total number of residues in box i      */
  double lx;					/* x box length			          */
  double ly;					/* y box length				  */
  double lz;					/* z box length				  */
  double least_boxl;
  double hx;					/* one-half x box length                  */
  double hy;					/* one-half y box length		  */
  double hz; 					/* one-half z box length 		  */
  double least_boxh;
  double vol;					/* box volume                             */
  double rc2;					/* cutoff radius (squared)                */
  double rc4;					/* cutoff radius ^4						  */	
  double weight;				/* total weight                           */
  double dens;					/* density                                */
  double temp;					/* instantaneous temperature              */
  double press;					/* instantaneous pressure - atomic        */
  double pressmol;				/* instantaneous pressure - molecular     */
  double cpress[6];				/* pressure tensor components - atomic    */
  double cpressmol[6];			/* pressure tensor components - molecular */
  double nfree;					/* number degrees of freedom              */
} box[NBOXES];


  /* ----------------------------------------- */
  /* BoxProp                                   */
  /* This structure contains information on    */
  /* number of different molecules in each     */
  /* box.                                      */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct BoxProp {
  int first;
  int last;
  int nbox;
} bp[NBOXES][NTYPES];


  /* ----------------------------------------- */
  /* results                                   */
  /* This structure contains the variables for */
  /* the accumulation of the calculated        */
  /* simulation properties. These values are   */
  /* used to obtain statistics for these       */
  /* properties.                               */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct results {
  int    count;
  double tempa, tempb, tempc;
  double pressa, pressb, pressc;
  double pressmola, pressmolb, pressmolc;
  double cpressa[6], cpressb[6], cpressc[6];
  double cpressmola[6], cpressmolb[6], cpressmolc[6];
  double densa, densb, densc;
  double ebonda, ebondb, ebondc;
  double ebenda, ebendb, ebendc;
#ifndef NEUTRAL
  double eureya, eureyb, eureyc;
#endif
  double etorsa, etorsb, etorsc;
  double eimpra, eimprb, eimprc;
  double ecoulomba,ecoulombb,ecoulombc;
  double enbonda, enbondb, enbondc;
  double enbondsa, enbondsb, enbondsc;
  double etpota, etpotb, etpotc;
  double etpotsa, etpotsb, etpotsc;
  double etkina, etkinb, etkinc;
  double etotala, etotalb, etotalc;
  double etotalsa, etotalsb, etotalsc;
  double Ha, Hb, Hc;
  double Cva, Cvb, Cvc;
#ifdef ELASTIC
  double Cijklb[21], Cijklc[21];
#endif
} res[NBOXES], resb[NBOXES];

  /* ----------------------------------------- */
  /* RMSD reference pdb                        */
  /* ----------------------------------------- */
#ifdef RMSD
#define MAXPOINTS     400
#ifndef MAIN
extern
#endif
struct atom_rmsd {
  char   res[5];                                /* symbol for residue    */
  char   name[5];
  double x;                                     /* x-position of site */
  double y;                                     /* y-position of site */
  double z;                                     /* z-position of site */
  int nres;
  int count;
} ref[MAXPOINTS];

#ifndef MAIN
extern
#endif
int refN;

#ifndef MAIN
extern
#endif
int res_excl[MAXPOINTS];        /* This will read in the residues number of the excluded residu
es*/

#ifndef MAIN
extern
#endif
int count_excl;             /* counter to keep count of number of residues to be excluded */

#ifndef MAIN
extern
#endif
struct rmsd_modes_struct{
  int file;
  int atoms;
}rmsd_modes;
#endif //RMSD

  /* ----------------------------------------- */
  /* Seed for random number generator          */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
long idum;                 

#ifndef MAIN
extern
#endif
long idum2;  

#ifndef MAIN
extern
#endif
long iy;                 

#ifndef MAIN
extern
#endif
long iv[32];                 

#ifndef MAIN
extern
#endif
long idum_seed;
  /* ----------------------------------------- */
  /* Character string for simulation title     */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
char title[MAXT];          


  /* ----------------------------------------- */
  /* Character string for output file name     */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
char name[50];             

  /* ----------------------------------------- */
  /* Another global string for a file name     */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
char name1[50];


  /* ----------------------------------------- */
  /* NH_therm                                  */
  /* Struture containing information on the    */
  /* thermostats in the Nose Hoover Chain      */
  /* method.                                   */
  /* ----------------------------------------- */
struct NH_therm{
	double r;					/* positions of positions of thermostats				*/
	double v;					/* velocities of thermostats							*/
	double G;					/* "force" to update velocities							*/
};

  /* ----------------------------------------- */
  /* NoseHooverChain                           */
  /* Struture containing information on the    */
  /* Nose Hoover Chain method.                 */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct NoseHooverChain{
	int M;						/* number of chains in the Nose Hoover Chain method		*/
	int Nc;						/* number of timesteps for multiple time step method	*/
	int Nys;					/* number for higher order terms						*/
	double Q;					/* masses of thermostats								*/
	double* w;					/* higher order coefficients							*/
	struct NH_therm *zeta;				/* structure containing info on each thermostat			*/
} nhc[NBOXES];

  /* ----------------------------------------- */
  /* Global variables needed for Parinello     */
  /* Rhaman NPT                                */
  /* ----------------------------------------- */
#ifdef PR_NPT

#ifndef MAIN
extern
#endif
	double axes[NBOXES][9];

#ifndef MAIN
extern
#endif
	double axesi[NBOXES][9];
	
#ifndef MAIN
extern
#endif
	double va[NBOXES][9];

#ifndef MAIN
extern
#endif
	double Ga[NBOXES][9];

#ifndef MAIN
extern
#endif
	double vvab[NBOXES][9], frab[NBOXES][9];

#ifndef MAIN
extern
#endif
	double delta[9];
#ifndef MAIN
extern
#endif
	double keaxes[NBOXES];

#ifndef MAIN
extern
#endif
	double sum_ktzeta_i[NBOXES]; //pe 2-M thermostats


#ifndef MAIN
extern
#endif
	double ke_thermostat[NBOXES]; //ke thermostats

#ifndef MAIN
extern
#endif
	double cnf9ktzeta1[NBOXES]; //pe 1st thermostat

#ifndef MAIN
extern
#endif
	double cka[NBOXES]; //kinetic energy box 

#ifndef MAIN
extern
#endif
	double pr_pv[NBOXES]; //pv term 
#endif //PR_NPT


#ifndef MAIN
extern
#endif
	double Wg[NBOXES];


  /* ----------------------------------------- */
  /* imp_sasa                                  */
  /* Structure containing information used in  */
  /* the SASA implicit solvent model.          */
  /* ----------------------------------------- */
#ifdef SASA
#ifndef MAIN
extern
#endif
#ifndef BGO
struct imp_sasa {
  double r;						/* Radius of atom										*/ 
  double p;						/* atom type paramter									*/
  double S;						/* SASA of isolated atom								*/
  double sigma;					/* atomic solvation parameter							*/
  double A;						/* SASA of atom											*/
  double dAx;
  double dAy;
  double dAz;
} **sasa;
#endif
#ifdef BGO
struct bgo_sasa {
  double r;						/* Radius of residue optimized as parameters in POPS						*/ 
  double r_m;						/* Radius of residue measured with residue-residue distances at protein-ligand interfaces	*/ 
  double p;						/* Residue type paramter									*/
  double S;						/* SASA of isolated atom								*/
  double A;						/* SASA of atom											*/
  double dAx;
  double dAy;
  double dAz;
} **sasa;

struct tfe_sasa {
  double m_BB;
  double m_SC;
  double b_BB;
  double b_SC;
  double S_BB;						/* Reference SASA of residue backbone from O'Brien's MTM of residues						*/
  double S_SC;						/* Reference SASA of residue sidechain from O'Brien's MTM of residues						*/
} **tfe;
double* concentration;
#ifdef MEMBRANE
   double *concent_local;
#endif
double scale_factor;
double scale_factor_b;
int solvent_type;
double** beta;
#ifdef MEMINSIDE
  int number_in;
  int *mem_in_list;
#endif
#endif

#endif

  /* ----------------------------------------- */
  /* Hybrid MD MC needs following structure    */
  /* Structure containing information like #   */
  /* of MD steps, moves accepted etc.          */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct hybrid{
	long cyc_hybmd;				/* Number of MD steps per MC step						*/
	long cyc_swap;				/* Number of MC steps after which swap move is proposed */
	long hyb_acc[NBOXES];		/* Number of hybrid MC moves accepted					*/
	long swap_acc[NBOXES];		/* Number of swap moves accepted						*/
	int flag;
} sim_hyb;
#ifndef MAIN
extern
#endif
double PIVOT[3][2][3][3];

#ifndef MAIN
extern
#endif
struct pivot{
	long pivot_acc[NBOXES];
	double PMS[NBOXES];
	double theta[NBOXES];
} mc_pivot;
#ifndef MAIN
extern
#endif
struct trans{
	long trans_acc[NBOXES];
	double PMS[NBOXES];
	double delta[NBOXES];
} mc_trans;
#ifndef MAIN
extern
#endif
struct axial{
	long axial_acc[NBOXES];
	double PMS[NBOXES];
	double strain[NBOXES];
} mc_axial;
#ifndef MAIN
extern
#endif
struct arandoml{
	long rand_acc[NBOXES];
	double PMS[NBOXES];
	double delta[NBOXES];
} mc_rand;
#ifndef MAIN
extern
#endif
long mc_hmc_acc[NBOXES];

// This move is valid only with solvent and only if SLIST is defined.
#ifndef MAIN
extern
#endif
struct solv_rand{
	long solv_acc[NBOXES];
	double PMS[NBOXES];
	double delta[NBOXES];
} mc_solv;

#ifdef DOS
  /* ----------------------------------------- */
  /* Density of states MD MC need following:   */
  /* Structure containing information like #   */
  /* of bins,energy range and other parameters */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct dos{
	double e_begin;
	double e_end;
	double e_width;
	double T_begin;
	double T_end;
	double mod_f;
	int e_bins;
	long dos_acc;
	#ifdef STATS
	int data_count;
	#endif
	double flat;
} sim_dos[NBOXES];

#ifndef MAIN
extern
#endif
struct hist{
	double e_mid;
	double g_of_e;
	int h_of_e;
} dos_hist[NBOXES][NBINS], dos_hist_temp[NBOXES][NBINS];
#ifndef MAIN
extern
#endif
struct averages{
	double ebond;
	double ebend;
#ifndef NEUTRAL
	double eurey;
	double eureya;
#endif
	double etors;
	double eimpr;
	double elj;
	double ecoul;
	double epotens;
	double d_nc;
	double ebonda;
	double ebenda;
	double etorsa;
	double eimpra;
	double elja;
	double ecoula;
	double epotensa;
	double d_nca;
	double hel;
	double hela;
	double con;
	double cona;
	double con_2;
	double con_2a;
	int count;
#ifdef SASA
	double esasa;
	double esasaa;
	double esasab;
#endif
	double gyr;
	double gyra;
	double rmsd;
	double rmsda;
}average[NBOXES][NBINS];

#endif // DOS

#ifdef EWALD

#ifndef SPME
struct Complex
{
   double r;  /*real part*/
   double i;  /*imag part*/
};

#ifndef MAIN
extern
#endif
	double **kvec;

#ifndef MAIN
extern
#endif
	int *total_k_vectors;

#ifndef MAIN
extern
#endif
	struct Complex **e_ikx, **e_iky, **e_ikz;

#ifndef MAIN
extern
#endif
	struct Complex *sum;

#endif // endif not SPME

#ifdef SPME

#ifndef MAIN
extern
#endif
struct pmedat {
  double x,y,z;
} **mm, **dm, *ss;

#ifndef MAIN
extern
#endif
fftw_complex *qgrida;

#ifndef MAIN
extern
#endif
fftw_plan qgrid_back, qgrid_forward;

#ifndef MAIN
extern
#endif
double *bsp_modx,*bsp_mody,*bsp_modz,*bsp_arr;

#ifndef MAIN
extern
#endif
double *mn, *dmn;

#ifndef MAIN
extern
#endif
int gridmx;

#endif   // endif SPME

#endif //endif EWALD

#ifdef MPI
  /* ----------------------------------------- */
  /* This is for running the code on parallel  */
  /* processors using MPI					   */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct mpi_def{
	int my_rank;
	int p;
	int source;
	int dest;
	int tag ;
	int flag;	
	char message[100];
	MPI_Status status;
} mpi;

struct msg_swap_request{
	double potens;
#ifdef MMDOS
	double totals;
#endif
	double e_begin;
	double e_end;
	double mod_f;
#ifdef XEDOS
	double l_begin;
	double l_end;
	double d_nc;
#endif
  double kT;
 
};

struct msg_swap_accept{
	double potens;
#ifdef MMDOS
	double totals;
#endif
#ifdef XEDOS
	double d_nc;
#endif
  double scale;
  double kT; 
	double accept;
};	
#endif

#ifdef STYPE
  /* ----------------------------------------- */
  /* This structure contains information about */
  /* the hash table where excluded interaction */
  /* are stored.                               */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct hashlist{
	int *a, *b;		/* Site numbers involved in the exluded interction           */
	int n;			/* Number of occupants at a specific index of the hash table */
	int *p;			/* Identifies which i of xint[k][i] is at this position      */
} **table;
#ifndef MAIN
extern
#endif
int hashsize;
#endif//STYPE

#ifdef HESSIAN
#ifndef MAIN
extern
#endif
double **hessian;
#endif

#ifdef MMDOS
  /* ----------------------------------------- */
  /* Density of states MD MC need following:   */
  /* Structure containing information like #   */
  /* of bins,energy range and other parameters */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct dos{
	double e_begin;
	double e_end;
	double e_width;
	double T_begin;
	double T_end;
	double mod_f;
	double p_dke;			// probability of proposing a delta KE move.
	double d_max;			// maximum delta KE displacement
	int e_bins;
	long dos_acc;
	long dke_acc;
	double flat;
	#ifdef STATS
	int data_count;
	#endif
} sim_dos[NBOXES];

#ifndef MAIN
extern
#endif
struct hist{
	double e_mid;
	double g_of_e;
	double k_of_e;			// kinetice energy for each Total E bin
	int h_of_e;
} dos_hist[NBOXES][NBINS], dos_hist_temp[NBOXES][NBINS];
#ifndef MAIN
extern
#endif
struct averages{
	double ebond;
	double ebend;
#ifndef NEUTRAL
	double eurey;
	double eureya;
#endif
	double etors;
	double eimpr;
	double elj;
	double ecoul;
	double epotens;
	double d_nc;
	double ebonda;
	double ebenda;
	double etorsa;
	double eimpra;
	double elja;
	double ecoula;
	double epotensa;
	double d_nca;
	double hel;
	double hela;
	double con;
	double cona;
	int count;
#ifdef SASA
	double esasa;
	double esasaa;
	double esasab;
#endif
	double gyr;
	double gyra;
	double rmsd;
	double rmsda;
}average[NBOXES][NBINS];
#endif // MMDOS
#ifdef NMA
#ifndef MAIN
extern
#endif
double *pos_nma;
#ifndef MAIN
extern
#endif
double *enr_nma;
#ifndef MAIN
extern
#endif
double *for_nma;
#endif
/* The following for animation counter	*/
#ifndef MAIN
extern
#endif
struct counter{
	int count;
} frames[NBOXES];
#ifdef CONFIGT
#ifndef MAIN
extern
#endif
struct hessn{
	double hesx;
	double hesy;
	double hesz;
	double hesr;
	double num;
	double den;
	double T;
}config[NBOXES],config_temp[NBOXES];
#endif

#ifdef CTDOS
  /* ----------------------------------------- */
  /* Density of states MD MC need following:   */
  /* Structure containing information like #   */
  /* of bins,energy range and other parameters */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct dos{
	double e_begin;
	double e_end;
	double e_width;
	double T_begin;
	double T_end;
	double mod_f;
	int e_bins;
	long dos_acc;
	double flat;
	#ifdef STATS
	int data_count;
	#endif
} sim_dos[NBOXES];

#ifndef MAIN
extern
#endif
struct hist{
	double e_mid;
	double g_of_e;
	double gT_of_e;
	int h_of_e;	
	double ct_den;
	double ct_num;
	double config_T;
} dos_hist[NBOXES][NBINS], dos_hist_temp[NBOXES][NBINS];
#ifndef MAIN
extern
#endif
struct averages{
	double ebond;
	double ebend;
#ifndef NEUTRAL
	double eurey;
	double eureya;
#endif
	double etors;
	double eimpr;
	double elj;
	double ecoul;
	double epotens;
	double d_nc;
	double ebonda;
	double ebenda;
	double etorsa;
	double eimpra;
	double elja;
	double ecoula;
	double epotensa;
	double d_nca;
	double hel;
	double hela;
	double con;
	double cona;
	double con_2;
	double con_2a;
	int count;
#ifdef SASA
	double esasa;
	double esasaa;
	double esasab;
#endif
	double gyr;
	double gyra;
	double rmsd;
	double rmsda;
}average[NBOXES][NBINS];

#endif // CTDOS

#ifdef ELASTIC
#ifndef MAIN
extern
#endif
struct elasticity_tensor{
	double kinetic_term[21];
	double born_term[21];
	double Pij[6];
	double born_term_av[21];
	double Pij_av[6];
	double Pijkl_av[6];
}elast[NBOXES];
#endif //ELASTIC

#ifdef KONS
#ifndef MAIN
extern
#endif
struct cons_struct{
	int n;	   //number of sites
	int *site; //array of sites
}*cons;		   //array of boxes
#endif//KONS



#ifdef REST
#ifndef MAIN
extern
#endif
struct rest_struct_1{  //this structure is used if the 
	int n;             //atom is restrained to its 
	int *site;         //original site.
	double *k;
	double *x0;
	double *y0;
	double *z0;
	double *req;
} *rest1;

#ifndef MAIN
extern
#endif
struct rest_struct_2{  //this structure is used if
	int n;             //the distance between two 
	int *site1;		   //atoms is restrained to 
	int *site2;        //a fixed distrance (req)
	double *k;
	double *req;
} *rest2;

#endif //

#ifdef XEDOS
  /* ----------------------------------------- */
  /* Density of states MD MC need following:   */
  /* Structure containing information like #   */
  /* of bins,energy range and other parameters */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct dos{
	double l_begin;
	double l_end;
	double l_width;
	double T_begin;
	double T_end;
	double mod_f;
	int l_bins;
	long dos_acc;
	double flat;
#ifdef STATS
	int data_count;
#endif
} sim_dos[NBOXES];

#ifndef MAIN
extern
#endif
struct hist{
	double l_mid;
	double g_of_l;
	int h_of_l;	
} dos_hist[NBOXES][NBINS], dos_hist_temp[NBOXES][NBINS];
#ifndef MAIN
extern
#endif
struct averages{
	double ebond;
	double ebend;
#ifndef NEUTRAL
	double eurey;
	double eureya;
#endif
	double etors;
	double eimpr;
	double elj;
	double ecoul;
	double epotens;
	double d_nc;
	double ebonda;
	double ebenda;
	double etorsa;
	double eimpra;
	double elja;
	double ecoula;
	double epotensa;
	double d_nca;
	double hel;
	double hela;
	double con;
	double cona;
    double con_2;
	double con_2a;
	double force_1;
	double force_2;
	double force_1a;
	double force_2a;
	int count;
#ifdef SASA
	double esasa;
	double esasaa;
	double esasab;
#endif
	double gyr;
	double gyra;
	double rmsd;
	double rmsda;
}average[NBOXES][NBINS];

#endif // XEDOS

#ifdef SMD
#ifndef MAIN
extern
#endif
struct smd_struct{
	int site1;			
	int site2;
	int count;
	double k1;
	double k2;
	double x1;
	double y1;
	double z1;
	double x2;
	double y2;
	double z2;
	double v; 
	double r0;
	double r1;
	double r2;
	double ex;
	double ey;
	double ez;
	double r1x;
	double r1y;
	double r1z;
	double r2x;
	double r2y;
	double r2z;
	double dw2;
	double dw1;
	double effr1;
	double effr2;
	double effr1_av;
	double effr2_av;
} *steermd;
#endif //SMD
#ifdef NSMD
#ifndef MAIN
extern
#endif
struct smd_struct{
	int site1;			
	int site2;
	int count;
	double k1;
	double k2;
	double x1;
	double y1;
	double z1;
	double x2;
	double y2;
	double z2;
	double v; 
	double r0;
	double ex;
	double ey;
	double ez;
	double r1x;
	double r1y;
	double r1z;
	double effr1;
	double effr1_av;
	double f2_av;
} *steermd;
#endif //NSMD

#ifndef MAIN
extern 
#endif
	int SITE1;

#ifndef MAIN
extern 
#endif
	int SITE2;

#ifndef MAIN
extern 
#endif
	int SITE3;

#ifdef MC

#ifndef MAIN
extern
#endif
double **ffox, **ffoy, **ffoz, **ffnx, **ffny, **ffnz;

#ifdef CONFIGT
#ifndef MAIN
extern
#endif
double hesox, hesoy, hesoz, hesor;
#endif

#ifndef MAIN
extern
#endif
struct energy_mc {
  double o_poten;	double n_poten;
  double o_potens;	double n_potens;
  double o_nbond;	double n_nbond;
  double o_nbonds;	double n_nbonds;
  double o_bond;	double n_bond;
  double o_bend;	double n_bend;
  double o_tors;	double n_tors;
  double o_impr;	double n_impr;
  double o_coulomb; double n_coulomb;
#ifdef SASA
  double o_esasa;	double n_esasa;
#endif
#ifdef EWALD
  double o_rewald;  double n_rewald; 
  double o_rewalds; double n_rewalds;
  double o_kewald;	double n_kewald;
  double o_sewald;	double n_sewald;
  double o_iewald;	double n_iewald; 
  double o_tewald;	double n_tewald;
#endif
#ifndef NEUTRAL
  double o_urey;	double n_urey;
#endif
} enmc[NBOXES], enmc_temp[NBOXES];

#endif
#ifndef MAIN
extern 
#endif
	int fileindex;

#ifdef FX_EDOS

  /* ----------------------------------------- */
  /* Density of states MD MC need following:   */
  /* Structure containing information like #   */
  /* of bins,energy range and other parameters */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct dos{
	double e_begin;
	double e_end;
	double e_width;
	double T_begin;
	double T_end;
	double r_value; double r_target; double tol_w; 
	int trips_num;  int trips_target; int trips_old;
	int e_bins;
	long dos_acc;
	int sign_walker;
	#ifdef STATS
	int data_count;
	#endif
} sim_dos[NBOXES];

#ifndef MAIN
extern
#endif
struct hist{
	double e_mid;
	double g_of_e;
	int nw_of_e;				// The three histograms
	int np_of_e;
	int nm_of_e;				
	double f_of_e;				// Fraction of n+/nw
	double w_of_e;				// weights
	double m_of_e;				// slope = df/de

} dos_hist[NBOXES][NBINS], dos_hist_temp[NBOXES][NBINS], dos_hist_old[NBOXES][NBINS];
#ifndef MAIN
extern
#endif
struct averages{
	double ebond;
	double ebend;
	double eurey;
	double eureya;
	double etors;
	double eimpr;
	double elj;
	double ecoul;
	double epotens;
	double d_nc;
	double ebonda;
	double ebenda;
	double etorsa;
	double eimpra;
	double elja;
	double ecoula;
	double epotensa;
	double d_nca;
	double hel;
	double hela;
	double con;
	double cona;
	double con_2;
	double con_2a;
	int count;
#ifdef SASA
	double esasa;
	double esasaa;
	double esasab;
#endif
	double gyr;
	double gyra;
	double rmsd;
	double rmsda;
}average[NBOXES][NBINS];
#endif//FX_EDOS
#ifndef MAIN
extern
#endif
unsigned long n_iter;				// This is a global variable to keep track of number of iterations

#ifdef BGO
struct bgo_struct{
double eps;
double req;
double d_native;
double d_repulsive;
double *hydro_index;
#ifdef SASA
double *vdw_radii;    //van der waals radii in diameter, added by Shuai Wei
double *vdw_radii_measured;    //van der waals radii in diameter measured through protein-ligand binding interfaces by Shanshan Cheng and Charlie Brooks, added by Shuai Wei
double *p_sasa;       //POPS parameters for calculating SASA for BGO model added by Shuai Wei
double *ref_sasa;       //POPS reference surface area for calculating SASA for BGO model added by Shuai Wei
double *tfe;          //transfer free energy in specific solvent, added by Shuai Wei
#endif
double mass;
int Psite;                            // Site index in crd where      solute starts for energy partition
}*bgo;
#endif

#ifdef GOLIK
#ifndef MAIN
extern
#endif
struct golik_struct{
  double eps;
  double req;
  double d_native;
  double d_repulsive;
#ifndef DNA_GOLIK
  int** contact;
#endif
  double *hydro_index;
  double mass;
  int Psite;				// Site index in crd where 	solute starts for energy partition
}*golik;
#ifndef MAIN
extern
#endif
struct gosol_struct{
  int N;					// Number of colute type
  
  double *eps;				// eps for each cosolute types
  double *sig;				// sigma for each cosolute types
  double lambda;			// Scaling parameter for protein-cosolute interaction
  double *mass;
  int *Nsol;
}gosol;
#endif

/*#ifdef DNA_GOLIK
#ifndef MAIN
extern
#endif
int goncN;

#ifndef MAIN
extern
#endif
struct gonc_struct{
  int a;
  int b;
  int sig;
} **gonc;
#endif
*/
#ifdef WALL
#ifndef MAIN
extern
#endif
double npr_1;
double npr_2;
double npr_w;
double npr_p;
double npr_3;
struct wall_struct{
  int n;
  double angle;
  int angle_site;
  double angle_k;
  double *z;
  int *nsites;
  double **num_density;
  double **eps;
  double **sig;
  double **hydro_index;
}*wall;
#ifdef SASA
// values will be assigned in vdw_radii_bgo.c
double   mBBbw; 
double   mSCbw; 
double   bBBbw; 
double   bSCbw; 
double   sBBbw;
double   sSCbw;
double   pbw; // The parameter p for a surface 
double   sigW; // Transfer free energy on surface by a residue
#endif
struct sphere{
  double x; 
  double y; 
  double z;
  double r; //In the unit of A
};
#endif


#ifdef REM
#ifndef MAIN
extern
#endif
struct replica{//This structure stores information about replica in each box
  int n;
  double sig;
  double E_m;
  double rho;
  int tag;		// This is the replica number that is running in box k
} rep_exc[NBOXES];

#ifndef MAIN
extern
#endif
struct junction{//This structure stores information about each junction
  int n;
  double sig;
  double E_m;
  double rho;
  double flux;
  int count;
}rep_jnc[NBOXES];
#ifndef MAIN
extern
#endif
struct replica_opt{
  double Jm;
  double Je;
}rep_opt[NBOXES];
#endif//REM


#ifdef TDXEDOS
  /* ----------------------------------------- */
  /* Density of states MD MC need following:   */
  /* Structure containing information like #   */
  /* of bins,energy range and other parameters */
  /* ----------------------------------------- */
#ifndef MAIN
extern
#endif
struct dos{
	double x1_begin;   //The begining and ending points of
	double x1_end;       //the first reaction coordinate.
	double x2_begin;   //The begining and ending points of
	double x2_end;       //the second reaction coordinate.
	double x1_width;
	double x2_width;
	double T_begin;
	double T_end;
	double mod_f;
	int x1_bins;
	int x2_bins;
	long dos_acc;
	double flat;
#ifdef STATS
	int data_count;
#endif
} sim_dos[NBOXES];

#ifndef MAIN
extern
#endif
struct hist{
	double x1_mid;
	double x2_mid;
	double g_of_l;
	int    h_of_l;	
} dos_hist[NBOXES][NBINS1][NBINS2], dos_hist_temp[NBOXES][NBINS1][NBINS2];

#ifndef MAIN
extern
#endif
struct averages{
	double ebond;
	double ebend;
#ifndef NEUTRAL
	double eurey;
	double eureya;
#endif
	double etors;
	double eimpr;
	double elj;
	double ecoul;
	double epotens;
	double d_nc;
	double ebonda;
	double ebenda;
	double etorsa;
	double eimpra;
	double elja;
	double ecoula;
	double epotensa;
	double d_nca;
	double hel;
	double hela;
	double con;
	double cona;
	double con_2;
	double con_2a;
	double force_1;
	double force_2;
	double force_1a;
	double force_2a;
	int count;
#ifdef SASA
	double esasa;
	double esasaa;
	double esasab;
#endif
	double gyr;
	double gyra;
	double rmsd;
	double rmsda;
	double x1;
	double x2;
	double x1a;
	double x2a;
}average[NBOXES][NBINS1][NBINS2];

#ifndef MAIN
extern
#endif
struct swap_neighbors_struct{
  int n_boxes;
  int *n_nabors;
  int **nabor;
  double *pick_box_prob;
} swapping;
#endif //TDXEDOS

/* ========================================== */
/* Buffers for I/O.  There is one for each    */
/* file that is appended.  Those which are    */
/* overwritten do not need buffers.           */
/* ========================================== */
#ifndef MAIN
extern
#endif
char **simul_buff;

#ifndef MAIN
extern
#endif
char **simul_ptr;

#ifndef MAIN
extern
#endif
time_t *simul_time;

#ifndef MAIN
extern
#endif
char **ener_buff;

#ifndef MAIN
extern
#endif
char **ener_ptr;

#ifndef MAIN
extern
#endif
time_t *ener_time;

#ifndef MAIN
extern
#endif
char **ord_buff;

#ifndef MAIN
extern
#endif
char **ord_ptr;

#ifndef MAIN
extern
#endif
time_t *ord_time;

#ifndef MAIN
extern
#endif
char **swap_buff;

#ifndef MAIN
extern
#endif
char **swap_ptr;

#ifndef MAIN
extern
#endif
time_t *swap_time;

#if defined(SMD) || defined(NSMD)
	#ifndef MAIN
	extern
	#endif
	char **smd_buff;

	#ifndef MAIN
	extern
	#endif
	char **smd_ptr;

	#ifndef MAIN
	extern
	#endif
	time_t *smd_time;
#endif

#if defined(TWHAM) || defined(XWHAM)
  #ifndef MAIN
  extern
  #endif
  time_t *wham_time;
#endif

/* ========================================== */
/* variable needed for wham.                  */
/* ========================================== */
#ifdef TWHAM
#ifndef MAIN
extern
#endif
double *pe_max;

#ifndef MAIN
extern
#endif
double *pe_min;
#endif //TWHAM

#ifdef CHECKP
//Checkpoint flags (don't checkpoint these, haha)

#ifndef MAIN
extern
#endif
int restart_equil;
#ifndef MAIN
extern
#endif
int restart_prod;
#ifndef MAIN
extern
#endif
int terminating;
#ifndef MAIN
extern
#endif
char equiname[2][100];
#ifndef MAIN
extern
#endif
char prodname[2][100];
#ifndef MAIN
extern
#endif
int checkptr;	// 0 or 1 based on which .chk# file is more recent
#endif


/* ========================================== */
/* Defines for stresses and invariants        */
/* ========================================== */
#ifdef ZHOU
        #ifndef MAIN
                extern
        #endif
        struct stress_type {

                double  I1,  I2,  I3;           //Stress invariants
                double  P;                      //Mean stress

                int    cpxI, cpxJ;             //True if eigenvalues are complex.
                double  Eval_I[3];             //Real part of eigenvalues of stress tensor
                double  Evec[3][3];             //Eigenvectors [vector][component]
                double  J1,  J2,  J3;           //Deviatoric stress invariants
                double  Eval_J[3];             //Eigenvalues of stress deviator tensor
                double  VMS;                    //von Mises stress

        } **stress_i, **stress_isum, **stress_ibarsum, **stress_ibarsumsq;    //Invariants

        #ifndef MAIN
                extern  
        #endif
        double ****stresses, ****stresses_sum;          //Stress tensors

        #ifndef MAIN
                extern  
        #endif
        int stressStepCount;

        #ifndef MAIN
                extern
        #endif
        int stressBlockCount;

        #ifndef MAIN
                extern
        #endif
        int stressOutputCount;

        #ifndef MAIN
                extern
        #endif
        int stressTotalCount;

        #ifdef TRR
                #ifndef MAIN
                        extern
                #endif
                int stressFrameCount;

                #ifndef MAIN
                        extern
                #endif
                double *stress_movie_psum;              //Allocate for box 0 only for now
        #endif


#endif


