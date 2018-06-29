
#include "defines.h"
#ifdef STATUS
void save_config(unsigned long, unsigned long);
void pbc_chk(int);
void forces(int);
void calcvalue(int);
char namestat2[40];
void curr_status (int k, int flag)
{
 
  FILE *stat2;

#ifndef MPI
   sprintf(namestat2, "./INPUT/status%d.txt",k);
#endif
#ifdef MPI
   sprintf(namestat2, "./INPUT/status%d.txt",mpi.my_rank);
#endif

  if( NULL != (stat2=fopen(namestat2,"r"))){
	  
	  char namestat[40];
#ifndef MPI
	  sprintf(namestat, "./OUTPUT/BOX%d/status%d.log",k,k);
#endif
#ifdef MPI
	  sprintf(namestat, "./OUTPUT/BOX%d/status%d.log",mpi.my_rank,mpi.my_rank);
#endif
	  FILE *stat1;
	  stat1= fopen(namestat,"a");

	  fprintf(stat1,"curr_status function called from box = %d with flag = %d: Current stored energies are : \n",k, flag);
	  #ifndef SASA
		fprintf(stat1, "bond(%lf)	bend(%lf)	tors(%lf)	impr(%lf)	nbonds(%lf)	coul(%lf)	potens(%lf)	kin(%lf)	totals(%lf)	 \n", 
		en[k].bond,en[k].bend,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals);
	  #endif
	  #ifdef SASA
		fprintf(stat1,"bond(%lf)	bend(%lf)	tors(%lf)	impr(%lf)	nbonds(%lf)	coul(%lf)	potens(%lf)	kin(%lf)	sasa(%lf)	totals(%lf)	 \n",
		en[k].bond,en[k].bend,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].esasa,en[k].totals);
	  #endif// SASA
	  
	  save_config(0,0);
	  fprintf(stat1,"Current configurations saved in .crd.sav file and now calling pbc_check\n");
	  
	  pbc_chk(0);// iteration =0
	  fprintf(stat1,"pbc_chk function performed and now calling forces (%d)\n",k);
	  
	  forces(k);
	  calcvalue(k);
	  fprintf(stat1,"forces (%d} called and updated energies are as follows \n",k);
	  #ifndef SASA
		fprintf(stat1, "bond(%lf)	bend(%lf)	tors(%lf)	impr(%lf)	nbonds(%lf)	coul(%lf)	potens(%lf)	kin(%lf)	totals(%lf)	 \n", 
		en[k].bond,en[k].bend,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals);
	  #endif
	  #ifdef SASA
		fprintf(stat1,"bond(%lf)	bend(%lf)	tors(%lf)	impr(%lf)	nbonds(%lf)	coul(%lf)	potens(%lf)	kin(%lf)	sasa(%lf)	totals(%lf)	 \n",
		en[k].bond,en[k].bend,en[k].tors,en[k].impr,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].esasa,en[k].totals);
	  #endif// SASA
	  fprintf(stat1,"\n \n");
      fclose(stat1);
	  fclose(stat2);
  }
}  

#endif

