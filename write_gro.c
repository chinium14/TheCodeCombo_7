/* ======================================================================== */
/* write_gro.cpp	                                                        */
/*      Nitin Rathore                                                       */
/*		This subroutine saves the molecular configuration.in Gromacs format */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void write_gro (unsigned long block, unsigned long cycle) 
{
  /* ================================================================== */
  /*                                                                    */
  /* Variables and pointers needed in the subroutine.                   */
  /*                                                                    */
  /* ================================================================== */
  char name[50];
char amino_name[47][4]={"ABC","GLY","ALA","SER","CYS","VAL","THR","ILE","PRO","MET","ASP","ASN","LEU","LYS","GLU",
		"GLN","ARG","HIS","PHE","TYR","TRP","CYT","GUA","ADE","THY","HOH","SOD","CLA","LIP","CHO","EAM","POT","TRE", "GOL","F00", "F01",
		"F02","F03","F04","F05","F06","F07","F08","F09","F10","F11","MEO"};
  //char da[6] = "ABCDE";
  time_t now;
  FILE *fas;
  
  for(int k=0; k<sim.NB; k++) {
    time(&now);
  /* ================================================================== */
  /*                                                                    */
  /* Open the file and write the output.                                */
  /*                                                                    */
  /* ================================================================== */
#ifdef MPI
	  sprintf(name,"./OUTPUT/GROMACS/simul%d.gro",mpi.my_rank);
#endif
#ifndef MPI
	  sprintf(name,"./OUTPUT/GROMACS/simul%d.gro",k);
#endif
    fas = fopen(name,"w");
  /* ----------------------------------------------- */
  /* Write heading and simulation parameters.        */
  /* ----------------------------------------------- */    
    fprintf(fas,"Gro file created by user: %lu blocks of %lu cycles @ %.24s\n",
	    block,cycle,ctime(&now));
  /* ----------------------------------------------- */
  /* Additional variables needed.                    */
  /* ----------------------------------------------- */    
	int rescount =0;//double pdb =1.00;
    int count = 0;
 //   int nc = 0;
  /* ----------------------------------------------- */
  /* The following loops around the sites in each    */
  /* residue and prints the x,y,z coordinates.       */
  /* ----------------------------------------------- */    
	fprintf(fas,"%d \n",box[k].boxns);
	for(int m=0; m<sim.NC; m++) {
		for(int i=0; i<bp[k][m].nbox; i++) {
			for (int kk=0; kk<mol[m].Nres; kk++) {			  
				for(int j=0; j<residue[k][rescount].Nsite; j++) {
					fprintf(fas,"%5d", rescount+1);
					fprintf(fas,"%5s",amino_name[residue[k][rescount].type]);
					fprintf(fas,"%5s", atom[k][count].name);
					fprintf(fas,"%5d", count+1);
					fprintf(fas,"%8.3lf", atom[k][count].x/10.0);
					fprintf(fas,"%8.3lf", atom[k][count].y/10.0);
					fprintf(fas,"%8.3lf", atom[k][count].z/10.0);
					fprintf(fas,"%8.4lf",  vv[k][count].x/10.0);
					fprintf(fas,"%8.4lf",  vv[k][count].y/10.0);
					fprintf(fas,"%8.4lf\n",vv[k][count].z/10.0);
					count++;
				}//for j
				rescount++;
			}//for kk
		}//for i
	}//for m
	fprintf(fas,"%10.5lf%10.5lf%10.5lf\n",box[k].lx/10.0,box[k].ly/10.0,box[k].lz/10.0);
	fclose(fas);
  }//for k
}
