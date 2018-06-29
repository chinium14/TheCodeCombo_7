/* ======================================================================== */
/* save_config.cpp                                                          */
/*                                                                          */
/*		This subroutine saves the molecular configuration. It is called     */
/* from main() with a frequency set by sim.blockc.  This frequency is       */
/* set in the simul.input file.  The file output is the same format as the  */
/* simul#.crd file so a simulation can be started from a previous           */
/* simulation.  The file saves the configurations in all of the boxes to    */
/* different files.  The output file is named simul#.crd.sav. To use it to  */
/* start a simulation, simply delete .sav from the name.                    */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						block:		The iteration number                    */
/*                      cycle:		The number of iterations in the         */
/*                                  simulation.                             */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"
void write_gro(unsigned long, unsigned long);

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void save_config (unsigned long block, unsigned long cycle) 
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
//  char  amino_id[21]={'X','G','A','S','C','V','T','I','P','M','D','N','L','K','E','Q','R','H','F','Y','W'};
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
	  sprintf(name,"./OUTPUT/BOX%d/simul%d.crd.sav",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	  sprintf(name,"./OUTPUT/BOX%d/simul%d.crd.sav",k,k);
#endif
    fas = fopen(name,"w");
  /* ----------------------------------------------- */
  /* Write heading and simulation parameters.        */
  /* ----------------------------------------------- */    
/*    fprintf(fas,"Last configuration: %d blocks %d of cycles @ %.24s\n",
	    block,cycle,ctime(&now));
    fprintf(fas,"\n");
    for(int i=0; i<sim.NC; i++) {
      fprintf(fas,"%d       Number of molecules %c in box %d\n",
	      bp[k][i].nbox,da[i],k);
    }
    fprintf(fas,"\n");
    for(int i=0; i<sim.NC; i++) {
      fprintf(fas,"%d       Number of atoms/sites in molecule %c\n",
	      mol[i].Nsite,da[i]);
    }
	  fprintf(fas,"\n");
	for (int i=0; i<sim.NC; i++) {
      fprintf(fas,"%d       Number of residues in molecule %c\n",
		  mol[i].Nres,da[i]); 
    }
    fprintf(fas,"\n");
    fprintf(fas,"%2.3f     Length of box %d [Angstroms]\n",box[k].boxl,k);
    fprintf(fas,"\n");
    fprintf(fas,"xyz coordinates of sites [Angstroms]\n");
    fprintf(fas,"\n");
*/
  /* ----------------------------------------------- */
  /* Additional variables needed.                    */
  /* ----------------------------------------------- */    
	int rescount =0;double pdb =1.00;
    int count = 0;
 //   int nc = 0;
  /* ----------------------------------------------- */
  /* The following loops around the sites in each    */
  /* residue and prints the x,y,z coordinates.       */
  /* ----------------------------------------------- */    
	fprintf(fas,"%d \n",box[k].boxns);
	for(int m=0; m<sim.NC; m++) {
//		fprintf(fas,"%d \n",mol[m].Nsite);
		for(int i=0; i<bp[k][m].nbox; i++) {
			for (int kk=0; kk<mol[m].Nres; kk++) {			  
				for(int j=0; j<residue[k][rescount].Nsite; j++) {
					fprintf(fas,"%5d  ", count+1);
					fprintf(fas,"%5d  ", rescount+1);
					fprintf(fas,"%s ",amino_name[residue[k][rescount].type]);
					
       /* *********************************************** */
       /* These if-else statements print the atom name    */
	   /* at the correct spacing which changes if the     */
	   /* name has one, two, or three characters in it.   */
       /* *********************************************** */
					if(	atom[k][count].name[1]==0) fprintf(fas,"%s   ", atom[k][count].name);
					else if(atom[k][count].name[2]==0) fprintf(fas,"%s  ", atom[k][count].name);
					else if(atom[k][count].name[3]==0) fprintf(fas,"%s ",  atom[k][count].name);
					else fprintf(fas,"%s", atom[k][count].name);
	  
					fprintf(fas,"%10.5lf ", atnopbc[k][count].x);
					fprintf(fas,"%7.5lf ",  atnopbc[k][count].y);
					fprintf(fas,"%7.5lf  ", atnopbc[k][count].z);
					fprintf(fas,"%4.2lf  ", pdb);
					fprintf(fas,"%4.2lf\n", pdb-1);
	  	  
					count++;
				}//for j
				rescount++;
			}//for kk
		}//for i
	}//for m

	fclose(fas);
  }//for k
  write_gro(block, cycle);
}
