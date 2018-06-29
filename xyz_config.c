/* ======================================================================== */
/* xyz_config.cpp                                                           */
/*                                                                          */
/*		This subroutine writes the atom/site coordinates in the format      */
/* readable by RASMOL and the pdb format.  This subroutine is called from   */
/* main() with frequency set by sim.blockc. There are two parts to the      */
/* subroutine.  The first writes the coordinates in the correct format to   */
/* ./OUTPUT/simul.crd#.xyz,./OUTPUT/simul.crd#.pbc.xyz, and                 */
/* ./OUTPUT/simulfinal%d.pdb. For these files, the new data overwrites the  */
/* old data.  If MOVIE is defined, the second part of the subroutine is     */
/* executed where frames for an animation are stored in the folder          */
/* ./MOVIE.                                                                 */
/* Passed Parameters:                                                       */
/*						flag:	if 0 save as simulinitial.pdb else pdb		*/
/*				                files saved on the fly as simulfinal.pdb	*/
/* ======================================================================== */

#include "defines.h"
void helix(int);
/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void xyz_config (int flag, unsigned long iter)
{
  /* ================================================================== */
  /*                                                                    */
  /* Variables and pointer needed in the subroutine.                    */
  /*                                                                    */
  /* ================================================================== */
  char name[80]; char name1[80]; //char da[6] = "ABCDE";
char amino_name[47][4]={"ABC","GLY","ALA","SER","CYS","VAL","THR","ILE","PRO","MET","ASP","ASN","LEU","LYS","GLU",
		"GLN","ARG","HIS","PHE","TYR","TRP","CYT","GUA","ADE","THY","HOH","SOD","CLA","LIP","CHO","EAM","POT","TRE", "GOL","F00", "F01",
		"F02","F03","F04","F05","F06","F07","F08","F09","F10","F11","MEO"};
//	char  amino_id[21]={'X','G','A','S','C','V','T','I','P','M','D','N','L','K','E','Q','R','H','F','Y','W'};

  FILE *io; FILE *fo; 
  /* ================================================================== */
  /*                                                                    */
  /* Write the coordinates to the ./OUTPUT/simul.crd#.xyz file.  These  */
  /* first coordinates is for the simulation with no pbc.               */
  /*                                                                    */
  /* ================================================================== */
  for(int k=0; k<sim.NB; k++) {
#ifdef MPI
	  sprintf(name,"./OUTPUT/BOX%d/simul.crd%d.xyz",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	  sprintf(name,"./OUTPUT/BOX%d/simul.crd%d.xyz",k,k);
#endif  
    io = fopen(name,"w");
    fprintf(io,"   %d\n",box[k].boxns);
    fprintf(io,"\n");
    for(int i=0; i<box[k].boxns; i++) {
	  fprintf(io," %2c   ",  atnopbc[k][i].name[0]);
      fprintf(io," %15.6f  ",atnopbc[k][i].x);
      fprintf(io," %15.6f  ",atnopbc[k][i].y);
      fprintf(io," %15.6f\n",atnopbc[k][i].z);
    }
    fclose(io);
  }
  /* ================================================================== */
  /*                                                                    */
  /* Write the coordinates to the ./OUTPUT/simul.crd#.pbc.xyz file.     */
  /* These second set of coordinates is for the simulation with pbc.    */
  /*                                                                    */
  /* ================================================================== */
  for(int k=0; k<sim.NB; k++) {
#ifdef MPI
	    sprintf(name,"./OUTPUT/BOX%d/simul.crd%d.pbc.xyz",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	  sprintf(name,"./OUTPUT/BOX%d/simul.crd%d.pbc.xyz",k,k);
#endif
    io = fopen(name,"w");
    fprintf(io,"   %d\n",box[k].boxns);
    fprintf(io,"\n");
    for(int i=0; i<box[k].boxns; i++) {
      fprintf(io," %2s   ",  atom[k][i].name);
      fprintf(io," %15.6f  ",atom[k][i].x);
      fprintf(io," %15.6f  ",atom[k][i].y);
      fprintf(io," %15.6f\n",atom[k][i].z);
    }
    fclose(io);
  }
  /* ================================================================== */
  /*                                                                    */
  /* Write the coordinates to the ./OUTPUT/simulinitial.pdb file.	    */
  /* or as simulfinal.pdb depending on whether flag is 0 or 1		    */
  /*                                                                    */
  /* ================================================================== */
for(int k=0; k<sim.NB; k++) {
    if (flag==0){ 
#ifdef MPI
	  sprintf(name,"./OUTPUT/BOX%d/simulinitial%d.pdb",mpi.my_rank,mpi.my_rank);
	  sprintf(name1,"./OUTPUT/BOX%d/phipsi_init%d.out",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	  sprintf(name,"./OUTPUT/BOX%d/simulinitial%d.pdb",k,k);
	  sprintf(name1,"./OUTPUT/BOX%d/phipsi_init%d.out",k,k);
#endif
	}
		
	else if (flag==1){
#ifdef MPI
	  sprintf(name, "./OUTPUT/BOX%d/simulfinal%d.pdb",mpi.my_rank,mpi.my_rank);
	  sprintf(name1,"./OUTPUT/BOX%d/phipsi_fin%d.out",mpi.my_rank,mpi.my_rank);
#endif
#ifndef MPI
	  sprintf(name, "./OUTPUT/BOX%d/simulfinal%d.pdb",k,k);
	  sprintf(name1,"./OUTPUT/BOX%d/phipsi_fin%d.out",k,k);
#endif
	}
    io = fopen(name,"w");
	fo = fopen(name1,"w");
#ifndef SASA
	fprintf(io,"REMARK\tITER		NBONDS	COUL		POTNS		KINE		TOT\n");
	fprintf(io, "REMARK\t%lu	%.3lf	%.3lf	%.3lf	%.3lf	%.3lf\n", iter,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].totals);
#endif
#ifdef SASA
	fprintf(io,"REMARK\tITER		NBONDS	COUL		POTNS		KINE	SOLV	TOT\n");
	fprintf(io,"REMARK\t%lu	%.3lf	%.3lf	%.3lf	%.3lf	%.3lf	%.3lf\n",iter,en[k].nbonds,en[k].coulomb,en[k].potens,en[k].kinet,en[k].esasa,en[k].totals);
#endif
	int rescount =0;double pdb =1.00;
    int count = 0;
//    int nc = 0;
    for(int m=0; m<sim.NC; m++) {
      for(int i=0; i<bp[k][m].nbox; i++) {
		  for (int kk=0; kk<mol[m].Nres; kk++) {
			  
			  
			for(int j=0; j<residue[k][rescount].Nsite; j++) {
				fprintf(io,"ATOM");
				if (atom[k][count].name[3]==0) fprintf(io,"%7d  ", count+1);
				else fprintf(io,"%7d ", count+1);
				if(	atom[k][count].name[1]==0) fprintf(io,"%s   ", atom[k][count].name);
				else if(	atom[k][count].name[2]==0) fprintf(io,"%s  ", atom[k][count].name);
				else fprintf(io,"%s ", atom[k][count].name);
				fprintf(io,"%s",amino_name[residue[k][rescount].type]);
				fprintf(io,"%6d  ", rescount+1);
#ifdef NEUTRAL  	  
				fprintf(io,"%10.3lf ", atnopbc[k][count].x);
				fprintf(io,"%7.3lf ",  atnopbc[k][count].y);
				fprintf(io,"%7.3lf  ", atnopbc[k][count].z);
#endif
#ifndef NEUTRAL  	  
				fprintf(io,"%10.3lf ", atom[k][count].x);
				fprintf(io,"%7.3lf ",  atom[k][count].y);
				fprintf(io,"%7.3lf  ", atom[k][count].z);
#endif
				fprintf(io,"%4.2lf  ", pdb);
				fprintf(io,"%4.2lf\n", pdb-1);
	  	  
				count++;
			}
			rescount++;
		  }fprintf(io,"END\n");
      }// loop i ends
	}// loop m ends
/*#ifdef RESTRAIN	// RESTRAIN is the wrong keyword, if anyone ever sees this
	fprintf(io,"ATOM");
	fprintf(io,"%7d  ",box[k].boxns+1);
	fprintf(io,"SPR SPR");
	fprintf(io,"%6d  ", rescount+1);
	fprintf(io,"%10.3lf ", restraint[k].x1);
	fprintf(io,"%7.3lf ",  restraint[k].y1);
	fprintf(io,"%7.3lf  ", restraint[k].z1);
	fprintf(io,"%4.2lf  ", pdb);
	fprintf(io,"%4.2lf\n", pdb-1);
	fprintf(io,"END\n");
	fprintf(io,"ATOM");
	fprintf(io,"%7d  ",box[k].boxns+2);
	fprintf(io,"SPR SPR");
	fprintf(io,"%6d  ", rescount+2);
	fprintf(io,"%10.3lf ", restraint[k].x2);
	fprintf(io,"%7.3lf ",  restraint[k].y2);
	fprintf(io,"%7.3lf  ", restraint[k].z2);
	fprintf(io,"%4.2lf  ", pdb);
	fprintf(io,"%4.2lf\n", pdb-1);
	fprintf(io,"END\n");
#endif
*/
fclose(io); // io for pdb completed 
/*-------------------------------------------------------------------------	*/
/*	Now print the pwi and psi distribution along with helicity				*/
/*-------------------------------------------------------------------------	*/
	int phic =1;
	double hel= ordparam[k].hel;
	for(int i=0; i<torsN[k]; i++) {
		if(tors[k][i].psitag==1){
			fprintf(fo,"%d	%d	%d	%d	%d	%.4lf	%.4lf	",phic,tors[k][i].a+1,
			tors[k][i].b+1,tors[k][i].c+1,tors[k][i].d+1,tors[k][i].psia/tors[k][i].count*180.0/PI,tors[k][i].psi*180.0/PI);
		}
		else if(tors[k][i].phitag==1){
			fprintf(fo,"%d	%d	%d	%d	%d	%.4lf	%.4lf	%.4lf\n",phic,tors[k][i].a+1,
			tors[k][i].b+1,tors[k][i].c+1,tors[k][i].d+1,tors[k][i].phia/tors[k][i].count*180.0/PI,tors[k][i].phi*180.0/PI,hel);
			phic++;
		}
	}
fclose(fo); // io for psi and osi completed
}
}
