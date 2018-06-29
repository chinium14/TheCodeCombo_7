/* ======================================================================== */
/* movie.cpp																*/
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
void movie (int flag, unsigned long iter)
{
  /* ================================================================== */
  /*                                                                    */
  /* Variables and pointer needed in the subroutine.                    */
  /*                                                                    */
  /* ================================================================== */
  #if defined(MOVIE) || defined(CRD) || defined(PHI)
  char name[80]; char name1[80]; char da[6] = "ABCDE";
  char amino_name[47][4]={"ABC","GLY","ALA","SER","CYS","VAL","THR","ILE","PRO","MET","ASP","ASN","LEU","LYS","GLU",
		"GLN","ARG","HIS","PHE","TYR","TRP","CYT","GUA","ADE","THY","HOH","SOD","CLA","LIP","CHO","EAM","POT","TRE", "GOL","F00", "F01",
		"F02","F03","F04","F05","F06","F07","F08","F09","F10","F11","MEO"};
  FILE *io; FILE *fo;
  #endif

  /* ================================================================== */
  /*                                                                    */
  /* If MOVIE/PHI/CRD is defined, save frames for the animation to		*/
  /* unique files each time.                                            */
  /*                                                                    */
  /*		BOX%d/MOVIE	 will contain pdb files if MOVIE is defined		*/
  /*		BOX%d/PHI	 will contain phi files	if PHI	 is defined		*/
  /*		BOX%d/CRD	 will contain crd files	if CRD	 is defined		*/
  /* ================================================================== */

  /* ----------------------------------------------- */
  /* Make the directories for Windows.               */
  /* ----------------------------------------------- */    

#ifdef WIN
  if (flag == 0){
	for(int k=0; k<sim.NB; k++){
		#ifdef MOVIE
			#ifdef MPI
			sprintf(name,"./OUTPUT/BOX%d/MOVIE",mpi.my_rank);
			#endif
			#ifndef MPI
			sprintf(name,"./OUTPUT/BOX%d/MOVIE",k);
			#endif
			_mkdir(name);
		#endif//MOVIE
		
		#ifdef PHI
			#ifdef MPI
			sprintf(name,"./OUTPUT/BOX%d/PHI",mpi.my_rank);
			#endif
			#ifndef MPI
			sprintf(name,"./OUTPUT/BOX%d/PHI",k);
			#endif
			_mkdir(name);
		#endif//PHI

		#ifdef CRD
			#ifdef MPI
			sprintf(name,"./OUTPUT/BOX%d/CRD",mpi.my_rank);
			#endif
			#ifndef MPI
			sprintf(name,"./OUTPUT/BOX%d/CRD",k);
			#endif
			_mkdir(name);
		#endif//CRD
	}
  }
#endif//WIN

  
  /* ----------------------------------------------- */
  /* Make the directories for Unix.                  */
  /* ----------------------------------------------- */
#ifndef WIN
  if (flag == 0){
	for(int k=0; k<sim.NB; k++){
		#ifdef MOVIE
			#ifdef MPI
   		    sprintf(name,"./OUTPUT/BOX%d/MOVIE",mpi.my_rank);
			#endif
			#ifndef MPI
			sprintf(name,"./OUTPUT/BOX%d/MOVIE",k);
			#endif		 
			mkdir(name,0755);
		#endif//MOVIE
		
		#ifdef PHI
			#ifdef MPI
			sprintf(name,"./OUTPUT/BOX%d/PHI",mpi.my_rank);
			#endif
			#ifndef MPI
			sprintf(name,"./OUTPUT/BOX%d/PHI",k);
			#endif		 
			mkdir(name,0755);
		#endif//PHI

		#ifdef CRD
			#ifdef MPI
		    sprintf(name,"./OUTPUT/BOX%d/CRD",mpi.my_rank);
			#endif
			#ifndef MPI
			sprintf(name,"./OUTPUT/BOX%d/CRD",k);
			#endif		 
			mkdir(name,0755);
		#endif//CRD
	}
  }
#endif//NWIN


#ifdef MOVIE
	for(int k=0; k<sim.NB; k++) {
		#ifdef MPI
		sprintf(name,"./OUTPUT/BOX%d/MOVIE/rasmol%.4d.pdb",mpi.my_rank,frames[k].count);
		#endif

		#ifndef MPI
		sprintf(name,"./OUTPUT/BOX%d/MOVIE/rasmol%.4d.pdb",k,frames[k].count);
		#endif		
		io = fopen(name,"w");

	  /* ----------------------------------------------- */
	  /* Make the file.                                  */
	  /* ----------------------------------------------- */    
	#ifdef PR_NPT
		fprintf(io,"REMARK Box Matrix for Parinello Rahman NPT\n");
		fprintf(io,"REMARK %lf %lf %lf\n",axes[k][0],axes[k][1],axes[k][2]);
		fprintf(io,"REMARK %lf %lf %lf\n",axes[k][3],axes[k][4],axes[k][5]);
		fprintf(io,"REMARK %lf %lf %lf\n",axes[k][6],axes[k][7],axes[k][8]);
	#endif

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
  		  
	#ifndef NEUTRAL
					fprintf(io,"%10.3lf ", atom[k][count].x);
					fprintf(io,"%7.3lf ",  atom[k][count].y);
					fprintf(io,"%7.3lf  ", atom[k][count].z);
	#endif
	#ifdef NEUTRAL
					fprintf(io,"%10.3lf ", atnopbc[k][count].x);
					fprintf(io,"%7.3lf ",  atnopbc[k][count].y);
					fprintf(io,"%7.3lf  ", atnopbc[k][count].z);
	#endif
					fprintf(io,"%4.2lf  ", pdb);
					fprintf(io,"%4.2lf\n", pdb-1);
	  		  
					count++;
				}
				rescount++;
			 }fprintf(io,"END\n");
		  }// loop i ends
		}// loop m ends

	fclose(io);
	}//kloop ends here
#endif //MOVIE

#ifdef PHI
/*-------------------------------------------------------------------------	*/
/*	Now print the pwi and psi distribution along with helicity				*/
/*-------------------------------------------------------------------------	*/
	for(int k=0; k<sim.NB; k++) {
		#ifdef MPI
		sprintf(name1,"./OUTPUT/BOX%d/PHI/phi%.4d.out",mpi.my_rank,frames[k].count);
		#endif
		#ifndef MPI
		sprintf(name1,"./OUTPUT/BOX%d/PHI/phi%.4d.out",k,frames[k].count);
		#endif		
		fo = fopen(name1,"w");
		int phic =1;
		double hel= ordparam[k].hel;
		for(int i=0; i<torsN[k]; i++) {
			if(tors[k][i].psitag==1){
				fprintf(fo,"%d	%d	%d	%d	%d	%.4lf	%.4lf	",phic,tors[k][i].a+1,
				tors[k][i].b+1,tors[k][i].c+1,tors[k][i].d+1,tors[k][i].psia/tors[k][i].count,tors[k][i].psi);
			}
			else if(tors[k][i].phitag==1){
				fprintf(fo,"%d	%d	%d	%d	%d	%.4lf	%.4lf	%.4lf\n",phic,tors[k][i].a+1,
				tors[k][i].b+1,tors[k][i].c+1,tors[k][i].d+1,tors[k][i].phia/tors[k][i].count,tors[k][i].phi,hel);
				phic++;
			}
		}
		fclose(fo);
	}
#endif //PHI 

#ifdef CRD
	FILE *fas;char name2[80];
/*-------------------------------------------------------------------------	*/
/*	Now print the crd files to save the trajectory							*/
/*-------------------------------------------------------------------------	*/
	for(int k=0; k<sim.NB; k++) {
          if(flag ==1){
		#ifdef MPI
		sprintf(name2,"./OUTPUT/BOX%d/CRD/simul%.4d.crd",mpi.my_rank,frames[k].count);
		#endif
		#ifndef MPI
		sprintf(name2,"./OUTPUT/BOX%d/CRD/simul%.4d.crd",k,frames[k].count);
		#endif
          }else if(flag == 2){
                #ifdef MPI
                sprintf(name2,"./OUTPUT/BOX%d/CRD/after_swap_simul%.4d.crd",mpi.my_rank,frames[k].count);
                #endif
                #ifndef MPI
                sprintf(name2,"./OUTPUT/BOX%d/CRD/after_swap_simul%.4d.crd",k,frames[k].count);
                #endif

          }else if(flag == 3){
                #ifdef MPI
                sprintf(name2,"./OUTPUT/BOX%d/CRD/before_swap_simul%.4d.crd",mpi.my_rank,frames[k].count);
                #endif
                #ifndef MPI
                sprintf(name2,"./OUTPUT/BOX%d/CRD/before_swap_simul%.4d.crd",k,frames[k].count);
                #endif

          }		
		fas = fopen(name2,"w");

  /* ----------------------------------------------- */
  /* Write heading and simulation parameters.        */
  /* ----------------------------------------------- */    
/*		fprintf(fas,"File saved as a part of trajectory\n");
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
	}// kloop ends here
#endif //CRD

	for(int k=0; k<sim.NB; k++) frames[k].count ++;
}
