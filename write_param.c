
/* ======================================================================== */
/* write_param.cpp                                                          */
/*		written by Nitin Rathore											*/
/* This subroutine writes the topology files for all the types of molecules */
/* in a format compatible with GROMACS										*/
/* ======================================================================== */

#include "defines.h"
char amino_name[47][4]={"ABC","GLY","ALA","SER","CYS","VAL","THR","ILE","PRO","MET","ASP","ASN","LEU","LYS","GLU",
		"GLN","ARG","HIS","PHE","TYR","TRP","CYT","GUA","ADE","THY","HOH","SOD","CLA","LIP","CHO","EAM","POT","TRE", "GOL","F00", "F01",
		"F02","F03","F04","F05","F06","F07","F08","F09","F10","F11","MEO"};

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void write_param (void)
{

  /* ================================================================== */
  /*                                                                    */
  /* Variables and pointer needed in the subroutine.                    */
  /*                                                                    */
  /* ================================================================== */
  char nameb[50];
  FILE *witp; FILE *fitp; FILE *fatp;

  /* ================================================================== */
  /*                                                                    */
  /* Write the atom types and L-J parameters                            */
  /*                                                                    */
  /* ================================================================== */

	sprintf(nameb,"./OUTPUT/GROMACS/ff_charmm.atp");
	fatp = fopen(nameb,"w");
	fprintf(fatp,"; Force -field File created by rathore to be used with gromacs\n");
	fprintf(fatp,"; name     mass\n");
	for (int kk =0; kk<NUM_ATOM; kk++) {
		fprintf(fatp,"%5s",atom_type[kk].name);
		fprintf(fatp,"%12.4lf\n",atom_type[kk].mass);
	}//kk

	fclose(fatp);
  
    double chg = 0.0;
	sprintf(nameb,"./OUTPUT/GROMACS/ff_charmmnb.itp");
	fitp = fopen(nameb,"w");
	fprintf(fitp,"[ atomtypes ]\n");
	fprintf(fitp,"; Force -field File created by rathore to be used with gromacs\n");
	fprintf(fitp,"; name     mass    charge   ptype          sigma      epsilon\n");
	for (int kk =0; kk<NUM_ATOM; kk++) {
		fprintf(fitp,"%5s",atom_type[kk].name);
		fprintf(fitp,"%12.4lf",atom_type[kk].mass);
		fprintf(fitp,"%9.3lf",chg);
		fprintf(fitp,"     A");
		/* Now for sigma, it has to be in nm and also scaled by 2^1/6	*/
		fprintf(fitp,"%16.6e%16.6e\n",sigtmp[0][atom_type[kk].atomid]/(10*pow(2.0,-(5.0/6.0))),-4.184*epstmp[0][atom_type[kk].atomid]);
	}//kk

//	fclose(fitp);


  /* ================================================================== */
  /*                                                                    */
  /* Write the parameters to file only for box 0                        */
  /*                                                                    */
  /* ================================================================== */
//  for(int k=0; k<sim.NB; k++) {
	int k=0;
	  int nmolec=0;	int rescount =0;int count = 0;
	  for(int cc=0; cc< sim.NC; cc++){
		int cgnr =1;

			sprintf(nameb,"./OUTPUT/GROMACS/molecule%d.itp",cc);
			witp = fopen(nameb,"w");
			fprintf(witp,";Molecular topology molecule type # %d in box %d\n",cc,k);
			fprintf(witp,"[ moleculetype ]\n");
			fprintf(witp,";  name  nrexcl\n");
			fprintf(witp,"TRE  3\n\n");
			fprintf(witp,"[ atoms ]\n");
			fprintf(witp,"; nr	type	resnr	resid	atom	cgnr	charge\n");

		  /* ----------------------------------------------- */
		  /* First write the atom types, charges etc.		 */
		  /* ----------------------------------------------- */ 
			int lb  = molec.fsite[nmolec];
			int ub  = molec.lsite[nmolec];
			count   = lb;
			int res_lb= molec.fres[nmolec];
			rescount= res_lb;

			for (int kk=0; kk<mol[cc].Nres; kk++) {	  
				for(int j=0; j<residue[k][rescount].Nsite; j++) {
					fprintf(witp,"%5d", count+1-lb);
					fprintf(witp,"%5s", atnopbc[k][count].t);
					fprintf(witp,"%6d", rescount+1-res_lb);
					fprintf(witp,"%5s",amino_name[residue[k][rescount].type]);
					fprintf(witp,"%5s", atom[k][count].name);
					fprintf(witp,"%6d", cgnr);
					fprintf(witp,"%8.3lf\n",pott[k][count].qq);
					count++;
					if ((count+1-lb)%10 ==0) cgnr++;
				}//for j
				rescount++;
			}//for kk

		  /* ----------------------------------------------- */
		  /* Bond parameters                                 */
		  /* ----------------------------------------------- */ 
			if(bondN[k]>0 && (ub-lb)>1){
				fprintf(witp,"\n[ bonds ]\n");
				fprintf(witp,";   ai    aj  funct	c0			c1\n");
				int c_bond =1;
				for(int i=0; i<bondN[k]; i++) {
					if(bond[k][i].a >=lb && bond[k][i].a <ub){
						fprintf(witp,"%5d%6d%6d%16.6e%16.6e	;%d\n",
							bond[k][i].a+1-lb,bond[k][i].b+1-lb,1,(bond[k][i].req/10.0),(bond[k][i].krbond*200),c_bond);
						c_bond++;
					}
				}
			}

		  /* ----------------------------------------------- */
		  /* Bend parameters                                 */
		  /* ----------------------------------------------- */
			if(bendN[k]>0 && (ub-lb)>2){
				fprintf(witp,"\n[ angles ]\n");
				fprintf(witp,";  ai	   aj    ak   funct	c0			c1\n");
				int c_bend =1;
				for(int i=0; i<bendN[k]; i++) {
					if(bend[k][i].a >=lb && bend[k][i].a <ub){
						fprintf(witp,"%5d%6d%6d%6d%16.6e%16.6e	;%d\n",
							bend[k][i].a+1-lb,bend[k][i].b+1-lb,bend[k][i].c+1-lb,1,bend[k][i].angeq*180.0/PI,bend[k][i].krbend*2.0,c_bend);
						c_bend++;
					}
				}
			}


		  /* ----------------------------------------------- */
		  /* Torsion parametes                               */
		  /* ----------------------------------------------- */
			if(torsN[k]>0 && (ub-lb)>3){
				fprintf(witp,"\n[ dihedrals ]\n");
				fprintf(witp,";   ai    aj    ak   al    funct		c0		c1\n");
				int c_tors =1;
				for(int i=0; i<torsN[k]; i++) {
					if(tors[k][i].a >=lb && tors[k][i].a <ub){
						for(int tt =0; tt<6; tt++){
							if(tors[k][i].kphi[tt] !=0.0){
								fprintf(witp,"%5d%6d%6d%6d%6d%16.6e%16.6e%16.6e	;%d\n",tors[k][i].a+1-lb,tors[k][i].b+1-lb,tors[k][i].c+1-lb,
									tors[k][i].d+1-lb,1,tors[k][i].delphi[tt]*180.0/PI,tors[k][i].kphi[tt],tors[k][i].nphi[tt],c_tors);
								c_tors++;
							}
							else break;
						}
					}
				}
			}

		  /* ----------------------------------------------- */
		  /* Improper parameters                               */
		  /* ----------------------------------------------- */ 

			if(imprN[k]>0 && (ub-lb)>3){
				fprintf(witp,"\n[ dihedrals ]\n");
				fprintf(witp,";   ai    aj    ak   al    funct		c0		c1\n");
				int c_impr=1;
				for(int i=0; i<imprN[k]; i++) {
					if(impr[k][i].a >=lb && impr[k][i].a <ub){
						fprintf(witp,"%5d%6d%6d%6d%6d%16.6e%16.6e	;%d\n",impr[k][i].a+1-lb,impr[k][i].b+1-lb,impr[k][i].c+1-lb,
									impr[k][i].d+1-lb,2,impr[k][i].angeq*180.0/PI,impr[k][i].kimpr*2.0,c_impr);
						c_impr++;
					}
				}
			}
		  /* ----------------------------------------------- */
		  /* 1-4  pairs	(S14)                                */
		  /* ----------------------------------------------- */ 

			#ifdef S14
			if(in14N[k]>0 && (ub-lb)>3 ){
				int c_in14=1;
				fprintf(witp,"\n[ pairs ]\n");
				if(nmolec==0) fprintf(fitp,"\n[ pairtypes ]\n");
				fprintf(witp,";   i    j   func          c6      c12\n");
				fprintf(fitp,";   i    j   func          c6      c12\n");
				for(int i=0; i<in14N[k]; i++){
					if(in14[k][i].a >=lb && in14[k][i].a <ub){
//						double C6	= 4.0*ljset14[atom[k][in14[k][i].a+1-lb].atomid][atom[k][in14[k][i].d+1-lb].atomid].eps*
//										   pow((ljset14[atom[k][in14[k][i].a+1-lb].atomid][atom[k][in14[k][i].d+1-lb].atomid].sig/(10.0*pow(2,1.0/6.0))),6.0);
//							double C12	= 4.0*ljset14[atom[k][in14[k][i].a+1-lb].atomid][atom[k][in14[k][i].d+1-lb].atomid].eps*
//										   pow((ljset14[atom[k][in14[k][i].a+1-lb].atomid][atom[k][in14[k][i].d+1-lb].atomid].sig/(10.0*pow(2,1.0/6.0))),12.0);
//						fprintf(witp,"%5d%6d%6d%16.6e%16.6e	;%d\n",in14[k][i].a+1-lb,in14[k][i].d+1-lb,1,
//							ljset14[atom[k][in14[k][i].a+1-lb].atomid][atom[k][in14[k][i].d+1-lb].atomid].sig/(10.0*pow(2,1.0/6.0)),
//							ljset14[atom[k][in14[k][i].a+1-lb].atomid][atom[k][in14[k][i].d+1-lb].atomid].eps,c_in14);
//							C6,C12,c_in14);
						fprintf(witp,"%5d%6d%6d	;%d\n",in14[k][i].a+1-lb,in14[k][i].d+1-lb,1,c_in14);
						int tag14=1;
						for(int ii=0; ii<i; ii++){
							if((atom[k][in14[k][ii].a-lb].atomid == atom[k][in14[k][i].a-lb].atomid && atom[k][in14[k][ii].d-lb].atomid == atom[k][in14[k][i].d-lb].atomid) ||
								(atom[k][in14[k][ii].a-lb].atomid == atom[k][in14[k][i].d-lb].atomid && atom[k][in14[k][ii].d-lb].atomid == atom[k][in14[k][i].a-lb].atomid)) tag14=0;
						}
						if (tag14!=0){
							fprintf(fitp,"%5s%5s%6d%16.6e%16.6e	;%d\n",atom[k][in14[k][i].a-lb].t,atom[k][in14[k][i].d-lb].t,1,
								ljset14[atom[k][in14[k][i].a-lb].atomid][atom[k][in14[k][i].d-lb].atomid].sig/(10.0*pow(2,1.0/6.0)),
								ljset14[atom[k][in14[k][i].a-lb].atomid][atom[k][in14[k][i].d-lb].atomid].eps,c_in14);
						}
						c_in14++;
					}
				}
			}
			#endif //S14
			fprintf(witp,"\n");
			fclose(witp);
			nmolec	+= bp[k][cc].nbox;
	  }//cc
//  }//k not needed as k=0 only.
fclose(fitp);
}
