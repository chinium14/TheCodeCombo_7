/* ======================================================================== */
/* read_topology.cpp                                                        */
/*                                                                          */
/*		This subroutine reads in the molecular topology for all molecules.  */
/* The topology contains the information telling the atoms involved in each */
/* bond, bend, torsion, and improper dihedral that the system contains.     */
/* The file, named simul#.psf, comes from a CHARMM output.                  */
/*                                                                          */  
/* Note: There is a separate simul#.psf file for each box.                  */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"
int duplicate14_check(int, int,int, int*, int*, int);
/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void read_topology (int ibox)
{

  /* ================================================================== */
  /*                                                                    */
  /* Variables and pointers needed in the subroutine.                   */
  /*                                                                    */
  /* ================================================================== */
char amino_name[47][4]={"ABC","GLY","ALA","SER","CYS","VAL","THR","ILE","PRO","MET","ASP","ASN","LEU","LYS","GLU",
		"GLN","ARG","HIS","PHE","TYR","TRP","CYT","GUA","ADE","THY","HOH","SOD","CLA","LIP","CHO","EAM","POT","TRE", "GOL","F00", "F01",
		"F02","F03","F04","F05","F06","F07","F08","F09","F10","F11","MEO"};
//	char  amino_id[21]={'X','G','A','S','C','V','T','I','P','M','D','N','L','K','E','Q','R','H','F','Y','W'};
	int k = ibox;
  int t1; int flag=0; int rescount=0; int temp;
  char namet[50], tt[150];char am[5]; char tempc[5];
  FILE *ind;
  


  /* ================================================================== */
  /*                                                                    */
  /* Begin reading in the data from file to store in the above          */
  /* structures.                                                        */
  /*                                                                    */
  /* ================================================================== */
//  for(int k=0; k<sim.NB; k++) {
  /* ----------------------------------------------- */
  /* Open the file.                                  */
  /* ----------------------------------------------- */    

#ifdef MPI
	//sprintf(namet,"./INPUT/simul%d.psf",mpi.my_rank);
	sprintf(namet,"./INPUT/simul.psf");
#endif
#ifndef MPI
	//sprintf(namet,"./INPUT/simul%d.psf",k);
	sprintf(namet,"./INPUT/simul.psf");
#endif
     if( NULL == (ind=fopen(namet,"r")) ) {
      fprintf(stdout,"ERROR: molecular topology file %s does not exist!!!\n",namet);
      exit(1);
    }

		fgets(tt,150,ind);
		fgets(tt,150,ind);
		fgets(tt,150,ind);
		fgets(tt,150,ind);
		fgets(tt,150,ind);
	
//  int i=0;
  int Nsite=0;
  int Nmolec = -1;
  int sites;
  for (int cc=0; cc<sim.NC; cc++){

  /* ----------------------------------------------- */
  /* Read in the information about each atom/sites.  */
  /* ----------------------------------------------- */    
		fscanf(ind,"%d",&sites);  fgets(tt,150,ind);
		if(sites!= (mol[cc].Nsite*bp[k][cc].nbox)){
			fprintf(stdout,"Number of total sites in crd (%i) and psf (%i) for molecule type %d  do not match!!! \n", mol[cc].Nsite*bp[k][cc].nbox, sites, cc);
			exit(102);
		}
		for (int i=0; i<(mol[cc].Nsite*bp[k][cc].nbox); i++){
			if(mol[cc].Nsite == 1) Nmolec++;
      else if(i%(mol[cc].Nsite+1)==0) Nmolec++;
//		while(i<box[k].boxns){
			fscanf(ind,"%d",&temp);				// the site number
			fscanf(ind,"%6s",tempc);			// dummy variable to read in pept name in file
			fscanf(ind,"%d",&rescount);			// the residue number

                        if(rescount > box[k].boxnres)
                        {       
                                fprintf(stdout, "Atom %i of molecule %i in simul.psf has residue number %i, but total number of residues in box %i in mol.input is %i !!!\n", i, cc, rescount, k, box[k].boxnres);
                                exit(102);
                        } 

			atnopbc[k][Nsite].n=rescount-1;			// subtract off one to make first residue 0
      atnopbc[k][Nsite].molec = Nmolec;
			residue[k][rescount-1].Nsite++;
			residue[k][rescount-1].chain=cc;	// setting the chain to the molecule type
			
			fscanf(ind, "%4s", am);			// read in the amino acid name
       /* *********************************************** */
       /* Compare the amino acid name to list above to    */
	   /* assign the appropriate number code in .type.    */
       /* *********************************************** */
			for(int ii=1;ii<47;ii++){
				if(strcmp(amino_name[ii],am)==0)
			 residue[k][rescount-1].type=ii;
			}
       /* *********************************************** */
       /* Read in the site/atom name and id number.       */
       /* *********************************************** */
			fscanf(ind,"%4s",  atnopbc[k][Nsite].name);
			fscanf(ind,"%d",   &atnopbc[k][Nsite].atomid);
			for (int kk =0; kk<NUM_ATOM; kk++) {
       /* *********************************************** */
       /* Compare the id number read to those read in     */
	   /* read_atom.  If the same, copy the name to       */
	   /* atnopbc[k][i].t.  Note: This name is more       */
	   /* descriptive than those read in read_atom.       */
       /* *********************************************** */
			if (atnopbc[k][Nsite].atomid== atom_type[kk].atomid){
				strcpy(atnopbc[k][Nsite].t,atom_type[kk].name);
				break;
			}
			}
		
       /* *********************************************** */
       /* Read in charge.                                 */
       /* *********************************************** */
			fscanf(ind,"%lf",   &atnopbc[k][Nsite].q);
			fgets(tt,150,ind);
			
//			i++;
			Nsite++;
//			if(Nsite==mol[cc].Nsite){
///				Nsite=0;
//				break;
//			}
		}//for i
	}//for cc
    



  /* ----------------------------------------------- */
  /* Read in the information about each bond.        */
  /* ----------------------------------------------- */    
       /* *********************************************** */
       /* Read in the number of bonds.                    */
       /* *********************************************** */
    fscanf(ind,"%d\n",&bondN[k]);fgets(tt,80,ind);
       /* *********************************************** */
       /* Allocate memory for this number of bonds.       */
       /* *********************************************** */
    int rb = bondN[k] + MEXTRA;
    bond[k] = (struct bonds*) calloc(rb,sizeof(struct bonds));
		if(bond[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for bond[k]\n"); exit(11);}
       /* *********************************************** */
       /* Variables needed.                               */
       /* *********************************************** */
    int count=0; t1=0;
	for(int i=0; i<bondN[k]; i++) {
      
       /* *********************************************** */
       /* Read in the two atoms/site numbers involved     */
	   /* in each bond. Subtract 1 to make the site       */
	   /* numbers begin at 0.                             */
       /* *********************************************** */
      fscanf(ind,"%d", &t1); bond[k][i].a = t1-1;
      fscanf(ind,"%d", &t1); bond[k][i].b = t1-1;
	  count++;
	  for (int nb =0; nb<NUM_BOND; nb++) {
       /* *********************************************** */
       /* Compare the site id numbers(type of site) to    */
	   /* the list read in read_atom.  When they match,   */
	   /* assign that bond the appropriate req and        */
	   /* krbond.  If no matches are found, there is an   */
	   /* error in the input files so exit.               */
       /* *********************************************** */
		  if ((bond_prop[nb].type1 == atnopbc[k][bond[k][i].a].atomid && bond_prop[nb].type2 == atnopbc[k][bond[k][i].b].atomid) || 
			  (bond_prop[nb].type1 == atnopbc[k][bond[k][i].b].atomid && bond_prop[nb].type2 == atnopbc[k][bond[k][i].a].atomid)) {
			  bond[k][i].req    = bond_prop[nb].eqbond;
			  bond[k][i].krbond = bond_prop[nb].kbond;
			  flag=1;break; 
		  }
	  }
		  if (flag==0) {
		  fprintf(stdout,"ERROR: bond parameters not defined between atoms %i and %i !!!\n", bond[k][i].a+1, bond[k][i].b+1);
	      exit(1);
		  }
	  

       /* *********************************************** */
       /* The input file contains four pair on a line,    */
	   /* so after four pairs are read, go to the next    */
	   /* line.                                           */
       /* *********************************************** */
      if (count==4) {
		  fscanf(ind,"\n");
		  count=0;
	  }
       /* *********************************************** */
       /* If the input file doesn't have four bonds on    */
	   /* the last line, skip a line in the end.          */
       /* *********************************************** */
	  else if (i== bondN[k]-1) fscanf(ind, "\n");
       /* *********************************************** */
       /* Convert from kcal/mol to kJ/mol.                */
       /* *********************************************** */
	  bond[k][i].krbond  *= 4.184;
	  flag=0;
#ifdef BEAD
	  bond[k][i].krbond = 100.0 * sim.T[k] * RG * 0.001;
	  bond[k][i].req    = 0.25 * 4.0;
#endif

    } 
    
  /* ----------------------------------------------- */
  /* Read in the information about each bend.        */
  /* ----------------------------------------------- */    

       /* *********************************************** */
       /* Read in the number of bends.                    */
       /* *********************************************** */
    fscanf(ind,"%d\n",&bendN[k]);fgets(tt,80,ind);
       /* *********************************************** */
       /* Allocate memory for this number of bends.       */
       /* *********************************************** */
    int rc = bendN[k] + MEXTRA;
    bend[k] = (struct bends*) calloc(rc,sizeof(struct bends));
		if(bend[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for bend[i]\n"); exit(11);}
	count =0; t1=0;
       /* *********************************************** */
       /* Variables needed.                               */
       /* *********************************************** */
    for(int i=0; i<bendN[k]; i++) {
      
       /* *********************************************** */
       /* Read in the three atoms/site numbers involved   */
	   /* in each bend. Subtract 1 to make the site       */
	   /* numbers begin at 0.                             */
       /* *********************************************** */
      fscanf(ind,"%d", &t1);  bend[k][i].a= t1-1;
      fscanf(ind,"%d", &t1);  bend[k][i].b= t1-1;
      fscanf(ind,"%d", &t1);  bend[k][i].c= t1-1;
	  count++;
       /* *********************************************** */
       /* Compare the site id numbers(type of site) to    */
	   /* the list read in read_atom.  When they match,   */
	   /* assign that bend the appropriate eqbend and     */
	   /* kbend.  If no matches are found, there is an   */
	   /* error in the input files so exit.               */
       /* *********************************************** */
	  for (int nb =0; nb< NUM_BEND; nb++) {
		  if((bend_prop[nb].type1==atnopbc[k][bend[k][i].a].atomid && bend_prop[nb].type2==atnopbc[k][bend[k][i].b].atomid 
			  && bend_prop[nb].type3==atnopbc[k][bend[k][i].c].atomid) ||  (bend_prop[nb].type1==atnopbc[k][bend[k][i].c].atomid 
			  && bend_prop[nb].type2==atnopbc[k][bend[k][i].b].atomid && bend_prop[nb].type3==atnopbc[k][bend[k][i].a].atomid)) {
			  bend[k][i].angeq  = bend_prop[nb].eqbend;
			  bend[k][i].krbend = bend_prop[nb].kbend;
			  flag=1;break; 
		  }
	  }
		  if (flag==0) {
			  fprintf(stdout,"ERROR: bend parameter not defined between atoms %d, %d, %d !!!\n",bend[k][i].a+1,bend[k][i].b+1,bend[k][i].c+1);
	          exit(1);
		  }
	  
       /* *********************************************** */
       /* The input file contains three triples on a line */
	   /* so after three triples are read, go to the next */
	   /* line.                                           */
       /* *********************************************** */
	  if (count==3) {
		  fscanf(ind, "\n");
		  count =0;
	  }
       /* *********************************************** */
       /* If the input file doesn't have three bends on   */
	   /* the last line, skip a line in the end.          */
       /* *********************************************** */
	  else if (i== bendN[k]-1) fscanf(ind, "\n");
       /* *********************************************** */
       /* Convert from kcal/mol to kJ/mol and degrees to  */
	   /* radians.                                        */
       /* *********************************************** */
	  bend[k][i].krbend *= 4.184;
      bend[k][i].angeq  *= PI/180.0;flag=0;
    }
	

  /* ----------------------------------------------- */
  /* Read in the information about each torsion.     */
  /* ----------------------------------------------- */    

       /* *********************************************** */
       /* Read in the number of torsions.                 */
       /* *********************************************** */
   fscanf(ind,"%d\n",&torsN[k]);fgets(tt,80,ind);
       /* *********************************************** */
       /* Variables needed.                               */
       /* *********************************************** */
    int rd = torsN[k] + MEXTRA;
       /* *********************************************** */
       /* Allocate memory for this number of torsions.    */
       /* *********************************************** */
    tors[k] = (struct torsions*) calloc(rd,sizeof(struct torsions));
		if(tors[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for tors[k]\n"); exit(11);}
    count=0; t1=0;
	for(int i=0; i<torsN[k]; i++) {
       /* *********************************************** */
       /* Read in the four atoms/site numbers involved    */
	   /* in each bend. Subtract 1 to make the site       */
	   /* numbers begin at 0.                             */
       /* *********************************************** */
      fscanf(ind,"%d", &t1);  tors[k][i].a=t1-1;
      fscanf(ind,"%d", &t1);  tors[k][i].b=t1-1;
      fscanf(ind,"%d", &t1);  tors[k][i].c=t1-1;
      fscanf(ind,"%d", &t1);  tors[k][i].d=t1-1;
	  count++;int tt=0;
	  for (int nb=0; nb< NUM_TORS; nb++) {
       /* *********************************************** */
       /* Compare the site id numbers(type of site) to    */
	   /* the list read in read_atom.  When they match,   */
	   /* assign that torsion the appropriate kphi,nphi   */
	   /* and delphi. If no matches are found, there is   */
	   /* an error in the input files so exit.            */
       /* *********************************************** */
		  
		  if ((dihedral_prop[nb].type1==atnopbc[k][tors[k][i].a].atomid && dihedral_prop[nb].type2==atnopbc[k][tors[k][i].b].atomid 
			  &&dihedral_prop[nb].type3==atnopbc[k][tors[k][i].c].atomid  && dihedral_prop[nb].type4==atnopbc[k][tors[k][i].d].atomid)
			  || (dihedral_prop[nb].type1==atnopbc[k][tors[k][i].d].atomid && dihedral_prop[nb].type2==atnopbc[k][tors[k][i].c].atomid 
			  &&dihedral_prop[nb].type3==atnopbc[k][tors[k][i].b].atomid  && dihedral_prop[nb].type4==atnopbc[k][tors[k][i].a].atomid)) {
			  
			  tors[k][i].kphi[tt]= 4.184*dihedral_prop[nb].kphi;
			  tors[k][i].nphi[tt]= dihedral_prop[nb].nphi;
			  tors[k][i].delphi[tt]= PI*dihedral_prop[nb].delphi/180; tt++;
			  flag=1;//break;
		  }
		  else if ((( dihedral_prop[nb].type1==100 && dihedral_prop[nb].type4 ==100 && dihedral_prop[nb].type2==atnopbc[k][tors[k][i].b].atomid 
			  && dihedral_prop[nb].type3==atnopbc[k][tors[k][i].c].atomid) || (dihedral_prop[nb].type1==100 && dihedral_prop[nb].type4 ==100 
			  && dihedral_prop[nb].type2==atnopbc[k][tors[k][i].c].atomid &&dihedral_prop[nb].type3==atnopbc[k][tors[k][i].b].atomid)) && flag==0) {
			  tors[k][i].kphi[0]= 4.184*dihedral_prop[nb].kphi;
			  tors[k][i].nphi[0]= dihedral_prop[nb].nphi;
			  tors[k][i].delphi[0]= PI*dihedral_prop[nb].delphi/180;
			  flag=1;break; 
		  }
	  }
		  if (flag==0) {
			  fprintf(stdout,"ERROR: dihedral parameters not defined for atom types	%d	%d	%d	%d!!!\n",atnopbc[k][tors[k][i].a].atomid,
				  atnopbc[k][tors[k][i].b].atomid,atnopbc[k][tors[k][i].c].atomid,atnopbc[k][tors[k][i].d].atomid);
	          exit(1);
		  }
	  flag=0;
       /* *********************************************** */
       /* The input file contains two atoms quads on a    */
	   /* line, so after two quads are read, go to the    */
	   /* next line.                                      */
       /* *********************************************** */
	  if (count==2) {
		  fscanf(ind, "\n");
		  count =0;
	  }
       /* *********************************************** */
       /* If the input file doesn't have two torsions on  */
	   /* the last line, skip a line in the end.          */
       /* *********************************************** */
	  else if (i== torsN[k]-1) fscanf(ind, "\n");
	//  tors[k][i].kphi   *= 4.184;
    // tors[k][i].delphi *= PI/180.0; // need correction
    }

		/* ------------------------------------------------	*/
		/*				CHARMM 19 Convention				*/
		/* Make a list of the main chain phi and psi angles	*/
		/* Phi defined as: C(i-1) - N(i) - CA(i) - N(i+1)	*/
		/* Psi defined as:	N(i) - CA(i) - C(i)	 - N(i+1)	*/
		/*													*/
		/*		CA	=	12 if CH1E	alpha carbon			*/ 
		/*			=	13 if CH2E	alpha carbon			*/
		/*		N	=	31 (N)/34 (NH1E)/35 (NH2E)/36 (NH3E)*/
		/*			=	38 (NH1)/39 (NH2)/40 (NH3)			*/
		/*		C	=	11	for carbonyl carbon C=O			*/
		/* ------------------------------------------------	*/

#ifdef NEUTRAL
	for(int i=0; i<torsN[k]; i++) {
		if (atnopbc[k][tors[k][i].a].atomid == 11 &&(atnopbc[k][tors[k][i].b].atomid == 38 ||
			atnopbc[k][tors[k][i].b].atomid == 31 || atnopbc[k][tors[k][i].b].atomid == 34)&&
		   (atnopbc[k][tors[k][i].c].atomid == 12 || atnopbc[k][tors[k][i].c].atomid == 13)&& 
		    atnopbc[k][tors[k][i].d].atomid == 11){
			tors[k][i].phitag=1;
			tors[k][i].psitag=0;
			tors[k][i].thetag=0;
		}
		else if((atnopbc[k][tors[k][i].a].atomid == 38 || atnopbc[k][tors[k][i].a].atomid == 40 ||
				 atnopbc[k][tors[k][i].a].atomid == 31 || atnopbc[k][tors[k][i].a].atomid == 34 ||
				 atnopbc[k][tors[k][i].a].atomid == 35 || atnopbc[k][tors[k][i].a].atomid == 36 ||
				 atnopbc[k][tors[k][i].a].atomid == 39)&& 
				(atnopbc[k][tors[k][i].b].atomid == 12 || atnopbc[k][tors[k][i].b].atomid == 13)&&
				 atnopbc[k][tors[k][i].c].atomid == 11 &&(atnopbc[k][tors[k][i].d].atomid == 38 ||
				 atnopbc[k][tors[k][i].d].atomid == 31 || atnopbc[k][tors[k][i].d].atomid == 34)){
			tors[k][i].phitag=0;
			tors[k][i].psitag=1;
			tors[k][i].thetag=0;
		}
		else{
			tors[k][i].phitag=0;
			tors[k][i].psitag=0;
			tors[k][i].thetag=0;
		}
	}
#endif

		/* ------------------------------------------------	*/
		/*				CHARMM 22 Convention				*/
		/* Make a list of the main chain phi and psi angles	*/
		/* Phi defined as: C(i-1) - N(i) - CA(i) - N(i+1)	*/
		/* Psi defined as:	N(i) - CA(i) - C(i)	 - N(i+1)	*/
		/*													*/
		/*		CA	=	22 if CH1E	alpha carbon			*/ 
		/*			=	23 if CH2E	alpha carbon			*/
		/*			=	29 if Proline ring carbon			*/
		/*		N	=	50 (Proline N)/54 (NH1)/55 (NH2)	*/
		/*			=	56 (NH3)/59 (terminal proline NH2+)	*/
		/*		C	=	20,32,33  for carbonyl carbon C=O	*/
		/* ------------------------------------------------	*/

#ifndef NEUTRAL
#ifndef BETAPEP
	for(int i=0; i<torsN[k]; i++) {
		if ((atnopbc[k][tors[k][i].a].atomid == 20 || atnopbc[k][tors[k][i].a].atomid == 32 ||
			 atnopbc[k][tors[k][i].a].atomid == 33)&&(atnopbc[k][tors[k][i].b].atomid == 50 ||
			 atnopbc[k][tors[k][i].b].atomid == 54)&&
		    (atnopbc[k][tors[k][i].c].atomid == 22 || atnopbc[k][tors[k][i].c].atomid == 23 ||
			 atnopbc[k][tors[k][i].c].atomid == 29)&&(atnopbc[k][tors[k][i].d].atomid == 20 ||
			 atnopbc[k][tors[k][i].d].atomid == 32 || atnopbc[k][tors[k][i].d].atomid == 33)){
			tors[k][i].phitag=1;
			tors[k][i].psitag=0;
			tors[k][i].thetag=0;
		}
		else if((atnopbc[k][tors[k][i].a].atomid == 50 || atnopbc[k][tors[k][i].a].atomid == 54 ||
				 atnopbc[k][tors[k][i].a].atomid == 55 || atnopbc[k][tors[k][i].a].atomid == 56 ||
				 atnopbc[k][tors[k][i].a].atomid == 59)&&(atnopbc[k][tors[k][i].b].atomid == 22 || 
				 atnopbc[k][tors[k][i].b].atomid == 23 || atnopbc[k][tors[k][i].b].atomid == 29)&&
				(atnopbc[k][tors[k][i].c].atomid == 20 || atnopbc[k][tors[k][i].c].atomid == 32 || 
				 atnopbc[k][tors[k][i].c].atomid == 33)&&(atnopbc[k][tors[k][i].d].atomid == 50||
				 atnopbc[k][tors[k][i].d].atomid == 54)){
			tors[k][i].phitag=0;
			tors[k][i].psitag=1;
			tors[k][i].thetag=0;
		}
		else{
			tors[k][i].phitag=0;
			tors[k][i].psitag=0;
			tors[k][i].thetag=0;
		}
	}
#endif//no betapep
#ifdef BETAPEP
	for(int i=0; i<torsN[k]; i++) {

		if ((strcmp(atnopbc[k][tors[k][i].a].name,"N")==0 && strcmp(atnopbc[k][tors[k][i].d].name,"C")==0)||
			(strcmp(atnopbc[k][tors[k][i].d].name,"N")==0 && strcmp(atnopbc[k][tors[k][i].a].name,"C")==0)){
			tors[k][i].thetag=1;
			tors[k][i].phitag=0;
			tors[k][i].psitag=0;
		}
		else if ((strcmp(atnopbc[k][tors[k][i].a].name,"CB")==0 && strcmp(atnopbc[k][tors[k][i].d].name,"N")==0)||
			(strcmp(atnopbc[k][tors[k][i].d].name,"CB")==0 && strcmp(atnopbc[k][tors[k][i].a].name,"N")==0)){
			tors[k][i].thetag=0;
			tors[k][i].psitag=1;
			tors[k][i].phitag=0;
		}
		else if ((strcmp(atnopbc[k][tors[k][i].a].name,"C")==0 && strcmp(atnopbc[k][tors[k][i].d].name,"CA")==0)||
			(strcmp(atnopbc[k][tors[k][i].d].name,"C")==0 && strcmp(atnopbc[k][tors[k][i].a].name,"CA")==0)){
			tors[k][i].thetag=0;
			tors[k][i].psitag=0;
			tors[k][i].phitag=1;
		}
		else{
			tors[k][i].thetag=0;
			tors[k][i].phitag=0;
			tors[k][i].psitag=0;
		}
	}
#endif //beta pep

#endif

  /* ----------------------------------------------- */
  /* Read in the information about each improper     */
  /* dihedral.                                       */
  /* ----------------------------------------------- */    

       /* *********************************************** */
       /* Read in the number of improper dihedrals.       */
       /* *********************************************** */
   fscanf(ind,"%d\n",&imprN[k]);fgets(tt,80,ind);
       /* *********************************************** */
       /* Variables needed.                               */
       /* *********************************************** */
    int re = imprN[k] + MEXTRA;
       /* *********************************************** */
       /* Allocate memory for this number of improper     */
	   /* dihedrals.                                      */
       /* *********************************************** */
    impr[k] = (struct impropers*) calloc(re,sizeof(struct impropers));
		if(impr[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for impr[k]\n"); exit(11);}
    count=0; t1=0;
	for(int i=0; i<imprN[k]; i++) {
       /* *********************************************** */
       /* Read in the four atoms/site numbers involved    */
	   /* in each improper angle. Subtract 1 to make      */
	   /* the site numbers begin at 0.                    */
       /* *********************************************** */
      fscanf(ind,"%d", &t1);  impr[k][i].a=t1-1;
      fscanf(ind,"%d", &t1);  impr[k][i].b=t1-1;
      fscanf(ind,"%d", &t1);  impr[k][i].c=t1-1;
      fscanf(ind,"%d", &t1);  impr[k][i].d=t1-1;
	  count++;
	  for (int nb=0; nb< NUM_IMPR; nb++) {
       /* *********************************************** */
       /* Compare the site id numbers(type of site) to    */
	   /* the list read in read_atom.  When they match,   */
	   /* assign that improper dihedral the appropriate   */
	   /* kimpr and angeq. If no matches are found, there */
	   /* is an error in the input files so exit.         */
       /* *********************************************** */
		  
		  if ((improper_prop[nb].type1==atnopbc[k][impr[k][i].a].atomid && improper_prop[nb].type2==atnopbc[k][impr[k][i].b].atomid 
			  &&improper_prop[nb].type3==atnopbc[k][impr[k][i].c].atomid  && improper_prop[nb].type4==atnopbc[k][impr[k][i].d].atomid)
			  || (improper_prop[nb].type1==atnopbc[k][impr[k][i].d].atomid && improper_prop[nb].type2==atnopbc[k][impr[k][i].c].atomid 
			  &&improper_prop[nb].type3==atnopbc[k][impr[k][i].b].atomid  && improper_prop[nb].type4==atnopbc[k][impr[k][i].a].atomid)) {
			  
			  impr[k][i].kimpr= 4.184*improper_prop[nb].kimpr;
			  impr[k][i].angeq= PI*improper_prop[nb].angeq/180; 
			  flag=1;//break;
		  }
		  else if ((( improper_prop[nb].type2==100 && improper_prop[nb].type3 ==100 && improper_prop[nb].type1==atnopbc[k][impr[k][i].a].atomid 
			  && improper_prop[nb].type4==atnopbc[k][impr[k][i].d].atomid) || (improper_prop[nb].type2==100 && improper_prop[nb].type3 ==100 
			  && improper_prop[nb].type1==atnopbc[k][impr[k][i].d].atomid &&improper_prop[nb].type4==atnopbc[k][impr[k][i].a].atomid)) && flag==0) {
			  
			  impr[k][i].kimpr= 4.184*improper_prop[nb].kimpr;
			  impr[k][i].angeq= PI*improper_prop[nb].angeq/180; 
			  flag=1;break; 
		  }
	  }
		  if (flag==0) {
			  fprintf(stdout,"ERROR: Improper dihedral parameters not defined for atoms %i, %i, %i, %i !!!\n", impr[k][i].a+1, impr[k][i].b+1, impr[k][i].c+1, impr[k][i].d+1);
	          exit(1);
		  }
	  flag=0;
       /* *********************************************** */
       /* The input file contains two atoms quads on a    */
	   /* line, so after two quads are read, go to the    */
	   /* next line.                                      */
       /* *********************************************** */
	  if (count==2) {
		  fscanf(ind, "\n");
		  count =0;
	  }
       /* *********************************************** */
       /* If the input file doesn't have two torsions on  */
	   /* the last line, skip a line in the end.          */
       /* *********************************************** */
	  else if (i== imprN[k]-1) fscanf(ind, "\n");

    }//for imprN[k]
	
  /* ----------------------------------------------- */
  /* Read in the information about each hydrogen     */
  /* bond donor.                                     */
  /* ----------------------------------------------- */    

       /* *********************************************** */
       /* Read in the number of hydrogen bond donors.     */
       /* *********************************************** */
   fscanf(ind,"%d\n",&hdonN[k]);fgets(tt,80,ind);
       /* *********************************************** */
       /* Variables needed.                               */
       /* *********************************************** */
    int ri = hdonN[k] + MEXTRA;
       /* *********************************************** */
       /* Allocate memory for this number of donors.      */
       /* *********************************************** */
    hdon[k] = (struct donors*) calloc(ri,sizeof(struct donors));
		if(hdon[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for hdor[k]\n"); exit(11);}
    count=0; t1=0;
	for(int i=0; i<hdonN[k]; i++) {
       /* *********************************************** */
       /* Read in the two atoms/site numbers involved     */
	   /* in each hydrogen bond donor. Subtract 1 to make */ 
	   /* the site numbers begin at 0.                    */
       /* *********************************************** */
      fscanf(ind,"%d", &t1); hdon[k][i].a = t1-1;
      fscanf(ind,"%d", &t1); hdon[k][i].b = t1-1;
	  count++;	  

       /* *********************************************** */
       /* The input file contains four pair on a line,    */
	   /* so after four pairs are read, go to the next    */
	   /* line.                                           */
       /* *********************************************** */
      if (count==4) {
		  fscanf(ind,"\n");
		  count=0;
	  }
       /* *********************************************** */
       /* If the input file doesn't have four bonds on    */
	   /* the last line, skip a line in the end.          */
       /* *********************************************** */
	  else if (i== hdonN[k]-1) fscanf(ind, "\n");

    }

  /* ----------------------------------------------- */
  /* Read in the information about each hydrogen     */
  /* bond acceptor.                                  */
  /* ----------------------------------------------- */    

       /* *********************************************** */
       /* Read in the number of hydrogen bond acceptors.  */
       /* *********************************************** */
   fscanf(ind,"%d\n",&haccN[k]);fgets(tt,80,ind);
       /* *********************************************** */
       /* Variables needed.                               */
       /* *********************************************** */
    int rj = haccN[k] + MEXTRA;
       /* *********************************************** */
       /* Allocate memory for this number of acceptors.   */
       /* *********************************************** */
    hacc[k] = (struct acceptors*) calloc(rj,sizeof(struct acceptors));
		if(hacc[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for hacc[k]\n"); exit(11);}
    count=0; t1=0;
	for(int i=0; i<haccN[k]; i++) {
       /* *********************************************** */
       /* Read in the two atoms/site numbers involved     */
	   /* in each hydrogen bond acceptor. Subtract 1 to   */ 
	   /* make the site numbers begin at 0.               */
       /* *********************************************** */
      fscanf(ind,"%d", &t1); hacc[k][i].a = t1-1;
      fscanf(ind,"%d", &t1); hacc[k][i].b = t1-1;
	  count++;	  

       /* *********************************************** */
       /* The input file contains four pair on a line,    */
	   /* so after four pairs are read, go to the next    */
	   /* line.                                           */
       /* *********************************************** */
      if (count==4) {
		  fscanf(ind,"\n");
		  count=0;
	  }
       /* *********************************************** */
       /* If the input file doesn't have four bonds on    */
	   /* the last line, skip a line in the end.          */
       /* *********************************************** */
	  else if (i== haccN[k]-1) fscanf(ind, "\n");

    }

  /* ----------------------------------------------- */
  /* Read in the information about each fixed        */
  /* exlusion.                                       */
  /* ----------------------------------------------- */    

	int *nnb_list;
   

       /* *********************************************** */
       /* Read in the number of fixed exlusions.          */
       /* *********************************************** */
   fscanf(ind,"%d\n",&exNBN[k]);fgets(tt,80,ind);
       /* *********************************************** */
       /* Variables needed.                               */
       /* *********************************************** */
   int rk = exNBN[k]+1;
       /* *********************************************** */
       /* Allocate memory for this number of fixed        */
	   /* exclusions.                                     */
       /* *********************************************** */
   exNB[k] = (struct NBexclusions*) calloc(rk+MEXTRA,sizeof(struct NBexclusions));
		if(exNB[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for exNB[k]\n"); exit(11);}
    nnb_list = (int*) calloc(rk,sizeof(int));
		if(nnb_list == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for nnb_list\n"); exit(11);}
    count=0; t1=0;
	for(int i=0; i<exNBN[k]; i++) {
       /* *********************************************** */
       /* Read in the list of the NNB's                   */
       /* *********************************************** */
      fscanf(ind,"%d", &t1); nnb_list[i] = t1-1;
	  count++;	  

       /* *********************************************** */
       /* The input file contains 8 atoms on a line.  Go  */
	   /* to the next line after 8 are read.              */
       /* *********************************************** */
      if (count==8) {
		  fscanf(ind,"\n");
		  count=0;
	  }
       /* *********************************************** */
       /* If the input file doesn't have 8 atoms on       */
	   /* the last line, skip a line in the end.          */
       /* *********************************************** */
	  else if (i== exNBN[k]-1) fscanf(ind, "\n");

    }
       /* *********************************************** */
       /* Now, loop over the number of sites and get the  */
	   /* number in the list just read that are excluded  */
	   /* with each site.                                 */
       /* *********************************************** */
	int lower_bound=0;
	int upper_bound=0;
	count = 0;
	for(int i=0; i<box[k].boxns; i++){
		fscanf(ind,"%d", &upper_bound); 
		for(int j=lower_bound; j<upper_bound; j++){
			exNB[k][count].a = i;
			exNB[k][count].b = nnb_list[j];
			count++;
		}
		lower_bound = upper_bound;
	}

	free(nnb_list);

  /* ----------------------------------------------- */
  /* Determine and list the total number of 1-4      */
  /* interactions.                                   */
  /* ----------------------------------------------- */    
       /* *********************************************** */
       /* First, loop around all bends and compare        */
       /* *********************************************** */
	int Nsolute = 0;
        for(int i=0; i<molec.Nsolute; i++) Nsolute += mol[i].Nsite;
	int rh;	
	if(torsN[k]>10)  rh = 100*Nsolute+100;
	else rh = 100;
	int *a14;
	int *b14;
	int *c14;
	int *d14;

			/* Since we don't know how many 1-4 interactions there are */
	        /* a priori, these will be temproray variables holding the */
	        /* atoms in involved and they will be freed after in14[k]  */
	        /* is assigned.                                            */

   a14 = (int*) calloc(rh,sizeof(int));
		if(a14 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for a\n"); exit(11);}
   b14 = (int*) calloc(rh,sizeof(int));
		if(b14 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for b\n"); exit(11);}
   c14 = (int*) calloc(rh,sizeof(int));
		if(c14 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for c\n"); exit(11);}
   d14 = (int*) calloc(rh,sizeof(int));
		if(d14 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for d\n"); exit(11);}

	
	count = 0;
	int duplicated_14=0;
	for(int i=0; i<bendN[k]-1; i++){
		for(int j=i+1; j<bendN[k]; j++){
			if(bend[k][i].c == bend[k][j].b){
				if(bend[k][i].b == bend[k][j].a){
					duplicated_14 = duplicate14_check(k,bend[k][i].a,bend[k][j].c,a14,d14,count);
					if(duplicated_14 == 0){
						a14[count] = bend[k][i].a;
						b14[count] = bend[k][i].b;
						c14[count] = bend[k][i].c;
						d14[count] = bend[k][j].c;
						count++;
					}
				}else if(bend[k][i].b == bend[k][j].c){
					duplicated_14 = duplicate14_check(k,bend[k][i].a,bend[k][j].a,a14,d14,count);
					if(duplicated_14 == 0){
						a14[count] = bend[k][i].a;
						b14[count] = bend[k][i].b;
						c14[count] = bend[k][i].c;
						d14[count] = bend[k][j].a;
						count++;
					}
				}

			}else if(bend[k][i].a == bend[k][j].b){
				if(bend[k][i].b == bend[k][j].c){
					duplicated_14 = duplicate14_check(k,bend[k][j].a,bend[k][i].c,a14,d14,count);
					if(duplicated_14 == 0){
						a14[count] = bend[k][j].a;
						b14[count] = bend[k][j].b;
						c14[count] = bend[k][j].c;
						d14[count] = bend[k][i].c;
						count++;
					}
				}else if(bend[k][i].b == bend[k][j].a){
					duplicated_14 = duplicate14_check(k,bend[k][j].c,bend[k][i].c,a14,d14,count);
					if(duplicated_14 == 0){
						a14[count] = bend[k][j].c;
						b14[count] = bend[k][j].b;
						c14[count] = bend[k][j].a;
						d14[count] = bend[k][i].c;
						count++;
					}
				}
			}

		}//for j
	}//for i
	if(count > rh){ fprintf(stdout, "ERROR: Guess (%i) for number of 1-4 interactions (%i) was not big enough!!!", rh, count); exit(11);}
	in14N[k]=count;
    rh = in14N[k] + MEXTRA;
       /* *********************************************** */
       /* Allocate memory for this number of 1-4          */
	   /* interactions.                                   */
       /* *********************************************** */
    in14[k] = (struct interactions1_4*) calloc(rh,sizeof(struct interactions1_4));
		if(in14[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for in14[k]\n"); exit(11);}
       /* *********************************************** */
       /* Now, assign the 1-4 interactions determined     */
	   /* above to the in14[k] structure.                 */
       /* *********************************************** */
		for(int i=0; i<in14N[k]; i++){
			in14[k][i].a = a14[i];
			in14[k][i].b = b14[i];
			in14[k][i].c = c14[i];
			in14[k][i].d = d14[i];
		}
       /* *********************************************** */
       /* Free the temporary variables.                   */
       /* *********************************************** */
		free(a14);
		free(b14);
		free(c14);
		free(d14);



	fclose(ind);

  /* ----------------------------------------------- */
  /* Allocate memory for the xint variable that      */
  /* contains the excluded interactions information. */
  /* Then make the list.                             */
  /* ----------------------------------------------- */    
#ifdef STYPE
#ifndef DNA_GOLIK
#ifndef S14
	int rg = bondN[k] + bendN[k] + exNBN[k];
#endif
#ifdef S14
	int rg = bondN[k] + bendN[k] + exNBN[k] + in14N[k];
#endif

	xintN[k] = rg;
	xint[k] = (struct excinter*) calloc(rg+MEXTRA, sizeof(struct excinter));
		if(xint[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for xint[k]\n"); exit(11);}
	count=0;

	for(int i=0; i<bondN[k]; i++){
		xint[k][count].a = bond[k][i].a;
		xint[k][count].b = bond[k][i].b;
		xint[k][count].xflag = 1;
		count++;
	}

	for(int i=0; i<bendN[k]; i++){
		xint[k][count].a = bend[k][i].a;
		xint[k][count].b = bend[k][i].c;
		xint[k][count].xflag = 2;
		count++;	
	}

	for(int i=0; i<exNBN[k]; i++){
		xint[k][count].a = exNB[k][i].a;
		xint[k][count].b = exNB[k][i].b;
		xint[k][count].xflag = 3;
		count++;	
	}
#ifdef S14
	for(int i=0; i<in14N[k]; i++){
		xint[k][count].a = in14[k][i].a;
		xint[k][count].b = in14[k][i].d;
		xint[k][count].xflag = 4;
		count++;	
	}
#endif //S14
#endif//DNA_GOLIK
#endif //STYPE
  /* ----------------------------------------------- */
  /* Assign the Lennard Jones parameters to each     */
  /* atom/site.  The parameters were read in         */
  /* read atom into the epstmp and sigtmp arrays.    */
  /* ----------------------------------------------- */    
       /* *********************************************** */
       /* Assign memory for LJ paramters.                 */
       /* *********************************************** */
    int ra = box[k].boxns + MEXTRA;
    pott[k] = (struct inter*) calloc(ra,sizeof(struct inter));
      if(pott[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for pott[k]\n"); exit(11);}
#ifdef S14
    pott14[k] = (struct inter*) calloc(ra,sizeof(struct inter));
      if(pott14[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for pott[k]\n"); exit(11);}
#endif
	for(int i=0; i<box[k].boxns; i++) {
       /* *********************************************** */
       /* Assign the potential id numbers to each site i. */
       /* *********************************************** */
      //pott[k][i].atomid = atnopbc[k][i].atomid;
      #ifdef S14
      //pott[k][i].atomid = atnopbc[k][i].atomid;
      #endif
    for (int mm =0; mm<NUM_ATOM;mm++){
       /* *********************************************** */
       /* Compare the id number to those read in          */
       /* read_atom (atom_type array). If they are the    */
       /* same, assign the correct mass to the site.      */
       /* *********************************************** */
      if(atom_type[mm].atomid == atnopbc[k][i].atomid){
        pott[k][i].mas = atom_type[mm].mass; 
        #ifdef S14
        pott14[k][i].mas = atom_type[mm].mass; 
        #endif
        break;
      }
    }

       /* *********************************************** */
       /* Assign the correct eps, sig, and qq to each     */
       /* site.  Note: The index of epstmp and sigtmp are */
       /* over id numbers, not site numbers.              */
       /* *********************************************** */
      pott[k][i].eps  = -4.184* epstmp[0][atnopbc[k][i].atomid];
      pott[k][i].sig  = 2*sigtmp[0][atnopbc[k][i].atomid];
      pott[k][i].qq   = atnopbc[k][i].q;

#ifdef BEAD
      pott[k][i].eps  = -2.494353;  // kJ/mol (gives T* = 1 if T = 300K)
      pott[k][i].sig  = 4.0; //angstoms
#endif

#ifdef S14
      pott14[k][i].eps  = -4.184* epstmp[1][atnopbc[k][i].atomid];
      pott14[k][i].sig  = 2*sigtmp[1][atnopbc[k][i].atomid];
      pott14[k][i].qq   = atnopbc[k][i].q;
#endif
      

//#endif
    }

       /* *********************************************** */
       /* Allocate memory for the ljset structure to be   */
	   /* used in setlj().                                */
	   /* *********************************************** */
#ifndef STYPE
    int rf = box[k].boxns + MEXTRA;   
    ljset[k].pot = (struct inter**) calloc(rf,sizeof(struct inter*));
		if(ljset[k].pot == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ljset[k].pot\n"); exit(11);}
	for(int n=0; n<rf; n++){ 
      ljset[k].pot[n] = (struct inter*) calloc(rf,sizeof(struct inter));
		if(ljset[k].pot[n] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ljset[k].pot[n]\n"); exit(11);}
	}
#endif

#ifdef STYPE
	int rf = NUM_LJSET;   
	ljset = (struct inter**) calloc(rf,sizeof(struct inter*));
		if(ljset == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ljset\n"); exit(11);}
#ifdef S14
	ljset14 = (struct inter**) calloc(rf,sizeof(struct inter*));
		if(ljset14 == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ljset14\n"); exit(11);}
#endif
	for(int n=0; n<rf; n++){ 
      ljset[n] = (struct inter*) calloc(rf,sizeof(struct inter));
	  if(ljset[n] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ljset[n]\n"); exit(11);}
#ifdef S14
     ljset14[n] = (struct inter*) calloc(rf,sizeof(struct inter));
	  if(ljset14[n] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for ljset14[n]\n"); exit(11);}
#endif
	}
#endif

  /* ----------------------------------------------- */
  /* Create the Urey-Bradley list of 1-3             */
  /* interactions.                                   */
  /*                                                 */
  /* The Urey-Bradley list is created from searching */
  /* the bends for ones that contain Urey-Bradley    */
  /* parameters.  Since we don't know a priori how   */
  /* many there are, first I will create a temporary */
  /* structure with NUM_BEND.  I will then search    */
  /* through each bend and compare it to the list of */
  /* Urey-Bradley terms.  When the match, I will add */
  /* them to the temporary urey list.  After the     */
  /* final number is known, I will create the urey[k]*/
  /* array used in simulation.                       */
  /* ----------------------------------------------- */
#ifndef NEUTRAL
	struct{
		int n;
		int *a;
		int *b;
		double *k;
		double *Seq;
	}urey_temp;

	urey_temp.a	= (int*)calloc(bendN[k]+1,sizeof(int));
		if(urey_temp.a == NULL){ fprintf(stdout,"ERROR: cannot allocate memory for urey_type.a\n"); exit(11);}
	urey_temp.b	= (int*)calloc(bendN[k]+1,sizeof(int));
		if(urey_temp.b == NULL){ fprintf(stdout,"ERROR: cannot allocate memory for urey_type.b\n"); exit(11);}
	urey_temp.k		= (double*)calloc(bendN[k]+1,sizeof(double));
		if(urey_temp.k== NULL){ fprintf(stdout,"ERROR: cannot allocate memory for urey_temp.k\n"); exit(11);}
	urey_temp.Seq	= (double*)calloc(bendN[k]+1,sizeof(double));
		if(urey_temp.Seq == NULL){ fprintf(stdout,"ERROR: cannot allocate memory for urey_temp.Seq\n"); exit(11);}

	urey_temp.n = 0;
		
	for(int i=0; i<bendN[k];i++){	
	  for (int nu =0; nu< NUM_UREY; nu++) {
		  if ((urey_prop[nu].type1 == atnopbc[k][bend[k][i].a].atomid && urey_prop[nu].type3 == atnopbc[k][bend[k][i].c].atomid && 
			   urey_prop[nu].type2 == atnopbc[k][bend[k][i].b].atomid) || 
			  (urey_prop[nu].type1 == atnopbc[k][bend[k][i].c].atomid && urey_prop[nu].type3 == atnopbc[k][bend[k][i].a].atomid && 
			   urey_prop[nu].type2 == atnopbc[k][bend[k][i].b].atomid)) {
			  urey_temp.a[urey_temp.n]		= bend[k][i].a;
			  urey_temp.b[urey_temp.n]		= bend[k][i].c;
			  urey_temp.Seq[urey_temp.n]	= urey_prop[nu].Seq;
			  urey_temp.k[urey_temp.n]		= urey_prop[nu].k;
			  urey_temp.n++;
			  break; 
		  }
	  }
	}

	ureyN[k] = urey_temp.n;

    urey[k] = (struct ureys*) calloc(urey_temp.n+1,sizeof(struct ureys));
		if(urey[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for urey[k]\n"); exit(11);}

		for(int i=0; i<urey_temp.n; i++){
			urey[k][i].a = urey_temp.a[i];
			urey[k][i].b = urey_temp.b[i];
			urey[k][i].Seq = urey_temp.Seq[i];
			urey[k][i].k = urey_temp.k[i] * 4.184;
		}

  free(urey_temp.Seq);
  free(urey_temp.a);
  free(urey_temp.b);
  free(urey_temp.k);
#endif




}

/* ======================================================================== */
/* duplicate14_check.cpp                                                    */
/*                                                                          */
/*		This subroutine determines if the 1-4 interaction found above by    */
/* searching through the bends is already part of a bond or bend or if it   */
/* was specifically exluded as a fixed NB exclusion.  All these usually     */
/* happen in ring structures.  If the 1-4 interaction in question is one of */
/* the above elements, the subroutine returns 1, else, it returns 0.        */
/*                                                                          */
/* ======================================================================== */

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
int duplicate14_check(int ibox, int a, int d, int *a14, int *d14, int count )
{
	int k=ibox;
  /* ----------------------------------------------- */
  /* Loop over the bonds.                            */
  /* ----------------------------------------------- */    
	for(int i=0; i<bondN[k]; i++){
		if( (bond[k][i].a == a && bond[k][i].b == d)||
			(bond[k][i].b == a && bond[k][i].a == d))
				return(1);
	}

  /* ----------------------------------------------- */
  /* Loop over the bends.                            */
  /* ----------------------------------------------- */   
	for(int i=0; i<bendN[k]; i++){
		if( (bend[k][i].a == a && bend[k][i].c == d)||
			(bend[k][i].c == a && bend[k][i].a == d))
				return(1);
	}

  /* ----------------------------------------------- */
  /* Loop over the fixed NB exclusions.              */
  /* ----------------------------------------------- */   
	for(int i=0; i<exNBN[k]; i++){
		if( (exNB[k][i].a == a && exNB[k][i].b == d)||
			(exNB[k][i].b == a && exNB[k][i].a == d))
				return(1);
	}
  /* ----------------------------------------------- */
  /* Remove any double counts for 1-4	             */
  /* ----------------------------------------------- */   
	for(int i=0; i<count; i++){
		if( (a14[i] == a && d14[i] == d)||
			(d14[i] == a && a14[i] == d))
				return(1);
	}
	return(0);

}
