/* ======================================================================== */
/* read_atom.cpp                                                            */
/*                                                                          */
/*   This subroutine reads the atom types and parameters from the           */
/* database.  It reads the information contained in simul19.param if        */
/* NEUTRAL is not defined and simul19_sasa.param if NEUTRAL is defined.     */
/* This subroutine gets the information about the atom/site id number,      */
/* name, and molecular weight/mass.  It then reads in the bond information  */
/* for all possible bonds, followed by the bending, torision, and LJ        */
/* information.                                                             */  
/* 									    */	
/*  	                                                                    */
/* Note:  This subroutine basically just reads in the CHARMM force field    */
/* parameters. All of them are not necessarily used in a simulations. The   */
/* atoms, bonds, bends, torsions, and LJ parameters used are read in        */
/* read_topology and setlj.  In these subroutines, the each atom/structure  */
/* present in the simulation is compared to those in this list to get the   */
/* appropriate value.                                                       */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void read_atom (void) 
{
  /* ================================================================== */
  /*                                                                    */
  /* Variables, flags, and pointers needed in the subroutine.           */
  /*                                                                    */
  /* ================================================================== */
  int qw=0; //set equal to zero just to prevent a warning.
  char tmp1[5], tmp2[5], tmp3[5],tmp4[5]; 
  int flag1 =0; 
  int flag2 =0; 
  int flag3 =0; 
  int flag4 =0;
  double tmp5;
  char tt[150];

  FILE *fptr1;
  /* ================================================================== */
  /*                                                                    */
  /* Set the file pointers to the correct input file.                   */
  /*                                                                    */
  /* ================================================================== */

  /* ----------------------------------------------- */
  /* Different files are needed for the SASA         */
  /* implicit solvent model.                         */
  /* ----------------------------------------------- */    
//#ifdef DNA_GOLIK
//  sprintf(fname,"./INPUT/coarse_grain.param");
//#else
//#ifndef NEUTRAL
//#ifndef SUG
//  sprintf(fname,"./INPUT/simul29_prot_na_lipid.param");
//#endif
//#ifdef SUG
//	sprintf(fname,"./INPUT/simul22_prot_na_lip_sug.param");
//#endif
//#endif
//#ifdef NEUTRAL
//  sprintf(fname,"./INPUT/simul19_sasa.param");
//#endif
//#endif //DNA_GOLIK
  char fname[100];
  sprintf(fname,"./INPUT/%s",param_file_name);
  if( NULL == (fptr1=fopen(fname,"r")) ) {
      fprintf(stdout,"ERROR: configuration file %s does not exist!!!\n",param_file_name);
      exit(1);
  }

  /* ================================================================== */
  /*                                                                    */
  /* Start reading the file information.                                */
  /*                                                                    */
  /* ================================================================== */
	fgets(tt,150,fptr1);//skip the title line

  /* ----------------------------------------------- */
  /* Read in the atom/site information.              */
  /* ----------------------------------------------- */
  fscanf(fptr1, "%d",   &NUM_ATOM); fgets(tt,150,fptr1);
  atom_type = (struct atomtype*)calloc(NUM_ATOM,sizeof(struct atomtype));
  if(atom_type == NULL){fprintf(stdout, "ERROR: cannot allocate memory for atom_type\n"); exit(11);}

	for (int i=0; i<NUM_ATOM; i++) {
		fscanf(fptr1, "%d",   &atom_type[i].atomid);
		fscanf(fptr1, "%4s",  atom_type[i].name);
		fscanf(fptr1, "%lf\n", &atom_type[i].mass);
	}

  	//fgets(tt,150,fptr1);
  /* ----------------------------------------------- */
  /* Read in the bond parameters.                    */
  /* ----------------------------------------------- */  
  fscanf(fptr1, "%d",   &NUM_BOND); fgets(tt,150,fptr1);
  bond_prop = (struct bondprop*)calloc(NUM_BOND,sizeof(struct bondprop));
  if(bond_prop == NULL){fprintf(stdout, "ERROR: cannot allocate memory for bond_prop\n"); exit(11);}

       /* *********************************************** */
       /* The first loop controls the reading of the      */
	   /* bonding atoms and parameters.                   */
       /* *********************************************** */
	for (int i=0; i<NUM_BOND; i++) {
		flag1=flag2=0;

       /* *********************************************** */
       /* Read in the name of the two bonding sites.      */
       /* *********************************************** */
		fscanf(fptr1, "%4s",   tmp1);
		fscanf(fptr1, "%4s",   tmp2);
		
	   /* *********************************************** */
       /* This loop compares the names of the two bonding */
	   /* sites with those read in above to assign the    */
	   /* correct atomid number.                          */
       /* *********************************************** */
		for (int ii =0; ii<NUM_ATOM; ii++) {
			if (strcmp(tmp1,atom_type[ii].name)==0 ) {
				bond_prop[i].type1=atom_type[ii].atomid;
				flag1 =1;
			}
			if (strcmp(tmp2,atom_type[ii].name)==0 ) {
				bond_prop[i].type2=atom_type[ii].atomid;
				flag2=1;
			}
			if ((flag1*flag2)!=0) break;
		}
       /* *********************************************** */
       /* Read in the harmonic constant and equilibrium   */
	   /* bond distance.                                  */
       /* *********************************************** */
		fscanf(fptr1, "%lf",   &bond_prop[i].kbond);
		fscanf(fptr1, "%lf\n", &bond_prop[i].eqbond);
	}

	//fgets(tt, 150, fptr1);


  /* ----------------------------------------------- */
  /* Read in bending(angle) parameters.              */
  /* ----------------------------------------------- */
  fscanf(fptr1, "%d",   &NUM_BEND); fgets(tt,150,fptr1);
  bend_prop = (struct bendprop*)calloc(NUM_BEND,sizeof(struct bendprop));
  if(bend_prop == NULL){fprintf(stdout, "ERROR: cannot allocate memory for bend_prop\n"); exit(11);}  
       /* *********************************************** */
       /* The first loop controls the reading of the      */
	   /* bending atoms and parameters.                   */
       /* *********************************************** */
	for (int i=0; i<NUM_BEND; i++) {
		flag1=flag2= flag3=0;

       /* *********************************************** */
       /* Read in the name of the three sites of the      */
	   /* bend.                                           */
       /* *********************************************** */
		fscanf(fptr1, "%4s",   tmp1);
		fscanf(fptr1, "%4s",   tmp2);
		fscanf(fptr1, "%4s",   tmp3);

	   /* *********************************************** */
       /* This loop compares the names of the three       */
	   /* sites involved in the bending with those read   */
	   /* in above to assigns the correct atomid number.  */
       /* *********************************************** */
		for (int ii =0; ii<NUM_ATOM; ii++) {
			if (strcmp(tmp1,atom_type[ii].name)==0 ) {
				bend_prop[i].type1=atom_type[ii].atomid;
				flag1 =1;
			}
			if (strcmp(tmp2,atom_type[ii].name)==0 ) {
				bend_prop[i].type2=atom_type[ii].atomid;
				flag2=1;
			}
			if (strcmp(tmp3,atom_type[ii].name)==0 ) {
				bend_prop[i].type3=atom_type[ii].atomid;
				flag3=1;
			}
			if ((flag1*flag2*flag3)!=0) break;
		}
       /* *********************************************** */
       /* Read in the harmonic constant and equilibrium   */
	   /* bend angle.                                     */
       /* *********************************************** */
		fscanf(fptr1, "%lf",   &bend_prop[i].kbend);
		fscanf(fptr1, "%lf\n", &bend_prop[i].eqbend);
	}

	//fgets(tt, 150, fptr1);


  /* ----------------------------------------------- */
  /* Read in the dihedral parameters.                */
  /* ----------------------------------------------- */
  fscanf(fptr1, "%d",   &NUM_TORS); fgets(tt,150,fptr1);
  dihedral_prop = (struct dihedralprop*)calloc(NUM_TORS,sizeof(struct dihedralprop));
  if(dihedral_prop == NULL){fprintf(stdout, "ERROR: cannot allocate memory for dihedral_prop\n"); exit(11);}  
       /* *********************************************** */
       /* This first loop controls the reading of the     */
	   /* atoms/sites involved in the dihedral angle      */
	   /* and the associated parameters.                  */
       /* *********************************************** */
	for (int i=0; i<NUM_TORS; i++) {
		flag1=flag2= flag3=flag4=0;

       /* *********************************************** */
       /* Read in the names of the four sites of the      */
	   /* dihedral.                                       */
       /* *********************************************** */
		fscanf(fptr1, "%4s",   tmp1);
		fscanf(fptr1, "%4s",   tmp2);
		fscanf(fptr1, "%4s",   tmp3);
		fscanf(fptr1, "%4s",   tmp4);

	   /* *********************************************** */
       /* This loop compares the names of the four        */
	   /* sites involved in the dihedral with those read  */
	   /* in above to assigns the correct atomid number.  */
       /* *********************************************** */
		for (int ii =0; ii<NUM_ATOM; ii++) {
			if (strcmp(tmp1,atom_type[ii].name)==0 ) {
				dihedral_prop[i].type1=atom_type[ii].atomid;
				flag1 =1;
			}
			if (strcmp(tmp2,atom_type[ii].name)==0 ) {
				dihedral_prop[i].type2=atom_type[ii].atomid;
				flag2=1;
			}
			if (strcmp(tmp3,atom_type[ii].name)==0 ) {
				dihedral_prop[i].type3=atom_type[ii].atomid;
				flag3=1;
			}
			if (strcmp(tmp4,atom_type[ii].name)==0 ) {
				dihedral_prop[i].type4=atom_type[ii].atomid;
				flag4=1;
			}
			if ((flag1*flag2*flag3*flag4)!=0) break;
		}

       /* *********************************************** */
       /* Read in the dihedral force field parameters.    */
       /* *********************************************** */
		fscanf(fptr1, "%lf",   &dihedral_prop[i].kphi);
		fscanf(fptr1, "%lf",   &dihedral_prop[i].nphi);
		fscanf(fptr1, "%lf\n", &dihedral_prop[i].delphi);
	}
	fgets(tt,150,fptr1);

  /* ----------------------------------------------- */
  /* Read in the Lennard Jones parameters.           */
  /* ----------------------------------------------- */    
       /* *********************************************** */
       /* This first loop controls the reading of the     */
	   /* atoms/sites.                                    */
       /* *********************************************** */
	for (int i=0; i<NUM_ATOM-1+1; i++) {

       /* *********************************************** */
       /* Read in the name of the site.                   */
       /* *********************************************** */
		fscanf(fptr1, "%4s",   tmp1);
	   /* *********************************************** */
       /* This loop compares the name of the site with    */
	   /* the names read in above to assign the proper    */
	   /* atomid number.                                  */
       /* *********************************************** */
		for (int ii =0; ii<NUM_ATOM; ii++) {
			if (strcmp(tmp1,atom_type[ii].name)==0 ) {
				qw=atom_type[ii].atomid;
				break;
			}
		}
       /* *********************************************** */
       /* Read in the LJ parameters.                      */
       /* *********************************************** */
		fscanf(fptr1, "%lf",   &epstmp[0][qw]);
		fscanf(fptr1, "%lf\n", &sigtmp[0][qw]);
	}
	//fgets(tt,150,fptr1);
  /* ----------------------------------------------- */
  /* Read in the improper torsion parameters.        */
  /* ----------------------------------------------- */
  fscanf(fptr1, "%d",   &NUM_IMPR); fgets(tt,150,fptr1);
  improper_prop = (struct improperprop*)calloc(NUM_IMPR,sizeof(struct improperprop));
  if(improper_prop == NULL){fprintf(stdout, "ERROR: cannot allocate memory for improper_prop\n"); exit(11);}  

       /* *********************************************** */
       /* This first loop controls the reading of the     */
	   /* atoms/sites involved in the improper angle      */
	   /* and the associated parameters.                  */
       /* *********************************************** */
	for (int i=0; i<NUM_IMPR; i++) {
		flag1=flag2= flag3=flag4=0;

       /* *********************************************** */
       /* Read in the names of the four sites of the      */
	   /* improper dihedral.                              */
       /* *********************************************** */
		fscanf(fptr1, "%4s",   tmp1);
		fscanf(fptr1, "%4s",   tmp2);
		fscanf(fptr1, "%4s",   tmp3);
		fscanf(fptr1, "%4s",   tmp4);

	   /* *********************************************** */
       /* This loop compares the names of the four        */
	   /* sites involved in the dihedral with those read  */
	   /* in above to assigns the correct atomid number.  */
       /* *********************************************** */
		for (int ii =0; ii<NUM_ATOM; ii++) {
			if (strcmp(tmp1,atom_type[ii].name)==0 ) {
				improper_prop[i].type1=atom_type[ii].atomid;
				flag1 =1;
			}
			if (strcmp(tmp2,atom_type[ii].name)==0 ) {
				improper_prop[i].type2=atom_type[ii].atomid;
				flag2=1;
			}
			if (strcmp(tmp3,atom_type[ii].name)==0 ) {
				improper_prop[i].type3=atom_type[ii].atomid;
				flag3=1;
			}
			if (strcmp(tmp4,atom_type[ii].name)==0 ) {
				improper_prop[i].type4=atom_type[ii].atomid;
				flag4=1;
			}
			if ((flag1*flag2*flag3*flag4)!=0) break;
		}

       /* *********************************************** */
       /* Read in the dihedral force field parameters.    */
       /* *********************************************** */
		fscanf(fptr1, "%lf",   &improper_prop[i].kimpr);
		fscanf(fptr1, "%lf",   &tmp5);
		fscanf(fptr1, "%lf\n", &improper_prop[i].angeq);
	}
	fgets(tt,150,fptr1);
  /* ----------------------------------------------- */
  /* Read in the Lennard Jones special 1-4           */
  /* parameters.                                     */
  /* ----------------------------------------------- */    
       /* *********************************************** */
       /* This first loop controls the reading of the     */
	   /* atoms/sites.                                    */
       /* *********************************************** */
	for (int i=0; i<NUM_ATOM-1; i++) {

       /* *********************************************** */
       /* Read in the name of the site.                   */
       /* *********************************************** */
		fscanf(fptr1, "%4s",   tmp1);
	   /* *********************************************** */
       /* This loop compares the name of the site with    */
	   /* the names read in above to assign the proper    */
	   /* atomid number.                                  */
       /* *********************************************** */
		for (int ii =0; ii<NUM_ATOM; ii++) {
			if (strcmp(tmp1,atom_type[ii].name)==0 ) {
				qw=atom_type[ii].atomid;
				break;
			}
		}
       /* *********************************************** */
       /* Read in the LJ parameters.                      */
       /* *********************************************** */
		fscanf(fptr1, "%lf",   &epstmp[1][qw]);
		fscanf(fptr1, "%lf\n", &sigtmp[1][qw]);
	}
	//fgets(tt,150,fptr1);

  /* ----------------------------------------------- */
  /* Read in the Urey-Bradley 1-3 parameters.        */
  /* ----------------------------------------------- */    
#ifndef NEUTRAL
  fscanf(fptr1, "%d",   &NUM_UREY); fgets(tt,150,fptr1);
  urey_prop = (struct ureyprop*)calloc(NUM_UREY,sizeof(struct ureyprop));
  if(urey_prop == NULL){fprintf(stdout, "ERROR: cannot allocate memory for urey_prop\n"); exit(11);}  

  for (int i=0; i<NUM_UREY; i++) {
		flag1=flag2=flag3=0;

       /* *********************************************** */
       /* Read in the name of the three sites of the      */
	   /* bend.                                           */
       /* *********************************************** */
		fscanf(fptr1, "%4s",   tmp1);
		fscanf(fptr1, "%4s",   tmp2);
		fscanf(fptr1, "%4s",   tmp3);

	   /* *********************************************** */
       /* This loop compares the names of the     two     */
	   /* sites involved in the Urey-Bradley 1-3 distance */
	   /* with those read in above to assigns the correct */
	   /* atomid number.                                  */
       /* *********************************************** */
		for (int ii =0; ii<NUM_ATOM; ii++) {
			if (strcmp(tmp1,atom_type[ii].name)==0 ) {
				urey_prop[i].type1=atom_type[ii].atomid;
				flag1 =1;
			}
			if (strcmp(tmp2,atom_type[ii].name)==0 ) {
				urey_prop[i].type2=atom_type[ii].atomid;
				flag2=1;
			}
			if (strcmp(tmp3,atom_type[ii].name)==0 ) {
				urey_prop[i].type3=atom_type[ii].atomid;
				flag3=1;
			}

			if ((flag1*flag2*flag3)!=0) break;
		}
       /* *********************************************** */
       /* Read in the harmonic constant and equilibrium   */
	   /* 1-3 distance.                                   */
       /* *********************************************** */
		fscanf(fptr1, "%lf",   &urey_prop[i].k);
		fscanf(fptr1, "%lf\n", &urey_prop[i].Seq);
	}

	//fgets(tt, 150, fptr1);
#endif

  /* ================================================== */
  /* Read in the NBFIX if they are present.             */
  /* ================================================== */
  char key[7] = "NBFIX\n";
  char key2[6] = "NBFIX";
  char test_str[150];
  int ReturnCode;
 
  ReturnCode = fscanf(fptr1, "%s", test_str);  fgets(tt, 150, fptr1);
 
  /* This while loops reads the file until the keywork NBFIX is found. */
  while(ReturnCode != EOF){
    if(strcmp(test_str,key)==0 || strcmp(test_str,key2)==0)  break;
    ReturnCode = fscanf(fptr1, "%s", test_str);  fgets(tt, 150, fptr1);;
  }
  
  /* Now, read the NBFIX                                               */
  if(ReturnCode == EOF)
    NUM_NBFIX = 0;
  else{
//    nbfix = (struct nbfixes**) calloc(NUM_ATOM,sizeof(struct nbfixes*));
//      if(nbfix == NULL){fprintf(stdout, "ERROR: cannot allocate memory for nbfix\n"); exit(11);}
//    for(int ii=0; ii<NUM_ATOM; ii++){
//      nbfix[ii] = (struct nbfixes*) calloc(NUM_ATOM,sizeof(struct nbfixes));
//        if(nbfix == NULL){fprintf(stdout, "ERROR: cannot allocate memory for nbfix[ii]\n"); exit(11);}
//    }    
    fscanf(fptr1, "%d",   &NUM_NBFIX); fgets(tt,150,fptr1);  
    nbfix_prop = (struct nbfixprop*)calloc(NUM_NBFIX,sizeof(struct nbfixprop));
      if(nbfix_prop == NULL){fprintf(stdout, "ERROR: cannot allocate memory for nbfix_prop\n"); exit(11);}
    for (int i=0; i<NUM_NBFIX; i++) {
      flag1=flag2=0;

      /* *********************************************** */
      /* Read in the name of the two types of the        */
      /* nbfix.                                          */
      /* *********************************************** */
      fscanf(fptr1, "%4s",   tmp1);
      fscanf(fptr1, "%4s",   tmp2);

      /* *********************************************** */
      /* This loop compares the names of the two         */
      /* sites involved in the NBFIX                     */
      /* with those read in above to assigns the correct */
      /* atomid number.                                  */
      /* *********************************************** */
      for (int ii =0; ii<NUM_ATOM; ii++){
        if(strcmp(tmp1,atom_type[ii].name)==0 ) {
          nbfix_prop[i].type1=atom_type[ii].atomid;
          flag1 =1;
        }
        if(strcmp(tmp2,atom_type[ii].name)==0 ) {
          nbfix_prop[i].type2=atom_type[ii].atomid;
          flag2=1;
        }
        if ((flag1*flag2*flag3)!=0) break;
      }
      
      /* *********************************************** */
      /* Read in the eps and sigma.                      */
      /* *********************************************** */
      fscanf(fptr1, "%lf",   &nbfix_prop[i].eps);
      fscanf(fptr1, "%lf\n", &nbfix_prop[i].sig);
      
//      #ifdef STYPE
      /* *********************************************** */
      /* Now, create a list indexed by atom_type rather  */
      /* than the nbfix_prop which is indexed by         */
      /* NBFIX number.                                   */
      /* *********************************************** */      
//      nbfix[nbfix_prop[i].type1][nbfix_prop[i].type2].eps = nbfix_prop[i].eps;
//      nbfix[nbfix_prop[i].type2][nbfix_prop[i].type1].eps = nbfix_prop[i].eps;
//      nbfix[nbfix_prop[i].type1][nbfix_prop[i].type2].eps = nbfix_prop[i].sig;
//      nbfix[nbfix_prop[i].type2][nbfix_prop[i].type1].eps = nbfix_prop[i].sig;      
//      #endif
      
    }

    //fgets(tt, 150, fptr1); 
  }

fclose(fptr1);
}

