#ifdef DNA_GOLIK
/* ======================================================================== */
/* read_dnagolik.cpp	                                                      */
/* This suboutine reads the information about the coarse-grain dna model.   */
/*																		                                      */
/* Written by Thomas Knotts 19 Mar 05                                       */
/*                                                    											*/
/* ======================================================================== */
#include "defines.h"
void set_dnago_params(int);
int is_bond(int, int,int);
int is_bend(int,int,int);
void construct_dna_hash(int, double*, double*, double*);
int exclusion_type(int, double, int, int, int, int, int,int,double,double);
int interaction_type(int,int);

/* ======================================================================== */
/* Begin read_golik()                                                       */
/* ======================================================================== */
void read_dnagolik(void){
  golik    = (struct golik_struct*) calloc(sim.NB,sizeof(struct golik_struct));
		if(golik == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for golik\n"); exit(11);}
  int boxnum;

  char tt[150],file_name[50];
  FILE *file_ptr;

  /* ====================================================== */
  /* Open the golik.input file and read in parameters.      */
  /* ====================================================== */
  sprintf(file_name,"./INPUT/golik.input");
  if( NULL == (file_ptr=fopen(file_name,"r")) ) {
	fprintf(stdout,"ERROR: Input file %s does not exist!!!\n",file_name);
	exit(1);
  }

  fgets(tt,150,file_ptr);
  fgets(tt,150,file_ptr);
  //double factor = pow(2.0,-(1.0/6.0));

  for(int k=0; k<sim.NB; k++){    
	  fscanf(file_ptr,"%d",&boxnum);			fgets(tt,150,file_ptr);
	  if(boxnum != k){ fprintf(stdout, "ERROR: Box numbers in golik.input (%i) are not correct! (should be %i) \n", boxnum, k); exit(010);}

	  fscanf(file_ptr,"%d %lf %lf",&golik[k].Psite,&golik[k].eps, &golik[k].mass);		fgets(tt,150,file_ptr);
	  fscanf(file_ptr,"%lf",&golik[k].req);		fgets(tt,150,file_ptr);
	  fscanf(file_ptr,"%lf",&golik[k].d_native);	fgets(tt,150,file_ptr);
	  fscanf(file_ptr,"%lf",&golik[k].d_repulsive);	fgets(tt,150,file_ptr);
	  golik[k].eps *= 4.184; //switch from kcal/mol to kJ/mol

	  /* ====================================================== */
	  /* Calculate the equilibrium bond distances, bend angles, */
	  /* torsion angles, and sigmas.  Also determine the type   */
    /* of non-bonded interactions.                            */
	  /* ====================================================== */
		set_dnago_params(k);
	 

  }//k

  fclose(file_ptr);
}



/* ======================================================================== */
/* Begin set_go_params()                                                    */
/* ======================================================================== */
void set_dnago_params(int ibox){
	
	int k=ibox;

	char tt[80];
	//double factor = pow(2.0,-(1.0/6.0));


	/* ====================================================== */
	/* Allocate posx,y,z and read in ref.pdb info.            */
	/* ====================================================== */
//	struct pos_struct{
//		double x;
//		double y;
//		double z;
//	} *pos;
  double *rx;
  double *ry;
  double *rz;

//	pos  = (struct pos_struct*)calloc(box[k].boxns,sizeof(struct pos_struct));
//		if(pos == NULL){fprintf(stdout, "ERROR: cannot allocate memory for pos\n"); exit(11);}
	rx  = (double*)calloc(box[k].boxns,sizeof(double));
		if(rx == NULL){fprintf(stdout, "ERROR: cannot allocate memory for rx\n"); exit(11);}
	ry  = (double*)calloc(box[k].boxns,sizeof(double));
		if(ry == NULL){fprintf(stdout, "ERROR: cannot allocate memory for ry\n"); exit(11);}
	rz  = (double*)calloc(box[k].boxns,sizeof(double));
		if(rz == NULL){fprintf(stdout, "ERROR: cannot allocate memory for rz\n"); exit(11);}


	FILE *fp;
	char temp[7];
	//char temp1[5]="END";
	//char temp2[5]="REMA";
	char temp1[5]="ATOM";
  char temp5[7]="HETATM";
    char temp3[10],temp4[10];
	int itemp1,itemp2;

	if( NULL == (fp=fopen("./INPUT/ref.pdb","r")) ) {			// open reference file
		fprintf(stdout,"input file ref.pdb does not exist!!!\n");
		exit(10);
	} 
	
	int atom_num = 0;  
	while (!feof(fp)){
		fscanf(fp, "%6s",   temp);
		//if (strcmp(temp,temp1)==0 || strcmp(temp,temp2)==0 ) {
		if(strcmp(temp,temp1)!=0 && strcmp(temp,temp5)!=0){
            fgets(tt,80,fp);
			continue;
		}
		fscanf(fp, "%d",   &itemp1);
		if(itemp1!=atom_num+1){
        	fprintf(stdout,"ERROR: ref.pdb site numbers wrong. atom_num = %d itemp1 = %d\n",atom_num,itemp1); 
			exit(10);
		}
		fscanf(fp,"%4s",  temp3);
		fscanf(fp,"%3s",  temp4);
		fscanf(fp,"%d",   &itemp2);
		fscanf(fp,"%lf",  &rx[atom_num]);
		fscanf(fp,"%lf",  &ry[atom_num]);
		fscanf(fp,"%lf",  &rz[atom_num]);
		fgets(tt,80,fp);
		atom_num++;
	}

	fclose(fp);




	/* ====================================================== */
	/* Set all of the bond parameters to the correct values.  */
	/* ====================================================== */
	  for(int i=0; i<bondN[k]; i++){
		  if(golik[k].req > 0.0001) bond[k][i].req = golik[k].req;
		  else{
			  double dx = rx[bond[k][i].a] - rx[bond[k][i].b];
			  double dy = ry[bond[k][i].a] - ry[bond[k][i].b];
			  double dz = rz[bond[k][i].a] - rz[bond[k][i].b];
			  double dr = sqrt(dx*dx + dy*dy + dz*dz);
			  bond[k][i].req = dr;	
		  }
	  }
  construct_dna_hash(k,rx,ry,rz);
#ifdef GOBT
  	/* ====================================================== */
	  /* Set all of the bend parameters to the correct values.  */
	  /* ====================================================== */
	  for(int i=0; i<bendN[k]; i++){
		/* ----------------------------------------- */
		/* Assign the atoms involved in the bending. */
		/* ----------------------------------------- */
			int ia = bend[k][i].a;
			int ib = bend[k][i].b;
			int ic = bend[k][i].c;

		  /* ----------------------------------------- */
		  /* Calculate the distances between the atoms */
		  /* involved in the bending.                  */
		  /* ----------------------------------------- */
			double dax = rx[ia] - rx[ib];
			double day = ry[ia] - ry[ib];
			double daz = rz[ia] - rz[ib];
			double dcx = rx[ic] - rx[ib];
			double dcy = ry[ic] - ry[ib];
			double dcz = rz[ic] - rz[ib];

		  /* ----------------------------------------- */
		  /* Apply minimum image convention.           */
		  /* ----------------------------------------- */
			double lx = box[k].lx;
			double ly = box[k].ly;
			double lz = box[k].lz;
			double hx = box[k].hx;
			double hy = box[k].hy;
			double hz = box[k].hz;

			if(dax >  hx) dax -= lx;
			else if(dax < -hx) dax += lx;
			if(day >  hy) day -= ly;
			else if(day < -hy) day += ly;
			if(daz >  hz) daz -= lz;
			else if(daz < -hz) daz += lz;


			if(dcx >  hx) dcx -= lx;
			else if(dcx < -hx) dcx += lx;
			if(dcy >  hy) dcy -= ly;
			else if(dcy < -hy) dcy += ly;
			if(dcz >  hz) dcz -= lz;
			else if(dcz < -hz) dcz += lz;

		  /* ----------------------------------------- */
		  /* Calculate the angle of the bend and then  */
		  /* the energy using the CHARMM force field.  */
		  /* ----------------------------------------- */
			double drai = 1.0 / sqrt(dax*dax + day*day + daz*daz);
			double erax = dax * drai; 
			double eray = day * drai; 
			double eraz = daz * drai; 
			double drci = 1.0 / sqrt(dcx*dcx + dcy*dcy + dcz*dcz);
			double ercx = dcx * drci; 
			double ercy = dcy * drci; 
			double ercz = dcz * drci; 

			double cosphi = erax*ercx + eray*ercy + eraz*ercz;
			bend[k][i].angeq = acos(cosphi);
	  }

	  /* ====================================================== */
	  /* Set all of the dihedral parameters to the correct		*/
	  /* values.                                                */
	  /* ====================================================== */

		for(int i=0; i<torsN[k]; i++) {
		  /* ----------------------------------------- */
		  /* Assign the atoms involved in the          */
		  /* dihedral angle.                           */
		  /* ----------------------------------------- */
			int ia = tors[k][i].a;
			int ib = tors[k][i].b;
			int ic = tors[k][i].c;
			int id = tors[k][i].d;

		  /* ----------------------------------------- */
		  /* Calculate the distances between the atoms */
		  /* involved in the dihedral angle.           */
		  /* ----------------------------------------- */
			double lx = box[k].lx;
			double ly = box[k].ly;
			double lz = box[k].lz;
			double hx = box[k].hx;
			double hy = box[k].hy;
			double hz = box[k].hz;

			double dabx = rx[ib] - rx[ia];
			double daby = ry[ib] - ry[ia];
			double dabz = rz[ib] - rz[ia];
			double dbcx = rx[ic] - rx[ib];
			double dbcy = ry[ic] - ry[ib];
			double dbcz = rz[ic] - rz[ib];
			double dcdx = rx[id] - rx[ic];
			double dcdy = ry[id] - ry[ic];
			double dcdz = rz[id] - rz[ic];

		  /* ----------------------------------------- */
		  /* Apply minimum image convention.           */
		  /* ----------------------------------------- */
			if(dabx >  hx) dabx -= lx;
			else if(dabx < -hx) dabx += lx;
			if(daby >  hy) daby -= ly;
			else if(daby < -hy) daby += ly;
			if(dabz >  hz) dabz -= lz;
			else if(dabz < -hz) dabz += lz;

			if(dbcx >  hx) dbcx -= lx;
			else if(dbcx < -hx) dbcx += lx;
			if(dbcy >  hy) dbcy -= ly;
			else if(dbcy < -hy) dbcy += ly;
			if(dbcz >  hz) dbcz -= lz;
			else if(dbcz < -hz) dbcz += lz;
    
			if(dcdx >  hx) dcdx -= lx;
			else if(dcdx < -hx) dcdx += lx;
			if(dcdy >  hy) dcdy -= ly;
			else if(dcdy < -hy) dcdy += ly;
			if(dcdz >  hz) dcdz -= lz;
			else if(dcdz < -hz) dcdz += lz;

		  /* ----------------------------------------- */
		  /* Calculate the bond lengths and unit       */
		  /* vectors involved in the dihedral angle.   */
		  /* ----------------------------------------- */
			double rab   = sqrt(dabx*dabx + daby*daby + dabz*dabz);
			double irab  = 1.0 / rab;
			double erabx = dabx * irab;
			double eraby = daby * irab;
			double erabz = dabz * irab;
			double rbc   = sqrt(dbcx*dbcx + dbcy*dbcy + dbcz*dbcz);
			double irbc  = 1.0 / rbc;
			double erbcx = dbcx * irbc;
			double erbcy = dbcy * irbc;
			double erbcz = dbcz * irbc;
			double rcd   = sqrt(dcdx*dcdx + dcdy*dcdy + dcdz*dcdz);
			double ircd  = 1.0 / rcd;
			double ercdx = dcdx * ircd;
			double ercdy = dcdy * ircd;
			double ercdz = dcdz * ircd;

		  /* ----------------------------------------- */
		  /* Calculate the cross and dot products      */
		  /* between unit vectors and the bond angles. */
		  /* ----------------------------------------- */
			double abbcx = eraby * erbcz - erabz * erbcy;
			double abbcy = erabz * erbcx - erabx * erbcz;
			double abbcz = erabx * erbcy - eraby * erbcx;
			double cosb  = -(erabx*erbcx + eraby*erbcy + erabz*erbcz);
			double isinb2 = 1.0 / (1.0 - cosb*cosb);
			double isinb  = sqrt(isinb2);

			double dccbx = ercdy * erbcz - ercdz * erbcy;
			double dccby = ercdz * erbcx - ercdx * erbcz;
			double dccbz = ercdx * erbcy - ercdy * erbcx;
			double cosc  = -(ercdx*erbcx + ercdy*erbcy + ercdz*erbcz);
			double isinc2 = 1.0 / (1.0 - cosc*cosc);
			double isinc  = sqrt(isinc2);


		  /* ----------------------------------------- */
		  /* Calculate the torsion/dihedral angle.     */
		  /* ----------------------------------------- */
			double abcdx = -(abbcy * dccbz - abbcz * dccby);
			double abcdy = -(abbcz * dccbx - abbcx * dccbz);
			double abcdz = -(abbcx * dccby - abbcy * dccbx);    
			double num = (abcdx*erbcx + abcdy*erbcy + abcdz*erbcz);
			double signum; 
			if (num >= 0.0) signum = 1.0;
			else signum = -1.0;
			double costau = -(abbcx*dccbx + abbcy*dccby + abbcz*dccbz) * isinb * isinc;
			if (costau > 1.0) costau = 1.0;
			if (costau < -1.0) costau = -1.0;
			tors[k][i].delphi[0] = signum * acos(costau) + PI; 
		}
#endif //GOBT
  free(rx);
  free(ry);
  free(rz);
}


/* ======================================================================== */
/* Begin construct_dna_hash()                                               */
/* ======================================================================== */

/* =============================================== */
/* This subroutine resets the hash table and to    */
/* delineate between different types of non-bonded */
/* interactions.                                   */
/* =============================================== */
void construct_dna_hash(int ibox, double *rx, double *ry, double *rz){
  int k=ibox;
  int ns = box[k].boxns;
  double factor = pow(2.0,-(1.0/6.0));
  
  /* --------------------------------------------- */
  /* Since I am using the same table variable that */
  /* is already allocated but contains useless     */
  /* information, I will reset the counters first. */
  /* --------------------------------------------- */
  for(int i=0; i<(NPRIME+1); i++) table[k][i].n=0;

  /* --------------------------------------------- */
  /* First loop through every pair and decide      */
  /* which type of interaction it is.  Since the   */
  /* number is not known a priori, allocate temp   */
  /* variables and after the number is known,      */
  /* assign the data to the xintN variable and     */
  /* free the temp variables.                      */
  /* --------------------------------------------- */
  #define TMPSIZE 68000
  struct nb_struct{
    int a;
    int b;
    int flag;
    double sig;
  }*nb;
  nb = (struct nb_struct*)calloc(TMPSIZE,sizeof(struct nb_struct));
  if(nb == 0 || nb == NULL){
    fprintf(stdout,"Cannot allocate nb in read_dnagolik.\n");
    exit(15);
  }
  int nbN=0;
  //goncN=0;
  double d_sum=0.0;
	int counter=0;
	
  /* --------------------------------------------- */
  /* Determine the maximum bond and bend distance  */
  /* to reduce search time.                        */
  /* --------------------------------------------- */
  double max_bond_dr = 0;
  double max_bend_dr = 0;
  for(int i=0; i<bondN[k]; i++){
    if(max_bond_dr < bond[k][i].req) max_bond_dr = bond[k][i].req;
  }
  for(int i=0; i<bendN[k]; i++){
    double dx = atom[k][bend[k][i].a].x - atom[k][bend[k][i].c].x;
    double dy = atom[k][bend[k][i].a].y - atom[k][bend[k][i].c].y;
    double dz = atom[k][bend[k][i].a].z - atom[k][bend[k][i].c].z;
    double dr = sqrt(dx*dx + dy*dy + dz*dz);
    if(max_bend_dr < dr) max_bend_dr = dr;
  }
  
  max_bond_dr = max_bond_dr + 3.0;//Add on a cushion.
  max_bend_dr = max_bend_dr + 3.0;	

  if(bendN[k]==0) max_bend_dr = golik[k].d_native+1.0;

  /* --------------------------------------------- */
  /* Flags                                         */
  /* 0 = non-native                                */
  /* 1 = native                                    */
  /* 2 = AT complimentary base pair                */
  /* 3 = CG complimentary base pair                */
  /* 4 = mismatch base pair                        */
  /* --------------------------------------------- */

  for(int a=0; a<ns-1; a++){
    for(int b=a+1; b<ns; b++){
      if(nbN >= TMPSIZE){fprintf(stdout,"Increase TMPSIZE in read_dnagolik.cpp. a=%d\n",a);exit(12);}
			double dx = rx[a] - rx[b];
			double dy = ry[a] - ry[b];
			double dz = rz[a] - rz[b];
			double dr = sqrt(dx*dx + dy*dy + dz*dz);
      int strand_a = atnopbc[k][a].molec;
      int strand_b = atnopbc[k][b].molec;

      int switch_flag = exclusion_type(k,dr,a,b,strand_a,strand_b,atnopbc[k][a].atomid,atnopbc[k][b].atomid,max_bond_dr, max_bend_dr);
      switch(switch_flag){
        /* ------------------------------------------------ */
        /* If the pair is a native contact.                 */
        /* ------------------------------------------------ */
        case 1:
          nb[nbN].a = a;
          nb[nbN].b = b;
          nb[nbN].sig = dr * factor;
          nb[nbN].flag = 1;
          nbN++;
          if(nbN > TMPSIZE){fprintf(stdout,"Increase TMPSIZE in read_dnagolik.cpp. a=%d\n",a);exit(12);}
          /* Keep track of average distance to calulate d_repulsive */
          d_sum = d_sum + dr;
				  counter ++;
          break;
        
        /* ------------------------------------------------ */
        /* If the pair is an AT complimentary base pair.    */
        /* ------------------------------------------------ */
/*        case 2:
          nb[nbN].a = a;
          nb[nbN].b = b;
          nb[nbN].flag = 2; 
          nbN++;
          if(nbN > TMPSIZE){fprintf(stdout,"Increase TMPSIZE in read_dnagolik.cpp. a=%d\n",a);exit(12);}
          break;*/
        /* ------------------------------------------------ */
        /* If the pair is a GC complimentary base pair.     */
        /* ------------------------------------------------ */
/*        case 3:
          nb[nbN].a = a;
          nb[nbN].b = b;
          nb[nbN].flag = 3;
          nbN++;
          if(nbN > TMPSIZE){fprintf(stdout,"Increase TMPSIZE in read_dnagolik.cpp. a=%d\n",a);exit(12);}
          break;*/
        /* ------------------------------------------------ */
        /* If the pair is a mismatched base pair.           */
        /* ------------------------------------------------ */
/*        case 4:
          nb[nbN].a = a;
          nb[nbN].b = b;
          nb[nbN].flag = 4;
          nbN++;
          if(nbN > TMPSIZE){fprintf(stdout,"Increase TMPSIZE in read_dnagolik.cpp. a=%d\n",a);exit(12);}
          break;*/
        /* ------------------------------------------------ */
        /* If the pair is bonded.                           */
        /* ------------------------------------------------ */
        case 5:
          nb[nbN].a = a;
          nb[nbN].b = b;
          nb[nbN].flag = 5;
          nbN++;
          if(nbN > TMPSIZE){fprintf(stdout,"Increase TMPSIZE in read_dnagolik.cpp. a=%d\n",a);exit(12);}
          break;
        /* ------------------------------------------------ */
        /* If the pair is a native contact and a bend.      */
        /* ------------------------------------------------ */
         case 6:
          nb[nbN].a = a;
          nb[nbN].b = b;
          nb[nbN].sig = dr * factor;
          nb[nbN].flag = 6;
          nbN++;
          if(nbN > TMPSIZE){fprintf(stdout,"Increase TMPSIZE in read_dnagolik.cpp. a=%d\n",a);exit(12);}
          /* Keep track of average distance to calulate d_repulsive */
          d_sum = d_sum + dr;
				  counter ++;
          break;

        /* ------------------------------------------------ */
        /* If the pair is an inter-strand sugar interaction */
        /* ------------------------------------------------ */
        #ifdef DRSAM
         case 8: 
          nb[nbN].a = a;
          nb[nbN].b = b;
          nb[nbN].flag = 8;
          nbN++;  
          if(nbN > TMPSIZE){fprintf(stdout,"Increase TMPSIZE in read_dnagolik.cpp. a=%d\n",a);exit(12);}
          break;  
        #endif

        /* ------------------------------------------------ */
        /* If the pair is a non-native contact.            */
        /* ------------------------------------------------ */
        default:
          break;
      }//switch

		}//b
	}//a

  	/* ====================================================== */
	  /* Now, if d_repulsive was input as zero, calculate the   */
    /* average native contact length and set d_repulsive      */
    /* equal to that.                                         */
	  /* ====================================================== */
    if(golik[k].d_repulsive == 0.0) golik[k].d_repulsive = d_sum/counter;

  	/* ====================================================== */
	  /* Now that d_repulsive is calculated, set the sigma of   */
    /* the sites equal to 2^(-1/6)*d_repulsive.  Then, reset  */
    /* ljset because setlj() has already been called.         */
    /* The variable, ljset, will then be the correct          */
    /* non-native sigmas which are used as the default pair   */
    /* identity in cnbnd_dnagolik.                            */
    /*                                                        */
    /* Note: If the param file contained the d_repulive value */
    /* for the LJ sigma, then you wouldn't have to do the     */
    /* next few lines.  However, that would require           */
    /* calculating the d_repulive by some outside script and  */
    /* setting up a different param file for each protein.    */
	  /* ====================================================== */
    if(k==0){
      for(int a=0; a<NUM_ATOM; a++){
        double sig = factor*golik[k].d_repulsive;
        int ta = atom_type[a].atomid;
        if(ta != 190 && ta != 195) sigtmp[0][ta] = sig/2.0;
        else sigtmp[0][ta] *= factor;
      }
    }

    for(int a=0; a<NUM_ATOM; a++) {
      for(int b=0; b<NUM_ATOM; b++) {
        double sig = 0.5 * (2.0 * sigtmp[0][atom_type[a].atomid] + 2.0 * sigtmp[0][atom_type[b].atomid]);//the 2.0 comes from the fact that CHARMM lists sig/2 
        ljset[atom_type[a].atomid][atom_type[b].atomid].sig = sig;
        ljset[atom_type[a].atomid][atom_type[b].atomid].rc2_non = sig*sig/factor/factor;
	    }
    }

  	/* ====================================================== */
	  /* Now that you know the interactions that are bonds,     */
    /* bends, native, bp, etc., put the infor into the xint   */
    /* variable.                                              */
	  /* ====================================================== */
    xintN[k] = nbN;
	  xint[k] = (struct excinter*) calloc(nbN+MEXTRA, sizeof(struct excinter));
		if(xint[k] == NULL){ fprintf(stdout, "ERROR: cannot allocate memory for xint[k]\n"); exit(11);}

    for(int i=0; i<xintN[k]; i++){
      xint[k][i].a = nb[i].a;
      xint[k][i].b = nb[i].b;
      xint[k][i].xflag = nb[i].flag;
      if(nb[i].flag==1 || nb[i].flag == 6) xint[k][i].sig = nb[i].sig; 
    }
  //}

  /* ================================================================== */
  /* Set the parameters for the types of iteractions.                   */
  /* ================================================================== */
/*  double epsilon = golik[k].eps;
  double eps_AT = epsilon*2.0/3.0*4.0;
  double eps_GC = epsilon*4.0;
  double eps_mm = epsilon;
  double factor = pow(2.0,-(1.0/6.0));
  double sig_AT = 2.9002;
  double sig_GC = 2.8694;
  double d_mm   = 1.0;
  double sig_non= factor*golik[k].d_repulsive;
  double sig_mm = factor*d_mm;
*/
  /* ================================================================== */
  /* Loop around number of site types interactions and apply the mixing */
  /* rules.                                                             */
  /* ================================================================== */
  for(int a=0; a<NUM_ATOM; a++){
    for(int b=0; b<NUM_ATOM; b++){
      int ta = atom_type[a].atomid;
      int tb = atom_type[b].atomid;
      int type = interaction_type(ta,tb);

/*      switch(type){
        case 2:
          sig = sig_AT;
          eps = eps_AT;
        case 3:
          sig = sig_GC;
          eps = eps_GC;
        case 4:
          sig = sig_mm;
          eps = eps_mm;
        default:
          sig = sig_non;
          eps = eps_non;
      }

      ljset[ta][tb].sig = sig;
      ljset[ta][tb].eps = eps;*/
      ljset[ta][tb].id = type;
	  }
  }

  /* ================================================================== */
  /* Set up hash table of excluded/native interactions.                 */
  /* ================================================================== */
 int maxn=0;

		for(int i=0; i<xintN[k]; i++){
			int a = xint[k][i].a;
			int b = xint[k][i].b;
			unsigned long int index = (a*b) % NPRIME;
			if(table[k][index].n == 0){
				table[k][index].a[0] = a;
				table[k][index].b[0] = b;
				table[k][index].p[0] = i;
			}
			else{
				int qa = table[k][index].n;
				table[k][index].a[qa] = a;
				table[k][index].b[qa] = b;
				table[k][index].p[qa]  = i;
			}
			table[k][index].n++;
			if( table[k][index].n>maxn)
      {
        maxn = table[k][index].n;
        if(maxn > hashsize){
          fprintf(stdout,"The hash table isn't big enough. Increase 'hashsize' in read_golik.cpp. i=%d, xintN[k]=%d\n",i,xintN[k]);
          exit(12);
        }
      }
		}//for i
  free(nb);
}

/* ======================================================================== */
/* Begin is_compliment()                                                    */
/* ======================================================================== */
int exclusion_type(int k, double dr, int a, int b, int strand_a, int strand_b,
int type_a, int type_b, double max_bond_dr, double max_bend_dr){
  /* ---------------------------------------------------------- */
  /* If the pair is bonded, there is no non-bonded interaction. */
  /* ---------------------------------------------------------- */
  if(dr < max_bond_dr){
    if(is_bond(k,a,b)!=0) return(5);
  }
  
  //{
    if((dr<golik[k].d_native) && (strand_a == strand_b))
    {
      if(dr < max_bend_dr){
        if(is_bend(k,a,b)==0) return(1);
      }
      return(6);
    }
    /* If the pair is an AT, return 1 */
    if((type_a == 3 && type_b ==5) || (type_a == 5 && type_b ==3))
      return(2);
    /* If the pair is an GC, return 2 */
    else if((type_a == 4 && type_b ==6) || (type_a == 6 && type_b ==4))
      return(3);
    else if((type_a >=3 && type_a <=6) && (type_b >=3 && type_b <=6))
      return(4);
    /* Tf the pair is a interstrand sugar */
        #ifdef DRSAM
    else if((type_a == 2 && type_b == 2) &&  (strand_a == strand_b))
      return(8);
        #endif
    else return(0);
  //}
}

int interaction_type(int type_a, int type_b){
    if((type_a == 3 && type_b ==5) || (type_a == 5 && type_b ==3)) //AT
      return(2);
    else if((type_a == 4 && type_b ==6) || (type_a == 6 && type_b ==4)) //CG
      return(3);
    else if((type_a >=3 && type_a <=6) && (type_b >=3 && type_b <=6)) //MM
      return(4);
       #ifdef DRSAM
    else if((type_a == 2 && type_b == 2) )  //SS
      return(7);
        #endif
    else return(0); //NON
}

/* ======================================================================== */
/* Begin is_bend()                                                          */
/* ======================================================================== */
int is_bend(int k, int a, int b){
  for (int n=0; n<bendN[k]; n++) {
	  if((a == bend[k][n].a && b== bend[k][n].c) || (a == bend[k][n].c && b== bend[k][n].a))
      return(1);
  }
  return(0);
} 
#endif

