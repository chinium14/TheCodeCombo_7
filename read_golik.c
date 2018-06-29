#ifdef GOLIK
#include "defines.h"
#ifndef DNA_GOLIK
/* ======================================================================== */
/* read_golik.cpp	                                                        */
/* This suboutine reads the information about the Go-like model of			*/
/* Hoang and Cieplak(JCP Vol. 112 Number 15 Page 6851 (2000)).              */
/*																		    */
/* Written by Thomas Knotts 8 Mar 04                                        */
/* Cosolutes added by Nitin Rathore											*/
/* ======================================================================== */
void set_go_params(int ibox);
int is_bond(int,int,int);
void read_golik(void){
  double hydro_index[47]={0.0,-0.4,1.8,-0.8,2.5,4.2,-0.7,4.5,1.6,1.9,-3.5,-3.5,3.8,-3.9,-3.5,
    -3.5,-4.5,-3.2,2.8,-1.3,-0.9,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
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
  double factor = pow(2.0,-(1.0/6.0));

  for(int k=0; k<sim.NB; k++){    
	  fscanf(file_ptr,"%d",&boxnum);			fgets(tt,150,file_ptr);
	  if(boxnum != k){ fprintf(stdout, "ERROR: Box numbers in golik.input (%i) are not correct! (should be %i)\n", boxnum, k); exit(010);}

	  fscanf(file_ptr,"%d %lf %lf",&golik[k].Psite,&golik[k].eps, &golik[k].mass);		fgets(tt,150,file_ptr);
	  fscanf(file_ptr,"%lf",&golik[k].req);		fgets(tt,150,file_ptr);
	  fscanf(file_ptr,"%lf",&golik[k].d_native);	fgets(tt,150,file_ptr);
	  fscanf(file_ptr,"%lf",&golik[k].d_repulsive);	fgets(tt,150,file_ptr);
//	  fscanf(file_ptr,"%d", &golik[k].N);		fgets(tt,150,file_ptr);
	  golik[k].contact  = (int**)calloc(golik[k].Psite,sizeof(int*));

	  if(golik[k].contact == NULL){fprintf(stdout, "ERROR: cannot allocate memory for golik[k].contact\n"); exit(11);}

	    if(golik[k].contact[i] == NULL){fprintf(stdout, "ERROR: cannot allocate memory for golik[k].contact[i]\n"); exit(11);}
	  }	

	  golik[k].eps *= 4.184; //switch from kcal/mol to kJ/mol


	  /* ====================================================== */
	  /* Calculate the equilibrium bond distances, bend angles, */
	  /* torsion angles, and sigmas.                            */
	  /* ====================================================== */
		set_go_params(k);
	  
	  /* ====================================================== */
	  /* Assign the masses to the mass in golik_input if given. */
	  /* ====================================================== */
	  if(golik[k].mass > 0.0001){
	  	for(int a=0; a< golik[k].Psite; a++){
		  	pott[k][a].mas	=	golik[k].mass;
	  	}
	  }
  
	 
	  /* ====================================================== */
	  /* Set all of the hydropathy parameters.                  */
	  /* ====================================================== */
	  golik[k].hydro_index  = (double*)calloc(box[k].boxns,sizeof(double));
		if(golik[k].hydro_index == NULL){fprintf(stdout, "ERROR: cannot allocate memory for golik[k].hydro_index\n"); exit(11);}    
	  int rescount = 0;
	  int atcount = 0;
	  for(int m=0; m<sim.NC; m++) {
		for(int i=0; i<bp[k][m].nbox; i++) {
		  for (int kk=0; kk<mol[m].Nres; kk++) {			  
			for(int j=0; j<residue[k][rescount].Nsite; j++) {
		  golik[k].hydro_index[atcount] = hydro_index[residue[k][rescount].type];
		  atcount++;
		}//for j
		rescount++;
		  }//for kk
		}//for i
	  }//for m


	  /* ====================================================== */
	  /* Search through the list and find 1-2 exclusions.       */
	  /* ====================================================== */

	  for(int a=0; a<golik[k].Psite; a++) {
	    for(int b=a+1; b<golik[k].Psite; b++) {
	      for (int n=0; n<bondN[k]; n++) {
			if ((a == bond[k][n].a && b== bond[k][n].b) || (a == bond[k][n].b && b== bond[k][n].a)) {
			  ljset[k].pot[a][b].eps = 0.0;
			  ljset[k].pot[a][b].sig = 0.0;
			  golik[k].contact[a][b] = 0;
			  break;
			}
		  }
		}
	  }


	  /* ====================================================== */
	  /* Now, if d_repulsive was input as zero, calculate the   */
      /* average native contact length and set d_repulsive      */
      /* equal to that.                                         */
	  /* ====================================================== */
	  if(golik[k].d_repulsive == 0.0){
		  double d_sum=0.0;
		  int counter=0;
		  for(int a=0; a<golik[k].Psite; a++){
			for(int b=a+1; b<golik[k].Psite; b++){
			  if(golik[k].contact[a][b] == 1){
				d_sum = d_sum + ljset[k].pot[a][b].sig/factor;
				counter ++;
			  }
			}
		  }

		  double d_repulsive = d_sum/counter;
		  golik[k].d_repulsive = d_repulsive;
		  for(int a=0; a<golik[k].Psite; a++){
			for(int b=a+1; b<golik[k].Psite; b++){
			  if(golik[k].contact[a][b] == 0){
				ljset[k].pot[a][b].sig = factor*golik[k].d_repulsive;
			  }
			}
		  }
	  }

		  

  }//k
fclose(file_ptr);
	/*---------------------------------------------	*/
	/* Now see if cosolutes are present or not and	*/
	/* if they are, allocate memory and assign		*/
	/* values to epsilon and sigma					*/
	/* --------------------------------------------	*/

  sprintf(file_name,"./INPUT/cosol.input");
  if( NULL != (file_ptr=fopen(file_name,"r")) ) {
	  fgets(tt,150,file_ptr);
	  fgets(tt,150,file_ptr);
	  fscanf(file_ptr,"%lf",&gosol.lambda);		fgets(tt,150,file_ptr);
	  fscanf(file_ptr,"%d", &gosol.N);			fgets(tt,150,file_ptr);

	  gosol.eps  = (double*)calloc(gosol.N,sizeof(double));
		if(gosol.eps == NULL){fprintf(stdout, "ERROR: cannot allocate memory for gosol.eps\n"); exit(11);}
	  gosol.sig  = (double*)calloc(gosol.N,sizeof(double));
		if(gosol.sig == NULL){fprintf(stdout, "ERROR: cannot allocate memory for gosol.sig\n"); exit(11);}
	  gosol.mass  = (double*)calloc(gosol.N,sizeof(double));
		if(gosol.mass == NULL){fprintf(stdout, "ERROR: cannot allocate memory for gosol.mass\n"); exit(11);}
	  gosol.Nsol  = (int*)calloc(gosol.N,sizeof(int));
		if(gosol.Nsol == NULL){fprintf(stdout, "ERROR: cannot allocate memory for gosol.Nsol\n"); exit(11);}

	  for (int i=0; i<gosol.N; i++){
		  fscanf(file_ptr,"%d",&gosol.Nsol[i]);
		  fscanf(file_ptr,"%lf",&gosol.eps[i]);
		  fscanf(file_ptr,"%lf",&gosol.sig[i]);
		  fscanf(file_ptr,"%lf",&gosol.mass[i]);fgets(tt,150,file_ptr);
		  gosol.eps[i] *= 4.184; //switch from kcal/mol to kJ/mol
	  }
	  fclose(file_ptr);
   
	   for(int k=0; k<sim.NB; k++){
		   int count=golik[k].Psite;
		   for(int j=0; j<gosol.N; j++){
			   for(int a=count; a< count+ gosol.Nsol[j]; a++){
				  pott[k][a].eps	=	gosol.eps[j];
				  pott[k][a].sig	=	gosol.sig[j];
				  pott[k][a].mas	=	gosol.mass[j];
			  }//a
			   count += gosol.Nsol[j] ;
		   }//j
	   }//k

		/* --------------------------------------------	*/
		/* Now assign LJ parameters to the PS and SS	*/
		/* interactions. Note PP interactions have been	*/
		/* already defined.								*/
		/* --------------------------------------------	*/
	   for(int k=0; k<sim.NB; k++){
		   for(int a=0; a<box[k].boxns-1; a++){
			  for(int b=golik[k].Psite; b<box[k].boxns; b++){
				  if(a<golik[k].Psite) { // PS interactions
					  ljset[k].pot[a][b].eps = sqrt(golik[k].eps*pott[k][b].eps); 
					  ljset[k].pot[a][b].sig = 0.5*(factor*golik[k].d_repulsive+pott[k][b].sig);
				  }
				  else{		// SS interactions
					  ljset[k].pot[a][b].eps = sqrt(pott[k][a].eps*pott[k][b].eps); 
					  ljset[k].pot[a][b].sig = 0.5*(pott[k][a].sig+pott[k][b].sig);
				  }
			  }//b
			}//a
	   }//k
  }// if cosol file exists
}

void set_go_params(int ibox){
	
	int k=ibox;
	char tt[80];
	double factor = pow(2.0,-(1.0/6.0));


	/* ====================================================== */
	/* Allocate posx,y,z and read in ref.pdb info.            */
	/* ====================================================== */
	struct pos_struct{
		double x;
		double y;
		double z;
	} *pos;

	pos  = (struct pos_struct*)calloc(box[k].boxns,sizeof(struct pos_struct));
		if(pos == NULL){fprintf(stdout, "ERROR: cannot allocate memory for pos\n"); exit(11);}


	FILE *fp;
	char temp[5];
	//char temp1[5]="END";
	//char temp2[5]="REMA";
	char temp1[5]="ATOM";
    char temp3[10],temp4[10];
	int itemp1,itemp2;

	if( NULL == (fp=fopen("./INPUT/ref.pdb","r")) ) {			// open reference file
		fprintf(stdout,"input file ref.pdb does not exist!!!\n");
		exit(10);
	} 
	
	int atom_num = 0;  
	while (!feof(fp)){
		fscanf(fp, "%4s",   temp);
		//if (strcmp(temp,temp1)==0 || strcmp(temp,temp2)==0 ) {
		if(strcmp(temp,temp1)!=0){
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
		fscanf(fp,"%lf",  &pos[atom_num].x);
		fscanf(fp,"%lf",  &pos[atom_num].y);
		fscanf(fp,"%lf",  &pos[atom_num].z);
		fgets(tt,80,fp);
		atom_num++;
	}

	fclose(fp);



	  /* ====================================================== */
	  /* Search for the native contacts and set lj parameters.  */
	  /* ====================================================== */  

	for(int a=0; a<golik[k].Psite-1; a++){
		for(int b=a+1; b<golik[k].Psite; b++){
			double dx = pos[a].x - pos[b].x;
			double dy = pos[a].y - pos[b].y;
			double dz = pos[a].z - pos[b].z;
			double dr = sqrt(dx*dx + dy*dy + dz*dz);
			ljset[k].pot[a][b].eps = golik[k].eps;
			if(dr<golik[k].d_native){  //if the distance is less than golick.contact, it is a native contact.
			 ljset[k].pot[a][b].sig = factor*dr;
			 golik[k].contact[a][b] = 1;
			}  
			else{
				ljset[k].pot[a][b].sig = factor*golik[k].d_repulsive;
				golik[k].contact[a][b] = 0;
			}
		}//b
	}//a


	/* ====================================================== */
	/* Set all of the bond parameters to the correct values.  */
	/* ====================================================== */
	  for(int i=0; i<bondN[k]; i++){
		if(golik[k].req > 0.0001) bond[k][i].req = golik[k].req;
		else{
			double dx = pos[bond[k][i].a].x - pos[bond[k][i].b].x;
			double dy = pos[bond[k][i].a].y - pos[bond[k][i].b].y;
			double dz = pos[bond[k][i].a].z - pos[bond[k][i].b].z;
			double dr = sqrt(dx*dx + dy*dy + dz*dz);
			bond[k][i].req = dr;	
		}
	  }

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
			double dax = pos[ia].x - pos[ib].x;
			double day = pos[ia].y - pos[ib].y;
			double daz = pos[ia].z - pos[ib].z;
			double dcx = pos[ic].x - pos[ib].x;
			double dcy = pos[ic].y - pos[ib].y;
			double dcz = pos[ic].z - pos[ib].z;

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

		for(int i=0; i<torsN[k]; i++) tors[k][i].delphi[0] = 0.0;

#endif //GOBT

}
#endif
/* ======================================================================== */
/* Begin is_bond()                                                          */
/* ======================================================================== */
int is_bond(int k, int a, int b){
  for (int n=0; n<bondN[k]; n++) {
	  if((a == bond[k][n].a && b== bond[k][n].b) || (a == bond[k][n].b && b== bond[k][n].a))
      return(1);
  }
  return(0);
} 
#endif
