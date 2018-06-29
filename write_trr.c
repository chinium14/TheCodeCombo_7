#ifdef TRR
/*	-----------------------------------------------------------------------	*/
/* write_trr.cpp															*/
/*		written by Rob Riggleman											*/
/* Based on Trio. This subroutine Writes positions, velocities, and forces	*/
/* to a traj.trr file which can be read and	analyzed using Gromacs			*/
/*  -----------------------------------------------------------------------	*/


#include "defines.h"

int reverse=1;									//set to 0 to turn off byte swapping
int precision=8;								//set to 4 to use float precision
FILE *das;
int write_int(long);
int write_float(float fp);
int write_bstring(char *str);
int write_vector(float *fp);
int strip_white(char *str);
void swap4(void *n);
void swap8(void *n);

int FLAG_x = 1;
int FLAG_v = 0;
int FLAG_f = 0;


void write_trr(int k, int cycle, int flag) {


char title[2]={""};
//char title[80];
//sprintf(title,"funk");

char nm[50];
/*
for(int i=0; i<100; i++) {
#ifdef MPI
	  sprintf(nm,"./OUTPUT/GROMACS/gromacs%d.trr%d",mpi.my_rank,i);
#endif
#ifndef MPI
	  sprintf(nm,"./OUTPUT/GROMACS/gromacs%d.trr%d",k,i);
#endif
	int fileindex = i;
	if(i == 99) {
	  fprintf(stdout,"There are too many simul.ouput files!!!\n");
	  exit(1);
	}
	else if( das=fopen(nm,"rb") ) {
		fclose(das);
		continue;
	}
	else {
	  das = fopen(nm,"ab");
	  break;
	}
}
*/
#ifdef MPI
	  sprintf(nm,"./OUTPUT/GROMACS/gromacs%d_%d.trr",mpi.my_rank,fileindex);
#endif
#ifndef MPI
	  sprintf(nm,"./OUTPUT/GROMACS/gromacs%d_%d.trr",k,fileindex);
#endif

das=fopen(nm,"ab");

//setting some variables gromacs will look for.  no clue what they do.
long ir_size=0,e_size=0,vir_size=0,pres_size=0,top_size=0,sym_size=0,nre=0;

//size of array to store vectors that define the simulation box
long box_size=precision*9;
long x_size =0; long v_size = 0; long f_size=0;
//size of position, velocity, and force sections
if (FLAG_x !=0) x_size=precision*3*box[k].boxns;
if (FLAG_v !=0) v_size=precision*3*box[k].boxns;
if (FLAG_f !=0) f_size=precision*3*box[k].boxns;

//no clue what these do
float lambda = 0.0;
long MAGIC=1993;

//version
long VERSION = 13;

long time_ind; float time_val;

if (flag==0) {
	time_val = (float) (sim.dt * 1.0*cycle);
	time_ind = cycle;
}
	
else if (flag==1 && (sim.ID2==3 || sim.ID2==5)){
	time_val= (float) (sim.dt*sim.cyc_eq+ sim.dtlong*1.0*cycle); // for multiple time step
	time_ind= sim.cyc_eq+cycle;
}
else {
	time_val = (float) (sim.dt * 1.0*(sim.cyc_eq+cycle));
	time_ind = sim.cyc_eq+cycle;
}

write_int(MAGIC);
write_int(VERSION);

//strip_white(title);
write_bstring(title);

//printf("%d %d %s\n",MAGIC,VERSION,title);

write_int(ir_size);
write_int(e_size);
write_int(box_size);
write_int(vir_size);
write_int(pres_size);
write_int(top_size);
write_int(sym_size);
write_int(x_size);
write_int(v_size);
write_int(f_size);
write_int(box[k].boxns);
write_int(time_ind);
write_int(nre);

//printf("%d %d %d %d %d %d %d %d %d %d %d %d %d\n",ir_size,e_size,box_size,vir_size,pres_size,
//	   top_size,sym_size,x_size,v_size,f_size,box[k].boxns,time_ind,nre);

write_float(time_val);
write_float(lambda);

//printf("%lf %lf\n",time_val,lambda);

//writes the array for  box dimensions
float bx[9];
int i;
for (i=0; i<9; i++)
	bx[i]=0.0;

bx[0]=(float)(box[k].lx/10.0);
bx[4]=(float)(box[k].ly/10.0);
bx[8]=(float)(box[k].lz/10.0);

write_vector(&bx[0]);
write_vector(&bx[3]);
write_vector(&bx[6]);

//printf("box:%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n\n",bx[0],bx[1],bx[2],bx[3],
//	   bx[4],bx[5],bx[6],bx[7],bx[8]);



//writes the array for positions
//converts units to nm
float pos[3];
//printf("atom posits\n");
if (FLAG_x !=0) {
	for (i=0; i<box[k].boxns; i++) {
		pos[0]=(float)(atom[k][i].x/10.0);				// figure out with or w/o pbc
		pos[1]=(float)(atom[k][i].y/10.0);
		pos[2]=(float)(atom[k][i].z/10.0);
		write_vector(pos);
	//	printf("%lf %lf %lf\n",pos[0],pos[1],pos[2]);
	}
}

//writes velocities
float vel[3];
//printf("atom velocities\n");
if (FLAG_v !=0) {
	for (i=0; i<box[k].boxns; i++) {
		vel[0]=(float)(vv[k][i].x/10.0);
		vel[1]=(float)(vv[k][i].y/10.0);
		vel[2]=(float)(vv[k][i].z/10.0);
		write_vector(vel);
	//	printf("%lf %lf %lf\n",vel[0],vel[1],vel[2]);
	}
}

//writes forces
float force[3];
//printf("atom forces\n");
if (FLAG_f !=0) {
	for (i=0; i<box[k].boxns; i++) {
		force[0]=(float)(ff[k][i].x/10.0);
		force[1]=(float)(ff[k][i].y/10.0);
		force[2]=(float)(ff[k][i].z/10.0);
		write_vector(force);
	//	printf("%lf %lf %lf\n",force[0],force[1],force[2]);
	}
}

fclose(das);
}



int write_int(long i) {
	if(!reverse) {
		if (fwrite(&i,4,1,das)!=1) 
			return 0;
		else
			return 1;
	}
	else {
		long si;
		si = i;
		swap4((void *) &si);
		if (fwrite(&si,4,1,das)!=1)
			return 0;
		else
			return 1;
	}
}


int write_float(float fp) {
	if (precision==4) {
		float sfp;
		sfp=fp;
		if (reverse)
			swap4((void *) &sfp);
		if(fwrite(&sfp,precision,1,das)!=1)
			return 0;
		else 
			return 1;
	}

	else {
		double sfp;
		sfp = (double) fp;
		if (reverse)
			swap8((void *) &sfp);

		if(fwrite(&sfp,precision,1,das)!=1)
			return 0;
		else
			return 1;
	}
}


int write_vector(float *fp) {
	if ( (!write_float(fp[0])) || (!write_float(fp[1])) ||
		(!write_float(fp[2])) ) {
		printf("failed to write vector!\n");
		exit(1);
	}
	else 
		return 1;
}


int write_bstring(char *str) {
//	printf("length: %d\nstring:\n%s\n",strlen(str),str);

	if(!write_int(strlen(str))) {
		fprintf(stdout,"failed to write string length!\n");
		exit(1);
	}
	
	if (fwrite(str,1,strlen(str),das)!=1) 
		return 0;
	else 
		return -1;
}


void swap4(void *n) {
	char temp[2];
	char *c = (char *) n;
	temp[0] = c[3];
	temp[1] = c[2];
	c[3] = c[0];
	c[2] = c[1];
	c[1] = temp[1];
	c[0] = temp[0];
	return;
}

void swap8(void *n) {
	char temp[4];
	char *c = (char *) n;
	temp[0]=c[7];
	temp[1]=c[6];
	temp[2]=c[5];
	temp[3]=c[4];
	c[7]=c[0];
	c[6]=c[1];
	c[5]=c[2];
	c[4]=c[3];
	c[3]=temp[3];
	c[2]=temp[2];
	c[1]=temp[1];
	c[0]=temp[0];
	return;
}



//strips leading and trailing white space from a string
int strip_white(char *buf) {
	int i,j,k;
	/* check arguments */
	if (!buf) return -1;
	if (!strlen(buf)) return -1;

	//ditches trailing white space
	for (i=strlen(buf)-1; buf[i]==' ' || buf[i]=='\t' || buf[i]=='\n' || buf[i]=='\r'; i--)
		buf[i]=0;

	//leading white space
	for (i=0;  buf[i]==' ' || buf[i]=='\t' || buf[i]=='\n' || buf[i]=='\r'; i++);
	if(i) {
		k=0;
		for (j=i; buf[j]; j++)
			buf[k++]=buf[j];
		
		buf[k]=0;
	}
	return strlen(buf);

}
#endif

