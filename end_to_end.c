/* ======================================================================== */
/* end_to_end.cpp															*/
/*                                                                          */
/*		This subroutine calculates the end to end distance between the		*/
/* site number "a" and site number "b" of box "k".                          */
/*                                                                          */
/* ======================================================================== */


#include "defines.h"
void min_image(int,double*,double*,double*);

void end_to_end(int k, int a, int b)
	{

//  double lx = box[k].lx;
//  double ly = box[k].ly;
//  double lz = box[k].lz;
//  double hx = box[k].hx;
//  double hy = box[k].hy;
//  double hz = box[k].hz;
//  double z_ref	=	29.0;			// This is based on ion channel settings.
  double dx, dy, dz, dr;

  if(a ==-1){
#ifdef IONC

	  /* ------------------------------------------ */
	  /* This is the case for ion channel when		*/
	  /* distance is defined based on cent. of mass	*/
	  /* ------------------------------------------	*/ 

	  int S[16]		=	{1,18,25,35,54,64,80,96,112,136,155,179,198,222,241,265};
	  double rcm	=	0.0;
	  int rcm_count	=	0;
	  	for(int i=0; i<16; i++){
			double rcmi	 =	atom[k][S[i]-1].z + atom[k][S[i]+275].z;
			rcm			+=	rcmi;
			rcm_count ++;	
		}
		rcm		/= rcm_count;
		dz	 = atnopbc[k][b].z - rcm;
		ordparam[k].d_nc = dz;

	  /* ------------------------------------------ */
	  /* This is the case for ion channel when		*/
	  /* distance is defined from a plane			*/
	  /* ------------------------------------------	*/
/*	  dx	 = atnopbc[k][b].x;
	  dy	 = atnopbc[k][b].y;
	  dz	 = atnopbc[k][b].z - z_ref;
	  dr  = fabs(dz);;
	  ordparam[k].d_nc = dr;
*/
#endif
  }
  else{
	  dx = atnopbc[k][a].x - atnopbc[k][b].x;
	  dy = atnopbc[k][a].y - atnopbc[k][b].y;
	  dz = atnopbc[k][a].z - atnopbc[k][b].z;

  /* ----------------------------------------------- */
  /* Calculate the distance between particles a & b. */
  /* Tale the square root and save the value as end  */
  /* to end distance for that box.				     */
  /* ----------------------------------------------- */
      double dr2 = dx*dx + dy*dy + dz*dz;
	  dr  = sqrt(dr2);
	  ordparam[k].d_nc = dr;
  }
//#ifdef XEDOS
  if(a ==-1){// force_2 should only be in z dir
	  ordparam[k].force_1 = ff[k][b].x;
	  ordparam[k].force_2 = ff[k][b].z;
  }
  else{
	  ordparam[k].force_1 = (ff[k][a].x*dx + ff[k][a].y*dy +ff[k][a].z*dz)/dr;
	  ordparam[k].force_2 = (ff[k][b].x*dx + ff[k][b].y*dy +ff[k][b].z*dz)/dr;
  }
//#endif
#ifdef SMD
	  ordparam[k].force_1 = (ff[k][a].x*steermd[k].ex + ff[k][a].y*steermd[k].ey +ff[k][a].z*steermd[k].ez);
	  ordparam[k].force_2 = (ff[k][b].x*steermd[k].ex + ff[k][b].y*steermd[k].ey +ff[k][b].z*steermd[k].ez);
	  steermd[k].r1x	  = -dx; // negative because dx is site 1- site 2
	  steermd[k].r1y	  = -dy;
	  steermd[k].r1z	  = -dz;
#endif
#ifdef NSMD
	  ordparam[k].force_1 = (ff[k][a].x*steermd[k].ex + ff[k][a].y*steermd[k].ey +ff[k][a].z*steermd[k].ez);
	  ordparam[k].force_2 = (ff[k][b].x*steermd[k].ex + ff[k][b].y*steermd[k].ey +ff[k][b].z*steermd[k].ez);
	  steermd[k].r1x	  = -dx; // negative because dx is site 1- site 2
	  steermd[k].r1y	  = -dy;
	  steermd[k].r1z	  = -dz;
#endif
	}

double calc_dist(int ibox, int a, int b){
	int k = ibox;
	double dx = atom[k][a].x - atom[k][b].x;
	double dy = atom[k][a].y - atom[k][b].y;
	double dz = atom[k][a].z - atom[k][b].z;
	min_image(k,&dx,&dy,&dz);
	return(sqrt(dx*dx + dy*dy + dz*dz));
}
