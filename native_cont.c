/* ======================================================================== */
/* native_cont.cpp                                                          */
/*                                                                          */
/* This subroutine determines the number of native contacts for the         */
/* current structure.  The native contacts are specified in the             */
/* contacts.inp file.                                                       */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"
void native_cont(int k){

	
  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;
  int con = 0;

	/* ----------------------------------------------- */
	/* Loop around the native contact list to          */
	/* determine if the contact exists.                */
	/* ----------------------------------------------- */
	if(contacts[k].n != 0){	
		
		for(int i=0; i<contacts[k].n; i++){
			contacts[k].code[i]=0;
			int a = contacts[k].a[i];
			int b = contacts[k].b[i];

			double dx = atom[k][a].x - atom[k][b].x;
			double dy = atom[k][a].y - atom[k][b].y;
			double dz = atom[k][a].z - atom[k][b].z;

	/* ----------------------------------------------- */
	/* Apply minimum image convention.                 */
	/* ----------------------------------------------- */
      if(dx >  hx) dx -= lx;
      else if(dx < -hx) dx += lx;
      if(dy >  hy) dy -= ly;
      else if(dy < -hy) dy += ly;
      if(dz >  hz) dz -= lz;
      else if(dz < -hz) dz += lz;

  /* ----------------------------------------------- */
  /* Calculate the distance between particles a & b. */
  /* Tale the square root and save the value as end  */
  /* to end distance for that box.				     */
  /* ----------------------------------------------- */
			double dr = sqrt(dx*dx + dy*dy + dz*dz);
			if(dr <= contacts[k].d[i]) {
				con++;
				contacts[k].code[i]=1;
			}
		
		}
	  
		ordparam[k].con = con;
	}//if

	int con_2=0;
	/* ----------------------------------------------- */
	/* Loop around the native contact 2 list to        */
	/* determine if the contact exists.                */
	/* ----------------------------------------------- */
	if(contacts_2[k].n != 0){	
		
		for(int i=0; i<contacts_2[k].n; i++){
			contacts_2[k].code[i]=0;
			int a = contacts_2[k].a[i];
			int b = contacts_2[k].b[i];

			double dx = atom[k][a].x - atom[k][b].x;
			double dy = atom[k][a].y - atom[k][b].y;
			double dz = atom[k][a].z - atom[k][b].z;

	/* ----------------------------------------------- */
	/* Apply minimum image convention.                 */
	/* ----------------------------------------------- */
      if(dx >  hx) dx -= lx;
      else if(dx < -hx) dx += lx;
      if(dy >  hy) dy -= ly;
      else if(dy < -hy) dy += ly;
      if(dz >  hz) dz -= lz;
      else if(dz < -hz) dz += lz;

  /* ----------------------------------------------- */
  /* Calculate the distance between particles a & b. */
  /* Tale the square root and save the value as end  */
  /* to end distance for that box.				     */
  /* ----------------------------------------------- */
			double dr = sqrt(dx*dx + dy*dy + dz*dz);
			if(dr <= contacts_2[k].d[i]) {
				con_2++;
				contacts_2[k].code[i]=1;
			}
		
		}
	  
		ordparam[k].con_2 = con_2;
	}//if
}//end subroutine
