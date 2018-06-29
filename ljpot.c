/* ======================================================================== */
/* ljpot.cpp                                                                */
/* This subroutine will set number of bonds, bends and tosions to zero      */
/* Also the epsilon and sigma will all be made identical to the desired		*/
/* value irrespective of each atom type. The product of charges for each	*/
/* atom pair is all set to 0. Therefore this ensures that energy due to		*/
/* bonds, bends, torsions, and coulombic is zero. The esasa is computed		*/
/* only if LJ is not defined. Thus if LJ is defines only lennard jones		*/
/* component is computed. Make sure that this subroutine is called after	*/
/* setlj , to ensure that no exclusions are defined.						*/
/* ======================================================================== */
#include "defines.h"
#ifdef LJ

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void ljpot ()
{
  /* ================================================================== */
  /*                                                                    */
  /* Variables needed in the subroutine.                                */
  /*                                                                    */
  /* ================================================================== */

  /* ------------------------------------------------------- */
  /* Parameters for Argon                                    */
  /* ------------------------------------------------------- */
	double sig  = 3.432;									/* put in dimensional sigma in A units			*/
	double eps  = 122.4;									/* put in dimensional epsilon in K units	    */
	double mass = 39.948; 									/* put in dimensional mass in AMU               */
  /* ------------------------------------------------------- */
  /* Parameters for Neon                                     */
  /* ------------------------------------------------------- */
	//double sig  = 2.789;									/* put in dimensional sigma in A units			*/
	//double eps  = 35.7;									/* put in dimensional epsilon in K units	    */
	//double mass = 20.183; 								/* put in dimensional mass in AMU               */

  /* ================================================================== */
  /*                                                                    */
  /* Reassign all site parameters to a noble gas value.                 */
  /*                                                                    */
  /* ================================================================== */
	
	for (int k=0; k<sim.NB; k++){
		for(int a=0; a<box[k].boxns; a++) {
//			if(a!=146 && a!=158){
				pott[k][a].eps = eps*RG*0.001;
				pott[k][a].sig = sig;
//			}
//			else{
//				pott[k][a].eps = eps*RG*0.001*0.3;
//				pott[k][a].sig = sig* 0.8;
//			}
#ifndef EWALD_DEBUG
				pott[k][a].qq  = 0.0;
#endif
#ifdef EWALD_DEBUG
				
				if(a%2==0){
					pott[k][a].qq  = -1.0;
					sprintf(atom[k][a].name,"NEG");
					if(a==box[k].boxns-1){
						pott[k][a].qq = 0.0;
						sprintf(atom[k][a].name,"NOQ");
					}
				}
				else{
					pott[k][a].qq = 1.0;
					sprintf(atom[k][a].name,"POS");
				}
#endif
				pott[k][a].mas	=mass;
		}// a ends
  /* ================================================================== */
  /*                                                                    */
  /* Stop calculation of all interactions from bonds, bends, torsions,  */
  /* and improper tosions.                                              */
  /*                                                                    */
  /* ================================================================== */
#ifdef EWALD_DEBUG
		double qqsum=0.0;
		for(int a=0; a<box[k].boxns; a++) qqsum += pott[k][a].qq;
		if(qqsum != 0.0){
			printf("sumqq = %lf LJ Charges not equal to zero\n", qqsum);
			exit(9);
		}
#endif
		bondN[k]=0;
		bendN[k]=0;
		torsN[k]=0;
		imprN[k]=0;

#ifdef SC
		int isp =0;
		for ( double ix=sig/2.0; ix<(box[k].lx]-sig/2); ix+=(1.002*sig)){
			for ( double iy=sig/2.0; iy<(box[k].ly-sig/2); iy+=(1.002*sig)){
				for ( double iz=sig/2.0; iz<(box[k].lz-sig/2); iz+=(1.002*sig)){
					if (isp<box[k].boxns){
						atom[k][isp].x	=	ix;atnopbc[k][isp].x=atom[k][isp].x;
						atom[k][isp].y	=	iy;atnopbc[k][isp].y=atom[k][isp].y;
						atom[k][isp].z	=	iz;atnopbc[k][isp].z=atom[k][isp].z;
						isp ++;
					}
				}
			}
		}
#endif
		
#ifdef FCC		
  /* ================================================================== */
  /*                                                                    */
  /* Assign the positions to an fcc lattice.                            */
  /*                                                                    */
  /* ================================================================== */
  /* -------------------------------------------------------- */
  /* Variables needed for position assignment.                */
  /* -------------------------------------------------------- */
	
		int nlin; //number of lines needed on each side of the box
		double a; //length of one unit cell
		int particle; //counter for the number of particles
		char ch;
		int i,xdir,ydir,zdir;

  /* -------------------------------------------------------- */
  /* Zero out the counter.                                    */
  /* -------------------------------------------------------- */
		particle = 0;
		ch = '0';


	
	
  /* -------------------------------------------------------- */
  /* Calculate out the number of particles needed on each     */
  /* side of the box.                                         */
  /* -------------------------------------------------------- */
		nlin = (int)pow((double)box[k].boxns/4.0,1.0/3.0);
		if ( ((double)nlin*(double)nlin*(double)nlin) < (double)box[k].boxns/4.0) nlin = nlin+1;//if N is not a cube root, add 1 to nlin
	
  /* -------------------------------------------------------- */
  /* Calculate the length of one unit cell.                   */
  /* -------------------------------------------------------- */
		a = box[k].lx/(double)nlin;


  /* -------------------------------------------------------- */
  /* Assign the positions to an fcc lattice.                  */
  /* -------------------------------------------------------- */
		for (zdir = 0; zdir<nlin; zdir++)
		{
			if (particle == box[k].boxns) break;
			for (ydir = 0; ydir<nlin; ydir++)
			{
				if (particle == box[k].boxns) break;
				for (xdir = 0; xdir<nlin; xdir++)
				{
					if (particle == box[k].boxns) break;
					for(i = 0; i<4; i++)
					{
						if (particle == box[k].boxns) break;
						switch(ch) {
							case '0':
								atom[k][particle].x	=	0.0 + (double)xdir*a; 
								atom[k][particle].y	=	0.0 + (double)ydir*a; 
								atom[k][particle].z	=	0.0 + (double)zdir*a;
								atnopbc[k][particle].x=atom[k][particle].x;
								atnopbc[k][particle].y=atom[k][particle].y;
								atnopbc[k][particle].z=atom[k][particle].z;
								ch = '1';
								particle++;
							
								break;
							case '1':
								atom[k][particle].x = 0.0   + (double)xdir*a;
								atom[k][particle].y = 0.5*a + (double)ydir*a;
								atom[k][particle].z = 0.5*a + (double)zdir*a;
								atnopbc[k][particle].x=atom[k][particle].x;
								atnopbc[k][particle].y=atom[k][particle].y;
								atnopbc[k][particle].z=atom[k][particle].z;
								ch = '2';
								particle++;
							
								break;
							case '2':
								atom[k][particle].x = 0.5*a + (double)xdir*a;
								atom[k][particle].y = 0.0 + (double)ydir*a;
								atom[k][particle].z = 0.5*a + (double)zdir*a;
								atnopbc[k][particle].x=atom[k][particle].x;
								atnopbc[k][particle].y=atom[k][particle].y;
								atnopbc[k][particle].z=atom[k][particle].z;
								ch = '3';
								particle++;
							
								break;
							case '3':
								atom[k][particle].x = 0.5*a + (double)xdir*a;
								atom[k][particle].y = 0.5*a + (double)ydir*a;
								atom[k][particle].z = 0.0 + (double)zdir*a;
								atnopbc[k][particle].x=atom[k][particle].x;
								atnopbc[k][particle].y=atom[k][particle].y;
								atnopbc[k][particle].z=atom[k][particle].z;
								ch = '0';
								particle++;
							
								break;
						}//switch		
					}// for i
				}// for zdir
			}// for ydir
		}// for xdir
#endif		
	}// k ends
}// subroutine ends
#endif


