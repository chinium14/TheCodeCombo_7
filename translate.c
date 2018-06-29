/* ======================================================================== */
/* translate.cpp                                                            */
/*            Translates the molecules in box "k" to the center on box	    */
/* ======================================================================== */
#ifdef TLATE
#include "defines.h"
void pbc_all(int);
void translate (int k)
{
	double xc = 0.0;
	double yc = 0.0;
	double zc = 0.0;
	double mc = 0.0;

//	for(int i=0; i<box[k].boxns; i++) {
	for(int i=0; i<molec.lsite[0]; i++) {		/* This is changed so that the system is centered based on center of mass of protein	*/
		xc += atnopbc[k][i].x * pott[k][i].mas;
		yc += atnopbc[k][i].y * pott[k][i].mas;
		zc += atnopbc[k][i].z * pott[k][i].mas;
		mc += pott[k][i].mas;
	}
	xc /= mc;
	yc /= mc;
	zc /= mc;

	for(int i=0; i<box[k].boxns; i++){
		atnopbc[k][i].x -= xc;
		atnopbc[k][i].y -= yc;
		atnopbc[k][i].z -= zc;
		atom[k][i].x     = atnopbc[k][i].x;
		atom[k][i].y     = atnopbc[k][i].y;
		atom[k][i].z     = atnopbc[k][i].z;
	}
	pbc_all(k);

	/* ---------------------------------------- */
	/* Check if the center of mass is at [000]	*/
	/* ----------------------------------------	*/
/*	xc =0.0; yc =0.0; zc =0.0;
	for(int i=0; i<box[k].boxns; i++) {
		xc += atnopbc[k][i].x * pott[k][i].mas;
		yc += atnopbc[k][i].y * pott[k][i].mas;
		zc += atnopbc[k][i].z * pott[k][i].mas;
		mc += pott[k][i].mas;
	}
*/

}

#endif

