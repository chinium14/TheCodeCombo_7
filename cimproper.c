/* ======================================================================== */
/* cimproper.cpp                                                            */
/*                                                                          */
/*		This subroutine calculates the improper torsion energies and        */
/* forces. It returns the value of the potential energy due to the improper */
/* torsions to be stored in en[k].eimpr. It is called from either forces()  */
/* or force_long().                                                         */
/*                                                                          */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

#ifdef PR_NPT
void min_image_npt_full(int, double*, double*, double*);
#endif

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

double cimproper (int ibox)
{

  int k = ibox;

  /* ================================================================== */
  /* Zero out the energy accumulator                                    */
  /* ================================================================== */
  double eimproper = 0.0;

  /* ================================================================== */
  /* Zero out the virial accumulators                                   */
  /* ================================================================== */
#ifdef PRESSURE
  double wxx = 0.0;
  double wyx = 0.0;
  double wyy = 0.0;
  double wzx = 0.0;
  double wzy = 0.0;
  double wzz = 0.0;
#endif

  /* ================================================================== */
  /* Get and assign information for pbc.                                */
  /* ================================================================== */
  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;


  /* ================================================================== */
  /* Calculate the energy and force associated with each dihedral       */
  /* angle. The loop is over the number of dihedral angles.             */
  /* ================================================================== */
  for(int i=0; i<imprN[k]; i++) {
  /* ----------------------------------------- */
  /* Assign the atoms involved in the          */
  /* dihedral angle.                           */
  /* ----------------------------------------- */
	int ia = impr[k][i].a;
    int ib = impr[k][i].b;
    int ic = impr[k][i].c;
    int id = impr[k][i].d;

  /* ----------------------------------------- */
  /* Calculate the distances between the atoms */
  /* involved in the dihedral angle.           */
  /* ----------------------------------------- */
    double dabx = atom[k][ib].x - atom[k][ia].x;
    double daby = atom[k][ib].y - atom[k][ia].y;
    double dabz = atom[k][ib].z - atom[k][ia].z;
    double dbcx = atom[k][ic].x - atom[k][ib].x;
    double dbcy = atom[k][ic].y - atom[k][ib].y;
    double dbcz = atom[k][ic].z - atom[k][ib].z;
    double dcdx = atom[k][id].x - atom[k][ic].x;
    double dcdy = atom[k][id].y - atom[k][ic].y;
    double dcdz = atom[k][id].z - atom[k][ic].z;

  /* ----------------------------------------- */
  /* Apply minimum image convention.           */
  /* ----------------------------------------- */
#ifndef PR_NPT
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
#endif
#ifdef PR_NPT
	min_image_npt_full(k, &dabx, &daby, &dabz);
	min_image_npt_full(k, &dbcx, &dbcy, &dbcz);
	min_image_npt_full(k, &dcdx, &dcdy, &dcdz);
#endif

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
  /* Calculate the improper torsion angle.     */
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
    double tau = signum * acos(costau);


  /* ----------------------------------------- */
  /* Calculate the force and energy.           */
  /* ----------------------------------------- */
    double fr = 0.0;
 //   double eimpr = 0.0;

    double rr = tau - impr[k][i].angeq;
	fr	=	-2.0*impr[k][i].kimpr*rr;
	eimproper += -fr*rr*0.5;


  /* ----------------------------------------- */
  /* Transform the forces to cartesian         */
  /* coordinates and accumulate.               */
  /* ----------------------------------------- */
    double fa    = -fr * irab * isinb2;
    double fax   = fa * abbcx;
    double fay   = fa * abbcy;
    double faz   = fa * abbcz;
    ff[k][ia].x += fax;
    ff[k][ia].y += fay;
    ff[k][ia].z += faz;

    double fb1   = fr * (rbc - rab*cosb) * irab *irbc * isinb2;
    double fb2   = fr * cosc * irbc * isinc2;
    ff[k][ib].x += (fb1 * abbcx + fb2 * dccbx);
    ff[k][ib].y += (fb1 * abbcy + fb2 * dccby);
    ff[k][ib].z += (fb1 * abbcz + fb2 * dccbz);

    double fc1   = fr * (rbc - rcd*cosc) * ircd *irbc *isinc2;
    double fc2   = fr * cosb * irbc * isinb2;
    double fcx   = fc1 * dccbx + fc2 * abbcx;
    double fcy   = fc1 * dccby + fc2 * abbcy;
    double fcz   = fc1 * dccbz + fc2 * abbcz;
    ff[k][ic].x += fcx;
    ff[k][ic].y += fcy;
    ff[k][ic].z += fcz;
    
    double fd    = -fr *ircd * isinc2;
    double fdx   = fd * dccbx;
    double fdy   = fd * dccby;
    double fdz   = fd * dccbz;
    ff[k][id].x += fdx;
    ff[k][id].y += fdy;
    ff[k][id].z += fdz;

  /* ----------------------------------------- */
  /* If PRESSURE is defined, accumulate the    */
  /* components of the virial tensor.          */
  /* ----------------------------------------- */
#ifdef PRESSURE
    double sxx = ( - fax * dabx + fcx * dbcx + fdx *(dbcx+dcdx) ); 
    double syx = ( - fay * dabx + fcy * dbcx + fdy *(dbcx+dcdx) ); 
    double syy = ( - fay * daby + fcy * dbcy + fdy *(dbcy+dcdy) ); 
    double szx = ( - faz * dabx + fcz * dbcx + fdz *(dbcx+dcdx) ); 
    double szy = ( - faz * daby + fcz * dbcy + fdz *(dbcy+dcdy) ); 
    double szz = ( - faz * dabz + fcz * dbcz + fdz *(dbcz+dcdz) ); 
	wxx += sxx;
	wyx += syx;
	wyy += syy;
	wzx += szx;
	wzy += szy;
	wzz += szz;
#ifdef ZHOU
	stresses[k][ia][0][0] += 0.25 * sxx;
	stresses[k][ia][1][0] += 0.25 * syx;
	stresses[k][ia][1][1] += 0.25 * syy;
	stresses[k][ia][2][0] += 0.25 * szx;
	stresses[k][ia][2][1] += 0.25 * szy;
	stresses[k][ia][2][2] += 0.25 * szz;

	stresses[k][ib][0][0] += 0.25 * sxx;
	stresses[k][ib][1][0] += 0.25 * syx;
	stresses[k][ib][1][1] += 0.25 * syy;
	stresses[k][ib][2][0] += 0.25 * szx;
	stresses[k][ib][2][1] += 0.25 * szy;
	stresses[k][ib][2][2] += 0.25 * szz;

	stresses[k][ic][0][0] += 0.25 * sxx;
	stresses[k][ic][1][0] += 0.25 * syx;
	stresses[k][ic][1][1] += 0.25 * syy;
	stresses[k][ic][2][0] += 0.25 * szx;
	stresses[k][ic][2][1] += 0.25 * szy;
	stresses[k][ic][2][2] += 0.25 * szz;

	stresses[k][id][0][0] += 0.25 * sxx;
	stresses[k][id][1][0] += 0.25 * syx;
	stresses[k][id][1][1] += 0.25 * syy;
	stresses[k][id][2][0] += 0.25 * szx;
	stresses[k][id][2][1] += 0.25 * szy;
	stresses[k][id][2][2] += 0.25 * szz;
#endif
#endif

  }

  /* ================================================================== */
  /*                                                                    */
  /* After the loop, if PRESSURE is defined, store the bonding          */
  /* contribution to the virial tensor to pvir[k].bond[i].              */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  pvir[k].impr[0] = wxx;
  pvir[k].impr[1] = wyx;
  pvir[k].impr[2] = wyy;
  pvir[k].impr[3] = wzx;
  pvir[k].impr[4] = wzy;
  pvir[k].impr[5] = wzz;
#endif

  return (eimproper);

}
