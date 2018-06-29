#ifdef GOLIK
/* ======================================================================== */
/* cgolik.cpp                                                               */
/*	This file contains subroutines that calculate the Go-type model of  */
/* Hoang and Cieplak(JCP Vol. 112 Number 15 Page 6851 (2000)).              */
/*									    */
/* Written by Thomas Knotts	8 Mar 04				    */
/* Cosolvent interaction added by Nitin Rathore				    */
/* ======================================================================== */


/* ======================================================================== */
/* cbb_golik(int)                                                           */
/*	This subroutine calculates the backbone potential for the Go-type   */
/* molel.								    */
/* ======================================================================== */

#include "defines.h"
#ifdef SPME
double erfc(double);
void pme_calc(int,double*);
#endif


double cbb_golik (int ibox)
{

  int k = ibox;
  /* ================================================================== */
  /*                                                                    */
  /* Zero out the energy accumulator                                    */
  /*                                                                    */
  /* ================================================================== */
  double ebond = 0.0;

  /* ================================================================== */
  /*                                                                    */
  /* Zero out the virial accumulators                                   */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  double wxx = 0.0;
  double wyx = 0.0;
  double wyy = 0.0;
  double wzx = 0.0;
  double wzy = 0.0;
  double wzz = 0.0;
#endif

  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;

  double k1 = -golik[k].eps;
  double k2 = k1 * 100;
  /* ================================================================== */
  /*                                                                    */
  /* Calculate the energy and force associated with each bond.  The     */
  /* loop is over the number of bonds.                                  */
  /*                                                                    */
  /* ================================================================== */
  for(int i=0; i<bondN[k]; i++) {
  /* ----------------------------------------- */
  /* Assign the atoms involved in the bonding. */
  /* ----------------------------------------- */
    int ia = bond[k][i].a;
    int ib = bond[k][i].b;
    double req = bond[k][i].req;
  /* ================================================================== */
  /* ----------------------------------------- */
  /* Calculate the distances between the atoms */
  /* involved in the bonding.                  */
  /* ----------------------------------------- */
    double dx = atom[k][ib].x - atom[k][ia].x;
    double dy = atom[k][ib].y - atom[k][ia].y;
    double dz = atom[k][ib].z - atom[k][ia].z;

  /* ----------------------------------------- */
  /* Apply minimum image convention.           */
  /* ----------------------------------------- */
    dx = MINIMG(dx,hx,lx);
    dy = MINIMG(dy,hy,ly);
    dz = MINIMG(dz,hz,lz);

  /* ----------------------------------------- */
  /* Calculate the length of the bond and then */
  /* the energy using the Go-type model.       */
  /* ----------------------------------------- */
    double drr = dx*dx + dy*dy + dz*dz;
    double dr  = sqrt(drr);
    double dri = 1.0 / dr;

    double rr = dr - req;
    double rr2= rr * rr;
    double fr1 = k1*rr*2.0;
    double fr2 = k2*rr*rr2*4.0;
    ebond += -rr*(fr1*0.5 + fr2*0.25);
    double fr = fr1 + fr2;
  /* ----------------------------------------- */
  /* Calculate the forces.                     */
  /* ----------------------------------------- */
    double fx = fr * dri * dx;
    double fy = fr * dri * dy;
    double fz = fr * dri * dz;

    ff[k][ia].x -= fx;
    ff[k][ia].y -= fy;
    ff[k][ia].z -= fz;
    ff[k][ib].x += fx;
    ff[k][ib].y += fy;
    ff[k][ib].z += fz;

  /* ----------------------------------------- */
  /* If PRESSURE is defined, accumulate the    */
  /* components of the virial tensor.          */
  /* ----------------------------------------- */
#ifdef PRESSURE
    double sxx = fx * dx;
    double syx = fy * dx;
    double syy = fy * dy;
    double szx = fz * dx;
    double szy = fz * dy;
    double szz = fz * dz;
	wxx += sxx;
	wyx += syx;
	wyy += syy;
	wzx += szx;
	wzy += szy;
	wzz += szz;
#ifdef ZHOU
	stresses[k][ia][0][0] += 0.5*sxx;
	stresses[k][ia][1][0] += 0.5*syx;
	stresses[k][ia][1][1] += 0.5*syy;
	stresses[k][ia][2][0] += 0.5*szx;
	stresses[k][ia][2][1] += 0.5*szy;
	stresses[k][ia][2][2] += 0.5*szz;

	stresses[k][ib][0][0] += 0.5*sxx;
	stresses[k][ib][1][0] += 0.5*syx;
	stresses[k][ib][1][1] += 0.5*syy;
	stresses[k][ib][2][0] += 0.5*szx;
	stresses[k][ib][2][1] += 0.5*szy;
	stresses[k][ib][2][2] += 0.5*szz;
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
  pvir[k].bond[0] = wxx;
  pvir[k].bond[1] = wyx;
  pvir[k].bond[2] = wyy;
  pvir[k].bond[3] = wzx;
  pvir[k].bond[4] = wzy;
  pvir[k].bond[5] = wzz;
#endif

  return (ebond);

}

/* ======================================================================== */
/* cba_golik(int)                                                           */
/*	This subroutine calculates the backbone potential for the Go-type   */
/* molel.								    */
/* ======================================================================== */
double cba_golik (int ibox)
{

  int k = ibox;
  /* ================================================================== */
  /* Zero out the energy accumulator                                    */
  /* ================================================================== */
  double ebend = 0.0;

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

  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;

#ifdef DNA_GOLIK
  double k_theta = -10.0*golik[k].eps*20;  //10.0=20/2
#else
  double k_theta = -10.0*golik[k].eps;  //10.0=20/2
#endif

#ifdef DRSAM
  k_theta = -1400*golik[k].eps;
#endif

  /* ================================================================== */
  /* Calculate the energy and force associated with each bend.  The     */
  /* loop is over the number of bends.                                  */
  /* ================================================================== */

  for(int i=0; i<bendN[k]; i++) {
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
    double dax = atom[k][ia].x - atom[k][ib].x;
    double day = atom[k][ia].y - atom[k][ib].y;
    double daz = atom[k][ia].z - atom[k][ib].z;
    double dcx = atom[k][ic].x - atom[k][ib].x;
    double dcy = atom[k][ic].y - atom[k][ib].y;
    double dcz = atom[k][ic].z - atom[k][ib].z;

  /* ----------------------------------------- */
  /* Apply minimum image convention.           */
  /* ----------------------------------------- */
    dax = MINIMG(dax,hx,lx);
    day = MINIMG(day,hy,ly);
    daz = MINIMG(daz,hz,lz);

    dcx = MINIMG(dcx,hx,lx);
    dcy = MINIMG(dcy,hy,ly);
    dcz = MINIMG(dcz,hz,lz);

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
    double phi = acos(cosphi);
    double isinphi = 1.0 / sqrt(1.0 - cosphi*cosphi);

    double rr = phi - bend[k][i].angeq;
    double fr = k_theta * rr * 2.0;
    ebend += -fr * rr * 0.5;

    double irasin = drai * isinphi * fr;
    double ircsin = drci * isinphi * fr;

  /* ----------------------------------------- */
  /* Calculate the forces.                     */
  /* ----------------------------------------- */
    double fxa = (erax * cosphi - ercx) * irasin;
    double fya = (eray * cosphi - ercy) * irasin;
    double fza = (eraz * cosphi - ercz) * irasin;
    double fxc = (ercx * cosphi - erax) * ircsin;
    double fyc = (ercy * cosphi - eray) * ircsin;
    double fzc = (ercz * cosphi - eraz) * ircsin;


    ff[k][ia].x += fxa;
    ff[k][ia].y += fya;
    ff[k][ia].z += fza;
    ff[k][ib].x += -(fxa + fxc);
    ff[k][ib].y += -(fya + fyc);
    ff[k][ib].z += -(fza + fzc);
    ff[k][ic].x += fxc;
    ff[k][ic].y += fyc;
    ff[k][ic].z += fzc;

  /* ----------------------------------------- */
  /* If PRESSURE is defined, accumulate the    */
  /* components of the virial tensor.          */
  /* ----------------------------------------- */
#ifdef PRESSURE
    double sxx = (fxa * dax + fxc * dcx);
    double syx = (fya * dax + fyc * dcx);
    double syy = (fya * day + fyc * dcy);
    double szx = (fza * dax + fzc * dcx);
    double szy = (fza * day + fzc * dcy);
    double szz = (fza * daz + fzc * dcz);
	wxx += sxx;
	wyx += syx;
	wyy += syy;
	wzx += szx;
	wzy += szy;
	wzz += szz;

#ifdef ZHOU
	double third = 1.0 / 3;
	stresses[k][ia][0][0] += third * sxx;
	stresses[k][ia][1][0] += third * syx;
	stresses[k][ia][1][1] += third * syy;
	stresses[k][ia][2][0] += third * szx;
	stresses[k][ia][2][1] += third * szy;
	stresses[k][ia][2][2] += third * szz;

	stresses[k][ib][0][0] += third * sxx;
	stresses[k][ib][1][0] += third * syx;
	stresses[k][ib][1][1] += third * syy;
	stresses[k][ib][2][0] += third * szx;
	stresses[k][ib][2][1] += third * szy;
	stresses[k][ib][2][2] += third * szz;

	stresses[k][ic][0][0] += third * sxx;
	stresses[k][ic][1][0] += third * syx;
	stresses[k][ic][1][1] += third * syy;
	stresses[k][ic][2][0] += third * szx;
	stresses[k][ic][2][1] += third * szy;
	stresses[k][ic][2][2] += third * szz;
#endif
#endif

  }

  /* ================================================================== */
  /*                                                                    */
  /* After the loop, if PRESSURE is defined, store the bending          */
  /* contribution to the virial tensor to pvir[k].bend[i].              */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  pvir[k].bend[0] = wxx;
  pvir[k].bend[1] = wyx;
  pvir[k].bend[2] = wyy;
  pvir[k].bend[3] = wzx;
  pvir[k].bend[4] = wzy;
  pvir[k].bend[5] = wzz;
#endif

  return (ebend);

}

/* ======================================================================== */
/* cda_golik(int)                                                           */
/*	This subroutine calculates the backbone potential for the Go-type   */
/* molel.  Dihedrals.							    */
/* ======================================================================== */
/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

double cda_golik (int ibox)
{

  int k = ibox;

  /* ================================================================== */
  /* Zero out the energy accumulator                                    */
  /* ================================================================== */
  double etorsion = 0.0;

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

#ifdef DNA_GOLIK
  double k_phi = 0.2 * golik[k].eps*20;
#else
  double k_phi = 0.2 * golik[k].eps;
#endif

#ifdef DRSAM
  k_phi = 28 * golik[k].eps;
#endif

  /* ================================================================== */
  /* Calculate the energy and force associated with each dihedral       */
  /* angle. The loop is over the number of dihedral angles.             */
  /* ================================================================== */

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
    dabx = MINIMG(dabx,hx,lx);
    daby = MINIMG(daby,hy,ly);
    dabz = MINIMG(dabz,hz,lz);

    dbcx = MINIMG(dbcx,hx,lx);
    dbcy = MINIMG(dbcy,hy,ly);
    dbcz = MINIMG(dbcz,hz,lz);

    dcdx = MINIMG(dcdx,hx,lx);
    dcdy = MINIMG(dcdy,hy,ly);
    dcdz = MINIMG(dcdz,hz,lz);

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
    double tau = signum * acos(costau);

  /* -----------------------------------------	*/
  /* Assign the values for phi and psi		*/
  /* -----------------------------------------	*/
	if(tors[k][i].phitag ==1)		tors[k][i].phi = tau;
	else if (tors[k][i].psitag ==1) tors[k][i].psi = tau;
	else if (tors[k][i].thetag ==1) tors[k][i].theta = tau;
  /* ----------------------------------------- */
  /* Calculate the force and energy.           */
  /* ----------------------------------------- */
#ifdef DNA_GOLIK
	double fr  =  k_phi*(sin(tau-tors[k][i].delphi[0]));
	etorsion  +=  k_phi*(1+cos(tau-tors[k][i].delphi[0]));
#else
	double fr  =  k_phi*3.0*(sin(3.0*tau-tors[k][i].delphi[0]));
	etorsion  +=  k_phi*(1+cos(3.0*tau-tors[k][i].delphi[0]));
#endif
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
  pvir[k].tors[0] = wxx;
  pvir[k].tors[1] = wyx;
  pvir[k].tors[2] = wyy;
  pvir[k].tors[3] = wzx;
  pvir[k].tors[4] = wzy;
  pvir[k].tors[5] = wzz;
#endif
  return (etorsion);
}

#ifndef DNA_GOLIK
/* ======================================================================== */
/* cnbnd_golik(int, double*)                                                */
/*                                                                          */
/*	This subroutine calculates the nonbonded native and non-native      */
/* energies of the Go-type model.                                           */
/*                                                                          */
/* Passed Parameters:                                                       */
/*		ibox:	The box number             		            */
/*              *energy:	    A pointer to an array defined in        */
/*                                  forces() or force_long().  It contains  */
/*                                  three entries.  Two for unshifted and   */
/*                                  shifted vdw energies and one for        */
/*                                  coulombic energies.                     */
/*                                                                          */
/* ======================================================================== */

void cnbnd_golik(int ibox, double *energy)
{

  int k = ibox;
  double fconv = (ELEC * ELEQ * ELEQ * NA * 0.001);   // convert atomic charges to units of kJ/mol
  /* ================================================================== */
  /*                                                                    */
  /* Zero out the energy accumulator                                    */
  /*                                                                    */
  /* ================================================================== */
  double enonbond	= 0.0;
  double enonbonds	= 0.0;
  double enonnative	= 0.0;
  double ecoulomb	= 0.0;
  double e_pp		= 0.0;
  double e_ps		= 0.0;
  double e_ss		= 0.0;

  /* ================================================================== */
  /*                                                                    */
  /* Zero out the virial accumulators                                   */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  double wxx = 0.0;
  double wyx = 0.0;
  double wyy = 0.0;
  double wzx = 0.0;
  double wzy = 0.0;
  double wzz = 0.0;
#endif//PRESSURE

  /* ================================================================== */
  /*                                                                    */
  /* Get and assign information for pbc and cutoffs                     */
  /*                                                                    */
  /* ================================================================== */
  double rc   = sim.rc;
  double rc2  = box[k].rc2;
  double rc4  = box[k].rc4;
  double rc_rep = golik[k].d_repulsive;
  double rc2_rep = rc_rep * rc_rep;


  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;

  /* ================================================================== */
  /*                                                                    */
  /* Calculate the energy and force on each molecule due to atoms in    */
  /* its neighborlist. The first loop is over sites/atoms 0 to N-1      */
  /* and the second is over the number of neighbors in atom a's         */
  /* neighborlist.                                                      */
  /*                                                                    */
  /* ================================================================== */

  for(int a=0; a<box[k].boxns-1; a++) {

  /* ----------------------------------------------- */
  /* Zero out the force accumulators for particle a. */
  /* ----------------------------------------------- */
    double fxa = 0.0;
    double fya = 0.0;
    double fza = 0.0;

  /* ----------------------------------------------- */
  /* Assign position of particle a.                  */
  /* ----------------------------------------------- */    
    double ax = atom[k][a].x;
    double ay = atom[k][a].y;
    double az = atom[k][a].z;

  /* ----------------------------------------------- */
  /* Loop around neighbors of particle a.            */
  /* ----------------------------------------------- */
    for(int bb=0; bb<nlist[k].count[a]; bb++) {
      int b = nlist[k].list[a][bb];
  /* ----------------------------------------------- */
  /* Calculate the x,y,z distances between           */
  /* particles a & b.                                */
  /* ----------------------------------------------- */
      double dx = ax - atom[k][b].x;
      double dy = ay - atom[k][b].y;
      double dz = az - atom[k][b].z;

  /* ----------------------------------------------- */
  /* Apply minimum image convention.                 */
  /* ----------------------------------------------- */
    dx = MINIMG(dx,hx,lx);
    dy = MINIMG(dy,hy,ly);
    dz = MINIMG(dz,hz,lz);

  /* ----------------------------------------------- */
  /* Calculate the distance between particles a & b. */
  /* If the distance is less than the cutoff,        */
  /* calculated the energy and the force. If it is   */
  /* greater, go to the next pair in the loop.       */
  /* ----------------------------------------------- */
      double dr2 = dx*dx + dy*dy + dz*dz;

      if(dr2 > rc2) continue;

  /* ----------------------------------------------- */
  /* Set the epsilon, sigma, and qq for the pair.    */
  /* ----------------------------------------------- */
	  double sig;
	  double eps;
	  double qq;
	  double fx, fy, fz;
	  int flag=0;
	  if (a<golik[k].Psite && b< golik[k].Psite){
		  flag = golik[k].contact[a][b];
	  }
	  else flag =1;				// need to put flag==1 for p-s and s-s interactions
	  if(flag == 0 && dr2 > rc2_rep) continue;
  /* ----------------------------------------------- */
  /* Set the site types of particles a and b.        */
  /* ----------------------------------------------- */  
      sig = ljset[k].pot[a][b].sig;
      eps = ljset[k].pot[a][b].eps;
      qq  = ljset[k].pot[a][b].qq;
      
  /* ----------------------------------------------- */
  /* Calculate the LJ force and energy.              */
  /* ----------------------------------------------- */
	  double force=0.0; double enbnd =0.0; double enbnds = 0.0;
	  if((eps*sig)!=0.0) {
	  
		  double sr2  = (sig * sig) / dr2;
		  double sr6  = sr2 * sr2 * sr2;
		  double sr12 = sr6 * sr6;

		  double sr2s  = (sig * sig) / rc2;
		  double sr6s  = sr2s * sr2s * sr2s;
		  double sr12s = sr6s * sr6s;

//		  double nu_nu	= (golik[k].hydro_index[a]*golik[k].hydro_index[b])/3.0;
//		  double lambda = gosol.lambda;
		  double nu_nu_ps ;
		  double nu_nu_ss ;
		  double lambdaps = gosol.lambda;

		  if (a< golik[k].Psite && b < golik[k].Psite){

			  force = (48.0*eps/dr2) * (sr12 - 0.5*sr6);
			  enbnd = 4.0*eps * (sr12 - sr6);
			  enbnds= 4.0*eps * (sr12s -sr6s);
			  
			  if(flag==0) {
				  e_pp += (enbnd +eps);
				  enonnative  += (enbnd + eps);
			  }
			  else {
				  e_pp += (enbnd -enbnds);
				  enonbond += enbnd;
				  enonbonds+= enbnds;
			  }
		  }
		  else if (a< golik[k].Psite && b > golik[k].Psite) {
			/* --------------------------------------------	*/
			/* Define the protein - solvent A interaction	*/
			/* --------------------------------------------	*/
			  if(b<(golik[k].Psite+gosol.Nsol[0])) 	nu_nu_ps = 0.0;
			/* --------------------------------------------	*/
			/* Define the protein - Solvent B interaction	*/
			/* --------------------------------------------	*/
			  else 			nu_nu_ps = 1.0;
			  
			  force = lambdaps*(48.0*eps/dr2) * (sr12 - 0.5*nu_nu_ps*sr6);
			  enbnd = lambdaps* 4.0 *eps * (sr12 - nu_nu_ps*sr6);
			  enbnds= lambdaps* 4.0 *eps * (sr12s -nu_nu_ps*sr6s);
			  e_ps += (enbnd -enbnds);
			  enonbond += enbnd;
			  enonbonds+= enbnds;
		  }
		  else if (a> golik[k].Psite && b > golik[k].Psite) {
			  if(a<(golik[k].Psite+gosol.Nsol[0]) && b > (golik[k].Psite+gosol.Nsol[0])) nu_nu_ss= 0.5;
			  else nu_nu_ss = 0.5;
			  force = (48.0*eps/dr2) * (sr12 - 0.5*nu_nu_ss*sr6);
			  enbnd = 4.0*eps * (sr12  - nu_nu_ss*sr6);
			  enbnds= 4.0*eps * (sr12s - nu_nu_ss*sr6s);
			  e_ss += (enbnd -enbnds);
			  enonbond += enbnd;
			  enonbonds+= enbnds;
		  }
/*			  
		  if(flag == 0){
			enonnative  += (enbnd + eps);
		  }
		  else{
			enonbond += enbnd;
			enonbonds+= enbnds;
		  }
*/
	  }
	  
  /* ----------------------------------------- */
  /* Accumulate the forces.                    */
  /* ----------------------------------------- */
	if(force != 0.0){
       /* *********************************************** */
       /* x-forces                                        */
       /* *********************************************** */
      fx          = force * dx;
      fxa        += fx;
      ff[k][b].x -= fx;

       /* *********************************************** */
       /* y-forces                                        */
       /* *********************************************** */
      fy          = force * dy;
      fya        += fy;
      ff[k][b].y -= fy;

       /* *********************************************** */
       /* z-forces                                        */
       /* *********************************************** */
      fz          = force * dz;
      fza        += fz;
      ff[k][b].z -= fz;
  /* ----------------------------------------- */
  /* If PRESSURE is defined, accumulate the    */
  /* components of the virial tensor.          */
  /* ----------------------------------------- */
#ifdef PRESSURE
      double sxx = fx * dx;
      double syx = fy * dx;
      double syy = fy * dy; 
      double szx = fz * dx;
      double szy = fz * dy;
      double szz = fz * dz;
	wxx += sxx;
	wyx += syx;
	wyy += syy;
	wzx += szx;
	wzy += szy;
	wzz += szz;
#ifdef ZHOU
	stresses[k][a][0][0] += 0.5*sxx;
	stresses[k][a][1][0] += 0.5*syx;
	stresses[k][a][1][1] += 0.5*syy;
	stresses[k][a][2][0] += 0.5*szx;
	stresses[k][a][2][1] += 0.5*szy;
	stresses[k][a][2][2] += 0.5*szz;

	stresses[k][b][0][0] += 0.5*sxx;
	stresses[k][b][1][0] += 0.5*syx;
	stresses[k][b][1][1] += 0.5*syy;
	stresses[k][b][2][0] += 0.5*szx;
	stresses[k][b][2][1] += 0.5*szy;
	stresses[k][b][2][2] += 0.5*szz;
#endif
#endif//PRESSURE
	}//if(force!=0.0)
    }//loop over number of neighbors of molecule a


  /* ----------------------------------------- */
  /* Assign the accumulated force from         */
  /* nonbonded interactions to the force on    */
  /* the atom.                                 */
  /* ----------------------------------------- */
    ff[k][a].x += fxa;
    ff[k][a].y += fya;
    ff[k][a].z += fza;

  

  }//loop over molecules a

  /* ================================================================== */
  /*                                                                    */
  /* After the loop, if PRESSURE is defined, store the bonding          */
  /* contribution to the virial tensor to pvir[k].bond[i].              */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  pvir[k].nbond[0] = wxx;
  pvir[k].nbond[1] = wyx;
  pvir[k].nbond[2] = wyy;
  pvir[k].nbond[3] = wzx;
  pvir[k].nbond[4] = wzy;
  pvir[k].nbond[5] = wzz;
#endif //PRESSURE

  /* ================================================================== */
  /*                                                                    */
  /* After the loop, assign the accumulated energies to "energy" array. */
  /*	energy[0] = enonbond						*/
  /*	energy[1] = enonbond shifted					*/
  /*	energy[2] = EPP protein-protein nonbond shifted			*/
  /*	energy[3] = EPS protein-solute nonbond shifted			*/
  /*	energy[4] = ESS solute-solute nonbond shifted			*/
  /*                                                                    */
  /* ================================================================== */
  energy[0] = enonbond+ enonnative;                 // non-shifted energy (LJ and coulomibic)
  energy[1] = enonbond - enonbonds+ enonnative;     // shifted energy (LJ and coulomibic) + enonnative
  energy[2] = e_pp;
  energy[3] = e_ps;
  energy[4] = e_ss;
}
#else //ifndef DNA_GOLIK
/* ======================================================================== */
/* cnbnd_dnagolik(int, double*)                                             */
/*                                                                          */
/* This subroutine calculates the nonbonded native and non-native           */
/* energies of the Go-type model.                                           */
/*                                                                          */
/* Passed Parameters:                                                       */
/*			ibox:	    The box number                          */
/*			*energy:    A pointer to an array defined in        */
/*                                  forces() or force_long().  It contains  */
/*                                  three entries.  Two for unshifted and   */
/*                                  shifted vdw energies and one for        */
/*                                  coulombic energies.                     */
/*                                                                          */
/* ======================================================================== */
void cnbnd_dnagolik(int ibox, double *energy)
{

  int k = ibox;

  #ifdef EWALD
  double kappa = sim.kappa[k];
  double kpi = 1.0 * kappa / sqrt(PI);
  double eewaldr  = 0.0;
  double eewaldrs = 0.0;
  double eewaldi  = 0.0;
  #endif

  //struct atoms *atomk;
  //atomk = atom[k];
  //double fconv = (ELEC * ELEQ * ELEQ * NA * 0.001);   // convert atomic charges to units of kJ/mol
  /* ================================================================== */
  /* Zero out the energy accumulator                                    */
  /* ================================================================== */
  double enonbond	= 0.0;
  double enonbonds	= 0.0;
  double enonnative	= 0.0;
  double ecoulomb	= 0.0;
  double ebp		= 0.0;
  double emm		= 0.0;
  double esolv          = 0.0;
  /* ================================================================== */
  /* Set the parameters for the types of iteractions.                   */
  /* ================================================================== */
  double epsilon = golik[k].eps;
  
  #ifdef DRSAM  
  double eps_AT = 2.0*epsilon;
  double eps_GC = 2.532*epsilon;
  #else
  double eps_AT = 4.0*epsilon*2.0/3.0;
  double eps_GC = 4.0*epsilon;
  #endif
  //double eps_AT = 4.116*epsilon*0.75;
  //double eps_GC = 4.116*epsilon*16.0/17.0;

  double factor = pow(2.0,-(1.0/6.0));
  double sig_AT = 2.9002;
  double sig_GC = 2.8694;
  //double sig_non= factor*golik[k].d_repulsive;
  #ifdef MM_ATTRACT
  double sig_mm = (sig_AT+sig_GC)/2.0;
  double eps_mm = epsilon;
  #else
  double d_mm   = 1.0;
  double sig_mm = factor*d_mm;
  double eps_mm = epsilon;
  #endif
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
  #ifdef EWALD_DEBUG
  double wxx_ewaldr = 0.0;
  double wyx_ewaldr = 0.0;
  double wyy_ewaldr = 0.0;
  double wzx_ewaldr = 0.0;
  double wzy_ewaldr = 0.0;
  double wzz_ewaldr = 0.0;

  double wxx_ewaldi = 0.0;
  double wyx_ewaldi = 0.0;
  double wyy_ewaldi = 0.0;
  double wzx_ewaldi = 0.0;
  double wzy_ewaldi = 0.0;
  double wzz_ewaldi = 0.0;
  #endif
#endif//PRESSURE

  /* ================================================================== */
  /*                                                                    */
  /* Get and assign information for pbc and cutoffs                     */
  /*                                                                    */
  /* ================================================================== */
  double rc   = sim.rc;
  double rci  = 1.0/rc;
  double rc2  = box[k].rc2;
  double rc2i = 1.0/box[k].rc2;
  //double rc4  = box[k].rc4;
  //double rc_non = golik[k].d_repulsive;
  //double rc2_non = rc_non * rc_non;
  #ifdef RFC
  double k_rf = (sim.epsRF[k] - 1.0) / (2.0*sim.epsRF[k] + 1.0)*rc2i*rci;
  double c_rf  = rci+k_rf*rc2;
  #endif

  #ifndef MM_ATTRACT
  double rc_mm = d_mm;
  double rc2_mm = rc_mm * rc_mm;
  #endif

  #ifdef DRSAM
  double eps_inf = 0.504982*epsilon;
  double A_I = 0.474876*(1.0 + 1.0 / (0.148378 + 10.9553*sim.Na_con[k]));
  #endif

  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;

  /* ================================================================== */
  /*                                                                    */
  /* Calculate the energy and force on each molecule due to atoms in    */
  /* its neighborlist. The first loop is over sites/atoms 0 to N-1      */
  /* and the second is over the number of neighbors in atom a's         */
  /* neighborlist.                                                      */
  /*                                                                    */
  /* ================================================================== */

  for(int a=0; a<box[k].boxns-1; a++) {

  /* ----------------------------------------------- */
  /* Zero out the force accumulators for particle a. */
  /* ----------------------------------------------- */
    //double fxa = 0.0;
    //double fya = 0.0;
    //double fza = 0.0;

  /* ----------------------------------------------- */
  /* Assign position of particle a.                  */
  /* ----------------------------------------------- */    
    //double ax = atom[k][a].x;
    //double ay = atom[k][a].y;
    //double az = atom[k][a].z;

  /* ----------------------------------------------- */
  /* Loop around neighbors of particle a.            */
  /* ----------------------------------------------- */
    for(int bb=0; bb<nlist[k].count[a]; bb++) {
      int b = nlist[k].list[a][bb];
  /* ----------------------------------------------- */
  /* Calculate the x,y,z distances between           */
  /* particles a & b.                                */
  /* ----------------------------------------------- */
      double dx = atom[k][a].x - atom[k][b].x;
      double dy = atom[k][a].y - atom[k][b].y;
      double dz = atom[k][a].z - atom[k][b].z;

  /* ----------------------------------------------- */
  /* Apply minimum image convention.                 */
  /* ----------------------------------------------- */
      dx = MINIMG(dx,hx,lx);
      dy = MINIMG(dy,hy,ly);
      dz = MINIMG(dz,hz,lz);

  /* ----------------------------------------------- */
  /* Calculate the distance between particles a & b. */
  /* If the distance is less than the cutoff,        */
  /* calculated the energy and the force. If it is   */
  /* greater, go to the next pair in the loop.       */
  /* ----------------------------------------------- */
      double dr2  = dx*dx + dy*dy + dz*dz;
      double dr2i = 1.0/dr2;

      if(dr2 > rc2) continue;

  /* ----------------------------------------------- */
  /* Set the epsilon, sigma, and qq for the pair.    */
  /* ----------------------------------------------- */
      double sig, eps, eps4;
      double sr2, sr4, sr6, sr10, sr12;
      double sr2s, sr6s, sr10s, sr12s;
      #ifdef DRSAM      
      int Nres_a;
      int Nres_b;
      int N_bp;
      double eps_N;
      double eps_s;
      double s_exp;
      double s_factor2;
      #endif

      #ifdef MM_ATTRACT
      double sr4s;
      #endif
      double fx, fy, fz;
      int xflag=nlist[k].xflag[a][bb];    
      int ta = atom[k][a].atomid;
      int tb = atom[k][b].atomid;

      sig     = ljset[ta][tb].sig;
      double rc2_non = ljset[ta][tb].rc2_non;
      //epsilon = ljset[ta][tb].eps;
      int flag = ljset[ta][tb].id;

      if(xflag){  //flag is 0 if there is no exclusion
        if(xflag == 1) flag = 1;
        else if(xflag == 5) flag = 5;
        else if(xflag == 6) flag = 6;
        else if(xflag == 8) flag = 0;
      }
      
  /* ----------------------------------------------- */
  /* Calculate the LJ force and energy.              */
  /* ----------------------------------------------- */
      double force = 0.0; //double enbnd =0.0; double enbnds = 0.0;
      #if defined(PRESSURE) && defined(EWALD_DEBUG)
      double force_rewald = 0.0;
      double force_iewald = 0.0;
      #endif
      switch(flag){
        case 1: //native
          sig = nlist[k].sig[a][bb];
          eps4 = 4.0*epsilon;
          sr2  = (sig * sig) * dr2i;
          sr6  = sr2 * sr2 * sr2;
	  sr12 = sr6 * sr6;

          sr2s  = (sig * sig) * rc2i;
          sr6s  = sr2s * sr2s * sr2s;
          sr12s = sr6s * sr6s;

          force      = (12.0*eps4*dr2i) * (sr12 - 0.5*sr6);
          enonbond  += eps4 * (sr12 - sr6);
          enonbonds += eps4 * (sr12s -sr6s);

          break;
        case 2: //AT compliment
          sig = sig_AT;
          eps4 = 4.0*eps_AT;
          sr2  = (sig * sig) * dr2i;
          sr6  = sr2 * sr2 * sr2;
          sr12 = sr6 * sr6;
          sr10 = sr12/sr2;

          sr2s  = (sig * sig) * rc2i;
          sr6s  = sr2s * sr2s * sr2s;
          sr12s = sr6s * sr6s;
          sr10s = sr12s/sr2s;
        
          force      = (60.0*eps4*dr2i) * (sr12 - sr10);
          ebp  += eps4 * (5.0*sr12 - 6.0*sr10);
          ebp  -= eps4 * (5.0*sr12s- 6.0*sr10s);
        
          break;
        case 3: //CG compliment
          sig = sig_GC;
          eps4 = 4.0*eps_GC;
          sr2  = (sig * sig) * dr2i;
          sr6  = sr2 * sr2 * sr2;
          sr12 = sr6 * sr6;
          sr10 = sr12/sr2;

          sr2s  = (sig * sig) * rc2i;
          sr6s  = sr2s * sr2s * sr2s;
          sr12s = sr6s * sr6s;
          sr10s = sr12s/sr2s;
        
          force      = (60.0*eps4*dr2i) * (sr12 - sr10);
          ebp  += eps4 * (5.0*sr12 - 6.0*sr10);
          ebp  -= eps4 * (5.0*sr12s- 6.0*sr10s);

          break;
        case 4: //mis-match
          #ifdef MM_ATTRACT
          sig = sig_mm;
          eps = eps_mm;
          sr2  = (sig * sig) * dr2i;
          sr4 = sr2 * sr2;
          sr12 = sr4 * sr4 * sr4;
      
          sr2s  = (sig * sig) * rc2i
	  sr4s  = sr2s * sr2s;
          sr12s = sr4s * sr4s * sr4s;

          force        = (96.0*eps*dr2i) * (sr12 - sr4);
          enonnative  += 4.0*eps*(2.0*sr12 - 6.0*sr4);
          enonnative  -= 4.0*eps*(2.0*sr12s - 6.0*sr4s);
          #else
          if(dr2<rc2_mm){
            sig = sig_mm;
            eps = eps_mm;
            sr2  = (sig * sig) * dr2i;
            sr6  = sr2 * sr2 * sr2;
            sr12 = sr6 * sr6;

            force        = (48.0*eps*dr2i) * (sr12 - 0.5*sr6);
            enonnative  += 4.0*eps*(sr12 - sr6) + eps;
          }
          #endif
          break;
        case 5://bonded pairs
          break;
        case 6://native that is also a bend
          sig = nlist[k].sig[a][bb];
          eps4 = 4.0*epsilon;
          sr2  = (sig * sig) * dr2i;
          sr6  = sr2 * sr2 * sr2;
          sr12 = sr6 * sr6;

          sr2s  = (sig * sig) * rc2i;
          sr6s  = sr2s * sr2s * sr2s;
          sr12s = sr6s * sr6s;

          force      = (12.0*eps4*dr2i) * (sr12 - 0.5*sr6);
          enonbond  += eps4 * (sr12 - sr6);
          enonbonds += eps4 * (sr12s -sr6s);

          break;
        #ifdef DRSAM
        case 7://attraction between sugars on different strands
          //For finding the number of residues on the shorter strand
          Nres_a = mol[residue[k][atnopbc[k][a].n].chain].Nres;
          Nres_b = mol[residue[k][atnopbc[k][b].n].chain].Nres;
          if (Nres_a < Nres_b) N_bp = Nres_a;
          else N_bp = Nres_b;

          eps_N = eps_inf * (1.0 - 1.0 / (1.40418 - 0.268231*N_bp));
          eps_s = A_I * eps_N;
          s_exp = exp(-0.1875*(sqrt(dr2) - 13.38));//alpha_s = 0.1875, r_s = 13.38
          s_factor2 = (1 - s_exp)*(1 - s_exp);

          force = -2.0 * eps_s * 0.1875 * s_exp * (1 / sqrt(dr2)) * (1.0 - s_exp);//alpha_s = 0.1875
          esolv += eps_s * s_factor2 - eps_s;

          break;
        #endif //DRSAM

        default://non-native
          if(dr2<rc2_non){
            //sig = sig_non;
            eps = epsilon;
            sr2  = (sig * sig) * dr2i;
            sr6  = sr2 * sr2 * sr2;
            sr12 = sr6 * sr6;

            force        = (48.0*eps*dr2i) * (sr12 - 0.5*sr6);
            enonnative  += 4.0*eps*(sr12 - sr6) + eps;
          }
          break;
      }	 

  /* ----------------------------------------------- */
  /* Calculate the Coulombic energy and forces.      */
  /* ----------------------------------------------- */
      #if defined(COULOMB) || defined(SPME)

      double qq = pott[k][a].qq * pott[k][b].qq * QFACTOR;

      if(qq != 0.0){
        #ifdef SPME 
        qq = qq/sim.epsRF[k];
        double dr = sqrt(dr2);
        double kr = kappa * dr;
        double erfckr = erfc(kr);
        double kpidrexpkr2 =2.0*kpi*dr*exp(-kr * kr);
        #endif      
        switch(flag){
          case 5: //bond
            #ifdef SPME
            //Intramolecular Part
            ecoulomb -= qq*(1.0-erfckr)/dr;
            force -= (qq/dr/dr2) * ((1.0-erfckr)-kpidrexpkr2);
            #if defined(PRESSURE) && defined(EWALD_DEBUG)
            force_iewald += (qq/dr/dr2) * ((1.0-erfckr)-kpidrexpkr2);
            eewaldi += qq*(1.0-erfckr)/dr;
            //printf("case=5, a=%d, b=%d\n",a,b);
            #endif
            #endif

            break;
          case 6: //bend
            #ifdef SPME
            //Intramolecular Part
            ecoulomb -= qq*(1.0-erfckr)/dr;
            force -= (qq/dr/dr2) * ((1.0-erfckr)-kpidrexpkr2);
            #if defined(PRESSURE) && defined(EWALD_DEBUG)
            force_iewald += (qq/dr/dr2) * ((1.0-erfckr)-kpidrexpkr2);
            eewaldi += qq*(1.0-erfckr)/dr;
            //printf("case=6, a=%d, b=%d\n",a,b);
            #endif
            #endif

            break;
          default:
            #ifdef FSH
            qq = qq/sim.epsRF[k]*dr2i;
            double dr = sqrt(dr2);
            double dr2_rc2 = dr2/rc2;
            double term1   = 1 - dr2_rc2;
            double term2   = 1 - 2*dr/rc + dr2_rc2;
            force	  += qq / dr * term1;
            ecoulomb	  += qq * dr * term2; 

            #elif defined(DEBHUCK)
            qq = qq / sim.epsRF[k];
            double lambdai = 1.0/sim.kappa[k];
            double dri  = sqrt(dr2i);
            double qqerlambda = qq*exp(-lambdai/dri)*dri;
            force += qqerlambda*(dri*lambdai+dr2i);
            ecoulomb += qqerlambda;

            #elif defined(RFC)
            double dri = sqrt(dr2i);
            double k_rfdr2=k_rf*dr2;
            force	  += (qq*dr2i) * (dri - 2.0 * k_rfdr2);
            double term1 = qq * (dri + k_rfdr2);
            double term2 = qq * c_rf;
            ecoulomb  += term1 - term2;

            #elif defined(SPME)
            //qq = qq / sim.epsRF[k];
            double dri = 1.0 / dr;
            //Real Part
            //ecoulomb += qq * (erfckr * dri - erfc(kappa * rc) / rc);
            ecoulomb += qq * erfckr * dri;
            force += (qq*dri/dr2) * (erfckr + kpidrexpkr2);
            #if defined(PRESSURE) && defined(EWALD_DEBUG)
            force_rewald += (qq*dri/dr2) * (erfckr + kpidrexpkr2);
            //eewaldr += qq * (erfckr * dri - erfc(kappa * rc) / rc);
            eewaldr += qq * erfckr * dri;
            #endif
            #endif
        }
      }
      #endif
 
  /* ----------------------------------------- */
  /* Accumulate the forces.                    */
  /* ----------------------------------------- */
      if(force != 0.0){
 	fx	    = force * dx;
 	ff[k][a].x += fx;
 	ff[k][b].x -= fx;

 	fy	    = force * dy;
 	ff[k][a].y += fy;
 	ff[k][b].y -= fy;

 	fz	    = force * dz;
 	ff[k][a].z += fz;
 	ff[k][b].z -= fz;

 	#ifdef PRESSURE
 	double sxx = fx * dx;
 	double syx = fy * dx;
 	double syy = fy * dy; 
 	double szx = fz * dx;
 	double szy = fz * dy;
 	double szz = fz * dz;
	wxx += sxx;
	wyx += syx;
	wyy += syy;
	wzx += szx;
	wzy += szy;
	wzz += szz;
#ifdef ZHOU
	stresses[k][a][0][0] += 0.5*sxx;
	stresses[k][a][1][0] += 0.5*syx;
	stresses[k][a][1][1] += 0.5*syy;
	stresses[k][a][2][0] += 0.5*szx;
	stresses[k][a][2][1] += 0.5*szy;
	stresses[k][a][2][2] += 0.5*szz;

	stresses[k][b][0][0] += 0.5*sxx;
	stresses[k][b][1][0] += 0.5*syx;
	stresses[k][b][1][1] += 0.5*syy;
	stresses[k][b][2][0] += 0.5*szx;
	stresses[k][b][2][1] += 0.5*szy;
	stresses[k][b][2][2] += 0.5*szz;
#endif
 	  #ifdef EWALD_DEBUG
 	  double fx_ewaldr = force_rewald * dx;
 	  double fy_ewaldr = force_rewald * dy;
 	  double fz_ewaldr = force_rewald * dz;
 	  double fx_ewaldi = -force_iewald * dx;
 	  double fy_ewaldi = -force_iewald * dy;
 	  double fz_ewaldi = -force_iewald * dz;

 	  /* ///////////////////////////////// */
 	  /* Ewald Real part contribution      */
 	  /* ///////////////////////////////// */
 	  wxx_ewaldr += fx_ewaldr * dx;
 	  wyx_ewaldr += fy_ewaldr * dx;
 	  wyy_ewaldr += fy_ewaldr * dy;
 	  wzx_ewaldr += fz_ewaldr * dx;
 	  wzy_ewaldr += fz_ewaldr * dy;
 	  wzz_ewaldr += fz_ewaldr * dz;

 	  /* ///////////////////////////////// */
 	  /* Ewald Intramolecular contribution */
 	  /* ///////////////////////////////// */
 	  wxx_ewaldi += fx_ewaldi * dx;
 	  wyx_ewaldi += fy_ewaldi * dx;
 	  wyy_ewaldi += fy_ewaldi * dy;
 	  wzx_ewaldi += fz_ewaldi * dx;
 	  wzy_ewaldi += fz_ewaldi * dy;
 	  wzz_ewaldi += fz_ewaldi * dz;
 	  #endif
 	#endif//PRESSURE
      }//if(force!=0.0)
      
    }//loop over number of neighbors of molecule a

    //ff[k][a].x += fxa;
    //ff[k][a].y += fya;
    //ff[k][a].z += fza;
  }//loop over molecules a
  /* ------------------------------------------ */
  /* Call the subroutines to calucate the       */
  /* k-space and self energy contributions      */
  /* to the electrostatic energy.               */
  /*                                            */
  /* enkewald[0] = k-space ewald energy         */
  /* enkewald[1] = self correction ewald energy */
  /* enkewald[2] = k-space - self energy        */
  /* ------------------------------------------ */
  #ifdef SPME
  double enkewald[3];
  pme_calc(k,enkewald);
  ecoulomb += enkewald[2];
  #endif
  /* ================================================================== */
  /*                                                                    */
  /* After the loop, if PRESSURE is defined, store the bonding          */
  /* contribution to the virial tensor to pvir[k].bond[i].              */
  /*                                                                    */
  /* ================================================================== */
  #ifdef PRESSURE
  pvir[k].nbond[0] = wxx;
  pvir[k].nbond[1] = wyx;
  pvir[k].nbond[2] = wyy;
  pvir[k].nbond[3] = wzx;
  pvir[k].nbond[4] = wzy;
  pvir[k].nbond[5] = wzz;

    #ifdef EWALD_DEBUG
    pvir[k].ewald_real[0] = wxx_ewaldr;
    pvir[k].ewald_real[1] = wyx_ewaldr;
    pvir[k].ewald_real[2] = wyy_ewaldr;
    pvir[k].ewald_real[3] = wzx_ewaldr;
    pvir[k].ewald_real[4] = wzy_ewaldr;
    pvir[k].ewald_real[5] = wzz_ewaldr;

    pvir[k].ewald_intra[0] = wxx_ewaldi;
    pvir[k].ewald_intra[1] = wyx_ewaldi;
    pvir[k].ewald_intra[2] = wyy_ewaldi;
    pvir[k].ewald_intra[3] = wzx_ewaldi;
    pvir[k].ewald_intra[4] = wzy_ewaldi;
    pvir[k].ewald_intra[5] = wzz_ewaldi;
    #endif

  #endif //PRESSURE

  /* ================================================================== */
  /*                                                                    */
  /* After the loop, assign the accumulated energies to "energy" array. */
  /*	energy[0] = enonbond											*/
  /*	energy[1] = enonbond shifted									*/
  /*	energy[2] = EPP protein-protein nonbond shifted					*/
  /*	energy[3] = EPS protein-solute nonbond shifted					*/
  /*	energy[4] = ESS solute-solute nonbond shifted					*/
  /*                                                                    */
  /* ================================================================== */
  energy[0] = enonbond+enonnative+ebp+emm+ecoulomb+esolv;               // non-shifted energy (LJ and coulomibic)
  energy[1] = enonbond-enonbonds+enonnative+ebp+emm+ecoulomb+esolv;     // shifted energy (LJ and coulomibic) + enonnative
  energy[2] = enonbond-enonbonds;
  energy[3] = enonnative+emm;
  energy[4] = ebp;
  energy[5] = ecoulomb;
  energy[11]= esolv;
  #ifdef EWALD
  energy[6] = eewaldr;
  energy[7] = eewaldrs;
  energy[8] = enkewald[0];
  energy[9] = enkewald[1];
  energy[10]= eewaldi;
  #endif
}

#endif
#endif //GOLIK
