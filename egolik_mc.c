#if defined(GOLIK) && defined(MC)
/* ======================================================================== */
/* egolik_mc.cpp                                                            */
/*																			*/
/* Written by Thomas A. Knotts IV, 1 Sep 2005.                              */
/* ======================================================================== */


/* ======================================================================== */
/* cbb_golik(int)                                                           */
/*		This subroutine calculates the backbone potential for the Go-type   */
/* molel.																	*/
/* ======================================================================== */

#include "defines.h"


double ebb_golik_mc (int ibox)
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
      if(dx >  hx) dx -= lx;
      else if(dx < -hx) dx += lx;
      if(dy >  hy) dy -= ly;
      else if(dy < -hy) dy += ly;
      if(dz >  hz) dz -= lz;
      else if(dz < -hz) dz += lz;

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

    ffox[k][ia] -= fx;
    ffoy[k][ia] -= fy;
    ffoz[k][ia] -= fz;
    ffox[k][ib] += fx;
    ffoy[k][ib] += fy;
    ffoz[k][ib] += fz;

  /* ----------------------------------------- */
  /* If PRESSURE is defined, accumulate the    */
  /* components of the virial tensor.          */
  /* ----------------------------------------- */
#ifdef PRESSURE
    wxx += fx * dx;
    wyx += fy * dx;
    wyy += fy * dy;
    wzx += fz * dx;
    wzy += fz * dy;
    wzz += fz * dz;
#endif

  }

  /* ================================================================== */
  /*                                                                    */
  /* After the loop, if PRESSURE is defined, store the bonding          */
  /* contribution to the virial tensor to pvir[k].bond[i].              */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  pviro[k].bond[0] = wxx;
  pviro[k].bond[1] = wyx;
  pviro[k].bond[2] = wyy;
  pviro[k].bond[3] = wzx;
  pviro[k].bond[4] = wzy;
  pviro[k].bond[5] = wzz;
#endif

  return (ebond);

}

/* ======================================================================== */
/* cba_golik(int)                                                           */
/*		This subroutine calculates the backbone potential for the Go-type   */
/* molel.																	*/
/* ======================================================================== */
double eba_golik_mc (int ibox)
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

#ifndef DNA_GOLIK
  double k_theta = -10.0*golik[k].eps;  //10.0=20/2
#else
  double k_theta = -10.0*golik[k].eps*20;  //10.0=20/2
  //double k_theta = -10.0*golik[k].eps;  //10.0=20/2
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
    double phi = acos(cosphi);
    double isinphi = 1.0 / sqrt(1.0 - cosphi*cosphi);

    double rr = phi - bend[k][i].angeq;
	double fr = k_theta * rr*2.0;
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


    ffox[k][ia] += fxa;
    ffoy[k][ia] += fya;
    ffoz[k][ia] += fza;
    ffox[k][ib] += -(fxa + fxc);
    ffoy[k][ib] += -(fya + fyc);
    ffoz[k][ib] += -(fza + fzc);
    ffox[k][ic] += fxc;
    ffoy[k][ic] += fyc;
    ffoz[k][ic] += fzc;

  /* ----------------------------------------- */
  /* If PRESSURE is defined, accumulate the    */
  /* components of the virial tensor.          */
  /* ----------------------------------------- */
#ifdef PRESSURE
    wxx += (fxa * dax + fxc * dcx);
    wyx += (fya * dax + fyc * dcx);
    wyy += (fya * day + fyc * dcy);
    wzx += (fza * dax + fzc * dcx);
    wzy += (fza * day + fzc * dcy);
    wzz += (fza * daz + fzc * dcz);
#endif

  }

  /* ================================================================== */
  /*                                                                    */
  /* After the loop, if PRESSURE is defined, store the bending          */
  /* contribution to the virial tensor to pvir[k].bend[i].              */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  pviro[k].bend[0] = wxx;
  pviro[k].bend[1] = wyx;
  pviro[k].bend[2] = wyy;
  pviro[k].bend[3] = wzx;
  pviro[k].bend[4] = wzy;
  pviro[k].bend[5] = wzz;
#endif

  return (ebend);

}

/* ======================================================================== */
/* cda_golik(int)                                                           */
/*		This subroutine calculates the backbone potential for the Go-type   */
/* molel.  Dihedrals.														*/
/* ======================================================================== */
/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

double eda_golik_mc (int ibox)
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

#ifndef DNA_GOLIK
  double k_phi = 0.2 * golik[k].eps;
#else
  double k_phi = 0.2 * golik[k].eps*20;
  //double k_phi = 0.2 * golik[k].eps;
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
    double tau = signum * acos(costau);

  /* -----------------------------------------	*/
  /* Assign the values for phi and psi			*/
  /* -----------------------------------------	*/
	if(tors[k][i].phitag ==1)		tors[k][i].phi = tau;
	else if (tors[k][i].psitag ==1) tors[k][i].psi = tau;
	else if (tors[k][i].thetag ==1) tors[k][i].theta = tau;
  /* ----------------------------------------- */
  /* Calculate the force and energy.           */
  /* ----------------------------------------- */
#ifdef DNA_GOLIK
	double fr	 =	k_phi*(sin(tau-tors[k][i].delphi[0]));
	etorsion	+=  k_phi*(1+cos(tau-tors[k][i].delphi[0]));
#else
	double fr    =   k_phi*3.0*(sin(3.0*tau-tors[k][i].delphi[0]));
    etorsion    +=  k_phi*(1+cos(3.0*tau-tors[k][i].delphi[0]));
#endif
  /* ----------------------------------------- */
  /* Transform the forces to cartesian         */
  /* coordinates and accumulate.               */
  /* ----------------------------------------- */
    double fa    = -fr * irab * isinb2;
    double fax   = fa * abbcx;
    double fay   = fa * abbcy;
    double faz   = fa * abbcz;
    ffox[k][ia] += fax;
    ffoy[k][ia] += fay;
    ffoz[k][ia] += faz;

    double fb1   = fr * (rbc - rab*cosb) * irab *irbc * isinb2;
    double fb2   = fr * cosc * irbc * isinc2;
    ffox[k][ib] += (fb1 * abbcx + fb2 * dccbx);
    ffoy[k][ib] += (fb1 * abbcy + fb2 * dccby);
    ffoz[k][ib] += (fb1 * abbcz + fb2 * dccbz);

    double fc1   = fr * (rbc - rcd*cosc) * ircd *irbc *isinc2;
    double fc2   = fr * cosb * irbc * isinb2;
    double fcx   = fc1 * dccbx + fc2 * abbcx;
    double fcy   = fc1 * dccby + fc2 * abbcy;
    double fcz   = fc1 * dccbz + fc2 * abbcz;
    ffox[k][ic] += fcx;
    ffoy[k][ic] += fcy;
    ffoz[k][ic] += fcz;
    
    double fd    = -fr *ircd * isinc2;
    double fdx   = fd * dccbx;
    double fdy   = fd * dccby;
    double fdz   = fd * dccbz;
    ffox[k][id] += fdx;
    ffoy[k][id] += fdy;
    ffoz[k][id] += fdz;

  /* ----------------------------------------- */
  /* If PRESSURE is defined, accumulate the    */
  /* components of the virial tensor.          */
  /* ----------------------------------------- */
#ifdef PRESSURE
    wxx += ( - fax * dabx + fcx * dbcx + fdx *(dbcx+dcdx) ); 
    wyx += ( - fay * dabx + fcy * dbcx + fdy *(dbcx+dcdx) ); 
    wyy += ( - fay * daby + fcy * dbcy + fdy *(dbcy+dcdy) ); 
    wzx += ( - faz * dabx + fcz * dbcx + fdz *(dbcx+dcdx) ); 
    wzy += ( - faz * daby + fcz * dbcy + fdz *(dbcy+dcdy) ); 
    wzz += ( - faz * dabz + fcz * dbcz + fdz *(dbcz+dcdz) ); 
#endif

  }
  /* ================================================================== */
  /*                                                                    */
  /* After the loop, if PRESSURE is defined, store the bonding          */
  /* contribution to the virial tensor to pvir[k].bond[i].              */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  pviro[k].tors[0] = wxx;
  pviro[k].tors[1] = wyx;
  pviro[k].tors[2] = wyy;
  pviro[k].tors[3] = wzx;
  pviro[k].tors[4] = wzy;
  pviro[k].tors[5] = wzz;
#endif
  return (etorsion);
}

#ifndef DNA_GOLIK
/* ======================================================================== */
/* cnbnd_golik(int, double*)                                                */
/*                                                                          */
/*		This subroutine calculates the nonbonded native and non-native      */
/* energies of the Go-type model.                                           */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                      *energy:	A pointer to an array defined in        */
/*                                  forces() or force_long().  It contains  */
/*                                  three entries.  Two for unshifted and   */
/*                                  shifted vdw energies and one for        */
/*                                  coulombic energies.                     */
/*                                                                          */
/* ======================================================================== */

void enbnd_golik_mc(int lb, int ub, int ibox, double *energy)
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

  for(int a=0; a<ub; a++) {

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
      if(dx >  hx) dx -= lx;
      else if(dx < -hx) dx += lx;
      if(dy >  hy) dy -= ly;
      else if(dy < -hy) dy += ly;
      if(dz >  hz) dz -= lz;
      else if(dz < -hz) dz += lz;

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
      ffox[k][b] -= fx;

       /* *********************************************** */
       /* y-forces                                        */
       /* *********************************************** */
      fy          = force * dy;
      fya        += fy;
      ffoy[k][b] -= fy;

       /* *********************************************** */
       /* z-forces                                        */
       /* *********************************************** */
      fz          = force * dz;
      fza        += fz;
      ffoz[k][b] -= fz;
  /* ----------------------------------------- */
  /* If PRESSURE is defined, accumulate the    */
  /* components of the virial tensor.          */
  /* ----------------------------------------- */
#ifdef PRESSURE
      wxx += fx * dx;
      wyx += fy * dx;
      wyy += fy * dy; 
      wzx += fz * dx;
      wzy += fz * dy;
      wzz += fz * dz;
#endif//PRESSURE
	}//if(force!=0.0)
    }//loop over number of neighbors of molecule a


  /* ----------------------------------------- */
  /* Assign the accumulated force from         */
  /* nonbonded interactions to the force on    */
  /* the atom.                                 */
  /* ----------------------------------------- */
    ffox[k][a] += fxa;
    ffoy[k][a] += fya;
    ffoz[k][a] += fza;

  

  }//loop over molecules a

  /* ================================================================== */
  /*                                                                    */
  /* After the loop, if PRESSURE is defined, store the bonding          */
  /* contribution to the virial tensor to pvir[k].bond[i].              */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  pviro[k].nbond[0] = wxx;
  pviro[k].nbond[1] = wyx;
  pviro[k].nbond[2] = wyy;
  pviro[k].nbond[3] = wzx;
  pviro[k].nbond[4] = wzy;
  pviro[k].nbond[5] = wzz;
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
/*		This subroutine calculates the nonbonded native and non-native        */
/* energies of the Go-type model.                                           */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                                        */
/*                      *energy:	A pointer to an array defined in          */
/*                                  forces() or force_long().  It contains  */
/*                                  three entries.  Two for unshifted and   */
/*                                  shifted vdw energies and one for        */
/*                                  coulombic energies.                     */
/*                                                                          */
/* ======================================================================== */
void enbnd_dnagolik_mc(int lb, int ub, int ibox, double *energy)
{

  int k = ibox;
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

  /* ================================================================== */
  /* Set the parameters for the types of iteractions.                   */
  /* ================================================================== */
  double epsilon = golik[k].eps;
  double eps_AT = epsilon*2.0/3.0*4.0;
  double eps_GC = epsilon*4.0;

  double factor = pow(2.0,-(1.0/6.0));
  double sig_AT = 2.9002;
  double sig_GC = 2.8694;
  double sig_non= factor*golik[k].d_repulsive;
  #ifdef MM_ATTRACT
  double sig_mm = (sig_AT+sig_GC)/2.0;
  double eps_mm = epsilon/16.0;
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
#endif//PRESSURE

  /* ================================================================== */
  /*                                                                    */
  /* Get and assign information for pbc and cutoffs                     */
  /*                                                                    */
  /* ================================================================== */
#ifdef RFC
  double rfcs = (sim.epsRF[k] - 1.0) / (2.0*sim.epsRF[k] + 1.0);
  double rfc  = rfcs / (rc2*rc);
#endif
  double rc   = sim.rc;
  double rc2  = box[k].rc2;
  double rc2i = 1.0/box[k].rc2;
  double rc4  = box[k].rc4;
  double rc_non = golik[k].d_repulsive;
  double rc2_non = rc_non * rc_non;
#ifndef MM_ATTRACT
  double rc_mm = d_mm;
  double rc2_mm = rc_mm * rc_mm;
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

  for(int a=0; a<ub; a++) {

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
      if(dx >  hx) dx -= lx;
      else if(dx < -hx) dx += lx;
      if(dy >  hy) dy -= ly;
      else if(dy < -hy) dy += ly;
      if(dz >  hz) dz -= lz;
      else if(dz < -hz) dz += lz;

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
	  double sig;
    double eps;
	  double eps4;
    double sr2;
    double sr6;
    double sr12;
    double sr10;
    double sr2s;
    #ifdef MM_ATTRACT
    double sr4s;
    #endif
    double sr6s;
    double sr12s;
    double sr10s;
	  double fx, fy, fz;
	  int xflag=nlist[k].xflag[a][bb];    
    int ta = atom[k][a].atomid;
	  int tb = atom[k][b].atomid;

    //sig     = ljset[ta][tb].sig;
		//epsilon = ljset[ta][tb].eps;
    int flag = ljset[ta][tb].id;

    if(xflag){  //flag is 0 if there is no exclusion
      if(xflag == 1) flag = 1;
      else if(xflag == 5) flag = 5;
      else if(xflag == 6) flag = 6;

    }

      
  /* ----------------------------------------------- */
  /* Calculate the LJ force and energy.              */
  /* ----------------------------------------------- */
	  double force=0.0; //double enbnd =0.0; double enbnds = 0.0;
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
      

		    sr2s  = (sig * sig) * rc2i;
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
      default://non-native
        if(dr2<rc2_non){
          sig = sig_non;
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
  /* This is a Debye-Huckel screened interaction.    */
  /* ----------------------------------------------- */
#ifdef COULOMB
  /* This was distance-dependent, force-shifted */
/*    switch(flag){
      case 5: //bond
        break;
      case 6: //bend
        break;
      default:
		    double qq = pott[k][a].qq * pott[k][b].qq * QFACTOR;
        if(qq != 0.0){
		      double dri4	= 1.0*dr2i*dr2i;
		      double terma	= qq/sim.epsRF[k];
		      double termb	= terma*dr2i;
		      double termc	= terma*(2.0/rc2-dr2/rc4);
		      double term3	= 2.0*terma*(dri4-1.0/rc4);
		      force		     += term3;
		      ecoulomb	   += termb-termc;
        }
*/
    switch(flag){
      case 5: //bond
        break;
      case 6: //bend
        break;
      default:
        double qq = pott[k][a].qq * pott[k][b].qq * QFACTOR;
        if(qq != 0.0){
	#ifdef DEBHUCK
	  qq = qq / sim.epsRF[k];
          double lambdai = 1.0/sim.kappa[k];
          double dri  = sqrt(dr2i);
          //double dri3 = dri*dr2i;
          double qqerlambda = qq*exp(-lambdai/dri)*dri;
          force += qqerlambda*(dri*lambdai+dr2i);
          ecoulomb += qqerlambda;
	#elif defined(RFC)
          double dri = sqrt(dr2i);
	  double rfcdr2=rfc*dr2;
          force	    += (qq/dr2) * (dri - 2.0 * rfcdr2);
	  double term1 = qq * (dri + rfcdr2);
	  double term2 = qq * (1.0 + rfcs) / rc;
          enonbond  += term1;
          enonbonds += term2;
          ecoulomb  += term1 - term2;
        #endif
        }
    }
#endif
 
  /* ----------------------------------------- */
  /* Accumulate the forces.                    */
  /* ----------------------------------------- */
	if(force != 0.0){
       /* *********************************************** */
       /* x-forces                                        */
       /* *********************************************** */
      fx          = force * dx;
      fxa        += fx;
      ffox[k][b] -= fx;

       /* *********************************************** */
       /* y-forces                                        */
       /* *********************************************** */
      fy          = force * dy;
      fya        += fy;
      ffoy[k][b] -= fy;

       /* *********************************************** */
       /* z-forces                                        */
       /* *********************************************** */
      fz          = force * dz;
      fza        += fz;
      ffoz[k][b] -= fz;
  /* ----------------------------------------- */
  /* If PRESSURE is defined, accumulate the    */
  /* components of the virial tensor.          */
  /* ----------------------------------------- */
#ifdef PRESSURE
      wxx += fx * dx;
      wyx += fy * dx;
      wyy += fy * dy; 
      wzx += fz * dx;
      wzy += fz * dy;
      wzz += fz * dz;
#endif//PRESSURE
	}//if(force!=0.0)
    }//loop over number of neighbors of molecule a


  /* ----------------------------------------- */
  /* Assign the accumulated force from         */
  /* nonbonded interactions to the force on    */
  /* the atom.                                 */
  /* ----------------------------------------- */
    ffox[k][a] += fxa;
    ffoy[k][a] += fya;
    ffoz[k][a] += fza;

  

  }//loop over molecules a

  /* ================================================================== */
  /*                                                                    */
  /* After the loop, if PRESSURE is defined, store the bonding          */
  /* contribution to the virial tensor to pvir[k].bond[i].              */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  pviro[k].nbond[0] = wxx;
  pviro[k].nbond[1] = wyx;
  pviro[k].nbond[2] = wyy;
  pviro[k].nbond[3] = wzx;
  pviro[k].nbond[4] = wzy;
  pviro[k].nbond[5] = wzz;
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
  energy[0] = enonbond+enonnative+ebp+emm+ecoulomb;               // non-shifted energy (LJ and coulomibic)
  energy[1] = enonbond-enonbonds+enonnative+ebp+emm+ecoulomb;     // shifted energy (LJ and coulomibic) + enonnative
  energy[2] = enonbond-enonbonds;
  energy[3] = enonnative+emm;
  energy[4] = ebp;
  energy[5] = ecoulomb;
}

#endif
#endif //GOLIK
