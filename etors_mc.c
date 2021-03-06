#ifdef MC

#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

double etors_mc (int ibox)
{

  int k = ibox;

  /* ================================================================== */
  /*                                                                    */
  /* Zero out the energy accumulator                                    */
  /*                                                                    */
  /* ================================================================== */
  double etorsion = 0.0;

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

  /* ================================================================== */
  /*                                                                    */
  /* Get and assign information for pbc.                                */
  /*                                                                    */
  /* ================================================================== */
  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;


  /* ================================================================== */
  /*                                                                    */
  /* Calculate the energy and force associated with each dihedral       */
  /* angle. The loop is over the number of dihedral angles.             */
  /*                                                                    */
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
    double fr = 0.0;
    double etors = 0.0;

	for(int tt =0; tt<6; tt++){
		if(tors[k][i].kphi[tt] !=0.0){
			fr		+=	tors[k][i].kphi[tt]*tors[k][i].nphi[tt]*(sin(tors[k][i].nphi[tt]*tau-tors[k][i].delphi[tt]));
			etors	+=  tors[k][i].kphi[tt]*(1+cos(tors[k][i].nphi[tt]*tau-tors[k][i].delphi[tt]));
		}
		else break;
	}
    etorsion += etors;

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
#endif
