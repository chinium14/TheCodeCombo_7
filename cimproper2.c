/* ======================================================================== */
/* cimproper2.cpp                                                           */
/*                                                                          */
/*		This subroutine calculates the analytical second derivative of the  */
/* improper torsion energies and											*/
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
#ifdef CONFIGT
#ifdef PR_NPT
void min_image_npt_full(int, double*, double*, double*);
#endif

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

double cimproper2 (int ibox)
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
    double xba = atom[k][ib].x - atom[k][ia].x;
    double yba = atom[k][ib].y - atom[k][ia].y;
    double zba = atom[k][ib].z - atom[k][ia].z;
    double xcb = atom[k][ic].x - atom[k][ib].x;
    double ycb = atom[k][ic].y - atom[k][ib].y;
    double zcb = atom[k][ic].z - atom[k][ib].z;
    double xdc = atom[k][id].x - atom[k][ic].x;
    double ydc = atom[k][id].y - atom[k][ic].y;
    double zdc = atom[k][id].z - atom[k][ic].z;

  /* ----------------------------------------- */
  /* Apply minimum image convention.           */
  /* ----------------------------------------- */
#ifndef PR_NPT
    if(xba >  hx) xba -= lx;
    else if(xba < -hx) xba += lx;
    if(yba >  hy) yba -= ly;
    else if(yba < -hy) yba += ly;
    if(zba >  hz) zba -= lz;
    else if(zba < -hz) zba += lz;

    if(xcb >  hx) xcb -= lx;
    else if(xcb < -hx) xcb += lx;
    if(ycb >  hy) ycb -= ly;
    else if(ycb < -hy) ycb += ly;
    if(zcb >  hz) zcb -= lz;
    else if(zcb < -hz) zcb += lz;
    
    if(xdc >  hx) xdc -= lx;
    else if(xdc < -hx) xdc += lx;
    if(ydc >  hy) ydc -= ly;
    else if(ydc < -hy) ydc += ly;
    if(zdc >  hz) zdc -= lz;
    else if(zdc < -hz) zdc += lz;
#endif
#ifdef PR_NPT
	min_image_npt_full(k, &xba, &yba, &zba);
	min_image_npt_full(k, &xcb, &ycb, &zcb);
	min_image_npt_full(k, &xdc, &ydc, &zdc);
#endif
	double xt = yba*zcb - ycb*zba;
    double yt = zba*xcb - zcb*xba;
    double zt = xba*ycb - xcb*yba;
    double xu = ycb*zdc - ydc*zcb;
    double yu = zcb*xdc - zdc*xcb;
    double zu = xcb*ydc - xdc*ycb;
    double xtu = yt*zu - yu*zt;
    double ytu = zt*xu - zu*xt;
    double ztu = xt*yu - xu*yt;
    double rt2 = xt*xt + yt*yt + zt*zt;
    double ru2 = xu*xu + yu*yu + zu*zu;
    double rtru = sqrt(rt2 * ru2);

  /* ----------------------------------------- */
  /* Calculate the bond lengths and unit       */
  /* vectors involved in the dihedral angle.   */
  /* ----------------------------------------- */
    double rba   = sqrt(xba*xba + yba*yba + zba*zba);
    double irab  = 1.0 / rba;
    double erabx = xba * irab;
    double eraby = yba * irab;
    double erabz = zba * irab;
    double rcb   = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
    double irbc  = 1.0 / rcb;
    double erbcx = xcb * irbc;
    double erbcy = ycb * irbc;
    double erbcz = zcb * irbc;
    double rdc   = sqrt(xdc*xdc + ydc*ydc + zdc*zdc);
    double ircd  = 1.0 / rdc;
    double ercdx = xdc * ircd;
    double ercdy = ydc * ircd;
    double ercdz = zdc * ircd;

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
	
	double dedphi   = -fr;
	double d2edphi2 = 2.0*impr[k][i].kimpr;
	 
//     abbreviations for first derivative chain rule terms

	double xca = atom[k][ic].x - atom[k][ia].x;
    double yca = atom[k][ic].y - atom[k][ia].y;
    double zca = atom[k][ic].z - atom[k][ia].z;
	double xdb = atom[k][id].x - atom[k][ib].x;
    double ydb = atom[k][id].y - atom[k][ib].y;
    double zdb = atom[k][id].z - atom[k][ib].z;
#ifndef PR_NPT
    if(xca >  hx) xca -= lx;
    else if(xca < -hx) xca += lx;
    if(yca >  hy) yca -= ly;
    else if(yca < -hy) yca += ly;
    if(zca >  hz) zca -= lz;
    else if(zca < -hz) zca += lz;

    if(xdb >  hx) xdb -= lx;
    else if(xdb < -hx) xdb += lx;
    if(ydb >  hy) ydb -= ly;
    else if(ydb < -hy) ydb += ly;
    if(zdb >  hz) zdb -= lz;
    else if(zdb < -hz) zdb += lz;
#endif
	
    double dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb);
    double dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb);
    double dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb);
    double dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb);
    double dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb);
    double dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb);

//     abbreviations for second derivative chain rule terms

    double xycb2 = xcb*xcb + ycb*ycb;
    double xzcb2 = xcb*xcb + zcb*zcb;
    double yzcb2 = ycb*ycb + zcb*zcb;
    double rcbxt = -2.0 * rcb * dphidxt;
    double rcbyt = -2.0 * rcb * dphidyt;
    double rcbzt = -2.0 * rcb * dphidzt;
    double rcbt2 = rcb * rt2;
    double rcbxu = 2.0 * rcb * dphidxu;
    double rcbyu = 2.0 * rcb * dphidyu;
    double rcbzu = 2.0 * rcb * dphidzu;
    double rcbu2 = rcb * ru2;
    double dphidxibt = yca*dphidzt - zca*dphidyt;
    double dphidxibu = zdc*dphidyu - ydc*dphidzu;
    double dphidyibt = zca*dphidxt - xca*dphidzt;
    double dphidyibu = xdc*dphidzu - zdc*dphidxu;
    double dphidzibt = xca*dphidyt - yca*dphidxt;
    double dphidzibu = ydc*dphidxu - xdc*dphidyu;
    double dphidxict = zba*dphidyt - yba*dphidzt;
    double dphidxicu = ydb*dphidzu - zdb*dphidyu;
    double dphidyict = xba*dphidzt - zba*dphidxt;
    double dphidyicu = zdb*dphidxu - xdb*dphidzu;
    double dphidzict = yba*dphidxt - xba*dphidyt;
    double dphidzicu = xdb*dphidyu - ydb*dphidxu;

//     chain rule terms for first derivative components
	
    double dphidxia = zcb*dphidyt - ycb*dphidzt;
    double dphidyia = xcb*dphidzt - zcb*dphidxt;
    double dphidzia = ycb*dphidxt - xcb*dphidyt;
    double dphidxib = dphidxibt + dphidxibu;
    double dphidyib = dphidyibt + dphidyibu;
    double dphidzib = dphidzibt + dphidzibu;
    double dphidxic = dphidxict + dphidxicu;
    double dphidyic = dphidyict + dphidyicu;
    double dphidzic = dphidzict + dphidzicu;
    double dphidxid = zcb*dphidyu - ycb*dphidzu;
    double dphidyid = xcb*dphidzu - zcb*dphidxu;
    double dphidzid = ycb*dphidxu - xcb*dphidyu;

//     chain rule terms for second derivative components

	double dxiaxia = rcbxt*dphidxia;
    double dxiaxic = rcbxt*dphidxict + xcb*xt/rcbt2;
    double dxiaxid = 0.0;
    double dyiayia = rcbyt*dphidyia;
    double dyiayic = rcbyt*dphidyict + ycb*yt/rcbt2;
    double dyiayid = 0.0;
    double dziazia = rcbzt*dphidzia;
    double dziazic = rcbzt*dphidzict + zcb*zt/rcbt2;
    double dziazid = 0.0;
    double dxibxic = -xcb*dphidxib/(rcb*rcb) - (yca*(zba*xcb+yt) - zca*(yba*xcb-zt))/rcbt2 - 2.0*(yt*zba-yba*zt)*dphidxibt/rt2
					 - (zdc*(ydb*xcb+zu)-ydc*(zdb*xcb-yu))/rcbu2 + 2.0*(yu*zdb-ydb*zu)*dphidxibu/ru2;
    double dxibxid = rcbxu*dphidxibu + xcb*xu/rcbu2;
    double dyibyic = -ycb*dphidyib/(rcb*rcb) - (zca*(xba*ycb+zt)-xca*(zba*ycb-xt))/rcbt2 - 2.0*(zt*xba-zba*xt)*dphidyibt/rt2
					 - (xdc*(zdb*ycb+xu)-zdc*(xdb*ycb-zu))/rcbu2 + 2.0*(zu*xdb-zdb*xu)*dphidyibu/ru2;
    double dyibyid = rcbyu*dphidyibu + ycb*yu/rcbu2;
    double dzibzic = -zcb*dphidzib/(rcb*rcb) - (xca*(yba*zcb+xt)-yca*(xba*zcb-yt))/rcbt2 - 2.0*(xt*yba-xba*yt)*dphidzibt/rt2
					 - (ydc*(xdb*zcb+yu)-xdc*(ydb*zcb-xu))/rcbu2 + 2.0*(xu*ydb-xdb*yu)*dphidzibu/ru2;
    double dzibzid = rcbzu*dphidzibu + zcb*zu/rcbu2;
    double dxicxid = rcbxu*dphidxicu - xcb*(zdb*ycb-ydb*zcb)/rcbu2;
    double dyicyid = rcbyu*dphidyicu - ycb*(xdb*zcb-zdb*xcb)/rcbu2;
    double dziczid = rcbzu*dphidzicu - zcb*(ydb*xcb-xdb*ycb)/rcbu2;
    double dxidxid = rcbxu*dphidxid;
    double dyidyid = rcbyu*dphidyid;
    double dzidzid = rcbzu*dphidzid;

//     get some second derivative chain rule terms by difference
    
    double dxiaxib = -dxiaxia - dxiaxic - dxiaxid;
    double dyiayib = -dyiayia - dyiayic - dyiayid;
    double dziazib = -dziazia - dziazic - dziazid;
    double dxibxib = -dxiaxib - dxibxic - dxibxid;
    double dyibyib = -dyiayib - dyibyic - dyibyid;
    double dzibzib = -dziazib - dzibzic - dzibzid;
    double dxicxic = -dxiaxic - dxibxic - dxicxid;
    double dyicyic = -dyiayic - dyibyic - dyicyid;
    double dziczic = -dziazic - dzibzic - dziczid;

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

    double fb1   = fr * (rcb - rba*cosb) * irab *irbc * isinb2;
    double fb2   = fr * cosc * irbc * isinc2;
    ff[k][ib].x += (fb1 * abbcx + fb2 * dccbx);
    ff[k][ib].y += (fb1 * abbcy + fb2 * dccby);
    ff[k][ib].z += (fb1 * abbcz + fb2 * dccbz);

    double fc1   = fr * (rcb - rdc*cosc) * ircd *irbc *isinc2;
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

	double hessxa  = dedphi*dxiaxia + d2edphi2*dphidxia*dphidxia;
	double hessxb  = dedphi*dxibxib + d2edphi2*dphidxib*dphidxib;
	double hessxc  = dedphi*dxicxic + d2edphi2*dphidxic*dphidxic;
	double hessxd  = dedphi*dxidxid + d2edphi2*dphidxid*dphidxid;
	double hessya  = dedphi*dyiayia + d2edphi2*dphidyia*dphidyia;
	double hessyb  = dedphi*dyibyib + d2edphi2*dphidyib*dphidyib;
	double hessyc  = dedphi*dyicyic + d2edphi2*dphidyic*dphidyic;
	double hessyd  = dedphi*dyidyid + d2edphi2*dphidyid*dphidyid;
	double hessza  = dedphi*dziazia + d2edphi2*dphidzia*dphidzia;
	double hesszb  = dedphi*dzibzib + d2edphi2*dphidzib*dphidzib;
	double hesszc  = dedphi*dziczic + d2edphi2*dphidzic*dphidzic;
	double hesszd  = dedphi*dzidzid + d2edphi2*dphidzid*dphidzid;
	

	config[k].hesx	+= hessxa + hessxb + hessxc + hessxd;
	config[k].hesy	+= hessya + hessyb + hessyc + hessyd;
	config[k].hesz	+= hessza + hesszb + hesszc + hesszd;


  /* ----------------------------------------- */
  /* If PRESSURE is defined, accumulate the    */
  /* components of the virial tensor.          */
  /* ----------------------------------------- */
#ifdef PRESSURE
    double sxx = ( - fax * xba + fcx * xcb + fdx *(xcb+xdc) ); 
    double syx = ( - fay * xba + fcy * xcb + fdy *(xcb+xdc) ); 
    double syy = ( - fay * yba + fcy * ycb + fdy *(ycb+ydc) ); 
    double szx = ( - faz * xba + fcz * xcb + fdz *(xcb+xdc) ); 
    double szy = ( - faz * yba + fcz * ycb + fdz *(ycb+ydc) ); 
    double szz = ( - faz * zba + fcz * zcb + fdz *(zcb+zdc) ); 
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
#endif

