#if defined(CONFIGT) && defined(GOLIK) && defined(MC)
/* ======================================================================== */
/* egolik2_mc.cpp                                                           */
/*                                                                          */
/* Written by Thomas A. Knotts IV, 1 Sep 2005                               */
/* ======================================================================== */


/* ======================================================================== */
/* cbb_golik(int)                                                           */
/*		This subroutine calculates the backbone potential for the Go-type     */
/* model.																	                                  */
/* ======================================================================== */

#include "defines.h"


double ebb_golik2_mc (int ibox)
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
    double dx2 = dx*dx;
    double dy2 = dy*dy;
    double dz2 = dz*dz;
    double drr = dx2 + dy2 + dz2;
    double dr  = sqrt(drr);
    double dri = 1.0 / dr;

    double rr = dr - req;
    double rr2= rr * rr;
    double fr1 = k1*rr*2.0;
    double fr2 = k2*rr*rr2*4.0;
    ebond += -rr*(fr1*0.5 + fr2*0.25);
    double fr = fr1 + fr2;

  /*------------------------------------------ */
  /* Calculate the contribution of restraints  */
  /* to diagonal of hessian		        			   */
  /*------------------------------------------ */
    double term1 = (-fr1-fr2)*dri;
    double term2 = -2.0*k1;
    double term3 = -12.0*k2*rr2;
    double term4 = term2 + term3 - term1;
    double term5 = term4/drr;
    double d2e11 = term1 + term5*dx2;
    double d2e22 = term1 + term5*dy2;
    double d2e33 = term1 + term5*dz2;
	  hesox	   += 2.0* d2e11;
	  hesoy	   += 2.0* d2e22;
	  hesoz	   += 2.0* d2e33;

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
double eba_golik2_mc (int ibox)
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
  //double k_theta = -10.0*golik[k].eps*20;  //10.0=20/2
  double k_theta = -10.0*golik[k].eps;  //10.0=20/2
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
    double xab = atom[k][ia].x - atom[k][ib].x;
    double yab = atom[k][ia].y - atom[k][ib].y;
    double zab = atom[k][ia].z - atom[k][ib].z;
    double xcb = atom[k][ic].x - atom[k][ib].x;
    double ycb = atom[k][ic].y - atom[k][ib].y;
    double zcb = atom[k][ic].z - atom[k][ib].z;

  /* ----------------------------------------- */
  /* Apply minimum image convention.           */
  /* ----------------------------------------- */
    if(xab >  hx) xab -= lx;
    else if(xab < -hx) xab += lx;
    if(yab >  hy) yab -= ly;
    else if(yab < -hy) yab += ly;
    if(zab >  hz) zab -= lz;
    else if(zab < -hz) zab += lz;


    if(xcb >  hx) xcb -= lx;
    else if(xcb < -hx) xcb += lx;
    if(ycb >  hy) ycb -= ly;
    else if(ycb < -hy) ycb += ly;
    if(zcb >  hz) zcb -= lz;
    else if(zcb < -hz) zcb += lz;

  /* ----------------------------------------- */
  /* Calculate the angle of the bend and then  */
  /* the energy using the CHARMM force field.  */
  /* ----------------------------------------- */
    double rab2 = xab*xab + yab*yab + zab*zab;
	  double rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
	  double xp	= ycb*zab - zcb*yab;
    double yp   = zcb*xab - xcb*zab;
    double zp	= xcb*yab - ycb*xab;
    double rp	= sqrt(xp*xp + yp*yp + zp*zp);
	  double dot  = xab*xcb + yab*ycb + zab*zcb;

  	double drai = 1.0 / sqrt(rab2);
    double erax = xab * drai; 
    double eray = yab * drai; 
    double eraz = zab * drai; 
    double drci = 1.0 / sqrt(rcb2);
    double ercx = xcb * drci; 
    double ercy = ycb * drci; 
    double ercz = zcb * drci; 

    double cosphi = erax*ercx + eray*ercy + eraz*ercz;
    double phi = acos(cosphi);
    double isinphi = 1.0 / sqrt(1.0 - cosphi*cosphi);

    double dt = phi - bend[k][i].angeq;
	  double fr = k_theta * dt * 2.0;
    ebend += -fr * dt * 0.5;

    double irasin = drai * isinphi * fr;
    double ircsin = drci * isinphi * fr;

  /* ----------------------------------------- */
  /* Calculate the 2nd derivatives for configT */
  /* ----------------------------------------- */
	double deddt	= -fr;
	double d2eddt2	= -2.0* k_theta;
//     first derivatives of bond angle with respect to coordinates

    double terma   = -1.0  / (rab2*rp);
    double termc   =  1.0  / (rcb2*rp);
	  double ddtdxia = terma * (yab*zp-zab*yp);
    double ddtdyia = terma * (zab*xp-xab*zp);
    double ddtdzia = terma * (xab*yp-yab*xp);
    double ddtdxic = termc * (ycb*zp-zcb*yp);
    double ddtdyic = termc * (zcb*xp-xcb*zp);
    double ddtdzic = termc * (xcb*yp-ycb*xp);
    double ddtdxib = -ddtdxia - ddtdxic;
    double ddtdyib = -ddtdyia - ddtdyic;
    double ddtdzib = -ddtdzia - ddtdzic;
//     abbreviations used in defining chain rule terms

    double xrab = 2.0 * xab / rab2;
    double yrab = 2.0 * yab / rab2;
	  double zrab = 2.0 * zab / rab2;
    double xrcb = 2.0 * xcb / rcb2;
    double yrcb = 2.0 * ycb / rcb2;
    double zrcb = 2.0 * zcb / rcb2;
    double rp2  = 1.0 / (rp*rp);
    double xabp = (yab*zp-zab*yp) * rp2;
    double yabp = (zab*xp-xab*zp) * rp2;
    double zabp = (xab*yp-yab*xp) * rp2;
    double xcbp = (ycb*zp-zcb*yp) * rp2;
    double ycbp = (zcb*xp-xcb*zp) * rp2;
    double zcbp = (xcb*yp-ycb*xp) * rp2;

//     chain rule terms for second derivative components

    double dxiaxia =  terma*(xab*xcb-dot) + ddtdxia*(xcbp-xrab);
//    double dxiayia =  terma*(zp+yab*xcb)  + ddtdxia*(ycbp-yrab);
//    double dxiazia =  terma*(zab*xcb-yp)  + ddtdxia*(zcbp-zrab);
    double dyiayia =  terma*(yab*ycb-dot) + ddtdyia*(ycbp-yrab);
//    double dyiazia =  terma*(xp+zab*ycb)  + ddtdyia*(zcbp-zrab);
    double dziazia =  terma*(zab*zcb-dot) + ddtdzia*(zcbp-zrab);
    double dxicxic =  termc*(dot-xab*xcb) - ddtdxic*(xabp+xrcb);
//    double dxicyic =  termc*(zp-ycb*xab)  - ddtdxic*(yabp+yrcb);
//    double dxiczic = -termc*(yp+zcb*xab)  - ddtdxic*(zabp+zrcb);
    double dyicyic =  termc*(dot-yab*ycb) - ddtdyic*(yabp+yrcb);
//    double dyiczic =  termc*(xp-zcb*yab)  - ddtdyic*(zabp+zrcb);
    double dziczic =  termc*(dot-zab*zcb) - ddtdzic*(zabp+zrcb);
    double dxiaxic =  terma*(yab*yab+zab*zab) - ddtdxia*xabp;
//    double dxiayic = -terma*xab*yab - ddtdxia*yabp;
//    double dxiazic = -terma*xab*zab - ddtdxia*zabp;
//    double dyiaxic = -terma*xab*yab - ddtdyia*xabp;
    double dyiayic =  terma*(xab*xab+zab*zab) - ddtdyia*yabp;
//    double dyiazic = -terma*yab*zab - ddtdyia*zabp;
//    double dziaxic = -terma*xab*zab - ddtdzia*xabp;
//    double dziayic = -terma*yab*zab - ddtdzia*yabp;
    double dziazic =  terma*(xab*xab+yab*yab) - ddtdzia*zabp;

//     get some second derivative chain rule terms by difference

	double dxibxia = -dxiaxia - dxiaxic;
//	double dxibyia = -dxiayia - dyiaxic;
//	double dxibzia = -dxiazia - dziaxic;
//	double dyibxia = -dxiayia - dxiayic;
	double dyibyia = -dyiayia - dyiayic;
//	double dyibzia = -dyiazia - dziayic;
//	double dzibxia = -dxiazia - dxiazic;
//	double dzibyia = -dyiazia - dyiazic;
	double dzibzia = -dziazia - dziazic;
	double dxibxic = -dxicxic - dxiaxic;
//	double dxibyic = -dxicyic - dxiayic;
//	double dxibzic = -dxiczic - dxiazic;
//	double dyibxic = -dxicyic - dyiaxic;
	double dyibyic = -dyicyic - dyiayic;
//	double dyibzic = -dyiczic - dyiazic;
//	double dzibxic = -dxiczic - dziaxic;
//	double dzibyic = -dyiczic - dziayic;
	double dzibzic = -dziczic - dziazic;
	double dxibxib = -dxibxia - dxibxic;
//	double dxibyib = -dxibyia - dxibyic;
//	double dxibzib = -dxibzia - dxibzic;
	double dyibyib = -dyibyia - dyibyic;
//	double dyibzib = -dyibzia - dyibzic;
	double dzibzib = -dzibzia - dzibzic;

	double hessxa  = deddt*dxiaxia + d2eddt2*ddtdxia*ddtdxia;
	double hessxb  = deddt*dxibxib + d2eddt2*ddtdxib*ddtdxib;
	double hessxc  = deddt*dxicxic + d2eddt2*ddtdxic*ddtdxic;
	double hessya  = deddt*dyiayia + d2eddt2*ddtdyia*ddtdyia;
	double hessyb  = deddt*dyibyib + d2eddt2*ddtdyib*ddtdyib;
	double hessyc  = deddt*dyicyic + d2eddt2*ddtdyic*ddtdyic;
	double hessza  = deddt*dziazia + d2eddt2*ddtdzia*ddtdzia;
	double hesszb  = deddt*dzibzib + d2eddt2*ddtdzib*ddtdzib;
	double hesszc  = deddt*dziczic + d2eddt2*ddtdzic*ddtdzic;

  hesox	+= hessxa + hessxb + hessxc;
	hesoy	+= hessya + hessyb + hessyc;
	hesoz	+= hessza + hesszb + hesszc;
    
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

double eda_golik2_mc (int ibox)
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
  //double k_phi = 0.2 * golik[k].eps*20;
  double k_phi = 0.2 * golik[k].eps;
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
  double k_costau_phi = k_phi*cos(tau-tors[k][i].delphi[0]);
	etorsion	+=  k_phi+k_costau_phi;
  double d2edphi2 = -k_costau_phi;
#else
	double fr    =   k_phi*3.0*(sin(3.0*tau-tors[k][i].delphi[0]));
  double k_costau_phi = k_phi*cos(3.0*tau-tors[k][i].delphi[0]);
  etorsion    +=  k_phi+k_costau_phi;
  double d2edphi2 = -9.0*k_costau_phi;
#endif
	double dedphi = -fr;
  /* ----------------------------------------- */
  /* Calculate the 2nd derivatives for configT */
  /* ----------------------------------------- */
    
    //     abbreviations for first derivative chain rule terms

	double xca = atom[k][ic].x - atom[k][ia].x;
  double yca = atom[k][ic].y - atom[k][ia].y;
  double zca = atom[k][ic].z - atom[k][ia].z;
	double xdb = atom[k][id].x - atom[k][ib].x;
  double ydb = atom[k][id].y - atom[k][ib].y;
  double zdb = atom[k][id].z - atom[k][ib].z;

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
	

	hesox	+= hessxa + hessxb + hessxc + hessxd;
	hesoy	+= hessya + hessyb + hessyc + hessyd;
	hesoz	+= hessza + hesszb + hesszc + hesszd;

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

    double fb1   = fr * (rcb - rba*cosb) * irab *irbc * isinb2;
    double fb2   = fr * cosc * irbc * isinc2;
    ffox[k][ib] += (fb1 * abbcx + fb2 * dccbx);
    ffoy[k][ib] += (fb1 * abbcy + fb2 * dccby);
    ffoz[k][ib] += (fb1 * abbcz + fb2 * dccbz);

    double fc1   = fr * (rcb - rdc*cosc) * ircd *irbc *isinc2;
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

void enbnd_golik2_mc(int lb, int ub, int ibox, double *energy)
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
      double dx2 = dx*dx;
      double dy2 = dy*dy;
      double dz2 = dz*dz;

      double dr2 = dx2 + dy2 + dz2;
      double dr2i= 1.0 / dr2;

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
      double eps4 = 4.0*eps;
      double term1 = 48.0*eps*dr2i;
      
  /* ----------------------------------------------- */
  /* Calculate the LJ force and energy.              */
  /* ----------------------------------------------- */
	  double force=0.0; double enbnd =0.0; double enbnds = 0.0;
	  if((eps*sig)!=0.0) {
	  
		  double sr2  = (sig * sig) * dr2i;
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

			  force = term1 * (sr12 - 0.5*sr6);
        double term2 = term1*dr2i*(14.0*sr12-4.0*sr6);

			  enbnd = eps4 * (sr12 - sr6);
			  enbnds= eps4 * (sr12s -sr6s);
	      
        hesox    += 2.0*(term2*dx2-force);
	      hesoy    += 2.0*(term2*dy2-force);
	      hesoz    += 2.0*(term2*dz2-force);
			  
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
        double lambda_term1 = lambdaps*term1;
			  force = lambda_term1 * (sr12 - 0.5*nu_nu_ps*sr6);
        double term3 = lambda_term1*dr2i*(14.0*sr12 - 4.0*nu_nu_ps*sr6);
			  
        hesox    += 2.0*(term3*dx2-force);
	      hesoy    += 2.0*(term3*dy2-force);
	      hesoz    += 2.0*(term3*dz2-force);

			  enbnd = lambdaps* eps4 * (sr12 - nu_nu_ps*sr6);
			  enbnds= lambdaps* eps4 * (sr12s -nu_nu_ps*sr6s);
			  e_ps += (enbnd -enbnds);
			  enonbond += enbnd;
			  enonbonds+= enbnds;
		  }
		  else if (a> golik[k].Psite && b > golik[k].Psite) {
			  if(a<(golik[k].Psite+gosol.Nsol[0]) && b > (golik[k].Psite+gosol.Nsol[0])) nu_nu_ss= 0.5;
			  else nu_nu_ss = 0.5;
			  force = term1 * (sr12 - 0.5*nu_nu_ss*sr6);
        double term4 = term1*dr2i*(14.0*sr12 - 4.0*nu_nu_ss*sr6);

        hesox    += 2.0*(term4*dx2-force);
	      hesoy    += 2.0*(term4*dy2-force);
	      hesoz    += 2.0*(term4*dz2-force);


			  enbnd = eps4 * (sr12  - nu_nu_ss*sr6);
			  enbnds= eps4 * (sr12s - nu_nu_ss*sr6s);
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

void enbnd_dnagolik2_mc(int lb, int ub, int ibox, double *energy)
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
  double eps_AT = epsilon*2.0/3.0*3.0;
  double eps_GC = epsilon*3.0;
  double eps_mm = epsilon;
  double factor = pow(2.0,-(1.0/6.0));
  double sig_AT = 2.9002;
  double sig_GC = 2.8694;
  double d_mm   = 1.0;
  double sig_non= factor*golik[k].d_repulsive;
  double sig_mm = factor*d_mm;

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
  double rc   = sim.rc;
  double rc2  = box[k].rc2;
  double rc2i = 1.0/box[k].rc2;
  double rc4  = box[k].rc4;
  double rc_non = golik[k].d_repulsive;
  double rc2_non = rc_non * rc_non;
  double rc_mm = d_mm;
  double rc2_mm = rc_mm * rc_mm;


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
      double dx2 = dx * dx;
      double dy2 = dy * dy;
      double dz2 = dz * dz;
      double dr2  = dx2 + dy2 + dz2;
      double dr2i = 1.0/dr2;

      if(dr2 > rc2) continue;

  /* ----------------------------------------------- */
  /* Set the epsilon, sigma, and qq for the pair.    */
  /* ----------------------------------------------- */
	  double sig;
    double eps;
	  double eps4;
    double term1,term2;
    double sr2;
    double sr6;
    double sr12;
    double sr10;
    double sr2s;
    double sr6s;
    double sr12s;
    double sr10s;
	  double fx, fy, fz;
	  int flag=nlist[k].xflag[a][bb]; 
      
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

        term1 = 12.0*eps4*dr2i;
			  force = term1 * (sr12 - 0.5*sr6);
        term2 = term1*dr2i*(14.0*sr12-4.0*sr6);

			  enonbond  += eps4 * (sr12 - sr6);
			  enonbonds += eps4 * (sr12s -sr6s);

        hesox    += 2.0*(term2*dx2-force);
	      hesoy    += 2.0*(term2*dy2-force);
	      hesoz    += 2.0*(term2*dz2-force);

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
        
        term1 = 60.0*eps4*dr2i;
			  force = term1 * (sr12 - sr10);
        term2 = term1*dr2i*(14.0*sr12-12.0*sr10);

			  ebp  += eps4 * (5.0*sr12 - 6.0*sr10);
			  ebp  -= eps4 * (5.0*sr12s- 6.0*sr10s);
        
        hesox    += 2.0*(term2*dx2-force);
	      hesoy    += 2.0*(term2*dy2-force);
	      hesoz    += 2.0*(term2*dz2-force);
       
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
        
        term1 = 60.0*eps4*dr2i;
			  force = term1 * (sr12 - sr10);
        term2 = term1*dr2i*(14.0*sr12-12.0*sr10);

			  ebp  += eps4 * (5.0*sr12 - 6.0*sr10);
			  ebp  -= eps4 * (5.0*sr12s- 6.0*sr10s);

        hesox    += 2.0*(term2*dx2-force);
	      hesoy    += 2.0*(term2*dy2-force);
	      hesoz    += 2.0*(term2*dz2-force);

        break;
      case 4: //mis-match
        if(dr2<rc2_mm){
          sig = sig_mm;
          eps = eps_mm;
		      sr2  = (sig * sig) * dr2i;
		      sr6  = sr2 * sr2 * sr2;
		      sr12 = sr6 * sr6;

          term1 = 48.0*eps*dr2i;
			    force = term1 * (sr12 - 0.5*sr6);
          term2 = term1*dr2i*(14.0*sr12-4.0*sr6);

			    enonnative  += 4.0*eps*(sr12 - sr6) + eps;

          hesox    += 2.0*(term2*dx2-force);
	        hesoy    += 2.0*(term2*dy2-force);
	        hesoz    += 2.0*(term2*dz2-force);
        }
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

        term1 = 12.0*eps4*dr2i;
			  force = term1 * (sr12 - 0.5*sr6);
        term2 = term1*dr2i*(14.0*sr12-4.0*sr6);

			  enonbond  += eps4 * (sr12 - sr6);
			  enonbonds += eps4 * (sr12s -sr6s);

        hesox    += 2.0*(term2*dx2-force);
	      hesoy    += 2.0*(term2*dy2-force);
	      hesoz    += 2.0*(term2*dz2-force);

        break;
      default://non-native
        if(dr2<rc2_non){
          sig = sig_non;
          eps = epsilon;
		      sr2  = (sig * sig) * dr2i;
		      sr6  = sr2 * sr2 * sr2;
		      sr12 = sr6 * sr6;

          term1 = 48.0*eps*dr2i;
			    force = term1 * (sr12 - 0.5*sr6);
          term2 = term1*dr2i*(14.0*sr12-4.0*sr6);

			    enonnative  += 4.0*eps*(sr12 - sr6) + eps;

          hesox    += 2.0*(term2*dx2-force);
	        hesoy    += 2.0*(term2*dy2-force);
	        hesoz    += 2.0*(term2*dz2-force);
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
    }*/
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
          double dr3i = dri*dr2i;
          double qqerlambda = qq*exp(-lambdai/dri)*dri;
          term1 = (dri*lambdai+dr2i);
	  double term_force = qqerlambda*term1;
          force += term_force;
          ecoulomb += qqerlambda;

          term2 = term1*lambdai+2.0*dr2i*lambdai + 3.0*dr3i;
          double term3 = term2 * dri * qqerlambda;

          hesox    += 2.0*(term3*dx2-term_force);
	  hesoy    += 2.0*(term3*dy2-term_force);
	  hesoz    += 2.0*(term3*dz2-term_force);
	  
	#elif defined(RFC)
	
          double dri = sqrt(dr2i);
	  double dr3i = dri * dr2i;
	  double term_force = qq * (dr3i - 2.0 * rfc);
          force	    += term_force;
	  term1 = qq * (dri + rfc * dr2);
	  term2 = qq * (1.0 + rfcs) / rc;
	  double term3 = qq * 3.0 * dr3i * dr2i;
          enonbond  += term1;
          enonbonds += term2;
          ecoulomb  += term1 - term2;
	
          hesox    += 2.0*(term3*dx2-term_force);
	  hesoy    += 2.0*(term3*dy2-term_force);
	  hesoz    += 2.0*(term3*dz2-term_force);
	
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
#endif //GOLIK
#endif

