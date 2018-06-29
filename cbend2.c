/* ======================================================================== */
/* cbend2.cpp                                                                */
/*                                                                          */
/*		This subroutine calculates the analytical second derivative of the  */
/* bending energies and forces.         */
/* It returns the value of the potential energy due to bending to be stored */
/* in en[k].ebend.  It is called from either forces() or force_short().     */
/*                                                                          */
/* Passed Parameters:                                                       */
/*						ibox:		The box number                          */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"
#ifdef CONFIGT
#ifdef PR_NPT
void min_image_npt_full(int, double *dx, double *dy, double *dz);
#endif
/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

double cbend2 (int ibox)
{

  int k = ibox;
  /* ================================================================== */
  /*                                                                    */
  /* Zero out the energy accumulator                                    */
  /*                                                                    */
  /* ================================================================== */
  double ebend = 0.0;

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
  /* ================================================================== */
  /*                                                                    */
  /* Calculate the energy and force associated with each bend.  The     */
  /* loop is over the number of bends.                                  */
  /*                                                                    */
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
#ifndef PR_NPT
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
#endif
#ifdef PR_NPT
	min_image_npt_full(k, &xab, &ycb, &zcb);
	min_image_npt_full(k, &xcb, &ycb, &zcb);
#endif

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

    double dt  = phi - bend[k][i].angeq;
	
	double fr  = -bend[k][i].krbend * dt*2.0;
    ebend += -fr * dt * 0.5;
    double irasin = drai * isinphi * fr;
    double ircsin = drci * isinphi * fr;

	double deddt	= -fr;
	double d2eddt2	= 2* bend[k][i].krbend;
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

	config[k].hesx	+= hessxa + hessxb + hessxc;
	config[k].hesy	+= hessya + hessyb + hessyc;
	config[k].hesz	+= hessza + hesszb + hesszc;



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
 
    double sxx = (fxa * xab + fxc * xcb);
    double syx = (fya * xab + fyc * xcb);
    double syy = (fya * yab + fyc * ycb);
    double szx = (fza * xab + fzc * xcb);
    double szy = (fza * yab + fzc * ycb);
    double szz = (fza * zab + fzc * zcb);
    wxx += sxx;
    wyx += syx;
    wyy += syy;
    wzx += szx;
    wzy += szy;
    wzz += szz;
#endif
 /* ----------------------------------------- */
/* If ZHOU is defined, accumulate the        */
/* components of the stress tensor.          */
/* ----------------------------------------- */
#ifdef ZHOU
  	double third = 1.0 / 3;
        stresses[k][ia][0][0] += sxx * third;
        stresses[k][ia][1][0] += syx * third;
        stresses[k][ia][1][1] += syy * third;
        stresses[k][ia][2][0] += szx * third;
        stresses[k][ia][2][1] += szy * third;
        stresses[k][ia][2][2] += szz * third;

        stresses[k][ib][0][0] += sxx * third;
        stresses[k][ib][1][0] += syx * third;
        stresses[k][ib][1][1] += syy * third;
        stresses[k][ib][2][0] += szx * third;
        stresses[k][ib][2][1] += szy * third;
        stresses[k][ib][2][2] += szz * third;

        stresses[k][ic][0][0] += sxx * third;
        stresses[k][ic][1][0] += syx * third;
        stresses[k][ic][1][1] += syy * third;
        stresses[k][ic][2][0] += szx * third;
        stresses[k][ic][2][1] += szy * third;
        stresses[k][ic][2][2] += szz * third;

#endif  //ZHOU

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
#endif
