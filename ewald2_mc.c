#if defined(CONFIGT) && defined(EWALD) && defined(MC)
/* ======================================================================== */
/* ewald2_mc.cpp                                                            */
/*                                                                          */
/* Written by Thomas A. Knotts IV, 1 Sep 2005                               */
/* ======================================================================== */

#include "defines.h"


  /* ----------------------------------------------- */
  /* Functions controlling complex numbers.  Found   */
  /* in cplx.c                                       */
  /* ----------------------------------------------- */    
struct Complex C_conj (struct Complex a);
struct Complex C_add  (struct Complex a, struct Complex b);
struct Complex C_sub  (struct Complex a, struct Complex b);
struct Complex C_mul  (struct Complex a, struct Complex b);
struct Complex C_div  (struct Complex a, struct Complex b);
struct Complex C_para (struct Complex a, struct Complex b);
struct Complex C_conj (struct Complex a);
struct Complex C_inv  (struct Complex a);
struct Complex C_sclmul (double, struct Complex a);
double C_mag (struct Complex a);
double C_ang (struct Complex a);





/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void ewald_kspace2_mc (int ibox, double *enkewald)
{
  /* ================================================================== */
  /*                                                                    */
  /* Variables needed for this subroutine.                              */
  /*                                                                    */
  /* ================================================================== */
  int k = ibox;
  double kappa = sim.kappa[k];
  double two_pi = 2.0 * PI;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;
  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double bv = lx * ly * lz;
  double iblx   = 1.0 / lx;
  double ibly   = 1.0 / ly;
  double iblz   = 1.0 / lz;

  double two_pi_sq   = two_pi * two_pi;
  double two_pi_iblx = two_pi*iblx;
  double two_pi_ibly = two_pi*ibly;
  double two_pi_iblz = two_pi*iblz;
  double two_pi_iblx_sq = two_pi_iblx*two_pi_iblx;
  double two_pi_ibly_sq = two_pi_ibly*two_pi_ibly;
  double two_pi_iblz_sq = two_pi_iblz*two_pi_iblz;
  int ns = box[k].boxns;
  double fconv = (ELEC * ELEQ * ELEQ * NA * 0.001);   // convert to units of kJ/mol
  int kmax_x = sim.kmax[k][0];
  int kmax_y = sim.kmax[k][1];
  int kmax_z = sim.kmax[k][2];
  int max_k_square = sim.max_k_square[k];
  double re1;
  double hes;

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
  /* Note:  Below we use Euler's Identity to do the mathmatics of       */
  /* complex expontials with sines and cosines.  The identities are:    */
  /* exp(ik.r)=cos(k.r)+i*sin(k.r) and exp(-k.r)=cos(k.r)-i*sin(k.r)    */
  /* Note 2: the exp(ik.r) term is needed for the calculation of the    */
  /* charge density, rho(k).                                            */
  /*                                                                    */
  /* ================================================================== */
 


  /* ================================================================== */
  /*                                                                    */
  /* Determine the e_ik's for (0,0,0), (1,1,1), and (-1,-1,-1)          */
  /* explicitly for each site.                                          */
  /*                                                                    */
  /* ================================================================== */
  for(int i=0; i<ns; i++) {
  /* ----------------------------------------- */
  /* Get the postion vector, r, for the        */
  /* exp(ik.r) part.                           */
  /* ----------------------------------------- */
    double dx = atom[k][i].x;
    double dy = atom[k][i].y;
    double dz = atom[k][i].z;
    

  /* ----------------------------------------- */
  /* Set the values of each e_ik for the three */
  /* k vectors (0,0,0) (1,1,1) and (-1,-1,-1). */
  /* ----------------------------------------- */
    e_ikx[i][0].r = 1.0;     e_ikx[i][0].i = 0.0;
    e_iky[i][0].r = 1.0;     e_iky[i][0].i = 0.0;
    e_ikz[i][0].r = 1.0;     e_ikz[i][0].i = 0.0;
    
    e_ikx[i][1].r = cos(two_pi_iblx * (double)dx);
    e_ikx[i][1].i = sin(two_pi_iblx * (double)dx);
    
    e_iky[i][1].r = cos(two_pi_ibly * (double)dy);
    e_iky[i][1].i = sin(two_pi_ibly * (double)dy);
   
    e_ikz[i][1].r = cos(two_pi_iblz * (double)dz);
    e_ikz[i][1].i = sin(two_pi_iblz * (double)dz);

    e_iky[i][-1] = C_conj(e_iky[i][1]);
    e_ikz[i][-1] = C_conj(e_ikz[i][1]);
  }

  /* ----------------------------------------- */
  /* Calculate the remaining e_ik's by         */
  /* recurrence.                               */
  /* ----------------------------------------- */
  for(int i=0; i<ns; i++){
	for(int kx=2; kx<=kmax_x; kx++){
		e_ikx[i][kx] = C_mul(e_ikx[i][kx-1],e_ikx[i][1]);
	}
	for(int ky=2; ky<=kmax_y; ky++){
		e_iky[i][ky] = C_mul(e_iky[i][ky-1],e_iky[i][1]);
		e_iky[i][-ky] = C_conj(e_iky[i][ky]);
	}
	for(int kz=2; kz<=kmax_z; kz++){
		e_ikz[i][kz] = C_mul(e_ikz[i][kz-1],e_ikz[i][1]);
		e_ikz[i][-kz] = C_conj(e_ikz[i][kz]);
	}
  }

  /* ================================================================== */
  /*                                                                    */
  /* Now, we can calculate the k-space contribution to the Coulomibic   */
  /* energy and forces according to the definition.                     */
  /*                                                                    */
  /* ================================================================== */
  /* ----------------------------------------- */
  /* Variable declarations.                     */
  /* ----------------------------------------- */
  double vk = 0.0;
  double vd = 0.0;
  int totk=0;
  struct Complex *e_ikr;

  e_ikr = (struct Complex*) calloc(ns,sizeof(struct Complex));

  double *fx;
  double *fy;
  double *fz;

  fx = (double*) calloc(ns,sizeof(double));
  fy = (double*) calloc(ns,sizeof(double));
  fz = (double*) calloc(ns,sizeof(double));
  struct Complex im1;
  double im2;
  double frc;
	

  for(int kk=0; kk<total_k_vectors[k]; kk++){
	  sum[kk].r = 0.0;
	  sum[kk].i = 0.0;
  }
  /* ----------------------------------------- */
  /* Loop first where kx=0.                    */
  /* ----------------------------------------- */  
  for(int ky=-kmax_y; ky<=kmax_y; ky++){
	for(int kz=-kmax_z; kz<=kmax_z; kz++){
		int ksq = ky*ky+kz*kz;
		if(ksq<=max_k_square && ksq!=0){
			sum[totk].r = 0.0;
			sum[totk].i = 0.0;
			for(int i=0; i<ns; i++){
				e_ikr[i]=C_mul(e_ikx[i][0],C_mul(e_iky[i][ky],e_ikz[i][kz]));
				sum[totk].r = sum[totk].r + pott[k][i].qq*e_ikr[i].r;
				sum[totk].i = sum[totk].i + pott[k][i].qq*e_ikr[i].i;
			}
#ifndef PRESSURE
			vk += kvec[k][totk]*(sum[totk].r*sum[totk].r+sum[totk].i*sum[totk].i);
#endif
			for(int i=0; i<ns; i++){
				im1 = C_mul(C_conj(e_ikr[i]),sum[totk]);
				im2 = 2.0*im1.i;
				re1 = im1.r;
				hes = 2.0*fconv*pott[k][i].qq*two_pi*two_pi*kvec[k][totk]*(pott[k][i].qq-re1);
				frc = fconv*pott[k][i].qq*two_pi*kvec[k][totk]*im2; //this is in kJ/ang/mol
				fy[i]  = fy[i] - frc*ibly*(double)ky;
				fz[i]  = fz[i] - frc*iblz*(double)kz;
				hesoy += hes*ibly*ibly*(double)ky*(double)ky;
				hesoz += hes*iblz*iblz*(double)kz*(double)kz;
			}
			
#ifdef PRESSURE
			double term1 = kvec[k][totk]*(sum[totk].r*sum[totk].r+sum[totk].i*sum[totk].i);	
			double term2 = 0.25/kappa/kappa;
			double term3 = 2.0/ksq;
			double term33 = 2.0 * term2 * two_pi_sq;
		        double term3yy = term33*ibly*ibly;
		        double term3zy = term33*iblz*ibly;
		        double term3zz = term33*iblz*iblz;

		        wyy += (1.0-(term3+term3yy)*(double)ky*(double)ky)*term1;
		        wzy += (0.0-(term3+term3zy)*(double)kz*(double)ky)*term1;
		        wzz += (1.0-(term3+term3zz)*(double)kz*(double)kz)*term1;

			vk += term1;
#endif	
			totk++;
		}//end if
	}
  }

  /* ----------------------------------------- */
  /* Now loop to get the rest of the vectors.  */
  /* ----------------------------------------- */
  for(int kx=1; kx<=kmax_x; kx++){
	  for(int ky=-kmax_y; ky<=kmax_y; ky++){
		  for(int kz=-kmax_z; kz<=kmax_z; kz++){
			int ksq = kx*kx+ky*ky+kz*kz;
			if(ksq<=max_k_square){
				sum[totk].r = 0.0;
				sum[totk].i = 0.0;
				for(int i=0; i<ns; i++){
					e_ikr[i]=C_mul(e_ikx[i][kx],C_mul(e_iky[i][ky],e_ikz[i][kz]));
					sum[totk].r = sum[totk].r + pott[k][i].qq*e_ikr[i].r;
					sum[totk].i = sum[totk].i + pott[k][i].qq*e_ikr[i].i;
				}
#ifndef PRESSURE
			vk += 2.0*kvec[k][totk]*(sum[totk].r*sum[totk].r+sum[totk].i*sum[totk].i);
#endif
			for(int i=0; i<ns; i++){
				im1 = C_mul(C_conj(e_ikr[i]),sum[totk]);
				im2 = 2.0*im1.i;
				re1 = im1.r;
				hes = 4.0*fconv*pott[k][i].qq*two_pi*two_pi*kvec[k][totk]*(pott[k][i].qq-re1);
				frc = 2.0*fconv*pott[k][i].qq*two_pi*kvec[k][totk]*im2; //this is in kJ/ang/mol
				fx[i] = fx[i] - frc*iblx*(double)kx;
				fy[i] = fy[i] - frc*ibly*(double)ky;
				fz[i] = fz[i] - frc*iblyz*(double)kz;
				hesox += hes*iblx*iblx*(double)kx*(double)kx;
				hesoy += hes*ibly*ibly*(double)ky*(double)ky;
				hesoz += hes*iblz*iblz*(double)kz*(double)kz;

			}

#ifdef PRESSURE
                        double term1 = 2.0*kvec[k][totk]*(sum[totk].r*sum[totk].r+sum[totk].i*sum[totk].i); 
                        double term2 = 0.25/kappa/kappa;
                        double term3 = 2.0/ksq;                                                            double term33 = 2.0 * term2 * two_pi_sq;
		        double term3xx = term33*iblx*iblx;
         		double term3yx = term33*ibly*iblx;
          		double term3yy = term33*ibly*ibly;
          		double term3zx = term33*iblz*iblx;
          		double term3zy = term33*iblz*ibly;
          		double term3zz = term33*iblz*iblz;

          		wxx += (1.0-(term3+term3xx)*(double)kx*(double)kx)*term1;
          		wyx += (0.0-(term3+term3yx)*(double)ky*(double)kx)*term1;
          		wyy += (1.0-(term3+term3yy)*(double)ky*(double)ky)*term1;
         		wzx += (0.0-(term3+term3zx)*(double)kz*(double)kx)*term1;
          		wzy += (0.0-(term3+term3zy)*(double)kz*(double)ky)*term1;
          		wzz += (1.0-(term3+term3zz)*(double)kz*(double)kz)*term1;
	
			vk += term1;
#endif				
				totk++;
			}//endif
		  }//for kz
	  }//for ky
  }//for kx
				

  /* ================================================================== */
  /*                                                                    */
  /* Convert the energy from charge units to kJ/mol and assign the      */
  /* values to enkewald[].												*/
  /*                                                                    */
  /* ================================================================== */
  enkewald[0] = vk * fconv;
  enkewald[1] = en[k].sewald;
  enkewald[2] = enkewald[0]-enkewald[1];
 
 
  for(int i=0; i<ns; i++){  
    ffox[k][i] += fx[i];
    ffoy[k][i] += fy[i];
    ffoz[k][i] += fz[i];
  }

//  }//for i

#ifdef PRESSURE

  pviro[k].ewald_kspace[0] = wxx*fconv;
  pviro[k].ewald_kspace[1] = wyx*fconv;
  pviro[k].ewald_kspace[2] = wyy*fconv;
  pviro[k].ewald_kspace[3] = wzx*fconv;
  pviro[k].ewald_kspace[4] = wzy*fconv;
  pviro[k].ewald_kspace[5] = wzz*fconv;

#endif
  free(e_ikr);
  free(fx);
  free(fy);
  free(fz);

}
 
/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
double erfc(double);
void ewald_rspace2_mc (int lb, int ub, int ibox, double *energy)
{

  int k = ibox;

  /* ================================================================== */
  /*                                                                    */
  /* Zero out the energy accumulator                                    */
  /*                                                                    */
  /* ================================================================== */
  double kappa = sim.kappa[k];
  double kpi  = 1.0 * kappa / sqrt(PI);
  double eewald_real  = 0.0;
  double eewald_reals = 0.0;
  double eewald_intra = 0.0;
  double hes1;
  double hes2;
  double fconv = (ELEC * ELEQ * ELEQ * NA * 0.001);
  double rc = sim.rc_ew;
  double rc2 = rc * rc;
 
  
  /* ================================================================== */
  /*                                                                    */
  /* Zero out the virial accumulators                                   */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
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

#endif//PRESSURE

  /* ================================================================== */
  /*                                                                    */
  /* Get and assign information for pbc and cutoffs                     */
  /*                                                                    */
  /* ================================================================== */
  double lx   = box[k].lx;
  double ly   = box[k].ly;
  double lz   = box[k].lz;
  double hx   = box[k].hx;
  double hy   = box[k].hy;
  double hz   = box[k].hz;
  double bh   = box[k].least_boxh;
  //double bh2  = bh * bh;

  /* ================================================================== */
  /*                                                                    */
  /* Loop around all the particles to get the real space and intra-     */
  /* molecular correction to k-space contributions to the ewald         */
  /* coulombic energy.                                                  */
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
    for(int bb=0; bb<nlist_ew[k].count[a]; bb++) {
      int b = nlist_ew[k].list[a][bb];

	  double fx_ewaldr;
	  double fy_ewaldr;
	  double fz_ewaldr;

	  double fx_ewaldi;
	  double fy_ewaldi;
	  double fz_ewaldi;

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
      if(dx >  hx) dx =- lx;
      if(dx < -hx) dx =+ lx;
      if(dy >  hy) dy =- ly;
      if(dy < -hy) dy =+ ly;
      if(dz >  hz) dz =- lz;
      if(dz < -hz) dz =+ lz;

  /* ----------------------------------------------- */
  /* Calculate the distance between particles a & b. */
  /* If the distance is less than half the box,      */
  /* length calculated the energy and the force. If  */
  /* it is  greater, go to the next pair in the loop.*/
  /* ----------------------------------------------- */
      double dr2 = dx*dx + dy*dy + dz*dz;
      if(dr2 > rc2) continue;
#ifndef STYPE
	  double qq  = ljset[k].pot[a][b].qq;
#endif

#ifdef STYPE
	  double qqab  = pott[k][a].qq * pott[k][b].qq * fconv;
	  int xflag = nlist_ew[k].xflag[a][bb];
	  double qq=qqab;
	  if(xflag){
		 if(xflag == 4) qq = E14FAC * qqab;
		 else qq = 0.0;
	  }
#endif //ifdef STYPE

	  double force_ewaldr = 0.0;	  
	  double force_ewaldi = 0.0;

  /* ----------------------------------------------- */
  /* Calculated Ewald coulombic interactions.        */
  /* ----------------------------------------------- */
       /* *********************************************** */
       /* Quantities needed.                              */
       /* *********************************************** */
      double dr = sqrt(dr2);
      double dri = 1.0 / dr;
	  double dri2 = dri * dri;
	  double dri3 = dri2 * dri;
	  double dri4 = dri2 * dri2;
	  double dri5 = dri3 * dri2;
	  double kr = kappa * dr;
	  double erfckr = erfc(kr);
	  double expkr2 = exp(-kr*kr);
   
       /* *********************************************** */
       /* Real and real-shifted part                      */
       /* *********************************************** */
	  if (qq!=0.0) {
		eewald_real += qq * erfckr * dri;
		eewald_reals += qq * erfc(kappa * rc) / rc;
		force_ewaldr = (qq*dri/dr2) * (erfckr + 2.0*kpi*dr*expkr2);
		hes1= erfckr*dri3 + 2.0*kpi*dri2*expkr2;
		hes2=-6.0*kpi*dri4*expkr2-4.0*kpi*kappa*kappa*dri2*expkr2-3.0*dri5*erfckr;
		hesox += -2.0*qq * (hes1 + dx*dx*hes2);
		hesoy += -2.0*qq * (hes1 + dy*dy*hes2);
		hesoz += -2.0*qq * (hes1 + dz*dz*hes2);
	  }
       /* *********************************************** */
       /* Intra-molecular part                            */
	   /* Note: If the pair is the 1-4 atoms in a         */
	   /* torsion, the coulomibic interaction is E14FAC of   */
	   /* the full interaction, so the intra-molecular    */
	   /* part of the Ewald summation should be           */
	   /* evaluated with (1.0-E14FAC) of the full interation.      */
       /* *********************************************** */
#ifndef STYPE
	  int test = ljset[k].pot[a][b].bond_flag + ljset[k].pot[a][b].bend_flag + ljset[k].pot[a][b].tors_flag + ljset[k].pot[a][b].exNB_flag;
#ifndef LJ
	  if (test > 0){
		double q1 = pott[k][a].qq;
		double q2 = pott[k][b].qq;
		double qq12 = q1*q2*fconv;
		if (ljset[k].pot[a][b].tors_flag == 1) qq12 *= (1.0-E14FAC);
		eewald_intra += qq12*(1.0-erfckr)*dri;
		force_ewaldi = (qq12*dri/dr2) * ((1-erfckr)-2.0*kpi*dr*exp(-kr*kr));
		hes1= (1-erfckr)*dri3 - 2.0*kpi*dri2*expkr2;
		hes2=-6.0*kpi*dri4*expkr2-4.0*kpi*kappa*kappa*dri2*expkr2+3.0*dri5*(1-erfckr);
		hesox += 2.0*qq12 * (hes1 - dx*dx*hes2);
		hesoy += 2.0*qq12 * (hes1 - dy*dy*hes2);
		hesoz += 2.0*qq12 * (hes1 - dz*dz*hes2);

	  }  
#endif
#endif //ifndef STYPE



#ifdef STYPE
#ifndef LJ
	  if (xflag){
		double qq12 = qqab;
		if (xflag == 4) qq12 *= (1.0-E14FAC);
		eewald_intra += qq12*(1.0-erfckr)*dri;
		force_ewaldi = (qq12*dri/dr2) * ((1-erfckr)-2.0*kpi*dr*exp(-kr*kr));
		hes1= (1-erfckr)*dri3 - 2.0*kpi*dri2*expkr2;
		hes2=-6.0*kpi*dri4*expkr2-4.0*kpi*kappa*kappa*dri2*expkr2+3.0*dri5*(1-erfckr);
		hesox += 2.0*qq12 * (hes1 - dx*dx*hes2);
		hesoy += 2.0*qq12 * (hes1 - dy*dy*hes2);
		hesoz += 2.0*qq12 * (hes1 - dz*dz*hes2);	  
	  }  
#endif
#endif

  /* ----------------------------------------- */
  /* Accumulate the forces.                    */
  /* ----------------------------------------- */
       /* *********************************************** */
       /* x-forces                                        */
       /* *********************************************** */
       /* ///////////////////////////////// */
       /* Ewald Real part contribution      */
       /* ///////////////////////////////// */
	  fx_ewaldr    = force_ewaldr*dx;
      fxa        += fx_ewaldr;
      ffox[k][b] -= fx_ewaldr;
       /* ///////////////////////////////// */
       /* Ewald Intramolecular contribution */
       /* ///////////////////////////////// */
	  fx_ewaldi    = -force_ewaldi*dx;
      fxa        += fx_ewaldi;
      ffox[k][b] -= fx_ewaldi;

       /* *********************************************** */
       /* y-forces                                        */
       /* *********************************************** */
       /* ///////////////////////////////// */
       /* Ewald Real part contribution      */
       /* ///////////////////////////////// */
	  fy_ewaldr    = force_ewaldr*dy;
      fya        += fy_ewaldr;
      ffoy[k][b] -= fy_ewaldr;
       /* ///////////////////////////////// */
       /* Ewald Intramolecular contribution */
       /* ///////////////////////////////// */
	  fy_ewaldi    = -force_ewaldi*dy;
      fya        += fy_ewaldi;
      ffoy[k][b] -= fy_ewaldi;

       /* *********************************************** */
       /* z-forces                                        */
       /* *********************************************** */
       /* ///////////////////////////////// */
       /* Ewald Real part contribution      */
       /* ///////////////////////////////// */
	  fz_ewaldr    = force_ewaldr*dz;
      fza        += fz_ewaldr;
      ffoz[k][b] -= fz_ewaldr;
       /* ///////////////////////////////// */
       /* Ewald Intramolecular contribution */
       /* ///////////////////////////////// */
	  fz_ewaldi    = -force_ewaldi*dz;
      fza        += fz_ewaldi;
      ffoz[k][b] -= fz_ewaldi;
  /* ----------------------------------------- */
  /* If PRESSURE is defined, accumulate the    */
  /* components of the virial tensor.          */
  /* ----------------------------------------- */
#ifdef PRESSURE
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

#endif//PRESSURE
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
  pviro[k].ewald_real[0] = wxx_ewaldr;
  pviro[k].ewald_real[1] = wyx_ewaldr;
  pviro[k].ewald_real[2] = wyy_ewaldr;
  pviro[k].ewald_real[3] = wzx_ewaldr;
  pviro[k].ewald_real[4] = wzy_ewaldr;
  pviro[k].ewald_real[5] = wzz_ewaldr;

  pviro[k].ewald_intra[0] = wxx_ewaldi;
  pviro[k].ewald_intra[1] = wyx_ewaldi;
  pviro[k].ewald_intra[2] = wyy_ewaldi;
  pviro[k].ewald_intra[3] = wzx_ewaldi;
  pviro[k].ewald_intra[4] = wzy_ewaldi;
  pviro[k].ewald_intra[5] = wzz_ewaldi;
#endif //PRESSURE

  /* ================================================================== */
  /*                                                                    */
  /* After the loop, assign the accumulated energies to "energy" array. */
  /*                                                                    */
  /* ================================================================== */  
  /* ------------------------------------------ */
  /* Total and assign all the different parts   */
  /* to energy[].                               */
  /* ------------------------------------------ */ 
  energy[0] = eewald_real;						// non-shifted real space energy
  energy[1] = eewald_real - eewald_reals;		// shifted real_space ewald energy
  energy[2] = eewald_intra;						// intra-molecular self energy 
}//end ewald_rspace




#endif //def EWALD

