#ifdef MC

#include "defines.h"

#ifdef NLIST

#ifdef EWALD
#ifdef SPME
void pme_calc_mc(int,double*);
#else
void ewald_kspace_mc (int,double*);
#endif
void ewald_rspace_mc (int, int, int,double*);
double erfc(double);
#endif




/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

void enbnd_mc(int lb, int ub, int ibox, double *energy)
{

  int k = ibox;
  double fconv = (ELEC * ELEQ * ELEQ * NA * 0.001);   // convert atomic charges to units of kJ/mol
  /* ================================================================== */
  /*                                                                    */
  /* Zero out the energy accumulator                                    */
  /*                                                                    */
  /* ================================================================== */
  double enonbond  = 0.0;
  double enonbonds = 0.0;
  double ecoulomb  = 0.0;

  #ifdef SPME
  /* ================================================================== */
  /* Initialize ewald parameters.  These are used if you do not call    */
  /* ewald_rspace.                                                      */
  /* ================================================================== */
  double kappa = sim.kappa[k];
  double kpi = 1.0 * kappa / sqrt(PI);
  double eewald_real = 0.0;
  double eewald_reals = 0.0;
  double eewald_intra = 0.0;
  #endif

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

#ifdef RFC
  double rfcs = (sim.epsRF[k] - 1.0) / (2.0*sim.epsRF[k] + 1.0);
  double rfc  = rfcs / (rc2*rc);
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
	  if ((a>=lb) || (b>=lb && b<ub)){
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

  /* ----------------------------------------------- */
  /* Set the site types of particles a and b.        */
  /* ----------------------------------------------- */
#ifdef STYPE
      int ta = atom[k][a].atomid;
	  int tb = atom[k][b].atomid;

	  int xflag = nlist[k].xflag[a][bb];
  /* ----------------------------------------------- */
  /* Check to see if a and b are a 1-4 interactions. */
  /* if they are, set sig and eps to special values. */
  /* ----------------------------------------------- */

#ifdef S14	  
	  
	  if(xflag == 4){
		sig = ljset14[ta][tb].sig;
		eps = ljset14[ta][tb].eps;  
	  }else{
		sig = ljset[ta][tb].sig;
		eps = ljset[ta][tb].eps;
	  }
#endif
#ifndef S14
		sig = ljset[ta][tb].sig;
		eps = ljset[ta][tb].eps;

#endif //#ifndef S14

#ifdef COULOMB
		qq = pott[k][a].qq * pott[k][b].qq * fconv;
#endif
#ifdef SPME
		qq = pott[k][a].qq * pott[k][b].qq * fconv;
		double qqab = qq;
#endif


  /* ----------------------------------------------- */
  /* Check to see if a and b are an excluded         */
  /* interaction.  If they are, change sig, eps,     */
  /* and qq.                                         */
  /* ----------------------------------------------- */
#ifndef LJ

			if(xflag){//xflag is zero if there is no exclusion
				if(xflag != 4){
					sig = 0.0;
					eps = 0.0;
#ifdef COULOMB		
					qq = 0.0;
#endif
#ifdef SPME
			 		qq = 0.0;
#endif							
				}
#ifdef COULOMB
				else qq *= E14FAC;
#endif
#ifdef SPME
				else qq *= E14FAC;
#endif	
			}
#endif//#ifndef LJ

#endif//#ifdef STYPE

#ifndef STYPE	  
	  sig = ljset[k].pot[a][b].sig;
      eps = ljset[k].pot[a][b].eps;
      qq  = ljset[k].pot[a][b].qq;
#endif
	  



  /* ----------------------------------------------- */
  /* Calculate the LJ force and energy.              */
  /* ----------------------------------------------- */
	  double force=0.0;
	  if((eps*sig)!=0.0) {
	  
      double sr2  = (sig * sig) / dr2;
      double sr6  = sr2 * sr2 * sr2;
      double sr12 = sr6 * sr6;

      double sr2s  = (sig * sig) / rc2;
      double sr6s  = sr2s * sr2s * sr2s;
      double sr12s = sr6s * sr6s;
#ifndef LJ
      force = (12.0*eps/dr2) * (sr12 - sr6);
	  enonbond  +=  eps * (sr12 - 2.0*sr6);
	  enonbonds +=  eps * (sr12s - 2.0*sr6s);
#endif
#ifdef LJ
	  force = (48.0*eps/dr2) * (sr12 - 0.5*sr6);
	  enonbond  +=  4.0*eps * (sr12  - sr6);
      enonbonds +=  4.0*eps * (sr12s - sr6s);
#endif
	  }
	  
  /* ----------------------------------------------- */
  /* Calculated coulombic interactions.              */
  /* ----------------------------------------------- */
#ifdef COULOMB
	if (qq!=0.0) {

		#ifdef RFC
			  double dri = 1.0 / sqrt(dr2);
			  force	    += (qq/dr2) * (dri - 2.0 * rfc * dr2);
			  enonbond  += qq * (dri + rfc * dr2);
			  enonbonds += qq * (1.0 + rfcs) / rc;
			  ecoulomb  += qq * (dri + rfc * dr2) - qq * (1.0 + rfcs) / rc;
		#endif

		#ifdef SHIFTF
				#ifdef RDIE
					double dri4	= 1.0/dr2/dr2;
					double terma	= qq/sim.epsRF[k];
					double termb	= terma/dr2;
					double termc	= terma*(2.0/rc2-dr2/rc4);
					double term3	= 2.0*terma*(dri4-1.0/rc4);
					force		   += term3;
					enonbond	   += termb;
					enonbonds	   += termc;
					ecoulomb	   += termb-termc;
				#endif	//RDIE
				#ifndef RDIE
					double dri		= 1.0 / sqrt(dr2);
					double dri4		= 1.0/dr2/dr2;
					double terma	= qq/sim.epsRF[k]*dri;
					double termc	= terma*(2.0*dr2/rc2-1.0/dri4/rc4);
					double term3	= terma*(1.0/dr2+2.0/rc2-3.0*dr2/rc4);
					
					force	  += term3;
					enonbond  += terma;
					enonbonds += termc;
					ecoulomb  += terma-termc;
				#endif //NO RDIE
		#endif //FSHIFT

		#ifdef SHIFTV
				#ifdef RDIE
					double terma = qq/sim.epsRF[k]/dr2;			//unlike before dielectric to be picked from input file;
					double termb = qq/sim.epsRF[k]/rc2;
					force	  += 2.0*terma/dr2;
					enonbond  += terma;				
					enonbonds += termb;
					ecoulomb  += terma-termb;
				#endif	//RDIE
				#ifndef RDIE
					double dri = 1.0 / sqrt(dr2);
					double terma = qq/sim.epsRF[k]*dri;			//unlike before dielectric to be picked from input file;
					double termb = qq/sim.epsRF[k]/rc;
					force	  += terma/dr2;
					enonbond  += terma;
					enonbonds += termb;
					ecoulomb  += terma-termb;
				#endif //NO RDIE
		#endif //VSHIFT

	}//if qq!=0
#endif	// COULOMB
#ifdef SPME
  /* ----------------------------------------------- */
  /* If SPME is defined, the real part and           */
  /* intra-molecular part of the coulombic           */
  /* interaction is calculated here.                 */
  /* ----------------------------------------------- */
       /* *********************************************** */
       /* Quantities needed.                              */
       /* *********************************************** */
      double dr = sqrt(dr2);
      double dri = 1.0 / dr;
      double kr = kappa * dr;
      double erfckr = erfc(kr);
	  double kpidrexpkr2 =2.0*kpi*dr*exp(-kr * kr);

       /* *********************************************** */
       /* Real and real-shifted part                      */
       /* *********************************************** */
      if (qq!=0.0) {
        eewald_real += qq * erfckr * dri;
        eewald_reals += qq * erfc(kappa * rc) / rc;
        force += (qq*dri/dr2) * (erfckr + kpidrexpkr2);
      }
       /* *********************************************** */
       /* Intra-molecular part                            */
       /* Note: If the pair is the 1-4 atoms in a         */
       /* torsion, the coulomibic interaction is 0.4 of   */
       /* the full interaction, so the intra-molecular    */
       /* part of the Ewald summation should be           */
       /* evaluated with 0.6 of the full interation.      */
       /* *********************************************** */
#ifndef STYPE
	int test = ljset[k].pot[a][b].bond_flag + ljset[k].pot[a][b].bend_flag + ljset[k].pot[a][b].tors_flag + ljset[k].pot[a][b].exNB_flag;
#ifndef LJ
      if (test > 0){
		double qqab = pott[k][a].qq * pott[k][b].qq * fconv;
        if (ljset[k].pot[a][b].tors_flag == 1) qqab *= (1.0-E14FAC);
        eewald_intra += qqab*(1.0-erfckr)*dri;
        force -= (qqab*dri/dr2) * ((1.0-erfckr)-kpidrexpkr2);
      }
#endif
#endif
	 
#ifdef STYPE
#ifndef LJ
      if (xflag){
        if (xflag == 4) qqab *= (1.0-E14FAC);
        eewald_intra += qqab*(1.0-erfckr)*dri;
        force -= (qqab*dri/dr2) * ((1.0-erfckr)-kpidrexpkr2);
      }       
#endif
#endif
#endif //SPME
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
		}// if valid interaction
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
  /*                                                                    */
  /* ================================================================== */
#ifndef EWALD
  energy[0] = enonbond;                 // non-shifted energy (LJ and coulomibic)
  energy[1] = enonbond - enonbonds;     // shifted energy (LJ and coulomibic)
#endif

#ifdef COULOMB
  energy[2] = ecoulomb; //(shifted coulombic energy only)
#endif

#ifdef EWALD
  /* ------------------------------------------ */
  /* Call the subroutines to calucate the       */
  /* k-space and self energy contributions      */
  /* to the electrostatic energy.               */
  /*                                            */
  /* enkewald[0] = k-space ewald energy         */
  /* enkewald[1] = self correction ewald energy */
  /* enkewald[2] = k-space - self energy        */
  /* ------------------------------------------ */
  double enkewald[3];
#ifdef SPME
    pme_calc_mc(k,enkewald);
#else
    ewald_kspace_mc(k,enkewald);
#endif
  /* ------------------------------------------ */
  /* Call the subroutine to calucate the        */
  /* real-space and intramolecular correction   */
  /* terms to the electrostatic energy.         */
  /*                                            */
  /* enrewald[0] = unshifted real-space energy  */
  /* enrewald[1] = shifted real-space energy    */
  /* enrewald[2] = intra-molecular correction   */
  /* ------------------------------------------ */
  double enrewald[3];
#ifdef SPME
  enrewald[0] = eewald_real;
  enrewald[1] = eewald_real - eewald_reals;
  enrewald[2] = eewald_intra;
#else
  ewald_rspace_mc(lb,ub,k,enrewald);
#endif  
  /* ------------------------------------------ */
  /* Total and assign all the different parts   */
  /* to energy[].                               */
  /* ------------------------------------------ */ 
  energy[0] = enonbond + enrewald[0] - enrewald[2] + enkewald[2];	// non-shifted energy (LJ and coulomibic)
  energy[1] = energy[0] - enonbonds;								// shifted energy (LJ and coulomibic)
  energy[2] = enrewald[0] + enkewald[2] - enrewald[2];              //total coulomic energy 
  energy[3] = enrewald[0] - enrewald[1];							// shifted real_space ewald energy
  energy[4] = enkewald[0];                                          // k-space ewald energy
  energy[5] = enkewald[1];                                          // self correction term
  energy[6] = enrewald[2];											// intra-molecular self energy
  energy[7] = enrewald[0];                                          // real-space ewald energy(unshifted)
#endif //EWALD
}

#endif
#endif//MC
