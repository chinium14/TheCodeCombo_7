
#ifdef MC
#ifdef SASA
#include "defines.h"



/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

double esasa_mc (int ibox)
{

  int k			= ibox;

  /* ================================================================== */
  /*                                                                    */
  /* Zero out the energy accumulator                                    */
  /*                                                                    */
  /* ================================================================== */
  double esasa  = 0.0;
#ifndef LJ
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
  /* Get and assign information for pbc and SASA parameters.            */
  /*                                                                    */
  /* ================================================================== */
  for(int a=0; a<box[k].boxns; a++)  sasa[k][a].A =sasa[k][a].S;
  double rp   = RPROBE;
  double lx = box[k].lx;
  double ly = box[k].ly;
  double lz = box[k].lz;
  double hx = box[k].hx;
  double hy = box[k].hy;
  double hz = box[k].hz;

 

  /* ================================================================== */
  /*                                                                    */
  /* Calculate the energy associated with the solvent.  The first loop  */
  /* is over atoms/sites 0 to N-1 and the second loop is over its       */
  /* neighbors.                                                         */
  /*                                                                    */
  /* Note: This loop cannot be combined with the force loop because the */
  /* force needs the atomic solvation parameter which can only be       */
  /* calculated after this loop is finished.                            */
  /*                                                                    */
  /* ================================================================== */
  for(int a=0; a<box[k].boxns-1; a++) {

   
  /* ----------------------------------------------- */
  /* Assign position of particle a.                  */
  /* ----------------------------------------------- */    
    double ax = atom[k][a].x;
    double ay = atom[k][a].y;
    double az = atom[k][a].z;

  /* ----------------------------------------------- */
  /* Loop around neighbors of particle a.            */
  /* ----------------------------------------------- */
    for(int bb=0; bb<nlist[k].count_sasa[a]; bb++) {
      int b = nlist[k].list_sasa[a][bb];

     
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
  /* accumulate the atomic solvation parameter, A.   */
  /* If the distance is greater, go to the next pair */ 
  /* in the loop.                                    */
  /* ----------------------------------------------- */
      double dr2	 = dx*dx + dy*dy + dz*dz;
	  double ra		 = sasa[k][a].r;
	  double rb		 = sasa[k][b].r;
	  double pa		 = sasa[k][a].p;
	  double pb		 = sasa[k][b].p;
	  double Sa		 = sasa[k][a].S;
	  double Sb		 = sasa[k][b].S;
	  double alpha_1 = PI*(ra+rp);			double alpha_b_1 = PI*(rb+rp);
	  double alpha_2 = (ra+rb+2.0*rp);
	  double alpha_3 = rb-ra;
	  double rc		 = alpha_2;
	  double rc2	 = rc*rc;

      if(dr2 > rc2) continue;
	  double rab = sqrt(dr2);
#ifndef STYPE
      double pab    = ljset[k].pot[a][b].pab;
#endif
#ifdef STYPE
	  double pab;
	  if(nlist[k].xflag_sasa[a][bb] == 1) pab = 0.8875;
	  else if(nlist[k].xflag_sasa[a][bb] == 3) pab = 0.0;
	  else pab = 0.3516;
#endif
	  /*-------------------------------------------------------	*/
	  /* Like the SASA implementation in CHARMM, the effect of	*/
	  /* hydrogen atoms is not considered while computing the	*/
	  /* the accessible surface area of a residue. This is done	*/
	  /* by setting pab = 0, if either a or b is a hydrogen		*/
	  /*-------------------------------------------------------	*/

	  if(atom[k][a].atomid ==1 || atom[k][a].atomid ==2 || atom[k][a].atomid ==3 || atom[k][a].atomid ==4 ||
		 atom[k][b].atomid ==1 || atom[k][b].atomid ==2 || atom[k][b].atomid ==3 || atom[k][b].atomid ==4){
		  pab=0.0;
	  }

	  double bab	= alpha_1*(alpha_2- rab)*(1.0 +alpha_3/rab);
	  double bba	= alpha_b_1*(alpha_2- rab)*(1.0 -alpha_3/rab);
	  double cab	= pa*pab/Sa;
	  double cba	= pb*pab/Sb;
	  sasa[k][a].A  *= (1.0-cab*bab);
	  sasa[k][b].A  *= (1.0-cba*bba);

	}// b loop ends here
  }// a loop ends here

  /* ================================================================== */
  /*                                                                    */
  /* Calculate the force associated with the solvent.  The first loop   */
  /* is over atoms/sites 0 to N-1 and the second loop is over its       */
  /* neighbors.                                                         */
  /*                                                                    */
  /* ================================================================== */
  for(int a=0; a<box[k].boxns-1; a++) {

		double sigA = sasa[k][a].sigma*sasa[k][a].A;

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
    for(int bb=0; bb<nlist[k].count_sasa[a]; bb++) {
      int b = nlist[k].list_sasa[a][bb];

      double fx, fy, fz;

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
  /* calculate the force. If the distance is         */
  /* greater, go to the next pair in the loop.       */ 
  /* ----------------------------------------------- */
      double dr2	 = dx*dx + dy*dy + dz*dz;
	  double ra		 = sasa[k][a].r;
	  double rb		 = sasa[k][b].r;
	  double pa		 = sasa[k][a].p;
	  double pb		 = sasa[k][b].p;
	  double Sa		 = sasa[k][a].S;
	  double Sb		 = sasa[k][b].S;
	  double alpha_1 = PI*(ra+rp);			double alpha_b_1 = PI*(rb+rp);
	  double alpha_2 = (ra+rb+2.0*rp);
	  double alpha_3 = rb-ra;
	  double rc		 = alpha_2;
	  double rc2	 = rc*rc;

      if(dr2 > rc2) continue;
	  double rab = sqrt(dr2);

#ifndef STYPE
      double pab    = ljset[k].pot[a][b].pab;
#endif
  /* ----------------------------------------------- */
  /* The value of pab changes if a and b are bonded. */ 
  /* If xflag_sasa == 1, they are bonded. If it is   */
  /* not 1, they are not bonded.                     */
  /* ----------------------------------------------- */
#ifdef STYPE
	  double pab;
	  if(nlist[k].xflag_sasa[a][bb] == 1) pab = 0.8875;
	  else if(nlist[k].xflag_sasa[a][bb] == 3) pab = 0.0;
	  else pab = 0.3516;
#endif


	  /*-------------------------------------------------------	*/
	  /* Like the SASA implementation in CHARMM, the effect of	*/
	  /* hydrogen atoms is not considered while computing the	*/
	  /* the accessible surface area of a residue. This is done	*/
	  /* by setting pab = 0, if either a or b is a hydrogen		*/
	  /*-------------------------------------------------------	*/

	  if(atom[k][a].atomid ==1 || atom[k][a].atomid ==2 || atom[k][a].atomid ==3 || atom[k][a].atomid ==4 ||
		 atom[k][b].atomid ==1 || atom[k][b].atomid ==2 || atom[k][b].atomid ==3 || atom[k][b].atomid ==4){
		  pab=0.0;
	  }

	  double bab	= alpha_1	*(alpha_2- rab)*(1.0 +alpha_3/rab);
	  double bba	= alpha_b_1	*(alpha_2- rab)*(1.0 -alpha_3/rab);
	  double cab	= pa*pab/Sa;
	  double cba	= pb*pab/Sb;
	  double dbab	= -alpha_1  *(1.0 + alpha_2*alpha_3/dr2);
	  double dbba	= -alpha_b_1*(1.0 - alpha_2*alpha_3/dr2);

	  double force	= (sigA*(cab/(1.0-cab*bab))*dbab + 
					   sasa[k][b].sigma*sasa[k][b].A*(cba/(1.0-cba*bba))*dbba)/rab;
	  

  /* ----------------------------------------- */
  /* Accumulate the forces.                    */
  /* ----------------------------------------- */
	  fx          = force * dx;
      fxa        += fx;
      ffox[k][b] -= fx;

      fy          = force * dy;
      fya        += fy;
      ffoy[k][b] -= fy;

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

#endif
    }// bloop ends here

  /* ----------------------------------------- */
  /* Assign the accumulated force from the     */
  /* implicit solvent interactions to the      */
  /* force on the atom.  Also accumulate       */
  /* the energy.                               */
  /* ----------------------------------------- */
    ffox[k][a] += fxa;
    ffoy[k][a] += fya;
    ffoz[k][a] += fza;
	esasa	   += sigA;
	

  }// aloop ends here
  esasa	   += sasa[k][box[k].boxns-1].A*sasa[k][box[k].boxns-1].sigma; 
  /* ================================================================== */
  /*                                                                    */
  /* After the loop, if PRESSURE is defined, store the bonding          */
  /* contribution to the virial tensor to pvir[k].bond[i].              */
  /*                                                                    */
  /* ================================================================== */
#ifdef PRESSURE
  pviro[k].sasa[0] = wxx;
  pviro[k].sasa[1] = wyx;
  pviro[k].sasa[2] = wyy;
  pviro[k].sasa[3] = wzx;
  pviro[k].sasa[4] = wzy;
  pviro[k].sasa[5] = wzz;

#endif
#endif // do all this only if LJ is not defined
return (esasa);

}

#endif
#endif//MC
