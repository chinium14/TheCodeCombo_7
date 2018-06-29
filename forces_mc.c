
#ifdef MC

#include "defines.h"


  /* ================================================================== */
  /*                                                                    */
  /* Function Prototypes                                                */
  /*                                                                    */
  /* ================================================================== */
double ebond_mc (int);
double ebend_mc (int);
#ifndef NEUTRAL
double eurey_mc (int);
#endif
double etors_mc (int);
double eimpr_mc (int);
#ifdef NLIST
void   enbnd_mc (int,int,int,double*);
#ifdef SASA								// SASA is only defined when neighbor list is defined
double esasa_mc(int);
#endif
#endif
#ifdef REST
double erest_mc(int);
#endif
#ifdef GOLIK
  double ebb_golik_mc(int);

  #ifdef GOBT
    double eba_golik_mc(int);
    double eda_golik_mc(int);
  #endif

  #ifndef DNA_GOLIK
    void enbnd_golik_mc(int,int,int, double*);
  #else
    void enbnd_dnagolik_mc(int,int,int,double*);
  #endif
#endif
#ifdef WALL
double ewall_mc(int);
#endif


#ifdef CONFIGT

double ebond2_mc (int);
double ebend2_mc (int);
#ifndef NEUTRAL
double eurey2_mc(int);
#endif
double etorsion2_mc(int);
double eimproper2_mc(int);
void   enbnd2_mc (int,int,int, double*);
#ifdef SASA
double esasa2_mc(int);
#endif//SASA
#ifdef REST
double erest2_mc(int);
#endif
#ifdef GOLIK
  double ebb_golik2_mc(int);

  #ifdef GOBT
    double eba_golik2_mc(int);
    double eda_golik2_mc(int);
  #endif

  #ifndef DNA_GOLIK
    void enbnd_golik2_mc(int,int,int, double*);
  #else
    void enbnd_dnagolik2_mc(int,int,int,double*);
  #endif
#endif//GOLIK
#ifdef WALL
double ewall2_mc(int);
#endif

#endif//CONFIGT


/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */

void forces_mc (int lb, int ub, int tag, int ibox)
{
  int k = ibox;

  /* ================================================================== */
  /* Zero out atomic forces arrays.                                     */
  /* ================================================================== */
  if(tag==0){
	  for(int i=0; i<box[k].boxns; i++) {
		  ffox[k][i] = 0.0;
		  ffoy[k][i] = 0.0;
		  ffoz[k][i] = 0.0;
	  }

    #ifdef CONFIGT
    hesox = 0.0;
    hesoy = 0.0;
    hesoz = 0.0;
    hesor = 0.0;
    #endif


	  /* ================================================================== */
	  /* Zero out the virial accumulators                                   */
	  /* ================================================================== */
	  #ifdef PRESSURE
	  for(int i=0; i<6; i++) {
		  pviro[k].nbond[i] = 0.0;
		  pviro[k].bond[i]  = 0.0;
		  pviro[k].bend[i]  = 0.0;
	    #ifndef NEUTRAL
		  pviro[k].urey[i]	 = 0.0;
	    #endif
		  pviro[k].tors[i]  = 0.0;
		  pviro[k].impr[i]  = 0.0;
	    #ifdef SASA
		  pviro[k].sasa[i]  = 0.0;
	    #endif//SASA
	    #ifdef EWALD
		  pviro[k].ewald_real[i] = 0.0;
		  pviro[k].ewald_kspace[i] = 0.0;
		  pviro[k].ewald_intra[i] = 0.0;
	    #endif//EWALD
	  }
	  #endif//PRESSURE
  
	  /* ================================================================== */
	  /* Call subroutines that calculate energy and forces from the         */
	  /* different contributions to the force field.                        */
	  /* ================================================================== */
	  #ifndef CONFIGT
      #ifndef GOLIK

	      enmc[k].o_bond = ebond_mc(k);
	      #ifdef REST
	      enmc[k].o_bond += erest_mc(k);
	      #endif
	      enmc[k].o_bend = ebend_mc(k);
	      #ifndef NEUTRAL
	      enmc[k].o_urey = eurey_mc(k);
	      #endif
	      enmc[k].o_tors = etors_mc(k);
	      enmc[k].o_impr = eimpr_mc(k);
	      double interen[8];
	      enbnd_mc(lb,ub,k,interen);
	      enmc[k].o_nbond  = interen[0];
	      enmc[k].o_nbonds = interen[1];
	      #ifdef COULOMB
	      enmc[k].o_coulomb = interen[2];
	      #endif
	      #ifdef EWALD
	      enmc[k].o_coulomb = interen[2];
	      enmc[k].o_tewald  = interen[2];
	      enmc[k].o_rewald  = interen[7];
	      enmc[k].o_rewalds = interen[3];
	      enmc[k].o_kewald  = interen[4];
	      enmc[k].o_sewald  = interen[5];
	      enmc[k].o_iewald  = interen[6];
	      #endif
	      #ifdef SASA
	      enmc[k].o_esasa = esasa_mc(k);
	      #endif

      #else // #ifndef GOLIK
    
        double interen[6];
        enmc[k].o_bond = ebb_golik_mc(k);	//bond
	#ifdef REST
	enmc[k].o_bond += erest_mc(k);
	#endif
        #ifdef GOBT
        enmc[k].o_bond += eba_golik_mc(k);	//bend
        enmc[k].o_bond += eda_golik_mc(k);	//tors
        #endif
        #ifndef DNA_GOLIK
        enbnd_golik_mc(lb,ub,k,interen);
        enmc[k].o_nbond	= interen[0];
        enmc[k].o_nbonds	= interen[1];
        enmc[k].o_bend	= interen[2];	// EPP
        enmc[k].o_urey	= interen[3];	// EPS
        enmc[k].o_impr	= interen[4];   // ESS
        #else
        enbnd_dnagolik_mc(lb,ub,k,interen);
        enmc[k].o_nbond	= interen[0];
        enmc[k].o_nbonds	= interen[1];
        enmc[k].o_bend	= interen[2];	// enative
        enmc[k].o_urey	= interen[3];	// enonnative+emm
        enmc[k].o_impr	= interen[4]; // ebp
        enmc[k].o_coulomb = interen[5];
        #endif
        #ifdef WALL
        enmc[k].o_tors	= ewall_mc(k);
        #endif

      #endif //#ifndef GOLIK


	  #else //#ifndef CONFIGT


      #ifndef GOLIK

	      enmc[k].o_bond = ebond2_mc(k);
        #ifdef REST
        enmc[k].o_bond += erest2_mc(k);
        #endif
	      enmc[k].o_bend = ebend2_mc(k);
	      #ifndef NEUTRAL
	      enmc[k].o_urey = eurey2_mc(k);
	      #endif 
	      enmc[k].o_tors = etorsion2_mc(k);
	      enmc[k].o_impr = eimproper2_mc(k);
	      double interen[8];
	      #ifdef NLIST
	      enbnd2_mc(lb,ub,k,interen);
	      #endif
	      enmc[k].o_nbond  = interen[0];
	      enmc[k].o_nbonds = interen[1];
	      #ifdef COULOMB
	      enmc[k].o_coulomb = interen[2];
	      #endif
	      #ifdef EWALD
	      enmc[k].o_coulomb = interen[2];
	      enmc[k].o_tewald  = interen[2];
	      enmc[k].o_rewald  = interen[7];
	      enmc[k].o_rewalds = interen[3];
	      enmc[k].o_kewald  = interen[4];
	      enmc[k].o_sewald  = interen[5];
	      enmc[k].o_iewald  = interen[6];
	      #endif
	      #ifdef SASA
	      enmc[k].o_esasa = esasa2_mc(k);
	      #endif

      #else //#ifndef GOLIK

        double interen[6];
        enmc[k].o_bond = ebb_golik2_mc(k);	//bond
        #ifdef REST
        enmc[k].o_bond += erest2_mc(k);
        #endif
        #ifdef GOBT
        enmc[k].o_bond += eba_golik2_mc(k);	//bend
        enmc[k].o_bond += eda_golik2_mc(k);	//tors
        #endif
        #ifndef DNA_GOLIK
        enbnd_golik2_mc(lb,ub,k,interen);
        enmc[k].o_nbond	= interen[0];
        enmc[k].o_nbonds	= interen[1];
        enmc[k].o_bend	= interen[2];	// EPP
        enmc[k].o_urey	= interen[3];	// EPS
        enmc[k].o_impr	= interen[4];   // ESS
        #else
        enbnd_dnagolik2_mc(lb,ub,k,interen);
        enmc[k].o_nbond	= interen[0];
        enmc[k].o_nbonds	= interen[1];
        enmc[k].o_bend	= interen[2];	// enative
        enmc[k].o_urey	= interen[3];	// enonnative+emm
        enmc[k].o_impr	= interen[4]; // ebp
        enmc[k].o_coulomb = interen[5];
        #endif
        #ifdef WALL
        enmc[k].o_tors	= ewall2_mc(k);
        #endif

      #endif //#ifndef GOLIK

	  #endif//ConfigT

    #ifdef GOLIK
      enmc[k].o_poten  = enmc[k].o_bond + enmc[k].o_tors + enmc[k].o_nbond;
      enmc[k].o_potens = enmc[k].o_bond + enmc[k].o_tors + enmc[k].o_nbonds;
    #else
	    #ifndef SASA
        enmc[k].o_poten  = enmc[k].o_bond + enmc[k].o_bend + enmc[k].o_tors + enmc[k].o_nbond + enmc[k].o_impr;
        enmc[k].o_potens = enmc[k].o_bond + enmc[k].o_bend + enmc[k].o_tors + enmc[k].o_nbonds + enmc[k].o_impr;
      #else
        enmc[k].o_poten  = enmc[k].o_bond + enmc[k].o_bend + enmc[k].o_tors + enmc[k].o_nbond+ enmc[k].o_esasa + enmc[k].o_impr;
        enmc[k].o_potens = enmc[k].o_bond + enmc[k].o_bend + enmc[k].o_tors + enmc[k].o_nbonds+ enmc[k].o_esasa + enmc[k].o_impr;
      #endif

      #ifndef NEUTRAL
        enmc[k].o_poten  += enmc[k].o_urey;
        enmc[k].o_potens += enmc[k].o_urey;
      #endif

    #endif

  }
  else{	//if tag ==1 i.e. fill new post move arrays

    for(int i=0; i<box[k].boxns; i++) {
		  ffox[k][i] = - ffox[k][i];
		  ffoy[k][i] = - ffoy[k][i];
		  ffoz[k][i] = - ffoz[k][i];
    }

	  #ifdef PRESSURE

    struct virial pvir_before;

    	pvir_before	= pviro[k];


/*	  for(int i=0; i<6; i++) {
		  pviro[k].nbond[i] = -pviro[k].nbond[i];
		  pviro[k].bond[i]  = -pviro[k].bond[i];
		  pviro[k].bend[i]  = -pviro[k].bend[i];
	    #ifndef NEUTRAL
		  pviro[k].urey[i]	 = -pviro[k].urey[i];
	    #endif
		  pviro[k].tors[i]  = -pviro[k].tors[i];
		  pviro[k].impr[i]  = -pviro[k].impr[i];
	    #ifdef SASA
		  pviro[k].sasa[i]  = -pviro[k].sasa[i];
	    #endif//SASA
	    #ifdef EWALD
		  pviro[k].ewald_real[i] = -pviro[k].ewald_real[i];
		  pviro[k].ewald_kspace[i] = -pviro[k].ewald_kspace[i];
		  pviro[k].ewald_intra[i] = -pviro[k].ewald_intra[i];
	    #endif//EWALD
	  }
    */
	  #endif//PRESSURE

    #ifdef CONFIGT
    hesox = -hesox;
    hesoy = -hesoy;
    hesoz = -hesoz;
    hesor = -hesor;
    #endif
    

  
	  /* ================================================================== */
	  /*                                                                    */
	  /* Call subroutines that calculate energy and forces from the         */
	  /* different contributions to the force field.                        */
	  /*                                                                    */
	  /* ================================================================== */
	  #ifndef CONFIGT
      #ifndef GOLIK

	      enmc[k].n_bond = ebond_mc(k);
	      #ifdef REST
	      enmc[k].n_bond += erest_mc(k);
	      #endif
	      enmc[k].n_bend = ebend_mc(k);
	      #ifndef NEUTRAL
	      enmc[k].n_urey = eurey_mc(k);
	      #endif
	      enmc[k].n_tors = etors_mc(k);
	      enmc[k].n_impr = eimpr_mc(k);
	      double interen[8];
	      enbnd_mc(lb,ub,k,interen);   
	      enmc[k].n_nbond  = interen[0];
	      enmc[k].n_nbonds = interen[1];
	      #ifdef COULOMB
	      enmc[k].n_coulomb = interen[2];
	      #endif
	      #ifdef EWALD
	      enmc[k].n_coulomb = interen[2];
	      enmc[k].n_tewald  = interen[2];
	      enmc[k].n_rewald  = interen[7];
	      enmc[k].n_rewalds = interen[3];
	      enmc[k].n_kewald  = interen[4];
	      enmc[k].n_sewald  = interen[5];
	      enmc[k].n_iewald  = interen[6];
	      #endif
	      #ifdef SASA
	      enmc[k].n_esasa = esasa_mc(k);
	      #endif
      #else //#ifndef GOLIK
        double interen[6];
        enmc[k].n_bond = ebb_golik_mc(k);	//bond
	#ifdef REST
	enmc[k].n_bond += erest_mc(k);
	#endif
        #ifdef GOBT
        enmc[k].n_bond += eba_golik_mc(k);	//bend
        enmc[k].n_bond += eda_golik_mc(k);	//tors
        #endif
        #ifndef DNA_GOLIK
        enbnd_golik_mc(lb,ub,k,interen);
        enmc[k].n_nbond	= interen[0];
        enmc[k].n_nbonds	= interen[1];
        enmc[k].n_bend	= interen[2];	// EPP
        enmc[k].n_urey	= interen[3];	// EPS
        enmc[k].n_impr	= interen[4];   // ESS
        #else
        enbnd_dnagolik_mc(lb,ub,k,interen);
        enmc[k].n_nbond	= interen[0];
        enmc[k].n_nbonds	= interen[1];
        enmc[k].n_bend	= interen[2];	// enative
        enmc[k].n_urey	= interen[3];	// enonnative+emm
        enmc[k].n_impr	= interen[4]; // ebp
        enmc[k].n_coulomb = interen[5];
        #endif
        #ifdef WALL
        enmc[k].n_tors	= ewall_mc(k);
        #endif

      #endif //#ifndef GOLIK


    #else //#ifndef CONFIGT


      #ifndef GOLIK

	      enmc[k].n_bond = ebond2_mc(k);
        #ifdef REST
        enmc[k].n_bond += erest2_mc(k);
        #endif
	      enmc[k].n_bend = ebend2_mc(k);
	      #ifndef NEUTRAL
	      enmc[k].n_urey = eurey2_mc(k);
	      #endif 
	      enmc[k].n_tors = etorsion2_mc(k);
	      enmc[k].n_impr = eimproper2_mc(k);
	      double interen[8];
	      #ifdef NLIST
	      enbnd2_mc(lb,ub,k,interen);
	      #endif
	      enmc[k].n_nbond  = interen[0];
	      enmc[k].n_nbonds = interen[1];
	      #ifdef COULOMB
	      enmc[k].n_coulomb = interen[2];
	      #endif
	      #ifdef EWALD
	      enmc[k].n_coulomb = interen[2];
	      enmc[k].n_tewald  = interen[2];
	      enmc[k].n_rewald  = interen[7];
	      enmc[k].n_rewalds = interen[3];
	      enmc[k].n_kewald  = interen[4];
	      enmc[k].n_sewald  = interen[5];
	      enmc[k].n_iewald  = interen[6];
	      #endif
	      #ifdef SASA
	      enmc[k].n_esasa = esasa2_mc(k);
	      #endif
      #else //#ifndef GOLIK

        double interen[6];
        enmc[k].n_bond = ebb_golik2_mc(k);	//bond
        #ifdef REST
        enmc[k].n_bond += erest2_mc(k);
        #endif
        #ifdef GOBT
        enmc[k].n_bond += eba_golik2_mc(k);	//bend
        enmc[k].n_bond += eda_golik2_mc(k);	//tors
        #endif
        #ifndef DNA_GOLIK
        enbnd_golik2_mc(lb,ub,k,interen);
        enmc[k].n_nbond	= interen[0];
        enmc[k].n_nbonds	= interen[1];
        enmc[k].n_bend	= interen[2];	// EPP
        enmc[k].n_urey	= interen[3];	// EPS
        enmc[k].n_impr	= interen[4];   // ESS
        #else
        enbnd_dnagolik2_mc(lb,ub,k,interen);
        enmc[k].n_nbond	= interen[0];
        enmc[k].n_nbonds	= interen[1];
        enmc[k].n_bend	= interen[2];	// enative
        enmc[k].n_urey	= interen[3];	// enonnative+emm
        enmc[k].n_impr	= interen[4]; // ebp
        enmc[k].n_coulomb = interen[5];
        #endif
        #ifdef WALL
        enmc[k].n_tors	= ewall2_mc(k);
        #endif

      #endif //#ifndef GOLIK

	  #endif//ConfigT


    #ifdef GOLIK
      enmc[k].n_poten  = enmc[k].n_bond + enmc[k].n_tors + enmc[k].n_nbond;
      enmc[k].n_potens = enmc[k].n_bond + enmc[k].n_tors + enmc[k].n_nbonds;
    #else
	    #ifndef SASA
		    enmc[k].n_poten  = enmc[k].n_bond + enmc[k].n_bend + enmc[k].n_tors + enmc[k].n_nbond + enmc[k].n_impr;
		    enmc[k].n_potens = enmc[k].n_bond + enmc[k].n_bend + enmc[k].n_tors + enmc[k].n_nbonds + enmc[k].n_impr;
	    #else
		    enmc[k].n_poten  = enmc[k].n_bond + enmc[k].n_bend + enmc[k].n_tors + enmc[k].n_nbond+ enmc[k].n_esasa + enmc[k].n_impr;
		    enmc[k].n_potens = enmc[k].n_bond + enmc[k].n_bend + enmc[k].n_tors + enmc[k].n_nbonds+ enmc[k].n_esasa + enmc[k].n_impr;
	    #endif

	    #ifndef NEUTRAL
		    enmc[k].n_poten  += enmc[k].n_urey;
		    enmc[k].n_potens += enmc[k].n_urey;
	    #endif

    #endif


		/* ---------------------------------------------------	*/
		/* Update the energies, forces, virials and hessians. 	*/
		/* ---------------------------------------------------	*/

		en[k].bond	+=enmc[k].n_bond	- enmc[k].o_bond;
		en[k].bend	+=enmc[k].n_bend	- enmc[k].o_bend;
		en[k].tors	+=enmc[k].n_tors	- enmc[k].o_tors;
		en[k].impr	+=enmc[k].n_impr	- enmc[k].o_impr;
		en[k].nbond	+=enmc[k].n_nbond	- enmc[k].o_nbond;
		en[k].nbonds	+=enmc[k].n_nbonds	- enmc[k].o_nbonds;
    #ifdef COULOMB
		en[k].coulomb += enmc[k].n_coulomb - enmc[k].o_coulomb;
    #endif
    #ifdef SASA
		en[k].esasa	+=enmc[k].n_esasa	- enmc[k].o_esasa;
    #endif
    #ifndef NEUTRAL
		en[k].urey	+=enmc[k].n_urey	- enmc[k].o_urey;
    #endif
    #ifdef EWALD
		en[k].coulomb += enmc[k].n_coulomb - enmc[k].o_coulomb;
		en[k].tewald += enmc[k].n_tewald - enmc[k].o_tewald;
		en[k].rewald += enmc[k].n_rewald - enmc[k].o_rewald;
		en[k].rewalds += enmc[k].n_rewalds - enmc[k].o_rewalds;
		en[k].kewald += enmc[k].n_kewald - enmc[k].o_kewald;
		en[k].sewald += enmc[k].n_sewald - enmc[k].o_sewald;
		en[k].iewald += enmc[k].n_iewald - enmc[k].o_iewald;
    #endif


	  for(int i=0; i<box[k].boxns; i++){
		  ff[k][i].x += ffox[k][i];
		  ff[k][i].y += ffoy[k][i];
		  ff[k][i].z += ffoz[k][i];
	  }
    #ifdef CONFIGT
    config[k].hesx += hesox;
    config[k].hesy += hesoy;
    config[k].hesz += hesoz;
    config[k].hesr += hesor;
    #endif

	  #ifdef PRESSURE
	  for(int i=0; i<6; i++) {
		  pvir[k].nbond[i]        += (pviro[k].nbond[i] - pvir_before.nbond[i]);
		  pvir[k].bond[i]         += (pviro[k].bond[i] - pvir_before.bond[i]);
		  pvir[k].bend[i]         += (pviro[k].bend[i] - pvir_before.bend[i]);
	    #ifndef NEUTRAL
		  pvir[k].urey[i]         += (pviro[k].urey[i] - pvir_before.urey[i]);
	    #endif
		  pvir[k].tors[i]         += (pviro[k].tors[i] - pvir_before.tors[i]);
		  pvir[k].impr[i]         += (pviro[k].impr[i] - pvir_before.impr[i]);
	    #ifdef SASA
		  pvir[k].sasa[i]         += (pviro[k].sasa[i] - pvir_before.sasa[i]);
	    #endif//SASA
	    #ifdef EWALD
		  pvir[k].ewald_real[i]   += (pviro[k].ewald_real[i] - pvir_before.ewald_real[i]);
		  pvir[k].ewald_kspace[i] += (pviro[k].ewald_kspace[i] - pvir_before.ewald_kspace[i]);
		  pvir[k].ewald_intra[i]  += (pviro[k].ewald_intra[i] - pvir_before.ewald_intra[i]);
	    #endif//EWALD
	  }
	  #endif//PRESSURE

  }// end else tag==1


}
#endif//MC

