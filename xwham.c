#ifdef XWHAM
/* ======================================================================== */
/* xwham.cpp                                                                */
/*   This file contains multiple subroutines.                               */
/* xwham_init                                                               */
/*   Controls the initialization of the histograms for each of the values   */
/*   in order#.output.                                                      */
/* xwham_accumulate                                                         */
/*   Driver to accumulate each of the histograms.                           */
/* xwham_average                                                            */
/*   Driver to calculate the average of the accumulated properties.         */
/* xwham_write                                                              */
/*   Driver to write the histograms to file.                                */
/*                                                                          */
/* Written by Thomas A. Knotts IV, 3 Jun 2008.                              */
/* ======================================================================== */

#include "defines.h"

void xwham_init(void){
  FILE *fp;
  char buf[150];
  char key1[3]="LB";
  char key2[3]="UB";
  char key3[3]="BN";
  int flag1 = 0;
  int flag2 = 0;
  int flag3 = 0;
  double lb, ub, bn;

  /* =============================================== */
  /* Read in the histogram ranges and bin size.      */
  /* =============================================== */
  if(NULL == (fp=fopen("./INPUT/xwham.inp","r"))){
    fprintf(stdout,"Input file xwham.inp does now exist!!!\n");
    exit(1);
  }

  while(fscanf(fp,"%s",buf)!=EOF){
    if(strcmp(buf,key1)==0){
      fscanf(fp,"%lf",&lb);  fgets(buf,150,fp);
      flag1 = 1;
    }
    else if(strcmp(buf,key2)==0){
      fscanf(fp,"%lf",&ub);  fgets(buf,150,fp);
      flag2 = 1;
    }
    else if(strcmp(buf,key3)==0){
      fscanf(fp,"%lf",&bn);  fgets(buf,150,fp);
      flag3 = 1;
    }
    else fgets(buf,150,fp);    
  }

  if (flag1 != 1){
    fprintf(stdout,"Could not find LB in xwham.inp!!!\n");
    exit(1);
  }

  if (flag2 != 1){
    fprintf(stdout,"Could not find UB in xwham.inp!!!\n");
    exit(1);
  }

  if (flag3 != 1){
    fprintf(stdout,"Could not find BN in xwham.inp!!!\n");
    exit(1);
  }

  /* =============================================== */
  /* Allocate the historams for each box with the    */
  /* range [lb,ub).                              */
  /* =============================================== */
  int n = (int)((ub-lb)/bn);
  
  for(int k=0; k<sim.NB; k++){
    h_ebond[k] = tak_histogram_calloc_uniform(n,lb,ub);
    h_ebend[k] = tak_histogram_calloc_uniform(n,lb,ub);
    h_eurey[k] = tak_histogram_calloc_uniform(n,lb,ub);
    h_etors[k] = tak_histogram_calloc_uniform(n,lb,ub);
    h_eimpr[k] = tak_histogram_calloc_uniform(n,lb,ub);
    h_elj[k]   = tak_histogram_calloc_uniform(n,lb,ub);
    h_eqq[k]   = tak_histogram_calloc_uniform(n,lb,ub);
    h_pe[k]    = tak_histogram_calloc_uniform(n,lb,ub);
    h_ke[k]    = tak_histogram_calloc_uniform(n,lb,ub);
    h_etot[k]  = tak_histogram_calloc_uniform(n,lb,ub);
    h_hel[k]   = tak_histogram_calloc_uniform(n,lb,ub);
    h_dnc[k]   = tak_histogram_calloc_uniform(n,lb,ub);
    h_con[k]   = tak_histogram_calloc_uniform(n,lb,ub);
    h_con_2[k] = tak_histogram_calloc_uniform(n,lb,ub);
    h_x1[k]    = tak_histogram_calloc_uniform(n,lb,ub);
    h_x2[k]    = tak_histogram_calloc_uniform(n,lb,ub);
    h_rmsd[k]  = tak_histogram_calloc_uniform(n,lb,ub);
    h_gyr[k]   = tak_histogram_calloc_uniform(n,lb,ub);
  }
}

static void xwham_file_write(int ibox, char file_name[], tak_histogram *h){
  const int k = ibox;
  FILE *fp;
  char path1[80];
  char path2[80];
  sprintf(path1,"./OUTPUT/BOX%d/XWHAM/",k);
  sprintf(path2,"_hist%d.dat",k);
  strcat(strcat(path1,file_name),path2);
  fp = fopen(path1, "w");
  tak_histogram_fwrite(fp,h);
  fclose(fp);
}

void xwham_write(int ibox, int flag){

  const int k = ibox;

  if(flag == 0){
    const double max_io_time_out = 1800;
    time_t now;
		time(&now);
    if(difftime(now,wham_time[k]) < max_io_time_out) return;
  }

  xwham_file_write(k,"pe",h_pe[k]);
  xwham_file_write(k,"dnc"  ,h_dnc[k]);
  xwham_file_write(k,"con"  ,h_con[k]);
  xwham_file_write(k,"con_2",h_con_2[k]);
  xwham_file_write(k,"x1"   ,h_x1[k]);
  xwham_file_write(k,"gry"  ,h_gyr[k]);
  time(&wham_time[k]);

}
void xwham_accumulate(int ibox){
  const int k = ibox;
  struct energy *enk = &en[k];
  struct struct_quant *ordk = &ordparam[k];
  const double xi = ordk->d_nc;
  tak_histogram_accumulate(h_ebond[k],xi,enk->bond);
  tak_histogram_accumulate(h_ebend[k],xi,enk->bend);
  tak_histogram_accumulate(h_eurey[k],xi,enk->urey);
  tak_histogram_accumulate(h_etors[k],xi,enk->tors);
  tak_histogram_accumulate(h_eimpr[k],xi,enk->impr);
  tak_histogram_accumulate(h_elj[k]  ,xi,(enk->nbonds-enk->coulomb));
  tak_histogram_accumulate(h_eqq[k]  ,xi,enk->coulomb);
  tak_histogram_accumulate(h_pe[k]   ,xi,enk->potens);
  tak_histogram_accumulate(h_ke[k]   ,xi,enk->kinet);
  tak_histogram_accumulate(h_etot[k] ,xi,enk->totals);
  tak_histogram_accumulate(h_hel[k]  ,xi,ordk->hel);
  tak_histogram_increment (h_dnc[k]  ,xi);  
  tak_histogram_accumulate(h_con[k]  ,xi,ordk->con);
  tak_histogram_accumulate(h_con_2[k],xi,ordk->con_2); 
  tak_histogram_accumulate(h_x1[k]   ,xi,ordk->x1);
  tak_histogram_accumulate(h_x2[k]   ,xi,ordk->x2);
  tak_histogram_accumulate(h_rmsd[k] ,xi,ordk->rmsd);
  tak_histogram_accumulate(h_gyr[k]  ,xi,ordk->gyr);
}



#endif
