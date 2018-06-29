/* ======================================================================== */
/* tak_histogram.h                                                          */
/*   This file contains the histogram structures and subroutine prototypes. */
/*                                                                          */
/* Written by Thomas A. Knotts IV,22 Sep 2005.                              */
/* ======================================================================== */
typedef struct{
  int n;            //Number of bins
  double xmin;      //Lower bound
  double xmax;      //Upper bound
  double bin_width; //Bin width
//  int    *count;    //Bin count
  double *vbin;     //Mid point value of bin
  double *bin;      //Accumulating value
}tak_histogram;

/*typedef struct{
  int n;          //Number of bins
  double *vbin;   //Mid point value of bin
  double *bin;    //Average
}tak_pdf;*/

tak_histogram* tak_histogram_alloc             (int n);
tak_histogram* tak_histogram_calloc            (int n);
tak_histogram* tak_histogram_calloc_uniform    (int n, double xmin, double xmax);
//int tak_histogram_set_ranges_uniform           (tak_histogram *h);
int tak_histogram_bin_find                     (tak_histogram *h,double x,int *index);
int tak_histogram_increment                    (tak_histogram *h, double x);
int tak_histogram_accumulate                   (tak_histogram *h, double x, double weight);
int tak_histogram_fwrite                       (FILE *fp, tak_histogram *h);
void tak_histogram_free                        (tak_histogram *h);
int tak_histogram_equal_bins_p                 (const tak_histogram *h1, const tak_histogram *h2);
tak_histogram* tak_histogram_clone             (const tak_histogram *src);
tak_histogram* tak_histogram_div               (const tak_histogram *h1, const tak_histogram *h2);

