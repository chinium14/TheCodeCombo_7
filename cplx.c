/*
** CPLX.C
** A collection of complex number routines.
**
**	struct Complex C_add (struct Complex a, struct Complex b)
**	struct Complex C_sub (struct Complex a, struct Complex b)
**	struct Complex C_mul (struct Complex a, struct Complex b)
**	struct Complex C_div (struct Complex a, struct Complex b)
**	struct Complex C_para (struct Complex a, struct Complex b)
**	   return sum, difference, product, quotient or parallel
**	   of complex quantities a and b.
**
**	struct Complex C_conj (struct Complex a)
**	struct Complex C_inv (struct Complex a)
**	   return complex conj or invert of complex quantity a
**
**	double C_mag (struct Complex a)
**	double C_ang (struct Complex a)
**	   return mag or angle (in rads) of complex quantity a
**
** copyright P. H. Anderson, Morgan State University, 5 May 96
*/
#include "defines.h"
#ifdef EWALD
#ifndef SPME
struct Complex C_conj (struct Complex a)
{
   struct Complex result;
   result.r = a.r;
   result.i = -a.i;
   return(result);
}

struct Complex C_add (struct Complex a, struct Complex b)
{
   struct Complex result;
   result.r = a.r + b.r;
   result.i = a.i + b.i;
   return(result);
}

struct Complex C_sub (struct Complex a, struct Complex b)
{
   struct Complex result;
   result.r = a.r - b.r;
   result.i = a.i - b.i;
   return(result);
}

struct Complex C_mul (struct Complex a, struct Complex b)
{
   struct Complex result;
   result.r = a.r*b.r - a.i*b.i;
   result.i = a.r*b.i + b.r*a.i;
   return(result);
}

struct Complex C_div (struct Complex a, struct Complex b)
{
   struct Complex c, result, num;
   double denom;

   c = C_conj(b);
   num = C_mul (a, c);
   denom = b.r*b.r + b.i*b.i + 1.2e-63;  /*to prevent division by zero*/

   result.r = num.r / denom;
   result.i = num.i / denom;
   return(result);
}

struct Complex C_para (struct Complex a, struct Complex b)
{
   struct Complex result, num, denom;

   num = C_mul (a, b);
   denom = C_add(a, b);

   result = C_div (num, denom);
   return(result);
}

struct Complex C_inv (struct Complex a)
{
   struct Complex result, num;
   num.r = 1.0;
   num.i = 0.0;
   result = C_div (num, a);
   return(result);
}

struct Complex C_sclmul (double s, struct Complex a){
	struct Complex result;
	result.r = s * a.r;
	result.i = s * a.i;
	return(result);
}
 double C_ang (struct Complex a)
 {
    double result;
    result = (double) atan2 ((double) a.i, (double) a.r + 1e-99);
    /* Note that 1e-99 is added to avoid computing the atan of
    ** 90 degrees
    */
    return (result);
 }

 double C_mag (struct Complex a)
 {
    double result;
    result = (double) sqrt ( (double) (a.r*a.r + a.i*a.i));
    return (result);
 }
#endif
#endif
