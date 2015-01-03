#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


typedef const double creal;
typedef const int cint;
typedef double real;

/* -------------------------------------------------------------------------- *
 * DESCRIPTION (reference phase where f(x,y) < 0):                            *
 * gaussian line in the square [0,1]x[0,1]                                    *
 * f(x,y) = y - yy0 - a0 exp[-ga (x - xx0)^2]                                 *
 * -------------------------------------------------------------------------- */

static double yy0 = 0.22;
static double a0 = 0.51;
static double xx0 = 0.541;
static double ga = 60.3;


/* -------------------------------------------------------------------------- */

void init()
{

  /* initialize random number generator. */
  srand(time(0));  
   
  double scalingfactor = 0.06;
  /* yy0 --> from 0.19 to 0.25 */
  yy0 = 0.19 + ((double)rand() / RAND_MAX)*scalingfactor;
  
  /* xx0 --> from 0.511 to 0.561 */
  xx0 = 0.511 + ((double)rand() / RAND_MAX)*scalingfactor;
  
  /* a0 --> from 0.48 to 0.54 */
  a0 = 0.48 + ((double)rand() / RAND_MAX)*scalingfactor;
  
  /* ga --> from 60.00 to 0.60.6 */
  ga = 60. + ((double)rand() / RAND_MAX)*scalingfactor;
  
  return;
}

/* -------------------------------------------------------------------------- */

real impl_func(creal xy[]) 
{
  double x,y,f0;

  x = xy[0];
  y = xy[1];

  f0 = y - yy0 - a0*exp(-ga*(x-xx0)*(x-xx0));

  return f0;
}

/* -------------------------------------------------------------------------- */

void check_area(creal area_n)
{
  real area_a;

  /* integration with MATHEMATICA with the given values of Y0,A0,GA,X0 */  
  /*area_a = 0.3364089454607542483401167;*/
  
  /* analytical integration with x in [0,1]  */
  area_a = yy0 + 0.5*a0*sqrt(M_PI/ga)*(erf(sqrt(ga)*(1.-xx0) )-erf(-sqrt(ga)*xx0));
  
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"-------------------- C: gaussian check --------------------\n");
  fprintf (stdout,"analytical area : %23.16e\n",area_a);
  fprintf (stdout,"numerical  area : %23.16e\n\n",area_n);
  fprintf (stdout,"absolute error  : %23.16e\n",fabs(area_a-area_n));
  fprintf (stdout,"relative error  : %23.16e\n",fabs(area_a-area_n)/area_a); 
  fprintf (stdout,"-----------------------------------------------------------\n");
  /*fprintf (stdout,"with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2\n");
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"analytical area :  3.3640894546075423e-01\n");
  fprintf (stdout,"numerical  area :  3.3640894546075722e-01\n\n");
  fprintf (stdout,"absolute error  :  2.9976021664879227e-15\n");
  fprintf (stdout,"relative error  :  8.9105899439812180e-15\n");
  fprintf (stdout,"------------------ C: end gaussian check ------------------\n");
  fprintf (stdout,"-----------------------------------------------------------\n");*/
  fprintf (stdout,"\n");

  return;
}
