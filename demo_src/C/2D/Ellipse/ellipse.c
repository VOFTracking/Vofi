#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ellipse.h"


typedef const double creal;
typedef const int cint;
typedef double real;

/* -------------------------------------------------------------------------- *
 * DESCRIPTION (reference phase where f(x,y) < 0):                            *
 * ellipse inside the square [0,1]x[0,1]                                      *
 * f(x,y) = c1*x^2 + c2*y^2 + c3*x*y + c4*x + c5*y - c6                       *
 * PARAMETERS:                                                                *
 * (xc,yc): center of the ellipse; ALPHA: angle between two axes x' and x;    *
 * (a1,b1): semiaxis along the two ellipse axes (local) x' and y'             *
 * -------------------------------------------------------------------------- */

static double a1 = 0.17;
static double b1 = 0.21;
static double alpha = 0.48;
static double xc = 0.523;
static double yc = 0.475;

/* -------------------------------------------------------------------------- */

void init(cint randominput)
{

  if(randominput) {
    
    /* initialize random number generator. */
    srand(time(0));  
  
    /* a1 --> from 0.15 to 0.19 */
    double axesscalingfactor = 0.04;
    a1 = 0.15 + ((double)rand() / RAND_MAX)*axesscalingfactor;
    /* a1 --> from 0.19 to 0.23 */
    b1 = 0.19 + ((double)rand() / RAND_MAX)*axesscalingfactor;
  
    /* alpha --> from zero to pi/2 */
    double alphascalingfactor = M_PI_2;
    alpha = ((double)rand() / RAND_MAX)*alphascalingfactor;
  
    /* xc & yc --> from 0.4 to 0.6 */
    double x0 = 0.45;
    double centerscalingfactor = 0.1;
    xc = x0 + ((double)rand() / RAND_MAX)*centerscalingfactor;
    yc = x0 + ((double)rand() / RAND_MAX)*centerscalingfactor;
  
    
  }
  
  return;
}

/* -------------------------------------------------------------------------- */

real impl_func(creal xy[]) 
{
  real x,y,a2,b2,ca,sa,c1,c2,c3,c4,c5,c6,f0;

  x = xy[0];
  y = xy[1];

  a2 = a1*a1;
  b2 = b1*b1;
  ca = cos(alpha);
  sa = sin(alpha);
  c1 = ca*ca/a2 + sa*sa/b2;
  c2 = sa*sa/a2 + ca*ca/b2;
  c3 = 2.*ca*sa*(b2-a2)/(a2*b2);
  c4 = -(2.*c1*xc + c3*yc);
  c5 = -(2.*c2*yc + c3*xc);
  c6 = 1.0 - (c1*xc*xc + c2*yc*yc + c3*xc*yc);
  
  f0 = c1*x*x + c2*y*y + c3*x*y + c4*x + c5*y - c6;

  return f0;
}

/* -------------------------------------------------------------------------- */

void check_area(creal area_n, cint randominput)
{
  real area_a;
  
  area_a = M_PI*a1*b1;
  
  
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"--------------------- C: ellipse check --------------------\n");
  fprintf (stdout," * ellipse inside the square [0,1]x[0,1] in a %dX%d grid  *\n", NMX, NMY);
  fprintf (stdout," * f(x,y) = c1*x^2 + c2*y^2 + c3*x*y + c4*x + c5*y - c6   *\n");
  fprintf (stdout," * PARAMETERS:                                            *\n");
  fprintf (stdout," * (xc,yc): center of the ellipse;                        *\n");
  fprintf (stdout," * ALPHA: angle between two axes x' and x;                *\n");                
  fprintf (stdout," * (a1,b1): semiaxis along the two ellipse axes x' and y' *\n");
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"a1:    %23.16e\n",a1);
  fprintf (stdout,"b1:    %23.16e\n",b1);
  fprintf (stdout,"alpha: %23.16e\n",alpha);
  fprintf (stdout,"xc:    %23.16e\n",xc);
  fprintf (stdout,"yc:    %23.16e\n",yc);
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"analytical area : %23.16e\n",area_a);
  fprintf (stdout,"numerical  area : %23.16e\n\n",area_n);
  fprintf (stdout,"absolute error  : %23.16e\n",fabs(area_a-area_n));
  fprintf (stdout,"relative error  : %23.16e\n",fabs(area_a-area_n)/area_a);
  fprintf (stdout,"-----------------------------------------------------------\n");
  if(!randominput) {
    fprintf (stdout,"with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2\n");
    fprintf (stdout,"-----------------------------------------------------------\n");
    fprintf (stdout,"analytical area :  1.1215485773315562e-01\n");
    fprintf (stdout,"numerical  area :  1.1215485773315680e-01\n\n");
    fprintf (stdout,"absolute error  :  1.1796119636642288e-15\n");
    fprintf (stdout,"relative error  :  1.0517707279971946e-14\n");
    fprintf (stdout,"------------------- C: end ellipse check ------------------\n");
    fprintf (stdout,"-----------------------------------------------------------\n");
  }
  fprintf (stdout,"\n");

  return;
}
