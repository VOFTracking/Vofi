#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ellipse.h"

typedef const double creal;
typedef const int cint;
typedef double real;

/* -------------------------------------------------------------------------- *
 * DESCRIPTION (reference phase where f(x,y) < 0):                            *
 * ellipse inside the square [0,1]x[0,1]                                      *
 * f(x,y) = c1*x^2 + c2*y^2 + c3*x*y + c4*x + c5*y - c6                       *
 * PARAMETERS:                                                                *
 * (XC,YC): center of the ellipse; ALPHA: angle between two axes x' and x;    *
 * (A1,B1): semiaxis along the two ellipse axes (local) x' and y'             *
 * -------------------------------------------------------------------------- */

real impl_func(creal xy[]) 
{
  real x,y,a2,b2,ca,sa,c1,c2,c3,c4,c5,c6,f0;

  x = xy[0];
  y = xy[1];

  a2 = A1*A1;
  b2 = B1*B1;
  ca = cos(ALPHA);
  sa = sin(ALPHA);
  c1 = ca*ca/a2 + sa*sa/b2;
  c2 = sa*sa/a2 + ca*ca/b2;
  c3 = 2.*ca*sa*(b2-a2)/(a2*b2);
  c4 = -(2.*c1*XC + c3*YC);
  c5 = -(2.*c2*YC + c3*XC);
  c6 = 1.0 - (c1*XC*XC + c2*YC*YC + c3*XC*YC);
  
  f0 = c1*x*x + c2*y*y + c3*x*y + c4*x + c5*y - c6;

  return f0;
}

/* -------------------------------------------------------------------------- */

void check_area(creal area_n)
{
  real area_a;
  
  area_a = MYPI*A1*B1;

  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"--------------------- C: ellipse check --------------------\n");
  fprintf (stdout,"analytical area : %23.16e\n",area_a);
  fprintf (stdout,"numerical  area : %23.16e\n\n",area_n);
  fprintf (stdout,"absolute error  : %23.16e\n",fabs(area_a-area_n));
  fprintf (stdout,"relative error  : %23.16e\n",fabs(area_a-area_n)/area_a);
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2\n");
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"analytical area :  1.1215485773315563e-01\n");
  fprintf (stdout,"numerical  area :  1.1215485773315678e-01\n\n");
  fprintf (stdout,"absolute error  :  1.1518563880485999e-15\n");
  fprintf (stdout,"relative error  :  1.0270231814560841e-14\n");
  fprintf (stdout,"------------------- C: end ellipse check ------------------\n");
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"\n");

  return;
}
