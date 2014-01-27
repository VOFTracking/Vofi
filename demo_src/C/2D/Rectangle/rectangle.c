#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rectangle.h"

typedef const double creal;
typedef const int cint;
typedef double real;

/* -------------------------------------------------------------------------- *
 * DESCRIPTION (reference phase where f(x,y) < 0):                            *
 * rectangle inside the square [0,1]x[0,1]                                    *
 * PARAMETERS:                                                                *
 * (XC,YC): center of the rectangle; ALPHA: angle between two axes x' and x;  *
 * (A1,B1): half sides along the local x' and y' axes                         *
 * -------------------------------------------------------------------------- */

real impl_func(creal xy[]) 
{
  real x,y,f0,f1,ca,sa;

  x = xy[0];
  y = xy[1];

  ca = cos(ALPHA);
  sa = sin(ALPHA);
  f1 = - A1 - ((x-XC)*ca+(y-YC)*sa);
  f0 = f1;
  /* f0 = MAX(f0,f1); */
  f1 = ((x-XC)*ca+(y-YC)*sa) - A1;
  f0 = MAX(f0,f1);
  f1 = -B1 - ((y-YC)*ca - (x-XC)*sa);
  f0 = MAX(f0,f1);
  f1 = ((y-YC)*ca - (x-XC)*sa) - B1;
  f0 = MAX(f0,f1);

  return f0;
}

/* -------------------------------------------------------------------------- */

void check_area(creal area_n)
{
  real area_a;
  
  area_a = 4.*A1*B1;

  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"-------------------- C: rectangle check -------------------\n");
  fprintf (stdout,"analytical area : %23.16e\n",area_a);
  fprintf (stdout,"numerical  area : %23.16e\n\n",area_n);
  fprintf (stdout,"absolute error  : %23.16e\n",fabs(area_a-area_n));
  fprintf (stdout,"relative error  : %23.16e\n",fabs(area_a-area_n)/area_a); 
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2\n");
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"analytical area :  2.3999999999999999e-01\n");
  fprintf (stdout,"numerical  area :  2.3999999999999996e-01\n\n");
  fprintf (stdout,"absolute error  :  2.7755575615628914e-17\n");
  fprintf (stdout,"relative error  :  1.1564823173178715e-16\n");
  fprintf (stdout,"----------------- C: end rectangle check ------------------\n");
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"\n");

  return;
}
