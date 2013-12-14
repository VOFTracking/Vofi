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

int cont_line(real *xx,real *yy,cint ntot)
{
  real ca,sa,xi,yi;

  ca = cos(ALPHA);
  sa = sin(ALPHA);
  xi = - A1;
  yi = - B1;
  xx[1] = XC + xi*ca - yi*sa;
  yy[1] = YC + xi*sa + yi*ca;
  xi = A1;
  xx[2] = XC + xi*ca - yi*sa;
  yy[2] = YC + xi*sa + yi*ca;
  yi = B1;
  xx[3] = XC + xi*ca - yi*sa;
  yy[3] = YC + xi*sa + yi*ca;
  xi = -A1;
  xx[4] = XC + xi*ca - yi*sa;
  yy[4] = YC + xi*sa + yi*ca;
  xx[5] = xx[1];
  yy[5] = yy[1];

  return 5;
}

/* -------------------------------------------------------------------------- */

void check_area(creal area_n)
{
  real area_a;
  
  area_a = 4.*A1*B1;

  fprintf (stderr,"-----------------------------------------------------\n");
  fprintf (stderr,"----------------- C: rectangle check ----------------\n\n");

  fprintf (stderr,"analytical area: %23.16e\n",area_a);
  fprintf (stderr,"numerical  area: %23.16e\n\n",area_n);
  fprintf (stderr,"absolute error  : %23.16e\n",fabs(area_a-area_n));
  fprintf (stderr,"relative error  : %23.16e\n",fabs(area_a-area_n)/area_a); 
  fprintf (stderr,"-------------- C: end rectangle check ---------------\n");
  fprintf (stderr,"-----------------------------------------------------\n\n");

  return;
}

/* -------------------------------------------------------------------------- */

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
