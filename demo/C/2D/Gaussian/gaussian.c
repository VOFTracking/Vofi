#include <stdio.h>
#include <math.h>
#include "gaussian.h"

typedef const double creal;
typedef const int cint;
typedef double real;

/* -------------------------------------------------------------------------- *
 * DESCRIPTION (reference phase where f(x,y) < 0):                            *
 * gaussian line in the square [0,1]x[0,1]                                    *
 * f(x,y) = y - YY0 - A0 exp[-GA (x - XX0)^2]                                 *
 * -------------------------------------------------------------------------- */

int cont_line(real *xx,real *yy,cint ntot)
{
  int i;
  double dx,xi;

  dx = 1./(ntot-1);
  for (i=1; i<=ntot; i++) {
    xi = (i-1)*dx;
    xx[i] = xi;
    yy[i] = YY0 + A0*exp(-GA*(xi-XX0)*(xi-XX0));
  }
  return ntot;
}

/* -------------------------------------------------------------------------- */

void check_area(creal area_n)
{
  real area_a;

  /* integration with MATHEMATICA with the given values of Y0,A0,GA,X0 */  
  area_a = 0.3364089454607542483401167;

  fprintf (stderr,"-------------------------------------------------------\n");
  fprintf (stderr,"------------------ C: gaussian check ------------------\n");
  fprintf (stderr,"analytical area : %23.16e\n",area_a);
  fprintf (stderr,"numerical  area : %23.16e\n\n",area_n);
  fprintf (stderr,"absolute error  : %23.16e\n",fabs(area_a-area_n));
  fprintf (stderr,"relative error  : %23.16e\n",fabs(area_a-area_n)/area_a); 
  fprintf (stderr,"-------------------------------------------------------\n");
  fprintf (stderr,"with Intel i7 3.4 GHz + Linux openSUSE 12.3 + gcc 4.7.2\n");
  fprintf (stderr,"-------------------------------------------------------\n");
  fprintf (stderr,"analytical area :  3.3640894546075423e-01\n");
  fprintf (stderr,"numerical  area :  3.3640894546075711e-01\n\n");
  fprintf (stderr,"absolute error  :  2.8865798640254070e-15\n");
  fprintf (stderr,"relative error  :  8.5805680942041364e-15\n");
  fprintf (stderr,"---------------- C: end gaussian check ----------------\n");
  fprintf (stderr,"-------------------------------------------------------\n");

  return;
}

/* -------------------------------------------------------------------------- */

real impl_func(creal xy[]) 
{
  double x,y,f0;

  x = xy[0];
  y = xy[1];

  f0 = y - YY0 - A0*exp(-GA*(x-XX0)*(x-XX0));

  return f0;
}
