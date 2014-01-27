#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sine_surf.h"

typedef const double creal;
typedef const int cint;
typedef double real;

/* -------------------------------------------------------------------------- *
 * DESCRIPTION (reference phase where f(x,y,z) < 0):                          *
 * sinusoidal surface inside the cube [0,1]x[0,1]x[0,1]                       *
 * f(x,y,z) = z - A0 - B0*sin(C1*pi*x + pi*D1)*sin(C1*pi*x + pi*E1)           *
 * -------------------------------------------------------------------------- */

real impl_func(creal xy[])
{
  double x,y,z,f0;

  x = xy[0];
  y = xy[1];
  z = xy[2];

  f0 = z - A0 - B0*sin(MYPI*(C1*x+D1))*sin(MYPI*(C1*y+E1));

  return f0;
}
/* -------------------------------------------------------------------------- */

void check_volume(creal vol_n)
{
  double vol_a;

  vol_a = 0.5;

  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"--------------- C: sinusoidal surface check ---------------\n");
  fprintf (stdout,"analytical volume: %23.16e\n",vol_a);
  fprintf (stdout,"numerical  volume: %23.16e\n\n",vol_n);
  fprintf (stdout,"absolute error   : %23.16e\n",fabs(vol_a-vol_n));
  fprintf (stdout,"relative error   : %23.16e\n",fabs(vol_a-vol_n)/vol_a); 
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2\n");
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"analytical volume:  5.0000000000000000e-01\n");
  fprintf (stdout,"numerical  volume:  5.0000000000000022e-01\n\n");
  fprintf (stdout,"absolute error   :  2.2204460492503131e-16\n");
  fprintf (stdout,"relative error   :  4.4408920985006262e-16\n");
  fprintf (stdout,"------------- C: end sinusoidal surface check -------------\n");
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"\n");

  return;
}
