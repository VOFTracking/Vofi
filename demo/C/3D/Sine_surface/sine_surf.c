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

  fprintf (stderr,"-------------------------------------------------------\n");
  fprintf (stderr,"------------- C: sinusoidal surface check -------------\n");
  fprintf (stderr,"analytical volume: %23.16e\n",vol_a);
  fprintf (stderr,"numerical  volume: %23.16e\n\n",vol_n);
  fprintf (stderr,"absolute error   : %23.16e\n",fabs(vol_a-vol_n));
  fprintf (stderr,"relative error   : %23.16e\n",fabs(vol_a-vol_n)/vol_a); 
  fprintf (stderr,"-------------------------------------------------------\n");
  fprintf (stderr,"with Intel i7 3.4 GHz + Linux openSUSE 12.3 + gcc 4.7.2\n");
  fprintf (stderr,"-------------------------------------------------------\n");
  fprintf (stderr,"analytical volume:  5.0000000000000000e-01\n");
  fprintf (stderr,"numerical  volume:  5.0000000000000022e-01\n\n");
  fprintf (stderr,"absolute error   :  2.2204460492503131e-16\n");
  fprintf (stderr,"relative error   :  4.4408920985006262e-16\n");
  fprintf (stderr,"----------- C: end sinusoidal surface check -----------\n");
  fprintf (stderr,"-------------------------------------------------------\n");

  return;
}
/* -------------------------------------------------------------------------- */

int cont_line_3D(double *xx,double *yy,double *zz, int ntot)
{
  int i,j;
  double mypi,dx,dy,x0,xi,y0,yi,z0,dz0,alx,aly;

  mypi = 2.*acos(0.0);
  z0 = 0.5;
  dz0 = 1./6.;

  dx = 1./(ntot);
  dy = dx;
  j=1;

  /* y=0, 0 <= x <= 1 */
  y0 = 0.; x0 = 0.;
  for (i=0; i<=ntot; i++) {
    xi = x0 + i*dx;
    yi = y0;
    xx[j] = xi;
    yy[j] = yi;
    alx = mypi*(1.6*xi + 1./7.);
    aly = mypi*(1.6*yi + 1./5.);
    zz[j] = z0 + dz0*sin(alx)*sin(aly);
    j++;
  }

  /* x=1, 0 <= y <= 1 */
  y0 = 0.; x0 = 1.;
  for (i=0; i<=ntot; i++) {
    xi = x0;
    yi = y0 + i*dy;
    xx[j] = xi;
    yy[j] = yi;
    alx = mypi*(1.6*xi + 1./7.);
    aly = mypi*(1.6*yi + 1./5.);
    zz[j] = z0 + dz0*sin(alx)*sin(aly);
    j++;
  }

  /* y=1, 0 <= x <= 1 */
  y0 = 1.; x0 = 1.;
  for (i=0; i<=ntot; i++) {
    xi = x0 - i*dx;
    yi = y0;
    xx[j] = xi;
    yy[j] = yi;
    alx = mypi*(1.6*xi + 1./7.);
    aly = mypi*(1.6*yi + 1./5.);
    zz[j] = z0 + dz0*sin(alx)*sin(aly);
    j++;
  }

  /* x=0, 0 <= y <= 1 */
  y0 = 1.; x0 = 0.;
  for (i=0; i<=ntot; i++) {
    xi = x0;
    yi = y0 - i*dy;
    xx[j] = xi;
    yy[j] = yi;
    alx = mypi*(1.6*xi + 1./7.);
    aly = mypi*(1.6*yi + 1./5.);
    zz[j] = z0 + dz0*sin(alx)*sin(aly);
    j++;
  }
  
  ntot = 4*(ntot+1);

  return ntot;
}
