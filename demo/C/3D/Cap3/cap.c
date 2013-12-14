#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cap.h"

typedef const double creal;
typedef const int cint;
typedef double real;


/* -------------------------------------------------------------------------- *
 * DESCRIPTION (reference phase where f(x,y,z) < 0):                          *
 * ellipsoidal cap inside the domain [-1,1]x[-1,1]x[0,1]                      *
 * f(x,y) = c1*x^2 + c2*y^2 + c3*x*y + c4*x + c5*y - c6                       *
 * f(x,y,z) = f(x,y) + (z-ZC)^2/C2                                            *
 * INPUT PARAMETERS:                                                          *
 * (XC,YC,ZC) center of the ellipsoid; ALPHA: angle between two axes x' and   *
 * x (in the x-y plane); (A1,B1,C1): semiaxis along the three ellipsoid       *
 * (local) axes: x',y',z'                                                     *
 * -------------------------------------------------------------------------- */

real impl_func(creal xy[])
{
  double x,y,z,A2,B2,C2,ca,sa,c1,c2,c3,c4,c5,c6,f0;

  x = xy[0];
  y = xy[1];
  z = xy[2];

  A2 = A1*A1;
  B2 = B1*B1;
  C2 = C1*C1;
  ca = cos(ALPHA);
  sa = sin(ALPHA);
  c1 = ca*ca/A2 + sa*sa/B2;
  c2 = sa*sa/A2 + ca*ca/B2;
  c3 = 2.*ca*sa*(B2-A2)/(A2*B2);
  c4 = -(2.*c1*XC + c3*YC);
  c5 = -(2.*c2*YC + c3*XC);
  c6 = 1.0 - (c1*XC*XC + c2*YC*YC + c3*XC*YC);
  
  f0 = c1*x*x + c2*y*y + c3*x*y + c4*x + c5*y - c6 + (z - ZC)*(z - ZC)/C2;

  return f0;
}

/* -------------------------------------------------------------------------- */

void check_volume(creal vol_n)
{
  double h0,vol_a;

  h0 = C1 + ZC;
  vol_a = MYPI*A1*B1*h0*h0*(1. - h0/(3.*C1))/C1;

  fprintf (stderr,"\n-----------------------------------------------------\n");
  fprintf (stderr,"--------------- C: cap check (3 cells) --------------\n\n");

  fprintf (stderr,"analytical volume: %23.16e\n",vol_a);
  fprintf (stderr,"numerical  volume: %23.16e\n\n",vol_n);
  fprintf (stderr,"absolute error   : %23.16e\n",fabs(vol_a-vol_n));
  fprintf (stderr,"relative error   : %23.16e\n\n",fabs(vol_a-vol_n)/vol_a); 
  fprintf (stderr,"------------- C: end cap check (3 cells) ------------\n");
  fprintf (stderr,"-----------------------------------------------------\n\n");

  return;
}

/* -------------------------------------------------------------------------- */

int cont_line_3D(double *xx,double *yy,double *zz, int ntot)
{
  int i;
  double pihalf,twopi,dphi,phi,r0,ca,sa,xi,yi,z0,abzc;
  double A2,B2,C2;

  A2 = A1*A1;
  B2 = B1*B1;
  C2 = C1*C1;
  ca = cos(ALPHA);
  sa = sin(ALPHA);
  pihalf = acos(0.0);

  twopi = 4.0*pihalf;
  dphi = twopi/(double)(ntot-1);
  z0 = 0.;
  abzc = A2*B2*(1. - (z0-ZC)*(z0-ZC)/C2);
  for (i=1; i<ntot; i++) {
    phi = (i-1)*dphi;
    r0 = abzc/(B2*cos(phi)*cos(phi) + A2*sin(phi)*sin(phi));
    r0 = sqrt(r0);
    xi = r0*cos(phi);
    yi = r0*sin(phi);
    xx[i] = XC + xi*ca - yi*sa;
    yy[i] = YC + xi*sa + yi*ca;
    zz[i] = z0;
  }

  xx[ntot] = xx[1];
  yy[ntot] = yy[1];
  zz[ntot] = zz[1];

  return ntot;
}
