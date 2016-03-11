/****************************************************************************
 * Copyright (C) 2015 by Simone Bnà(a), Sandro Manservisi(a),               *
 * Ruben Scardovelli(a), Philip Yecko(b) and Stephane Zaleski(c,d)          *
 * (a) DIN–Lab. di Montecuccolino, Università di Bologna,                   *
 *     Via dei Colli 16, 40136 Bologna, Italy                               *
 * (b) Physics Department, Cooper Union, New York, NY, USA                  *
 * (c) Sorbonne Universités, UPMC Univ Paris 06, UMR 7190,                  *
 *     Institut Jean Le Rond d’Alembert, F-75005, Paris, France             *
 * (d) CNRS, UMR 7190, Institut Jean Le Rond d’Alembert, F-75005,           *
 *     Paris, France                                                        *
 *                                                                          *
 * You should have received a copy of the CPC license along with Vofi.      *
 * If not, see http://cpc.cs.qub.ac.uk/licence/licence.html.                *
 *                                                                          *
 * e-mail: ruben.scardovelli@unibo.it                                       *
 *                                                                          *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sphere.h"

typedef const double creal;
typedef const int cint;
typedef double real;

/* -------------------------------------------------------------------------- *
 * DESCRIPTION (reference phase where f(x,y,z) < 0):                          *
 * ellipsoid/sphere inside the cube [0,1]x[0,1]x[0,1]                         *
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
  double vol_a,inv_frac;

  inv_frac = 8.;
  vol_a = 4.*MYPI*A1*B1*C1/(3.*inv_frac);

  fprintf (stdout,"------------------------------------------------------------------------------------\n");
  fprintf (stdout,"--------------------- C: 1/8 sphere check ------------------------------------------\n");
  fprintf (stdout," * sphere inside the cube [%.1f,%.1f]x[%.1f,%.1f]x[%.1f,%.1f] in a %dX%dX%d grid            *\n", X0, X0+H, Y0, Y0+H, Z0, Z0+H, NMX, NMY, NMZ);
  fprintf (stdout," * f(x,y) = c1*x^2 + c2*y^2 + c3*x*y + c4*x + c5*y - c6                            *\n");
  fprintf (stdout," * f(x,y,z) = f(x,y) + (z-ZC)^2/C2                                                 *\n");
  fprintf (stdout," * PARAMETERS:                                                                     *\n");
  fprintf (stdout," * (A1,B1,C1): semiaxis along the three ellipsoid (local) axes: x',y',z'           *\n");
  fprintf (stdout," * ALPHA: angle between two axes x' and x (in the x-y plane);                      *\n");                
  fprintf (stdout," * (XC,YC,ZC) center of the ellipsoid;                                             *\n");
  fprintf (stdout,"------------------------------------------------------------------------------------\n");
  fprintf (stdout,"a1:    %23.16e\n",A1);
  fprintf (stdout,"b1:    %23.16e\n",B1);
  fprintf (stdout,"b1:    %23.16e\n",C1);
  fprintf (stdout,"alpha: %23.16e\n",ALPHA);
  fprintf (stdout,"xc:    %23.16e\n",XC);
  fprintf (stdout,"yc:    %23.16e\n",YC);
  fprintf (stdout,"zc:    %23.16e\n",ZC);
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"analytical volume: %23.16e\n",vol_a);
  fprintf (stdout,"numerical  volume: %23.16e\n\n",vol_n);
  fprintf (stdout,"absolute error   : %23.16e\n",fabs(vol_a-vol_n));
  fprintf (stdout,"relative error   : %23.16e\n",fabs(vol_a-vol_n)/vol_a); 
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2\n");
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"analytical volume:  5.2359877559829882e-01\n");
  fprintf (stdout,"numerical  volume:  5.2359877559829937e-01\n\n");
  fprintf (stdout,"absolute error   :  5.5511151231257827e-16\n");
  fprintf (stdout,"relative error   :  1.0601848938211723e-15\n");
  fprintf (stdout,"----------------- C: end 1/8 sphere check -----------------\n");
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"\n");

  return;
}
