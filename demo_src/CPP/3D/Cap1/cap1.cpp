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

#include <cstdio>
#include <iostream> 
#include <iomanip> 
#include <cstdlib>
#include <cmath>
#include "cap1.h"

typedef const double creal;
typedef const int cint;
typedef double real;
using namespace std;

//* -------------------------------------------------------------------------- *
//* DESCRIPTION (reference phase where f(x,y,z) < 0):                          *
//* ellipsoidal cap inside the domain [0,1]x[0,1]x[0,1]                        *
//* f(x,y) = c1*x^2 + c2*y^2 + c3*x*y + c4*x + c5*y - c6                       *
//* f(x,y,z) = f(x,y) + (z-ZC)^2/C2                                            *
//* INPUT PARAMETERS:                                                          *
//* (XC,YC,ZC) center of the ellipsoid; ALPHA: angle between two axes x' and   *
//* x (in the x-y plane); (A1,B1,C1): semiaxis along the three ellipsoid       *
//* (local) axes: x',y',z'                                                     *
//* -------------------------------------------------------------------------- *

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

//* -------------------------------------------------------------------------- *

void check_volume(creal vol_n)
{
  double h0,vol_a;

  h0 = C1 + ZC;
  vol_a = MYPI*A1*B1*h0*h0*(1. - h0/(3.*C1))/C1;

  cout << "-----------------------------------------------------------" << endl;
  cout << "----------------- CPP: cap check (1 cell) -----------------" << endl;
  cout << "analytical volume: " << scientific << setw(23) << setprecision(16) 
       << vol_a << endl;
  cout << "numerical  volume: " << scientific << setw(23) << setprecision(16) 
       << vol_n << endl << endl;
  cout << "absolute error   : " << scientific << setw(23) << setprecision(16) 
       << fabs(vol_a-vol_n) << endl;
  cout << "relative error   : " << scientific << setw(23) << setprecision(16) 
       << fabs(vol_a-vol_n)/vol_a << endl;
  cout << "-----------------------------------------------------------" << endl;
  cout << "with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2" << endl;
  cout << "-----------------------------------------------------------" << endl;
  cout << "analytical volume:  9.4090699975015856e-03" << endl;	    
  cout << "numerical  volume:  9.4090641039290823e-03" << endl << endl;
  cout << "absolute error   :  5.8935725032877029e-09" << endl;	    
  cout << "relative error   :  6.2637141660681008e-07" << endl;         
  cout << "--------------- CPP: end cap check (1 cell) ---------------" << endl;
  cout << "-----------------------------------------------------------" << endl 
       << endl;
  
  return;
}
