#include <cstdio>
#include <iostream> 
#include <iomanip> 
#include <cstdlib>
#include <cmath>
#include "sphere.h"

typedef const double creal;
typedef const int cint;
typedef double real;
using namespace std;

//* -------------------------------------------------------------------------- *
//* DESCRIPTION (reference phase where f(x,y,z) < 0):                          *
//* ellipsoid/sphere inside the cube [0,1]x[0,1]x[0,1]                         *
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
  double vol_a,inv_frac;

  inv_frac = 8.;
  vol_a = 4.*MYPI*A1*B1*C1/(3.*inv_frac);

  cout << "-----------------------------------------------------------" << endl;
  cout << "------------------- CPP: 1/8 sphere check -----------------" << endl;
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
  cout << "analytical volume:  5.2359877559829882e-01" << endl;	     
  cout << "numerical  volume:  5.2359877559829937e-01" << endl << endl; 
  cout << "absolute error   :  5.5511151231257827e-16" << endl;	     
  cout << "relative error   :  1.0601848938211723e-15" << endl;         
  cout << "----------------- CPP: end 1/8 sphere check ---------------" << endl;
  cout << "-----------------------------------------------------------" << endl 
       << endl;
  
  return;
}
