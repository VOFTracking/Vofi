#include <cstdio>
#include <iostream> 
#include <iomanip> 
#include <cstdlib>
#include <cmath>
#include "sine_surf.h"

typedef const double creal;
typedef const int cint;
typedef double real;
using namespace std;

//* -------------------------------------------------------------------------- *
//* DESCRIPTION (reference phase where f(x,y,z) < 0):                          *
//* sinusoidal surface inside the cube [0,1]x[0,1]x[0,1]                       *
//* f(x,y,z) = z - A0 - B0*sin(C1*pi*x + pi*D1)*sin(C1*pi*x + pi*E1)           *
//* -------------------------------------------------------------------------- *

real impl_func(creal xy[])
{
  double x,y,z,f0;

  x = xy[0];
  y = xy[1];
  z = xy[2];

  f0 = z - A0 - B0*sin(MYPI*(C1*x+D1))*sin(MYPI*(C1*y+E1));

  return f0;
}

//* -------------------------------------------------------------------------- *

void check_volume(creal vol_n)
{
  double vol_a;

  vol_a = 0.5;

  cout << "-----------------------------------------------------------" << endl;
  cout << "--------------- CPP: sinusoidal surface check -------------" << endl;
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
  cout << "analytical volume:  5.0000000000000000e-01" << endl;	     
  cout << "numerical  volume:  5.0000000000000011e-01" << endl << endl; 
  cout << "absolute error   :  1.1102230246251565e-16" << endl;	     
  cout << "relative error   :  2.2204460492503131e-16" << endl;         
  cout << "--------------- CPP: sinusoidal surface check -------------" << endl;
  cout << "-----------------------------------------------------------" << endl 
       << endl;
  
  return;
}
