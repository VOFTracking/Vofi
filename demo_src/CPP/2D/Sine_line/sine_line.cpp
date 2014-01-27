#include <cstdio>
#include <iostream> 
#include <iomanip> 
#include <cstdlib>
#include <cmath>
#include "sine_line.h"

typedef const double creal;
typedef const int cint;
typedef double real;
using namespace std;

//* -------------------------------------------------------------------------- *
//* DESCRIPTION (reference phase where f(x,y) < 0):                            *
//* sinusoidal line in the square [0,1]x[0,1]                                  *
//* f(x,y) = y - B0*sin(C0 pi x+ pi/D0) - A0                                   *
//* -------------------------------------------------------------------------- *

real impl_func(creal xy[]) 
{
  double x,y,f0;

  x = xy[0];
  y = xy[1];

  f0 = y - A0 - B0*sin(C0*MYPI*x + MYPI/D0);

  return f0;
}

//* -------------------------------------------------------------------------- *

void check_area(creal area_n)
{
  real area_a;

  //* analytical integration with x in [0,1]  *
  area_a = A0 + B0*(-cos((C0 + 1./D0)*MYPI) + cos(MYPI/D0))/(C0*MYPI);
  
  cout << "-----------------------------------------------------------" << endl;
  cout << "------------------- CPP: sine line check ------------------" << endl;
  cout << "analytical area : " << scientific << setw(23) << setprecision(16) 
       << area_a << endl;
  cout << "numerical  area : " << scientific << setw(23) << setprecision(16) 
       << area_n << endl << endl;
  cout << "absolute error  : " << scientific << setw(23) << setprecision(16) 
       << fabs(area_a-area_n) << endl;
  cout << "relative error  : " << scientific << setw(23) << setprecision(16) 
       << fabs(area_a-area_n)/area_a << endl;
  cout << "-----------------------------------------------------------" << endl;
  cout << "with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2" << endl;
  cout << "-----------------------------------------------------------" << endl;
  cout << "analytical area :  5.0000000000000000e-01" << endl;
  cout << "numerical  area :  4.9999999999993749e-01" << endl << endl;
  cout << "absolute error  :  6.2505556286396313e-14" << endl;
  cout << "relative error  :  1.2501111257279263e-13" << endl;
  cout << "----------------- CPP: end sine line check ----------------" << endl;
  cout << "-----------------------------------------------------------" << endl 
       << endl;

  return;
}
