#include <cstdio>
#include <iostream> 
#include <iomanip> 
#include <cstdlib>
#include <cmath>
#include "gaussian.h"

typedef const double creal;
typedef const int cint;
typedef double real;
using namespace std;

//* -------------------------------------------------------------------------- *
//* DESCRIPTION (reference phase where f(x,y) < 0):                            *
//* gaussian line in the square [0,1]x[0,1]                                    *
//* f(x,y) = y - YY0 - A0 exp[-GA (x - XX0)^2]                                 *
//* -------------------------------------------------------------------------- *

real impl_func(creal xy[]) 
{
  double x,y,f0;

  x = xy[0];
  y = xy[1];

  f0 = y - YY0 - A0*exp(-GA*(x-XX0)*(x-XX0));

  return f0;
}

//* -------------------------------------------------------------------------- *

void check_area(creal area_n)
{
  real area_a;

  //* integration with MATHEMATICA with the given values of Y0,A0,GA,X0 *  
  area_a = 0.3364089454607542483401167;
  
  cout << "-----------------------------------------------------------" << endl;
  cout << "------------------- CPP: gaussian check -------------------" << endl;
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
  cout << "analytical area :  3.3640894546075423e-01" << endl;
  cout << "numerical  area :  3.3640894546075717e-01" << endl << endl;
  cout << "absolute error  :  2.9420910152566648e-15" << endl;
  cout << "relative error  :  8.7455790190926772e-15" << endl;
  cout << "----------------- CPP: end gaussian check -----------------" << endl;
  cout << "-----------------------------------------------------------" << endl 
       << endl;

  return;
}
