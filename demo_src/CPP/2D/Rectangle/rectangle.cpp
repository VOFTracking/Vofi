#include <cstdio>
#include <iostream> 
#include <iomanip> 
#include <cstdlib>
#include <cmath>
#include "rectangle.h"

typedef const double creal;
typedef const int cint;
typedef double real;
using namespace std;

//* -------------------------------------------------------------------------- *
//* DESCRIPTION (reference phase where f(x,y) < 0):                            *
//* rectangle inside the square [0,1]x[0,1]                                    *
//* PARAMETERS:                                                                *
//* (XC,YC): center of the rectangle; ALPHA: angle between two axes x' and x;  *
//* (A1,B1): half sides along the local x' and y' axes                         *
//* -------------------------------------------------------------------------- *

real impl_func(creal xy[]) 
{
  real x,y,f0,f1,ca,sa;

  x = xy[0];
  y = xy[1];

  ca = cos(ALPHA);
  sa = sin(ALPHA);
  f1 = - A1 - ((x-XC)*ca+(y-YC)*sa);
  f0 = f1;
  /* f0 = MAX(f0,f1); */
  f1 = ((x-XC)*ca+(y-YC)*sa) - A1;
  f0 = MAX(f0,f1);
  f1 = -B1 - ((y-YC)*ca - (x-XC)*sa);
  f0 = MAX(f0,f1);
  f1 = ((y-YC)*ca - (x-XC)*sa) - B1;
  f0 = MAX(f0,f1);

  return f0;
}

//* -------------------------------------------------------------------------- *

void check_area(creal area_n)
{
  real area_a;
  
  area_a = 4.*A1*B1;

  cout << "-----------------------------------------------------------" << endl;
  cout << "------------------- CPP: rectangle check ------------------" << endl;
  cout << "analytical  area : " << scientific << setw(23) << setprecision(16) 
       << area_a << endl;
  cout << "numerical  area  : " << scientific << setw(23) << setprecision(16) 
       << area_n << endl << endl;
  cout << "absolute error   : " << scientific << setw(23) << setprecision(16) 
       << fabs(area_a-area_n) << endl;
  cout << "relative error   : " << scientific << setw(23) << setprecision(16) 
       << fabs(area_a-area_n)/area_a << endl;
  cout << "-----------------------------------------------------------" << endl;
  cout << "with Intel i7 3.4 GHz + Linux openSUSE 12.3 + gcc 4.7.2 -O3" << endl;
  cout << "-----------------------------------------------------------" << endl;
  cout << "analytical area :  2.3999999999999999e-01" << endl;
  cout << "numerical  area :  2.3999999999999996e-01" << endl << endl;
  cout << "absolute error  :  2.7755575615628914e-17" << endl;
  cout << "relative error  :  1.1564823173178715e-16" << endl;
  cout << "----------------- CPP: end rectangle check ----------------" << endl;
  cout << "-----------------------------------------------------------" << endl 
       << endl;
  
  return;
}
