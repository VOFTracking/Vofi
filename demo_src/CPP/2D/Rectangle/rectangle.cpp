#include <cstdio>
#include <iostream> 
#include <iomanip> 
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "rectangle.h"

typedef const double creal;
typedef const int cint;
typedef double real;
using namespace std;

//* -------------------------------------------------------------------------- *
//* DESCRIPTION (reference phase where f(x,y) < 0):                            *
//* rectangle inside the square [0,1]x[0,1]                                    *
//* PARAMETERS:                                                                *
//* (xc,yc): center of the rectangle; alpha: angle between two axes x' and x;  *
//* (a1,b1): half sides along the local x' and y' axes                         *
//* -------------------------------------------------------------------------- *

static double a1 = 0.2;
static double b1 = 0.3;
static double alpha = 0.;
static double xc = 0.52;
static double yc = 0.44;

/* -------------------------------------------------------------------------- */

void init()
{

  /* initialize random number generator. */
  srand(time(0));  
  
  /* a1 --> from 0.17 to 0.23 */
  double axesscalingfactor = 0.06;
  a1 = 0.17 + ((double)rand() / RAND_MAX)*axesscalingfactor;
  /* a1 --> from 0.27 to 0.33 */
  b1 = 0.27 + ((double)rand() / RAND_MAX)*axesscalingfactor;
  
  /* alpha --> from zero to pi/2 */
  double alphascalingfactor = M_PI_2;
  alpha = ((double)rand() / RAND_MAX)*alphascalingfactor;
  
  /* xc & yc --> from 0.4 to 0.6 */
  double x0 = 0.45;
  double centerscalingfactor = 0.1;
  xc = x0 + ((double)rand() / RAND_MAX)*centerscalingfactor;
  yc = x0 + ((double)rand() / RAND_MAX)*centerscalingfactor;
  
  return;
}

/* -------------------------------------------------------------------------- */

real impl_func(creal xy[]) 
{
  real x,y,f0,f1,ca,sa;

  x = xy[0];
  y = xy[1];

  ca = cos(alpha);
  sa = sin(alpha);
  f1 = - a1 - ((x-xc)*ca+(y-yc)*sa);
  f0 = f1;
  /* f0 = MAX(f0,f1); */
  f1 = ((x-xc)*ca+(y-yc)*sa) - a1;
  f0 = MAX(f0,f1);
  f1 = -b1 - ((y-yc)*ca - (x-xc)*sa);
  f0 = MAX(f0,f1);
  f1 = ((y-yc)*ca - (x-xc)*sa) - b1;
  f0 = MAX(f0,f1);

  return f0;
}

//* -------------------------------------------------------------------------- *

void check_area(creal area_n)
{
  real area_a;
  
  area_a = 4.*a1*b1;

  cout << "-----------------------------------------------------------" << endl;
  cout << "------------------- CPP: rectangle check ------------------" << endl;
  cout << "analytical area : " << scientific << setw(23) << setprecision(16) 
       << area_a << endl;
  cout << "numerical  area : " << scientific << setw(23) << setprecision(16) 
       << area_n << endl << endl;
  cout << "absolute error  : " << scientific << setw(23) << setprecision(16) 
       << fabs(area_a-area_n) << endl;
  cout << "relative error  : " << scientific << setw(23) << setprecision(16) 
       << fabs(area_a-area_n)/area_a << endl;
  cout << "-----------------------------------------------------------" << endl;
//   cout << "with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2" << endl;
//   cout << "-----------------------------------------------------------" << endl;
//   cout << "analytical area :  2.3999999999999999e-01" << endl;
//   cout << "numerical  area :  2.3999999999999996e-01" << endl << endl;
//   cout << "absolute error  :  2.7755575615628914e-17" << endl;
//   cout << "relative error  :  1.1564823173178715e-16" << endl;
//   cout << "----------------- CPP: end rectangle check ----------------" << endl;
//   cout << "-----------------------------------------------------------" << endl; 
  cout << endl;
  
  return;
}
