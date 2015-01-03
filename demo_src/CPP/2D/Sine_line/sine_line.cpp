#include <cstdio>
#include <iostream> 
#include <iomanip> 
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "sine_line.h"

typedef const double creal;
typedef const int cint;
typedef double real;
using namespace std;

//* -------------------------------------------------------------------------- *
//* DESCRIPTION (reference phase where f(x,y) < 0):                            *
//* sinusoidal line in the square [0,1]x[0,1]                                  *
//* f(x,y) = y - b0*sin(c0 pi x+ pi/d0) - a0                                   *
//* -------------------------------------------------------------------------- *

static double a0 = 0.5;
static double b0 = 0.25;
static double c0 = 4.0;
static double d0 = 14.0;


/* -------------------------------------------------------------------------- */

void init()
{

  /* initialize random number generator. */
  srand(time(0));  
  
  /* a0 --> from 0.45 to 0.55 */
  double scalingfactor = 0.1;
  a0 = 0.45 + ((double)rand() / RAND_MAX)*scalingfactor;
  /* b0 --> from 0.20 to 0.30 */
  b0 = 0.20 + ((double)rand() / RAND_MAX)*scalingfactor;
  /* c0 --> from 0.35 to 0.45 */
  c0 = 0.35 + ((double)rand() / RAND_MAX)*scalingfactor;  
  /* d0 --> from 13.95 to 14.05 */
  d0 = 13.95 + ((double)rand() / RAND_MAX)*scalingfactor;
  
  return;
}

/* -------------------------------------------------------------------------- */

real impl_func(creal xy[]) 
{
  double x,y,f0;

  x = xy[0];
  y = xy[1];

  f0 = y - a0 - b0*sin(c0*M_PI*x + M_PI/d0);

  return f0;
}

//* -------------------------------------------------------------------------- *

void check_area(creal area_n)
{
  real area_a;

  //* analytical integration with x in [0,1]  *
  area_a = a0 + b0*(-cos((c0 + 1./d0)*M_PI) + cos(M_PI/d0))/(c0*M_PI);
  
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
//   cout << "with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2" << endl;
//   cout << "-----------------------------------------------------------" << endl;
//   cout << "analytical area :  5.0000000000000000e-01" << endl;
//   cout << "numerical  area :  4.9999999999993749e-01" << endl << endl;
//   cout << "absolute error  :  6.2505556286396313e-14" << endl;
//   cout << "relative error  :  1.2501111257279263e-13" << endl;
//   cout << "----------------- CPP: end sine line check ----------------" << endl;
//   cout << "-----------------------------------------------------------" << endl; 
  cout << endl;

  return;
}
