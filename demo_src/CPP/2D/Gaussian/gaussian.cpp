#include <cstdio>
#include <iostream> 
#include <iomanip> 
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "gaussian.h"

typedef const double creal;
typedef const int cint;
typedef double real;
using namespace std;

//* -------------------------------------------------------------------------- *
//* DESCRIPTION (reference phase where f(x,y) < 0):                            *
//* gaussian line in the square [0,1]x[0,1]                                    *
//* f(x,y) = y - yy0 - a0 exp[-ga (x - xx0)^2]                                 *
//* -------------------------------------------------------------------------- *

static double yy0 = 0.22;
static double a0 = 0.51;
static double xx0 = 0.541;
static double ga = 60.3;


/* -------------------------------------------------------------------------- */

void init(cint randominput)
{

  if(randominput) {
  
    /* initialize random number generator. */
    srand(time(0));  
   
    double scalingfactor = 0.06;
    /* yy0 --> from 0.19 to 0.25 */
    yy0 = 0.19 + ((double)rand() / RAND_MAX)*scalingfactor;
  
    /* xx0 --> from 0.511 to 0.561 */
    xx0 = 0.511 + ((double)rand() / RAND_MAX)*scalingfactor;
  
    /* a0 --> from 0.48 to 0.54 */
    a0 = 0.48 + ((double)rand() / RAND_MAX)*scalingfactor;
  
    /* ga --> from 60.00 to 0.60.6 */
    ga = 60. + ((double)rand() / RAND_MAX)*scalingfactor;
  
  }
  return;
}

/* -------------------------------------------------------------------------- */

real impl_func(creal xy[]) 
{
  double x,y,f0;

  x = xy[0];
  y = xy[1];

  f0 = y - yy0 - a0*exp(-ga*(x-xx0)*(x-xx0));

  return f0;
}

//* -------------------------------------------------------------------------- *

void check_area(creal area_n, cint randominput)
{
  real area_a;

  /* analytical integration with x in [0,1]  */
  area_a = yy0 + 0.5*a0*sqrt(M_PI/ga)*(erf(sqrt(ga)*(1.-xx0) )-erf(-sqrt(ga)*xx0));
  
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
  if(!randominput) {
    cout << "with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2" << endl;
    cout << "-----------------------------------------------------------" << endl;
    cout << "analytical area :  3.3640894546075428e-01" << endl;
    cout << "numerical  area :  3.3640894546075717e-01" << endl << endl;
    cout << "absolute error  :  2.8865798640254070e-15" << endl;
    cout << "relative error  :  8.5805680942041349e-15" << endl;
    cout << "----------------- CPP: end gaussian check -----------------" << endl;
    cout << "-----------------------------------------------------------" << endl;
  }
  cout << endl;

  return;
}
