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
//* f(x,y,z) = z - a0 - b0*sin(c1*pi*x + pi*d1)*sin(c1*pi*x + pi*e1)           *
//* -------------------------------------------------------------------------- *

static double a0 = 0.5;
static double b0 = 1./6.;
static double c1 = 1.6;
static double d1 = 1./7.;
static double e1 = 1./5.;


/* -------------------------------------------------------------------------- */

void init()
{

  /* initialize random number generator. */
  srand(time(0));  

  double scalingfactor = 0.05;
 
  /* b0 --> from 1.15 to 1.20 */
  b0 = 1.15 + ((double)rand() / RAND_MAX)*scalingfactor;
  
  /* c1 --> from 1.6 to 1.65 */
  c1 = 1.6 + ((double)rand() / RAND_MAX)*scalingfactor;
   
  /* d1 --> from 0.12 to 0.17 */
  d1 = 0.12 + ((double)rand() / RAND_MAX)*scalingfactor;
 
  /* e1 --> from 0.15 to 0.25 */
  e1 = 0.15 + ((double)rand() / RAND_MAX)*scalingfactor;
  
  return;
}


real impl_func(creal xy[])
{
  double x,y,z,f0;

  x = xy[0];
  y = xy[1];
  z = xy[2];

  f0 = z - a0 - b0*sin(M_PI*(c1*x+d1))*sin(M_PI*(c1*y+e1));

  return f0;
}

//* -------------------------------------------------------------------------- *

void check_volume(creal vol_n)
{
  double vol_a;

  //vol_a = 0.5;
  vol_a = a0 + (b0/(c1*M_PI*c1*M_PI))*(cos(d1*M_PI) - cos((d1+c1)*M_PI))*(cos(e1*M_PI) - cos((e1+c1)*M_PI)); 

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
//   cout << "with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2" << endl;
//   cout << "-----------------------------------------------------------" << endl;
//   cout << "analytical volume:  5.0000000000000000e-01" << endl;	     
//   cout << "numerical  volume:  5.0000000000000011e-01" << endl << endl; 
//   cout << "absolute error   :  1.1102230246251565e-16" << endl;	     
//   cout << "relative error   :  2.2204460492503131e-16" << endl;         
//   cout << "------------- CPP: end sinusoidal surface check -----------" << endl;
//   cout << "-----------------------------------------------------------" << endl 
  cout << endl;
  
  return;
}
