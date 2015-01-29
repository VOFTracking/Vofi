#include <cstdio>
#include <iostream> 
#include <iomanip> 
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "ellipse.h"

typedef const double creal;
typedef const int cint;
typedef double real;
using namespace std;

//* -------------------------------------------------------------------------- *
//* DESCRIPTION (reference phase where f(x,y) < 0):                            *
//* ellipse inside the square [0,1]x[0,1]                                      *
//* f(x,y) = c1*x^2 + c2*y^2 + c3*x*y + c4*x + c5*y - c6                       *
//* PARAMETERS:                                                                *
//* (xc,yc): center of the ellipse; ALPHA: angle between two axes x' and x;    *
//* (a1,b1): semiaxis along the two ellipse axes (local) x' and y'             *
//* -------------------------------------------------------------------------- *

static double a1 = 0.17;
static double b1 = 0.21;
static double alpha = 0.48;
static double xc = 0.523;
static double yc = 0.475;

/* -------------------------------------------------------------------------- */

void init(cint randominput)
{
  
  if(randominput) {
    /* initialize random number generator. */
    srand(time(0));  
  
    /* a1 --> from 0.15 to 0.19 */
    double axesscalingfactor = 0.04;
    a1 = 0.15 + ((double)rand() / RAND_MAX)*axesscalingfactor;
    /* a1 --> from 0.19 to 0.23 */
    b1 = 0.19 + ((double)rand() / RAND_MAX)*axesscalingfactor;
  
    /* alpha --> from zero to pi/2 */
    double alphascalingfactor = M_PI_2;
    alpha = ((double)rand() / RAND_MAX)*alphascalingfactor;
  
    /* xc & yc --> from 0.4 to 0.6 */
    double x0 = 0.45;
    double centerscalingfactor = 0.1;
    xc = x0 + ((double)rand() / RAND_MAX)*centerscalingfactor;
    yc = x0 + ((double)rand() / RAND_MAX)*centerscalingfactor;
  
  }
  
  return;
}

/* -------------------------------------------------------------------------- */

real impl_func(creal xy[]) 
{
  real x,y,a2,b2,ca,sa,c1,c2,c3,c4,c5,c6,f0;

  x = xy[0];
  y = xy[1];

  a2 = a1*a1;
  b2 = b1*b1;
  ca = cos(alpha);
  sa = sin(alpha);
  c1 = ca*ca/a2 + sa*sa/b2;
  c2 = sa*sa/a2 + ca*ca/b2;
  c3 = 2.*ca*sa*(b2-a2)/(a2*b2);
  c4 = -(2.*c1*xc + c3*yc);
  c5 = -(2.*c2*yc + c3*xc);
  c6 = 1.0 - (c1*xc*xc + c2*yc*yc + c3*xc*yc);
  
  f0 = c1*x*x + c2*y*y + c3*x*y + c4*x + c5*y - c6;

  return f0;
}

//* -------------------------------------------------------------------------- *

void check_area(creal area_n, cint randominput)
{
  real area_a;
  
  area_a = M_PI*a1*b1;

  cout << "-----------------------------------------------------------" << endl;
  cout << "-------------------- CPP: ellipse check -------------------" << endl;
  cout << "analytical area : " << scientific << setw(23) << setprecision(16) 
       << area_a << endl;
  cout << "numerical  area : " << scientific << setw(23) << setprecision(16) 
       << area_n << endl << endl;
  cout << "absolute error  : " << scientific << setw(23) << setprecision(16) 
       << fabs(area_a-area_n) << endl;
  cout << "relative error  : " << scientific << setw(23) << setprecision(16) 
       << fabs(area_a-area_n)/area_a << endl;
  if(!randominput) {
    cout << "-----------------------------------------------------------" << endl;
    cout << "with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2" << endl;
    cout << "-----------------------------------------------------------" << endl;
    cout << "analytical area :  1.1215485773315563e-01" << endl;
    cout << "numerical  area :  1.1215485773315677e-01" << endl << endl;
    cout << "absolute error  :  1.1379786002407855e-15" << endl;
    cout << "relative error  :  1.0146494081855289e-14" << endl;
    cout << "------------------ CPP: end ellipse check -----------------" << endl;
    cout << "-----------------------------------------------------------" << endl; 
  }
  cout << endl;
  
  return;
}
