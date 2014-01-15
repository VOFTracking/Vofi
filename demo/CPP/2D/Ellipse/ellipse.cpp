#include <cstdio>
#include <iostream> 
#include <iomanip> 
#include <cstdlib>
#include <cmath>
#include "ellipse.h"
typedef const double creal;
typedef const int cint;
typedef double real;
using namespace std;

/* -------------------------------------------------------------------------- *
 * DESCRIPTION (reference phase where f(x,y) < 0):                            *
 * ellipse inside the square [0,1]x[0,1]                                      *
 * f(x,y) = c1*x^2 + c2*y^2 + c3*x*y + c4*x + c5*y - c6                       *
 * PARAMETERS:                                                                *
 * (XC,YC): center of the ellipse; ALPHA: angle between two axes x' and x;    *
 * (A1,B1): semiaxis along the two ellipse axes (local) x' and y'             *
 * -------------------------------------------------------------------------- */

int cont_line(real *xx,real *yy,cint ntot)
{
  int i;
  double dphi,phi,r0,ca,sa,xi,yi,a2,b2;

  ca = cos(ALPHA);
  sa = sin(ALPHA);
  a2 = A1*A1;
  b2 = B1*B1;

  dphi = 2.0*MYPI/(double)(ntot-1);

  for (i=1; i<ntot; i++) {
    phi = (i-1)*dphi;
    r0 = b2*a2/(b2*cos(phi)*cos(phi) + a2*sin(phi)*sin(phi));
    r0 = sqrt(r0);
    xi = r0*cos(phi);
    yi = r0*sin(phi);
    xx[i] = XC + xi*ca - yi*sa;
    yy[i] = YC + xi*sa + yi*ca;
  }

  xx[ntot] = xx[1];
  yy[ntot] = yy[1];

  return ntot;
}

/* -------------------------------------------------------------------------- */
void check_area(creal area_n)
{
  real area_a;
  
  area_a = MYPI*A1*B1;

  cerr << "-------------------------------------------------------" << endl;
  cerr << "------------------- C: ellipse check ------------------\n" << endl;
  cerr << "analytical  area : " << scientific << setw(23) << setprecision(16) << area_a << endl;
  cerr << "numerical  area  : " << scientific << setw(23) << setprecision(16) << area_n << endl << endl;
  cerr << "absolute error   : " << scientific << setw(23) << setprecision(16) << fabs(area_a-area_n) << endl;
  cerr << "relative error   : " << scientific << setw(23) << setprecision(16) << fabs(area_a-area_n)/area_a << endl;
  cerr << "-------------------------------------------------------" << endl;
  cerr << "with Intel i7 3.4 GHz + Linux openSUSE 12.3 + gcc 4.7.2" << endl;
  cerr << "-------------------------------------------------------" << endl;
  cerr << "analytical area :  1.1215485773315563e-01" << endl;
  cerr << "numerical  area :  1.1215485773315680e-01" << endl << endl;
  cerr << "absolute error  :  1.1657341758564144e-15" << endl;
  cerr << "relative error  :  1.0393969547266394e-14" << endl;
  cerr << "----------------- C: end ellipse check ----------------" << endl;
  cerr << "-------------------------------------------------------" << endl;
  
  return;
}

/* -------------------------------------------------------------------------- */
real impl_func(creal xy[]) 
{
  real x,y,a2,b2,ca,sa,c1,c2,c3,c4,c5,c6,f0;

  x = xy[0];
  y = xy[1];

  a2 = A1*A1;
  b2 = B1*B1;
  ca = cos(ALPHA);
  sa = sin(ALPHA);
  c1 = ca*ca/a2 + sa*sa/b2;
  c2 = sa*sa/a2 + ca*ca/b2;
  c3 = 2.*ca*sa*(b2-a2)/(a2*b2);
  c4 = -(2.*c1*XC + c3*YC);
  c5 = -(2.*c2*YC + c3*XC);
  c6 = 1.0 - (c1*XC*XC + c2*YC*YC + c3*XC*YC);
  
  f0 = c1*x*x + c2*y*y + c3*x*y + c4*x + c5*y - c6;

  return f0;
}
