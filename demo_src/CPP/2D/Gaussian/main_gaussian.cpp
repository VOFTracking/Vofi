#include <cstdio>
#include <cmath>
#include "vofi.h"
#include "gaussian.h"

#define NDIM  3
#define N2D   2

extern void check_area(creal);
extern real impl_func(creal []);

//* -------------------------------------------------------------------------- *
//* PROGRAM TO INITIALIZE THE COLOR FUNCTION SCALAR FIELD CC                   *
//* [X0,X0+H] X [Y0,Y0+H]                                                      *
//* -------------------------------------------------------------------------- *

int main()
{
  cint nx=NMX, ny=NMY, ndim0=N2D;
  int i,j,itrue;
  real cc[NMX][NMY],x0[NDIM],xloc[NDIM];
  double h0,fh,area_n;
  
//* -------------------------------------------------------------------------- *
//* initialization of the color function with local Gauss integration          * 
//* -------------------------------------------------------------------------- *
  
  h0 = H/nx;                                                  //* grid spacing *
  itrue = 1;

  //* put starting point in the center of the square *
  x0[0] = 0.5;
  x0[1] = 0.5; 
  x0[2] = 0.;
  
  //* get the characteristic value fh of the implicit function *
  fh = Get_fh(impl_func,x0,h0,ndim0,itrue);

  //* put now starting point in (X0,Y0) to initialize the color function *
  x0[0] = X0; 
  x0[1] = Y0; 

  xloc[2] = 0.;                //* xloc: minor vertex of each cell of the grid *
  for (i=0; i<nx; i++)
    for (j=0; j<ny; j++) {
      xloc[0] = x0[0] + i*h0;
      xloc[1] = x0[1] + j*h0;
      
      cc[i][j] = Get_cc(impl_func,xloc,h0,fh,ndim0);
   }

  //* final global check *
  area_n = 0.0;

  for (i=0;i<nx; i++)
    for (j=0;j<ny; j++) 
      area_n += cc[i][j];
  
  area_n = area_n*h0*h0;

  check_area(area_n);

  return 0;
}
