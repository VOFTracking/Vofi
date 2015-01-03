#include <stdio.h>
#include <math.h>
#include "vofi.h"
#include "gaussian.h"

#define NDIM  3
#define N2D   2

extern void check_area(creal);
extern real impl_func(creal []);
extern void init();

/* -------------------------------------------------------------------------- *
 * PROGRAM TO INITIALIZE THE COLOR FUNCTION SCALAR FIELD CC                   *
 * -------------------------------------------------------------------------- */

int main()
{
  cint nx=NMX, ny=NMY, ndim0=N2D;
  int i,j,itrue;
  real cc[NMX][NMY],x0[NDIM],xloc[NDIM];
  double h0,fh,area_n;
  
/* -------------------------------------------------------------------------- *
 * initialization of the color function with local Gauss integration          * 
 * -------------------------------------------------------------------------- */
  
  /* grid spacing */
  h0 = H/nx;                                                  
  itrue = 1;

  /* put starting point in the center of the square */
  x0[0] = 0.5;
  x0[1] = 0.5; 
  x0[2] = 0.;
  
  init();
  
  /* get the characteristic value fh of the implicit function */
  fh = vofi_get_fh(impl_func,x0,h0,ndim0,itrue);

  /* put now starting point in (X0,Y0) to initialize the color function */
  x0[0] = X0; 
  x0[1] = Y0; 

  /* xloc: minor vertex of each cell of the grid */
  xloc[2] = 0.;                
  for (i=0; i<nx; i++)
    for (j=0; j<ny; j++) {
      xloc[0] = x0[0] + i*h0;
      xloc[1] = x0[1] + j*h0;

      cc[i][j] = vofi_get_cc(impl_func,xloc,h0,fh,ndim0);
   }

  /* final global check */     
  area_n = 0.0;

  for (i=0;i<nx; i++)
    for (j=0;j<ny; j++) 
      area_n += cc[i][j];
  
  area_n = area_n*h0*h0;

  check_area(area_n);

  return 0;
}
