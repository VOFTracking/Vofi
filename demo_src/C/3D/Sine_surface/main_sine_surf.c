#include <stdio.h>
#include <math.h>
#include <string.h>
#include "vofi.h"
#include "sine_surf.h"

#define NDIM  3
#define N2D   2
#define N3D   3

extern void check_volume(vofi_creal, vofi_cint);
extern vofi_real impl_func(vofi_creal []);
extern void init(vofi_cint);

/* -------------------------------------------------------------------------- *
 * PROGRAM TO INITIALIZE THE COLOR FUNCTION SCALAR FIELD                      *
 * -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  vofi_cint nx=NMX,ny=NMY,nz=NMZ,ndim0=N3D;
  int i,j,k,itrue;
  vofi_real cc[NMX][NMY][NMZ],x0[NDIM],xloc[NDIM];
  double h0,fh,vol_n;
  int randominput = 0;  
  int count;
  
  for (count = 1; count < argc; ++count)
  {
    if (!strcmp(argv[count], "-r") || !strcmp(argv[count], "--randominput")) 
    {
      randominput = 1;
    }
  }

/* -------------------------------------------------------------------------- *
 * initialization of the color function with local Gauss integration          * 
 * -------------------------------------------------------------------------- */
   
  h0 = H/nx;                                                  /* grid spacing */
  itrue = 1;

  /* starting point to get fh */
  x0[0] = 0.5;
  x0[1] = 0.5; 
  x0[2] = 0.5; 
  
  init(randominput);

  /* get the characteristic value fh of the implicit function */
  fh = vofi_Get_fh(impl_func,x0,h0,ndim0,itrue);
 
  /* put now starting point in (X0,Y0,Z0) to initialize the color function */
  x0[0] = X0; 
  x0[1] = Y0; 
  x0[2] = Z0; 

  /* xloc: minor vertex of each cell of the grid */
  for (i=0; i<nx; i++)
    for (j=0; j<ny; j++) 
      for (k=0; k<nz; k++) {
      xloc[0] = x0[0] + i*h0;
      xloc[1] = x0[1] + j*h0;
      xloc[2] = x0[2] + k*h0;
      
      cc[i][j][k] = vofi_Get_cc(impl_func,xloc,h0,fh,ndim0);
   }

  /* final global check */     
  vol_n = 0.0;

  for (i=0;i<nx; i++)
    for (j=0;j<ny; j++) 
      for (k=0;k<nz; k++) 
	 vol_n += cc[i][j][k];
  
  vol_n = vol_n*h0*h0*h0;

  check_volume(vol_n, randominput);

  return 0;
}
