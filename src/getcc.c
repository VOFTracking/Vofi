#include "stddecl.h"

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * Driver to compute the volume fraction value in a given cell in two and     *
 * three dimensions                                                           *
 * INPUT:  pointer to the implicit function, starting point x0, grid          * 
 * spacing h0, characteristic function value fh, space dimension ndim0        *
 * OUTPUT: cc: volume fraction value                                          *
 * -------------------------------------------------------------------------- */

real Get_cc(integrand impl_func,creal x0[],creal h0,creal fh,cint ndim0) 
{
  int nsub;
  real pdir[NDIM],sdir[NDIM],tdir[NDIM],side[NSEG];
  real cc;
  dir_data icps; 

  icps = get_dirs(impl_func,x0,pdir,sdir,tdir,h0,fh,ndim0);
  if (icps.icc >= 0)
    cc = (real) icps.icc;
  else {
    nsub = get_limits(impl_func,x0,side,pdir,sdir,tdir,h0,ndim0);
    if (ndim0 == 2) 
      cc = get_area(impl_func,x0,side,pdir,sdir,h0,nsub,icps.ipt);
    else 
      cc =  get_volume(impl_func,x0,side,pdir,sdir,tdir,h0,nsub,icps.ipt);
  }
  
  return cc;
}
