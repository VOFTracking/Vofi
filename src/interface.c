#include "vofi_stddecl.h"
#include "vofi.h"

/* ------------------------------------------------------------------- *
 * DESCRIPTION:                                                        *
 * FORTRAN to C interface for the function vofi_Get_fh                      *
 * INPUT and OUTPUT: see vofi_Get_fh                                        *
 * ------------------------------------------------------------------- */

real EXPORT(vofi_get_fh)(integrand impl_func,creal x0[],creal *H0,cint *Ndim0,cint *iX0)
{
  creal h0 = *H0;
  cint ndim0 = *Ndim0, ix0 = *iX0;
  real Fh;

  Fh = vofi_Get_fh(impl_func,x0,h0,ndim0,ix0); 

  return Fh;
}

/* ------------------------------------------------------------------- *
 * DESCRIPTION:                                                        *
 * FORTRAN to C interface for the function vofi_Get_cc                      *
 * INPUT and OUTPUT: see vofi_Get_cc                                        *
 * ------------------------------------------------------------------- */

real EXPORT(vofi_get_cc)(integrand impl_func,creal x0[],creal *H0,creal *Fh,cint *Ndim0)
{
  creal h0= *H0, fh = *Fh;
  cint ndim0 = *Ndim0;
  real CC;
  
  CC = vofi_Get_cc(impl_func,x0,h0,fh,ndim0);

  return CC;
}
