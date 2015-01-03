/*
        vofi.h
                Prototypes for the vofi library
*/

#ifndef VOFI_H
#define VOFI_H

typedef const double creal;
typedef double real;
typedef const int cint;
typedef double (*integrand) (creal []);

#ifdef __cplusplus
extern "C" {
#endif
  
/** C/C++ APIs */
real vofi_get_fh(integrand,creal [],creal,cint,cint); 
real vofi_get_cc(integrand,creal [],creal,creal,cint);

#ifdef __cplusplus
}
#endif

#endif