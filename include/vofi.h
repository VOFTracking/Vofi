/*
        vofi.h
                Prototypes for the vofi library
*/

#ifndef VOFI_H
#define VOFI_H

typedef const double  vofi_creal;
typedef double  vofi_real;
typedef const int  vofi_cint;
typedef double (*integrand) ( vofi_creal []);

#ifdef __cplusplus
extern "C" {
#endif
  
/** C/C++ APIs */
vofi_real vofi_Get_fh(integrand,vofi_creal [],vofi_creal,vofi_cint,vofi_cint); 
vofi_real vofi_Get_cc(integrand,vofi_creal [],vofi_creal,vofi_creal,vofi_cint);

#ifdef __cplusplus
}
#endif

#endif