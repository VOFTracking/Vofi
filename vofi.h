/*
        vofi.h
                Prototypes for the vofi library
*/

#ifdef __cplusplus
extern "C" {
#endif

typedef const double creal;
typedef double real;
typedef const int cint;
typedef double (*integrand) (creal []);

double Get_fh(integrand,creal [],creal,cint,cint); 
double Get_cc(integrand,creal [],creal,creal,cint); 

#ifdef __cplusplus
}
#endif
