#ifndef VOFI_STDDECL_H
#define VOFI_STDDECL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
#define Extern extern "C"
#else
#define Extern extern
#endif

#define PREFIX(s) s

#if NOUNDERSCORE
#define SUFFIX(s) s
#else
#define SUFFIX(s) s##_
#endif

#define EXPORT(s) EXPORT_(PREFIX(s))
#define EXPORT_(s) SUFFIX(s)

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SGN0P(a) ((a<0) ? -1 : 1)
#define Sq(a) ((a)*(a))
#define Sq3(a) (a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
#define Sqd3(a,b) ((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))
#define SHFT4(a,b,c,d)  (a)=(b); (b)=(c); (c)=(d)
#define CPSF(s,t,f,g) (s)=(t); (f)=(g)

#define EPS_M    1.5e-07
#define EPS_LOC  1.5e-07
#define EPS_E    5.0e-07
#define EPS_R    1.0e-14
#define EPS_NOT0 1.0e-50
#define NDIM     3
#define NVER     4
#define NEND     2
#define NLSX     3
#define NLSY     3
#define NLSZ     3
#define NSEG    10

typedef double real;
typedef const double creal;
typedef const int cint;
typedef int * const int_cpt;
typedef double (*integrand) (creal []);

/* xval: coordinates of the minimum or where the sign has changed, fval: local
   function value, sval: distance from the starting point, if applicable,
   iat: sign to have f>0, or if = 0 there is no minimum or no sign change */  
typedef struct {
  real xval[NDIM]; 
  real fval; 
  real sval; 
  int iat; 
} min_data;

/* ivs,ivt: indices to locate the vertex in the face, igs,igt: on/off indices
   for the gradient components, iat: sign to have f>0, or if = 0 there is no 
   minimum or no sign change */ 
typedef struct {
  int ivs; int ivt;               
  int igs; int igt; 
  int iat;    
} chk_data;

/* icc: full/empty/cut cell (1/0/-1); ipt: tentative number of integration 
   points; isb: number of subdivisions, not yet implemented (hence: 0/1) */
typedef struct {
  int icc; int ipt; int isb;
} dir_data;

/* function prototypes */
/* */

/** Fortran APIs */
real EXPORT(vofi_get_fh)(integrand,creal [],creal *,cint *,cint *); 
real EXPORT(vofi_get_cc)(integrand,creal [],creal *,creal *,cint *);

real get_segment_zero(integrand,creal [],creal [],creal [],creal,cint); 
int check_side_consistency(integrand,creal [],creal [],creal [],creal);
chk_data check_face_consistency(integrand,creal [],creal [],creal [],
                                creal [],creal);

dir_data get_dirs(integrand,creal [],real [],real [],real [],creal,creal,cint);
int get_limits(integrand,creal [],real [],creal [],creal [],creal [],creal,cint);

void get_side_intersections(integrand,real [],creal [],real [],creal [],
			    creal,int_cpt);
void get_face_intersections(integrand,min_data,creal [],real [],creal [],
			    creal [],creal,int_cpt);
min_data get_segment_min(integrand,creal [],creal [],creal [],creal,cint,cint);
min_data get_face_min(integrand,creal [],creal [],creal [],chk_data,creal);

real get_area(integrand,creal [],creal [],creal [],creal [],creal,cint,cint);
real get_volume(integrand,creal [],creal [],creal [],creal [],creal [],creal,
		  cint,cint);

#endif


