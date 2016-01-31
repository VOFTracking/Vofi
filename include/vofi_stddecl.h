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

typedef double vofi_real;
typedef const double vofi_creal;
typedef const int vofi_cint;
typedef int * const vofi_int_cpt;
typedef double (*integrand) (vofi_creal []);

/* xval: coordinates of the minimum or where the sign has changed, fval: local
   function value, sval: distance from the starting point, if applicable,
   iat: sign to have f>0, or if = 0 there is no minimum or no sign change */  
typedef struct {
  vofi_real xval[NDIM]; 
  vofi_real fval; 
  vofi_real sval; 
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

/* Fortran APIs */
vofi_real EXPORT(vofi_get_fh)(integrand,vofi_creal [],vofi_creal *,vofi_cint *,vofi_cint *); 
vofi_real EXPORT(vofi_get_cc)(integrand,vofi_creal [],vofi_creal *,vofi_creal *,vofi_cint *);

vofi_real vofi_get_segment_zero(integrand,vofi_creal [],vofi_creal [],vofi_creal [],vofi_creal,vofi_cint); 
int vofi_check_side_consistency(integrand,vofi_creal [],vofi_creal [],vofi_creal [],vofi_creal);
chk_data vofi_check_face_consistency(integrand,vofi_creal [],vofi_creal [],vofi_creal [],
                                vofi_creal [],vofi_creal);

dir_data vofi_get_dirs(integrand,vofi_creal [],vofi_real [],vofi_real [],vofi_real [],vofi_creal,vofi_creal,vofi_cint);
int vofi_get_limits(integrand,vofi_creal [],vofi_real [],vofi_creal [],vofi_creal [],vofi_creal [],vofi_creal,vofi_cint);

void vofi_get_side_intersections(integrand,vofi_real [],vofi_creal [],vofi_real [],vofi_creal [],
			    vofi_creal,vofi_int_cpt);
void vofi_get_face_intersections(integrand,min_data,vofi_creal [],vofi_real [],vofi_creal [],
			   vofi_creal [],vofi_creal,vofi_int_cpt);
min_data vofi_get_segment_min(integrand,vofi_creal [],vofi_creal [],vofi_creal [],vofi_creal,vofi_cint,vofi_cint);
min_data vofi_get_face_min(integrand,vofi_creal [],vofi_creal [],vofi_creal [],chk_data,vofi_creal);

vofi_real vofi_get_area(integrand,vofi_creal [],vofi_creal [],vofi_creal [],vofi_creal [],vofi_creal,vofi_cint,vofi_cint);
vofi_real vofi_get_volume(integrand,vofi_creal [],vofi_creal [],vofi_creal [],vofi_creal [],vofi_creal [],vofi_creal,
		  vofi_cint,vofi_cint);

#endif


