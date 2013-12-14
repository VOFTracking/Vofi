#include "stddecl.h"

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * subdivide the side along the secondary/tertiary (2/3) direction to define  *
 * rectangles/rectangular hexahedra with or without the interface             *
 * INPUT: pointer to the implicit function, starting point x0, primary,       *
 * secondary, tertiary directions pdir, sdir, tdir, grid spacing h0,          *
 * subdivision direction stdir (2/3)                                          *
 * OUTPUT: nsub: total number of subdivisions; array lim_intg: start/end of   *
 * each subdivision (lim_intg[0] = 0, lim_intg[nsub] = h0)                    *
 * -------------------------------------------------------------------------- */

int get_limits(integrand impl_func,creal x0[],real lim_intg[],creal pdir[],
	       creal sdir[],creal tdir[],creal h0,cint stdir)
{
  int i,j,k,iv,nsub,nvp,nvn;
  real fv[NVER],x1[NDIM],x2[NDIM],fe[NEND],ds,ls;
  chk_data fvga; 
  min_data xfsa;
  
  lim_intg[0] = 0.;
  nsub = 1;  
  if (stdir == 2) {                     /* get the internal limits along sdir */
    for (j=0;j<2;j++) {                                          /* two sides */
      for (i=0;i<NDIM;i++) {
	x1[i] = x0[i] + j*pdir[i]*h0;
	x2[i] = x1[i] + sdir[i]*h0;
      }
      fe[0] = impl_func(x1);
      fe[1] = impl_func(x2);
      get_side_intersections(impl_func,fe,x1,lim_intg,sdir,h0,&nsub);
    }
  }
  else {                                /* get the external limits along tdir */
    for (k=0;k<2;k++) {                                         /* two planes */
      /* DEBUG 1 */
      /* DEBUG 1 */
      fprintf(stderr,"\nPlane: %2d \n",k+1);

      nvp = nvn = iv = 0;
      for (j=0;j<2;j++) {                                        /* two sides */
	/* DEBUG 2 */
	/* DEBUG 2 */
	fprintf(stderr,"Check side: %2d \n",j+1);

	for (i=0;i<NDIM;i++) { 
	  x1[i] = x0[i] + k*pdir[i]*h0+j*sdir[i]*h0;
	  x2[i] = x1[i] + tdir[i]*h0;
	}
	fe[0] = impl_func(x1);
	fv[iv++] = fe[0];
	fe[1] = impl_func(x2);
	fv[iv++] = fe[1];
	if (fe[0]*fe[1] >= 0.0) {
	  if ((fe[0]+fe[1]) > 0.)
	    nvp += 2;
	  else
	    nvn += 2;
	}
	get_side_intersections(impl_func,fe,x1,lim_intg,tdir,h0,&nsub);
	/* DEBUG 3 */
	/* DEBUG 3 */
	fprintf(stderr,"nsub: %2d \n\n",nsub);

      }
      if (nvp == 4 || nvn == 4) {         /* get the extra limits in the face */
	/* DEBUG 4 */
	/* DEBUG 4 */
	fprintf(stderr,"\nCheck whole face: %2d,  nvp,nvn: %2d %2d \n",k+1,nvp,nvn);

	xfsa.iat = 0;
	for (i=0;i<NDIM;i++) 
	  x1[i] = x0[i] + k*pdir[i]*h0;
	fvga = check_face_consistency(impl_func,fv,x1,sdir,tdir,h0); 
	if (fvga. iat != 0)
	  xfsa = get_face_min(impl_func,x1,sdir,tdir,fvga,h0);
	if (xfsa.iat != 0)
	  get_face_intersections(impl_func,xfsa,x1,lim_intg,sdir,tdir,h0,&nsub);
	/* DEBUG 5 */
	/* DEBUG 5 */
	fprintf(stderr,"nsub: %2d \n\n",nsub);

      } 
    }
  }
  lim_intg[nsub] = h0;

  /* DEBUG 6 */
  /* DEBUG 6 */
  fprintf(stderr,"\n before ordering and removal \n");
  for (i=0;i<=nsub;i++)
    fprintf(stderr,"ns:%2d dx,dx/h: %.12f %.12f \n",i,lim_intg[i],lim_intg[i]/h0);

  for (j=2;j<nsub;j++) {                         /* order limits from 0 to h0 */
    ls = lim_intg[j];
    i = j-1;
    while (i > 0 && lim_intg[i] > ls) {
      lim_intg[i+1] = lim_intg[i];
      i--;
    }
    lim_intg[i+1] = ls;
  }

  i = 0;
  while (i<nsub) {                             /*remove zero-length intervals */
    ds = lim_intg[i+1] - lim_intg[i];
    if (ds < EPS_R) {
      for (j=i;j<nsub;j++)
	lim_intg[j] = lim_intg[j+1];
      nsub--;
    }
    i++;
  }
  lim_intg[0] = 0.;                                        /* just for safety */
  lim_intg[nsub] = h0;
  /* DEBUG 7 */
  /* DEBUG 7 */
  fprintf(stderr,"\n final \n");
  for (i=0;i<=nsub;i++)
    fprintf(stderr,"ns:%2d dx,dx/h: %.12f %.12f \n",i,lim_intg[i],lim_intg[i]/h0);


  return nsub;    
}
  /* DEBUG 7 
  fprintf(stderr,"\n final \n");
  for (i=0;i<=nsub;i++)
    fprintf(stderr,"ns:%2d dx,dx/h: %f %f \n",i,lim_intg[i],lim_intg[i]/h0);
*/
