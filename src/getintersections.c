#include "stddecl.h"

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * compute the interface intersections, if any, with a given cell side;       *
 * these are new internal/external limits of integration                      *
 * INPUT: pointer to the implicit function, function value at the endpoints   *
 * fe, starting point x0, direction stdir, grid spacing h0                    *
 * OUTPUT: nsub: updated number of subdivisions; array lim_intg: updated      *
 * start of new subdivisions                                                  *
 * -------------------------------------------------------------------------- */

void get_side_intersections(integrand impl_func,real fe[],creal x0[],real 
		       lim_intg[],creal stdir[],creal h0,int_cpt nsub)
{
  int f_iat; 
  cint true_sign=1,max_iter=50;   
  real dh0,fh0,ss;    
  min_data xfsa;
  
  if (fe[0]*fe[1] < 0.0) {
    dh0 = get_segment_zero(impl_func,fe,x0,stdir,h0,true_sign);
    if (fe[0] > 0.0)
      dh0 = h0 - dh0;
    lim_intg[*nsub] = dh0;
    (*nsub)++;
  }
  else {
    f_iat = check_side_consistency(impl_func,fe,x0,stdir,h0); 
    if (f_iat != 0) {
      xfsa = get_segment_min(impl_func,fe,x0,stdir,h0,f_iat,max_iter);
      if (xfsa.iat != 0) {
	fh0 = fe[1];
	fe[1] = xfsa.fval;
	dh0 = get_segment_zero(impl_func,fe,x0,stdir,xfsa.sval,true_sign);
	if (fe[0] > 0.0 || fe[1] < 0.0)
	  dh0 = xfsa.sval - dh0;
	lim_intg[*nsub] = dh0;
	(*nsub)++;
	ss = h0 - xfsa.sval;
	fe[0] = fe[1]; 
	fe[1] = fh0;
	dh0 = get_segment_zero(impl_func,fe,xfsa.xval,stdir,ss,true_sign);
	if (fe[0] > 0.0 || fe[1] < 0.0)
	  dh0 = ss - dh0;
	lim_intg[*nsub] = xfsa.sval + dh0;
	(*nsub)++;
      }
    }
  }
    
  return;
}    

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * get the external limits of integration that are inside the face            *
 * INPUT: pointer to the implicit function, structure with point position     *
 * with negative f value and function sign attribute, starting point x0,      * 
 * secondary and tertiary directions sdir and tdir, grid spacing h0           *
 * OUTPUT: nsub: updated number of subdivisions; array lim_intg: updated      *
 * start of new subdivisions                                                  *
 * -------------------------------------------------------------------------- */

void get_face_intersections(integrand impl_func,min_data xfsa,creal x0[],real 
	      lim_intg[],creal sdir[],creal tdir[],creal h0,int_cpt nsub)
{
  int i,k,iter,js,jt,not_conv,ipt,ist,f_iat;
  cint max_iter = 50;
  real pt0[NDIM],pt1[NDIM],pt2[NDIM],ptt[NDIM],mp0[NDIM],mp1[NDIM];
  real ss[NDIM],exdir[NDIM],indir[NDIM],fe[NEND];
  real ss0,ds0,fpt0,sss,sst,ssx,ssy,tol2,normdir,d1,d2,a1,a2;
  creal tol = EPS_M; 
 
  /* GRAPHICS I */
  tol2 = 2.*tol;
  for (i=0;i<NDIM;i++)
    pt0[i] = xfsa.xval[i];
  fe[0] = xfsa.fval;
  f_iat = xfsa.iat;

  /* initialize internal direction and get js and jt indices */
  for (i=0;i<NDIM;i++) {
    indir[i] = sdir[i];
    if (sdir[i] > 0.5)
      js = i;
    if (tdir[i] > 0.5)
      jt = i;
  }

  for (i=0;i<NDIM;i++)
    pt2[i] = pt1[i] = pt0[i];

  /* get zero or boundary point pt2 along secondary direction with ss -> h0 */
  ss0 = x0[js] + h0 - pt0[js];         
  pt2[js] = x0[js] + h0;            
  fe[1] = f_iat*impl_func(pt2);
  if (fe[1] > 0.) {
    ds0 = get_segment_zero(impl_func,fe,pt0,indir,ss0,f_iat);
    pt2[js] = pt0[js] + ds0;
  } 
  /* DEBUG 1 */
  /* DEBUG 1
  fprintf(stderr,"ini. pt. right : %17.9e %17.9e %17.9e \n", pt2[0],pt2[1],pt2[2]);
 */
  fprintf(stderr,"sec/right : %17.9e %17.9e %17.9e \n", pt2[0],pt2[1],pt2[2]);


  /* get zero or boundary point pt1 along secondary direction with ss -> 0 */
  ss0 = pt0[js] - x0[js];    
  pt1[js] = x0[js];                  
  indir[js] = -1.;
  fe[1] = f_iat*impl_func(pt1);
  if (fe[1] > 0.) {
    ds0 = get_segment_zero(impl_func,fe,pt0,indir,ss0,f_iat);
    pt1[js] = pt0[js] - ds0;
  }
  /* DEBUG 2 */
  /* DEBUG 2
  fprintf(stderr,"ini. pt. left  : %17.9e %17.9e %17.9e \n", pt1[0],pt1[1],pt1[2]);
 */
  fprintf(stderr,"sec/left  : %17.9e %17.9e %17.9e \n", pt1[0],pt1[1],pt1[2]);

  for (i=0;i<NDIM;i++)
    pt0[i] = 0.5*(pt1[i] + pt2[i]);                  /* starting midpoint pt0 */
  fpt0 = f_iat*impl_func(pt0);
  ss0 = pt2[js]-pt1[js];
  /* DEBUG 3 */
  /* DEBUG 3
  fprintf(stderr,"ini. pt. centre: %17.9e %17.9e %17.9e \n", pt0[0],pt0[1],pt0[2]);
 */
  fprintf(stderr,"sec/centre: %17.9e %17.9e %17.9e \n", pt0[0],pt0[1],pt0[2]);

  /* now get the two external limits along the tertiary direction */
  for (k=-1;k<=1;k=k+2) {
    /* DEBUG 4 */
    /* DEBUG 4
    fprintf(stderr,"-------------------------------------------------\n");
    fprintf(stderr,"down/up (-1/1): %2d \n",k);
 */
    fprintf(stderr,"-------------------------------------------------\n");
    fprintf(stderr,"TOP/BOTTOM: %2d \n",k);

    iter = 0;
    not_conv = 1;
    for (i=0;i<NDIM;i++)                     /* initialize external direction */
      exdir[i] = tdir[i];
    exdir[jt] = k;
    if (k < 0)
      sst = pt0[jt] - x0[jt];   
    else
      sst = x0[jt] + h0 - pt0[jt];   
    for (i=0;i<NDIM;i++) {
      mp1[i] = pt0[i];
      pt1[i] = mp1[i] + sst*exdir[i];
    }
    fe[0] = fpt0;
    fe[1] = f_iat*impl_func(pt1);
    sss = ss0;
    while (not_conv  && iter < max_iter) {    /* iterative loop for the limit */
      /* DEBUG 5 */
      /* DEBUG 5
      fprintf(stderr,"iter: %2d, sss: %e sst: %e \n",iter+1,sss,sst);
 */
      fprintf(stderr,"iter: %2d, sss: %e sst: %e \n",iter+1,sss,sst);
      if (fe[1] > 0.) {
	ds0 = get_segment_zero(impl_func,fe,mp1,exdir,sst,f_iat);
	sst = ds0;
      }
      /* DEBUG 6 */
      /* DEBUG 6
      fprintf(stderr,"sst (real): %e \n",sst);
 */
      fprintf(stderr,"sst (real): %e \n",sst);
      for (i=0;i<NDIM;i++) {
	pt1[i] = mp1[i] + sst*exdir[i];         /* zero along the secant line */
	mp0[i] = mp1[i];
	mp1[i] = pt2[i] = ptt[i] = pt1[i];
      }
      /* DEBUG 7 */
      /* DEBUG 7
      fprintf(stderr,"tp/bot: %17.9e %17.9e %17.9e \n", pt1[0],pt1[1],pt1[2]);
 */
      fprintf(stderr,"tp/bot: %17.9e %17.9e %17.9e \n", pt1[0],pt1[1],pt1[2]);

      /* try to get other zero along the secondary direction */ 
      ipt = ist = 0;             
      ptt[js] += tol;            
      fe[0] = f_iat*impl_func(ptt);
      if (fe[0] < 0.) {
	ipt = 1;
	ssx = x0[js] + h0 - ptt[js];
	indir[js] = 1.;
      }
      else {
        ptt[js] -= tol2;
	fe[0] = f_iat*impl_func(ptt);
	if (fe[0] < 0.) {
	  ipt = 1;
	  ssx = ptt[js] - x0[js];
	  indir[js] = -1.;
	}
      }
      if (ipt) {             
	sss = MIN(1.2*sss,ssx); /* get the segment length along secondary dir */
	for (i=0;i<NDIM;i++) 
	  pt2[i] = ptt[i] + sss*indir[i];
	fe[1] = f_iat*impl_func(pt2);
	while (fe[1] < 0. && ist < 3 && sss < ssx) {
	  sss = MIN(3.*sss,ssx);
	  if (ist == 2)
	    sss = ssx;
	  for (i=0;i<NDIM;i++) 
	    pt2[i] = ptt[i] + sss*indir[i];
	  fe[1] = f_iat*impl_func(pt2);
	  ist++;
	}
        if (fe[0]*fe[1] < 0.) {         /* get other zero along secondary dir */
	  ds0 = get_segment_zero(impl_func,fe,ptt,indir,sss,f_iat);
	  for (i=0;i<NDIM;i++) 
	    pt2[i] = ptt[i] + ds0*indir[i];
	}
	/* DEBUG 8 */
	/* DEBUG 8
	fprintf(stderr,"zero 2: %17.9e %17.9e %17.9e \n", pt2[0],pt2[1],pt2[2]);
 */
	fprintf(stderr,"zero 2: %17.9e %17.9e %17.9e \n", pt2[0],pt2[1],pt2[2]);
	for (i=0;i<NDIM;i++)                        /* get midpoint and width */ 
	  mp1[i] = 0.5*(pt1[i] + pt2[i]);
	fe[0] = f_iat*impl_func(mp1);
	sss = fabs(pt1[js]-pt2[js]);
      }
      /* DEBUG 9 */
      /* DEBUG 9
      fprintf(stderr,"midpt : %17.9e %17.9e %17.9e \n", mp1[0],mp1[1],mp1[2]);
 */
      fprintf(stderr,"midpt : %17.9e %17.9e %17.9e \n", mp1[0],mp1[1],mp1[2]);

      normdir = 0.;                             
      for (i=0;i<NDIM;i++) {
	exdir[i] =  mp1[i] - mp0[i];                      /* secant direction */
        normdir += exdir[i]*exdir[i];
      }  
      normdir = sqrt(normdir) + EPS_NOT0; 
      for (i=0;i<NDIM;i++) {              
	exdir[i] = exdir[i]/normdir;                 /* unit secant direction */
	d1 = SGN0P(exdir[i]);
	d2 = fabs(exdir[i]) + EPS_NOT0;
	a1 = (x0[i] - mp1[i])/(d1*d2);
	a2 = (x0[i] + h0 - mp1[i])/(d1*d2);
	ss[i] = MAX(a1,a2);
      }
      ssy = MIN(ss[0],ss[1]);      
      ssy = MIN(ssy,ss[2]);             
      sst = MIN(1.2*sst,ssy);
      
      if (!ipt || sss < tol2 || sst < EPS_R) {       /* convergence criterion */ 
        not_conv = 0;                       
	/* DEBUG 10 */
	/* DEBUG 10
        fprintf(stderr,"Converged! k: %2d, ipt:%2d sss,sst: %e %e \n",k,ipt,sss,sst);
 */
        fprintf(stderr,"Converged! k: %2d, ipt:%2d sss,sst: %e %e \n",
		k,ipt,sss,sst);
      }
      else {
        for (i=0;i<NDIM;i++) 
          pt1[i] = mp1[i] + sst*exdir[i];
        fe[1] = f_iat*impl_func(pt1);
	ist = 0;                   /* get the segment length along secant dir */
        while (fe[1] < 0. && ist < 3 && sst < ssy) { 
          sst = MIN(3.*sst,ssy);
	  if (ist == 2)
	    sst = ssy;
          for (i=0;i<NDIM;i++) 
            pt1[i] = mp1[i] + sst*exdir[i];
          fe[1] = f_iat*impl_func(pt1);
	  ist++;
        }
      }  
      iter++;
    }
    lim_intg[*nsub] = mp1[jt] - x0[jt];
    (*nsub)++;
  }

  return;
}

/* -------------------------------------------------------------------------- */
  /* GRAPHICS I 
  double xx[3],yy[3],zz[3];
  xx[0] = yy[0] = zz[0] = 0.;
  if (sdir[0] < 0.5 && tdir[0] < 0.5) {
    xx[1] = xx[2] = x0[0];
  }
  else if (sdir[1] < 0.5 && tdir[1] < 0.5) {
    yy[1] = yy[2] = x0[1];
  }
  else {
    zz[1] = zz[2] = x0[2];
  }
  
  if (sdir[0] > 0.5) {
    xx[1] = x0[0]; xx[2] = xx[1] + h0;
  }
  else if (sdir[1] > 0.5) {
    yy[1] = x0[1]; yy[2] = yy[1] + h0;
  }
  else {
    zz[1] = x0[2]; zz[2] = zz[1] + h0;
  }
   END GRAPHICS I */     


  /* GRAPHICS II
  if (extdir[0] > 0.5) {
    xx[1] = xx[2] = x[0] + side[nseg-1];
  }
  else if (extdir[1] > 0.5) {
    yy[1] = yy[2] = x[1] + side[nseg-1];
  }
  else {
     zz[1] = zz[2] = x[2] + side[nseg-1];
  }   
  line_3D(xx,yy,zz,2,1,1,4,1,0,0);

  if (extdir[0] > 0.5) {
    xx[1] = xx[2] = x[0] + side[nseg-2];
  }
  else if (extdir[1] > 0.5) {
    yy[1] = yy[2] = x[1] + side[nseg-2];
  }
  else {
     zz[1] = zz[2] = x[2] + side[nseg-2];
  }       
  line_3D(xx,yy,zz,2,1,1,4,1,0,0);
 */
  /* END GRAPHICS II */


/*
	  fprintf(stderr,"iteration for sss: fs1,fs0: %e %e sss, ssx %e %e \n",
		  fs[1],fs[0],ssx,sss);
          fprintf(stderr,"iteration for sst: fs1,fs0: %e %e sst, ss1 %e %e \n",
                  fs[1],fs[0],sst,ss1);
    ixs = MAX(0,k);
    xgfs[ixs] = mp1[jt];
 */ 
