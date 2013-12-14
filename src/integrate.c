#include "stddecl.h"
#include "GL.h"

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * compute the normalized cut area with a Gauss-Legendre quadrature           *
 * INPUT: pointer to the implicit function, starting point x0, internal       *
 * limits of integration int_lim_intg, primary and secondary directions pdir  *
 * and  sdir, grid spacing h0, number of internal subdivisions nintsub,       *
 * tentative number of internal integration points nintpt                     *
 * OUTPUT: area: normalized value of the cut area or 2D volume fraction       *
 * -------------------------------------------------------------------------- */

double get_area(integrand impl_func,creal x0[],creal int_lim_intg[],creal 
		pdir[],creal sdir[],creal h0,cint nintsub,cint nintpt)
{
  int i,ns,k,npt,cut_rect;
  cint true_sign = 1;
  real x1[NDIM],x20[NDIM],x21[NDIM],fe[NEND];
  real area,ds,cs,xis,ht,GL_1D;
  creal *ptinw, *ptinx;

  /* GRAPHICS I */
  area = 0.;
  for (i=0;i<NDIM;i++) 
    x1[i] = x0[i] + pdir[i]*h0;

  /* DEBUG 1 */
  /* DEBUG 1 */
  fprintf(stderr,"\nIntegration over %2d rectangles \n",nintsub);

  for (ns=1;ns<=nintsub;ns++) {                   /* loop over the rectangles */
    ds = int_lim_intg[ns] - int_lim_intg[ns-1];
    cs = 0.5*(int_lim_intg[ns] + int_lim_intg[ns-1]);
    cut_rect = 0;
    for (i=0;i<NDIM;i++) {
      x20[i] = x0[i] + sdir[i]*cs;
      x21[i] = x1[i] + sdir[i]*cs;
    }    
    fe[0] = impl_func(x20);
    fe[1] = impl_func(x21);
    if (fe[0]*fe[1] <= 0.)
      cut_rect = 1;        
    
    if (!cut_rect) {                     /* no interface: full/empty rectangle */
      if (fe[0] < 0.0)
	area += ds*h0; 
      /* DEBUG 2 */
    }
    else {                   /* cut rectangle: internal numerical integration */
      if (ds < 0.1*h0) 
	npt = 4;
      else if (ds < 0.2*h0)
	npt = 8;
      else if (ds < 0.4*h0)
	npt = MIN(nintpt,12);
      else if (ds < 0.6*h0) 
	npt = MIN(nintpt,16);
      else
	npt = MIN(nintpt,20);

      switch (npt) {
      case 4:
	ptinx = csi04;
        ptinw = wgt04;
        break;
      case 8:
	ptinx = csi08;
        ptinw = wgt08;
        break;
      case 12:
	ptinx = csi12;
	ptinw = wgt12;
        break;
      case 16:
	ptinx = csi16;
	ptinw = wgt16;
	break;
      default:
	ptinx = csi20;
	ptinw = wgt20;
	break;
      }

      GL_1D = 0.;
      /* DEBUG 3 */
      /* DEBUG 3 */
      fprintf(stderr,"band: %2d; internal GL integration: %2d pts \n",ns,npt);

      for (k=0;k<npt;k++) {
	xis = cs + 0.5*ds*(*ptinx);
	for (i=0;i<NDIM;i++) {
	  x20[i] = x0[i] + sdir[i]*xis;
	  x21[i] = x1[i] + sdir[i]*xis;
	}
	fe[0] = impl_func(x20);
	fe[1] = impl_func(x21);
	ht = get_segment_zero(impl_func,fe,x20,pdir,h0,true_sign);
	/* DEBUG 4 */
	/* DEBUG 4 */
 	fprintf(stderr,"k:%2d %17.10f %17.10f %17.10f \n",k,*ptinx,xis/h0,ht/h0);

	GL_1D += (*ptinw)*ht;
	ptinx++;
	ptinw++;
	/* GRAPHICS II */
      }
      /* GRAPHICS III */
      area += 0.5*ds*GL_1D;
    }
  }

  area = area/(h0*h0);                               /* normalized area value */
  return area;
}                       

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * compute the normalized cut volume with a double Gauss-Legendre quadrature  *
 * INPUT: pointer to the implicit function, starting point x0, external       *
 * limits of integration ext_lim_intg, primary, secondary and tertiary        *
 * directions pdir, sdir and tdir, grid spacing h0, number of external        *
 * subdivisions nextsub, tentative number of internal integration points      *
 * nintpt                                                                     *  
 * OUTPUT: vol: normalized value of the cut volume or 3D volume fraction      *
 * -------------------------------------------------------------------------- */

double get_volume(integrand impl_func,creal x0[],creal ext_lim_intg[],creal 
		  pdir[],creal sdir[],creal tdir[],creal h0,cint nextsub,
		  cint nintpt)
{
  int i,ns,k,nexpt,cut_hexa,f_iat,nintsub;
  cint stdir=2,max_iter=50;
  real x1[NDIM],x2[NDIM],x3[NDIM],fe[NEND],int_lim_intg[NSEG];
  real vol,ds,cs,xis,f1,f2,area_n,GL_1D;
  creal *ptexw, *ptexx;
  min_data xfsa;
  
  vol = 0.;
  
  /* DEBUG 1 */
  /* DEBUG 1 */
  fprintf(stderr,"\nIntegration over %2d hexahedra \n",nextsub);

  for (ns=1;ns<=nextsub;ns++) {        /* loop over the rectangular hexahedra */
    ds = ext_lim_intg[ns] - ext_lim_intg[ns-1];        
    cs = 0.5*(ext_lim_intg[ns] + ext_lim_intg[ns-1]);
    cut_hexa = 0;
    for (i=0;i<NDIM;i++) { 
      x1[i] = x0[i] + tdir[i]*cs;
      x2[i] = x1[i] + pdir[i]*h0;
    }
    f1 = impl_func(x1);
    f2 = impl_func(x2);
    if (f1*f2 <= 0.)
      cut_hexa = 1;                        
    if (!cut_hexa) {            /* check lower side along secondary direction */
      fe[0] = f1;
      for (i=0;i<NDIM;i++)  
	x3[i] = x1[i] + sdir[i]*h0;
      fe[1] = impl_func(x3);
      if (fe[0]*fe[1] <= 0.)
	cut_hexa = 1;        
      else {
	f_iat = check_side_consistency(impl_func,fe,x1,sdir,h0); 
	if (f_iat != 0) { 
	  xfsa = get_segment_min(impl_func,fe,x1,sdir,h0,f_iat,max_iter);
	  cut_hexa = xfsa.iat;        
	}
      }
    }
    if (!cut_hexa) {            /* check upper side along secondary direction */
      fe[0] = f2;
      for (i=0;i<NDIM;i++)  
	x3[i] = x2[i] + sdir[i]*h0;
      fe[1] = impl_func(x3);
      if (fe[0]*fe[1] <= 0.)
	cut_hexa = 1;        
      else {
	f_iat = check_side_consistency(impl_func,fe,x2,sdir,h0); 
	if (f_iat != 0) {
	  xfsa = get_segment_min(impl_func,fe,x2,sdir,h0,f_iat,max_iter);
	  cut_hexa = xfsa.iat;        
	}
      }
    }

    if (!cut_hexa) {                   /* no interface: full/empty hexahedron */ 
      if (f1 < 0.)
	vol += ds;
      /* DEBUG 2 */
      /* DEBUG 2 */
      fprintf(stderr,"band: %2d, empty/full: %e \n",ns,f1);

    }
    else {                  /* cut hexahedron: external numerical integration */
      if (ds < 0.1*h0) {
	nexpt = 8;
	ptexx = csi08;
        ptexw = wgt08;
      }
      else if (ds < 0.3*h0) {
	nexpt = 12;
	ptexx = csi12;
        ptexw = wgt12;
      }
      else if (ds < 0.5*h0) {
	nexpt = 16;
	ptexx = csi16;
        ptexw = wgt16;
      }
      else {
	nexpt = 20;
	ptexx = csi20;
        ptexw = wgt20;
      }
      GL_1D = 0.;
      /* DEBUG 3 */
      /* DEBUG 3 */
      fprintf(stderr,"band: %2d, external GL integration, %2d pts \n",ns,nexpt);

      for (k=0;k<nexpt;k++) {
	xis = cs + 0.5*ds*(*ptexx);
	for (i=0;i<NDIM;i++) 
	  x1[i] = x0[i] + tdir[i]*xis;
	nintsub = get_limits(impl_func,x1,int_lim_intg,pdir,sdir,tdir,h0,stdir);
	area_n = get_area(impl_func,x1,int_lim_intg,pdir,sdir,h0,nintsub,nintpt);
	/* DEBUG 4 */
	/* DEBUG 4 */
        fprintf(stderr,"k: %2d glx,y,area: %22.15e %.15e %.15e \n",k,*ptexx,xis,area_n);

 	GL_1D += (*ptexw)*area_n;
	ptexx++;
	ptexw++;
      }
      vol += 0.5*ds*GL_1D;
    }
  }

  vol = vol/h0;                                    /* normalized volume value */
  return vol;
}                       

/* -------------------------------------------------------------------------- *
 * Graphical and debug sections for get_area                                  *
 * -------------------------------------------------------------------------- */
  /* GRAPHICS I
  double xx[3],yy[3],zz[3],x3[21],y3[21],z3[21]; 
  xx[0] = yy[0] = zz[0] = 0.;
  x3[0] = y3[0] = z3[0] = 0.;
  */

      /* DEBUG 2 */
      /* plot full/empty circles */
      /*
      if (fabs(ds-h0) < EPS_M && ndim == 2) {      
	xx[1] = x0[0] + 0.5*h0;
	yy[1] = x0[1] + 0.5*h0;
	if(fs[0] < 0.0)
	  line_2D(xx,yy,1,0,0,0,5,13,4);
	else
	  line_2D(xx,yy,1,0,0,0,5,12,4);
      }  */                                                       /* end plot */
	/* GRAPHICS II */
        /* plot segment
	if (fs[0] < 0.0) {                                    
	  xx[1] = x20[0]; xx[2] = xx[1] + pdir[0]*ht;
	  yy[1] = x20[1]; yy[2] = yy[1] + pdir[1]*ht;
	  zz[1] = x20[2]; zz[2] = zz[1] + pdir[2]*ht;
	}
	else {
	  xx[1] = x21[0]; xx[2] = xx[1] - pdir[0]*ht;
	  yy[1] = x21[1]; yy[2] = yy[1] - pdir[1]*ht;
	  zz[1] = x21[2]; zz[2] = zz[1] - pdir[2]*ht;
	}
	if (ndim == 3)
	  x3[k+1] = xx[2]; y3[k+1] = yy[2]; z3[k+1] = zz[2]; 
 */

      /* GRAPHICS III */
      /* plot line in 3D     
      if (ndim == 3)
	line_3D(x3,y3,z3,npt,1,1,5,1,13,0);
*/

/* -------------------------------------------------------------------------- *
 * Graphical and debug sections for get_volume                                *
 * -------------------------------------------------------------------------- */

