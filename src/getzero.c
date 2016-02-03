/****************************************************************************
 * Copyright (C) 2015 by Simone Bnà(a), Sandro Manservisi(a),               *
 * Ruben Scardovelli(a), Philip Yecko(b) and Stephane Zaleski(c,d)          *
 * (a) DIN–Lab. di Montecuccolino, Università di Bologna,                   *
 *     Via dei Colli 16, 40136 Bologna, Italy                               *
 * (b) Physics Department, Cooper Union, New York, NY, USA                  *
 * (c) Sorbonne Universités, UPMC Univ Paris 06, UMR 7190,                  *
 *     Institut Jean Le Rond d’Alembert, F-75005, Paris, France             *
 * (d) CNRS, UMR 7190, Institut Jean Le Rond d’Alembert, F-75005,           *
 *     Paris, France                                                        *
 *                                                                          *
 * You should have received a copy of the CPC license along with Vofi.      *
 * If not, see http://cpc.cs.qub.ac.uk/licence/licence.html.                *
 *                                                                          *
 * e-mail: ruben.scardovelli@unibo.it                                       *
 *                                                                          *
 ****************************************************************************/

/**
 * @file getzero.c
 * @authors Simone Bnà, Sandro Manservisi, Ruben Scardovelli, 
 *          Philip Yecko and Stephane Zaleski 
 * @date  12 November 2015
 * @brief it computes the zero in a given segment.
 */


#include "vofi_stddecl.h"

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * compute the zero in a given segment of length s0, the zero is strictly     *
 * bounded, i.e. f(0)*f(s0) < 0                                               *
 * METHOD: standard hybrid method with a combination of secant and bisection  *
 * INPUT: pointer to the implicit function, function value at the             *
 * endpoints fe, starting point x0, direction dir, segment length s0, sign    *
 * attribute f_sign                                                           *
 * OUTPUT: sz: length of the segment where f is negative                      *
 * -------------------------------------------------------------------------- */

vofi_real vofi_get_segment_zero(integrand impl_func,vofi_creal fe[],vofi_creal x0[],
                           vofi_creal dir[],vofi_creal s0,vofi_cint f_sign)
{
  int not_conv,iss,i,iter;
  vofi_cint max_iter=25;
  vofi_real xs[NDIM],sl,sr,ss,sold,fl,fr,fs,fold,dss,dsold,dfs,sz;


  if (fe[0] > 0.0) {                                    /* sl where f(sl) < 0 */
    sl = s0;
    sr = 0.0;
    fl = fe[1];
    fr = fe[0];
    iss = 1;
  }
  else {
    sl = 0.0;
    sr = s0;
    fl = fe[0];
    fr = fe[1];
    iss = 0;
  }

  if (fabs(fl) <= fabs(fr)) {               /* ss where |fs| = MIN(|fl|,|fr|) */
    ss = sl;
    fs = fl;
  }
  else {
    ss = sr;
    fs = fr;
  }
  xs[2] = 0.;
  dsold = dss = s0;
  not_conv = 1;
  dfs = (fr-fl)/(sr-sl);
  iter = 0;
  
  while (not_conv  && iter < max_iter) {                    /* iterative loop */

    if ( ((ss-sr)*dfs-fs)*((ss-sl)*dfs-fs) > 0.0 ||
	 fabs(2.*fs) > fabs(dsold*dfs) ) {                  /* bisection step */
      dsold = dss;
      sold = ss;
      dss = 0.5*(sr - sl);
      ss = sl + dss;
    }
    else {                                                     /* secant step */
      dsold = dss;
      sold = ss;
      dss = fs/dfs;
      ss = ss - dss;
    }
    iter++;
    if (fabs(dss) < EPS_R)                           /* convergence criterion */
      not_conv = 0;

    if (not_conv) {                       /* new fs and dfs, bracket the zero */
      for (i=0; i<NDIM; i++)
        xs[i] = x0[i] + ss*dir[i];
      fold = fs;
      fs = f_sign*impl_func(xs);
      dfs = (fs-fold)/(ss-sold);
      if (fs < 0.0)
	sl = ss;
      else
	sr = ss;
    }
  }

  if (!not_conv)                                /* segment length where f < 0 */
    sz = (1-iss)*ss + iss*(s0-ss);
  else {                                               /* too many iterations */
    fprintf(stderr,"Root finding: too many iterations! \n");
    sz = -1.;
  }

  return sz;
}
