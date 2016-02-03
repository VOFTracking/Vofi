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
 * @file getcc.h
 * @authors Simone Bnà, Sandro Manservisi, Ruben Scardovelli, 
 *          Philip Yecko and Stephane Zaleski 
 * @date  12 November 2015
 * @brief Driver to compute the integration limits and the volume fraction 
 *        in two and three dimensions.
 */

#include "vofi_stddecl.h"

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * Driver to compute the volume fraction value in a given cell in two and     *
 * three dimensions                                                           *
 * INPUT:  pointer to the implicit function, starting point x0, grid          * 
 * spacing h0, characteristic function value fh, space dimension ndim0        *
 * OUTPUT: cc: volume fraction value                                          *
 * -------------------------------------------------------------------------- */

vofi_real vofi_Get_cc(integrand impl_func,vofi_creal x0[],vofi_creal h0,vofi_creal fh,vofi_cint ndim0) 
{
  int nsub;
  vofi_real pdir[NDIM],sdir[NDIM],tdir[NDIM],side[NSEG];
  vofi_real cc;
  dir_data icps; 

  icps = vofi_get_dirs(impl_func,x0,pdir,sdir,tdir,h0,fh,ndim0);
  if (icps.icc >= 0)
    cc = (vofi_real) icps.icc;
  else {
    nsub = vofi_get_limits(impl_func,x0,side,pdir,sdir,tdir,h0,ndim0);
    if (ndim0 == 2) 
      cc = vofi_get_area(impl_func,x0,side,pdir,sdir,h0,nsub,icps.ipt);
    else 
      cc =  vofi_get_volume(impl_func,x0,side,pdir,sdir,tdir,h0,nsub,icps.ipt);
  }
  
  return cc;
}
