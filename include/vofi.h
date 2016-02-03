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
 * This file is part of Vofi.                                               *
 * This is supplement to the papers:                                        *
 * [1] S Bnà, S Manservisi, R Scardovelli, P Yecko, S Zaleski,              *
 *     "Numerical integration of implicit functions for the initialization  *
 *     of the VOF function", Computers & Fluids 113, 42-52,                 *
 *     doi:10.1016/j.compfluid.2014.04.010                                  *
 * [2] S Bnà, S Manservisi, R Scardovelli, P Yecko, S Zaleski,              *
 *     "VOFI -- A library to initialize the volume fraction scalar field",  *
 *     Computer Physics Communications, Computer Physics Communications,    *
 *     2015, doi:10.1016/j.cpc.2015.10.026                                  *
 *                                                                          *
 * You should have received a copy of the CPC license along with Vofi.      *
 * If not, see http://cpc.cs.qub.ac.uk/licence/licence.html.                *
 *                                                                          *
 * e-mail: ruben.scardovelli@unibo.it                                       *
 *                                                                          *
 ****************************************************************************/

/**
 * @file vofi.h
 * @authors Simone Bnà, Sandro Manservisi, Ruben Scardovelli, 
 *          Philip Yecko and Stephane Zaleski 
 * @date  12 November 2015
 * @brief Prototypes for the vofi library.
 *
 * For any further information on the algorithm implemented in the code
 * @see http://www.sciencedirect.com/science/article/pii/S0010465515004087
 * @see http://www.sciencedirect.com/science/article/pii/S0045793014001480
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

/**
 * @brief Starting from point x0 get a zero of the implicit function given by the
 * user, using gradient ascent/descent, then compute its absolute value at a
 * distance hb along the normal direction.
 * @param impl_func pointer to the implicit function
 * @param x0 starting point
 * @param h0 grid spacing
 * @param ndim0 space dimension
 * @param ix0 switch for @p x0 (ix0=1: point x0 is given; ix0=0: use the default value for x0)     
 * @param fh "characteristic" function value
 * @note C/C++ API
 */
vofi_real vofi_Get_fh(integrand,vofi_creal [],vofi_creal,vofi_cint,vofi_cint); 

/**
 * @brief Driver to compute the volume fraction value in a given cell in two 
 * and three dimensions.
 * @param impl_func pointer to the implicit function
 * @param x0 starting point
 * @param h0 grid spacing
 * @param fh characteristic function value
 * @param ndim0 space dimension
 * @param cc volume fraction value
 * @note C/C++ API
 */
vofi_real vofi_Get_cc(integrand,vofi_creal [],vofi_creal,vofi_creal,vofi_cint);

#ifdef __cplusplus
}
#endif

#endif