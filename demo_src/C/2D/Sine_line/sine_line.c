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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "sine_line.h"

typedef const double creal;
typedef const int cint;
typedef double real;

/* -------------------------------------------------------------------------- *
 * DESCRIPTION (reference phase where f(x,y) < 0):                            *
 * sinusoidal line in the square [0,1]x[0,1]                                  *
 * f(x,y) = y - b0*sin(c0 pi x+ pi/d0) - a0                                   *
 * -------------------------------------------------------------------------- */

static double a0 = 0.5;
static double b0 = 0.25;
static double c0 = 4.0;
static double d0 = 14.0;


/* -------------------------------------------------------------------------- */

void init(cint randominput)
{
  
  if(randominput) {

    /* initialize random number generator. */
    srand(time(0));  
  
    /* a0 --> from 0.45 to 0.55 */
    double scalingfactor = 0.1;
    a0 = 0.45 + ((double)rand() / RAND_MAX)*scalingfactor;
    /* b0 --> from 0.20 to 0.30 */
    b0 = 0.20 + ((double)rand() / RAND_MAX)*scalingfactor;
    /* c0 --> from 0.35 to 0.45 */
    c0 = 0.35 + ((double)rand() / RAND_MAX)*scalingfactor;  
    /* d0 --> from 13.95 to 14.05 */
    d0 = 13.95 + ((double)rand() / RAND_MAX)*scalingfactor;
  
  }
  
  return;
}

/* -------------------------------------------------------------------------- */

real impl_func(creal xy[]) 
{
  double x,y,f0;

  x = xy[0];
  y = xy[1];

  f0 = y - a0 - b0*sin(c0*M_PI*x + M_PI/d0);

  return f0;
}

/* -------------------------------------------------------------------------- */

void check_area(creal area_n, cint randominput)
{
  real area_a;

  /* analytical integration with x in [0,1]  */
  area_a = a0 + b0*(-cos((c0 + 1./d0)*M_PI) + cos(M_PI/d0))/(c0*M_PI);

  fprintf (stdout,"------------------------------------------------------------------------\n");
  fprintf (stdout,"--------------------- C: sine line check -------------------------------\n");
  fprintf (stdout," * sinusoidal line in the square [%.1f,%.1f]x[%.1f,%.1f] in a %dX%d grid   *\n", X0, X0+H, Y0, Y0+H, NMX, NMY);
  fprintf (stdout," * f(x,y) = y - b0*sin(c0 pi x+ pi/d0) - a0                            *\n");
  fprintf (stdout,"------------------------------------------------------------------------\n");
  fprintf (stdout,"a0:    %23.16e\n",a0);
  fprintf (stdout,"b0:    %23.16e\n",b0);
  fprintf (stdout,"c0:    %23.16e\n",c0);
  fprintf (stdout,"d0:    %23.16e\n",d0);
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"analytical area : %23.16e\n",area_a);
  fprintf (stdout,"numerical  area : %23.16e\n\n",area_n);
  fprintf (stdout,"absolute error  : %23.16e\n",fabs(area_a-area_n));
  fprintf (stdout,"relative error  : %23.16e\n",fabs(area_a-area_n)/area_a); 
  fprintf (stdout,"-----------------------------------------------------------\n");
  if(!randominput) {
    fprintf (stdout,"with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2\n");
    fprintf (stdout,"-----------------------------------------------------------\n");
    fprintf (stdout,"analytical area :  5.0000000000000000e-01\n");
    fprintf (stdout,"numerical  area :  4.9999999999993749e-01\n\n");
    fprintf (stdout,"absolute error  :  6.2505556286396313e-14\n");
    fprintf (stdout,"relative error  :  1.2501111257279263e-13\n");
    fprintf (stdout,"----------------- C: end sine line check ------------------\n");
    fprintf (stdout,"-----------------------------------------------------------\n");
  }
  fprintf (stdout,"\n");

  return;
}
