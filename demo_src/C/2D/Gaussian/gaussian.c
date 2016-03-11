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
#include "gaussian.h"

typedef const double creal;
typedef const int cint;
typedef double real;

/* -------------------------------------------------------------------------- *
 * DESCRIPTION (reference phase where f(x,y) < 0):                            *
 * gaussian line in the square [0,1]x[0,1]                                    *
 * f(x,y) = y - yy0 - a0 exp[-ga (x - xx0)^2]                                 *
 * -------------------------------------------------------------------------- */

static double yy0 = 0.22;
static double a0 = 0.51;
static double xx0 = 0.541;
static double ga = 60.3;


/* -------------------------------------------------------------------------- */

void init(cint randominput)
{

  if(randominput) {
    
    /* initialize random number generator. */
    srand(time(0));  
   
    double scalingfactor = 0.06;
    /* yy0 --> from 0.19 to 0.25 */
    yy0 = 0.19 + ((double)rand() / RAND_MAX)*scalingfactor;
  
    /* xx0 --> from 0.511 to 0.561 */
    xx0 = 0.511 + ((double)rand() / RAND_MAX)*scalingfactor;
  
    /* a0 --> from 0.48 to 0.54 */
    a0 = 0.48 + ((double)rand() / RAND_MAX)*scalingfactor;
  
    /* ga --> from 60.00 to 0.60.6 */
    ga = 60. + ((double)rand() / RAND_MAX)*scalingfactor;
  
  }
  
  return;
}

/* -------------------------------------------------------------------------- */

real impl_func(creal xy[]) 
{
  double x,y,f0;

  x = xy[0];
  y = xy[1];

  f0 = y - yy0 - a0*exp(-ga*(x-xx0)*(x-xx0));

  return f0;
}

/* -------------------------------------------------------------------------- */

void check_area(creal area_n, cint randominput)
{
  real area_a;

  /* analytical integration with x in [0,1]  */
  area_a = yy0 + 0.5*a0*sqrt(M_PI/ga)*(erf(sqrt(ga)*(1.-xx0) )-erf(-sqrt(ga)*xx0));
  
  fprintf (stdout,"----------------------------------------------------------------------\n");
  fprintf (stdout,"--------------------- C: gaussian check ------------------------------\n");
  fprintf (stdout," * gaussian line in the square [%.1f,%.1f]x[%.1f,%.1f] in a %dX%d grid   *\n", X0, X0+H, Y0, Y0+H, NMX, NMY);
  fprintf (stdout," * f(x,y) = y - yy0 - a0 exp[-ga (x - xx0)^2]                        *\n");
  fprintf (stdout,"----------------------------------------------------------------------\n");
  fprintf (stdout,"yy0:   %23.16e\n",yy0);
  fprintf (stdout,"xx0:   %23.16e\n",xx0);
  fprintf (stdout,"a0:    %23.16e\n",a0);
  fprintf (stdout,"ga:    %23.16e\n",ga);
  fprintf (stdout,"--------------------------------------------------------------------\n");
  fprintf (stdout,"analytical area : %23.16e\n",area_a);
  fprintf (stdout,"numerical  area : %23.16e\n\n",area_n);
  fprintf (stdout,"absolute error  : %23.16e\n",fabs(area_a-area_n));
  fprintf (stdout,"relative error  : %23.16e\n",fabs(area_a-area_n)/area_a); 
  fprintf (stdout,"-----------------------------------------------------------\n");
  if(!randominput) {
    fprintf (stdout,"with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2\n");
    fprintf (stdout,"-----------------------------------------------------------\n");
    fprintf (stdout,"analytical area :  3.3640894546075428e-01\n");
    fprintf (stdout,"numerical  area :  3.3640894546075722e-01\n\n");
    fprintf (stdout,"absolute error  :  2.9420910152566648e-15\n");
    fprintf (stdout,"relative error  :  8.7455790190926756e-15\n");
    fprintf (stdout,"------------------ C: end gaussian check ------------------\n");
    fprintf (stdout,"-----------------------------------------------------------\n");
  }
  fprintf (stdout,"\n");
      
  return;
}
