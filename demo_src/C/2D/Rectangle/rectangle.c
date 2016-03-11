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
#include "rectangle.h"

typedef const double creal;
typedef const int cint;
typedef double real;

/* -------------------------------------------------------------------------- *
 * DESCRIPTION (reference phase where f(x,y) < 0):                            *
 * rectangle inside the square [0,1]x[0,1]                                    *
 * PARAMETERS:                                                                *
 * (xc,yc): center of the rectangle; alpha: angle between two axes x' and x;  *
 * (a1,b1): half sides along the local x' and y' axes                         *
 * -------------------------------------------------------------------------- */

static double a1 = 0.2;
static double b1 = 0.3;
static double alpha = 0.;
static double xc = 0.52;
static double yc = 0.44;

/* -------------------------------------------------------------------------- */

void init(cint randominput)
{

  if(randominput) {
    
    /* initialize random number generator. */
    srand(time(0));  
  
    /* a1 --> from 0.17 to 0.23 */
    double axesscalingfactor = 0.06;
    a1 = 0.17 + ((double)rand() / RAND_MAX)*axesscalingfactor;
    /* a1 --> from 0.27 to 0.33 */
    b1 = 0.27 + ((double)rand() / RAND_MAX)*axesscalingfactor;
  
    /* alpha --> from zero to pi/2 */
    double alphascalingfactor = M_PI_2;
    alpha = ((double)rand() / RAND_MAX)*alphascalingfactor;
  
    /* xc & yc --> from 0.4 to 0.6 */
    double x0 = 0.45;
    double centerscalingfactor = 0.1;
    xc = x0 + ((double)rand() / RAND_MAX)*centerscalingfactor;
    yc = x0 + ((double)rand() / RAND_MAX)*centerscalingfactor;
  
  }
  
  return;
}

/* -------------------------------------------------------------------------- */

real impl_func(creal xy[]) 
{
  real x,y,f0,f1,ca,sa;

  x = xy[0];
  y = xy[1];

  ca = cos(alpha);
  sa = sin(alpha);
  f1 = - a1 - ((x-xc)*ca+(y-yc)*sa);
  f0 = f1;
  /* f0 = MAX(f0,f1); */
  f1 = ((x-xc)*ca+(y-yc)*sa) - a1;
  f0 = MAX(f0,f1);
  f1 = -b1 - ((y-yc)*ca - (x-xc)*sa);
  f0 = MAX(f0,f1);
  f1 = ((y-yc)*ca - (x-xc)*sa) - b1;
  f0 = MAX(f0,f1);

  return f0;
}

/* -------------------------------------------------------------------------- */

void check_area(creal area_n, cint randominput)
{
  real area_a;
  
  area_a = 4.*a1*b1;

  fprintf (stdout,"----------------------------------------------------------------------\n");
  fprintf (stdout,"--------------------- C: rectangle check -----------------------------\n");
  fprintf (stdout," * rectangle inside the square [%.1f,%.1f]x[%.1f,%.1f] in a %dX%d grid   *\n", X0, X0+H, Y0, Y0+H, NMX, NMY);
  fprintf (stdout," * f(x,y) = c1*x^2 + c2*y^2 + c3*x*y + c4*x + c5*y - c6              *\n");
  fprintf (stdout," * PARAMETERS:                                                       *\n");
  fprintf (stdout," * (a1,b1): half sides along the local x' and y' axes                *\n");
  fprintf (stdout," * ALPHA: angle between two axes x' and x;                           *\n");                
  fprintf (stdout," * (xc,yc): center of the rectangle;                                 *\n");
  fprintf (stdout,"----------------------------------------------------------------------\n");
  fprintf (stdout,"a1:    %23.16e\n",a1);
  fprintf (stdout,"b1:    %23.16e\n",b1);
  fprintf (stdout,"alpha: %23.16e\n",alpha);
  fprintf (stdout,"xc:    %23.16e\n",xc);
  fprintf (stdout,"yc:    %23.16e\n",yc);
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"analytical area : %23.16e\n",area_a);
  fprintf (stdout,"numerical  area : %23.16e\n\n",area_n);
  fprintf (stdout,"absolute error  : %23.16e\n",fabs(area_a-area_n));
  fprintf (stdout,"relative error  : %23.16e\n",fabs(area_a-area_n)/area_a); 
  fprintf (stdout,"-----------------------------------------------------------\n");
  if(!randominput) {
    fprintf (stdout,"with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2\n");
    fprintf (stdout,"-----------------------------------------------------------\n");
    fprintf (stdout,"analytical area :  2.3999999999999999e-01\n");
    fprintf (stdout,"numerical  area :  2.3999999999999996e-01\n\n");
    fprintf (stdout,"absolute error  :  2.7755575615628914e-17\n");
    fprintf (stdout,"relative error  :  1.1564823173178715e-16\n");
    fprintf (stdout,"----------------- C: end rectangle check ------------------\n");
    fprintf (stdout,"-----------------------------------------------------------\n");
  }
  fprintf (stdout,"\n");

  return;
}
