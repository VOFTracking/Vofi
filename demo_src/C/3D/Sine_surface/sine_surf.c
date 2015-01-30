#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "sine_surf.h"

typedef const double creal;
typedef const int cint;
typedef double real;

/* -------------------------------------------------------------------------- *
 * DESCRIPTION (reference phase where f(x,y,z) < 0):                          *
 * sinusoidal surface inside the cube [0,1]x[0,1]x[0,1]                       *
 * f(x,y,z) = z - a0 - b0*sin(c1*pi*x + pi*d1)*sin(c1*pi*x + pi*e1)           *
 * -------------------------------------------------------------------------- */

static double a0 = 0.5;
static double b0 = 1./6.;
static double c1 = 1.6;
static double d1 = 1./7.;
static double e1 = 1./5.;


/* -------------------------------------------------------------------------- */

void init(cint randominput)
{

  if(randominput) {
  
    /* initialize random number generator. */
    srand(time(0));  
  

    double scalingfactor = 0.05;
 
    /* b0 --> from 1.15 to 1.20 */
    b0 = 1.15 + ((double)rand() / RAND_MAX)*scalingfactor;
  
    /* c1 --> from 1.6 to 1.65 */
    c1 = 1.6 + ((double)rand() / RAND_MAX)*scalingfactor;
   
    /* d1 --> from 0.12 to 0.17 */
    d1 = 0.12 + ((double)rand() / RAND_MAX)*scalingfactor;
 
    /* e1 --> from 0.15 to 0.25 */
    e1 = 0.15 + ((double)rand() / RAND_MAX)*scalingfactor;
    
  }
  
  return;
}

/* -------------------------------------------------------------------------- */

real impl_func(creal xy[])
{
  double x,y,z,f0;

  x = xy[0];
  y = xy[1];
  z = xy[2];

  f0 = z - a0 - b0*sin(M_PI*(c1*x+d1))*sin(M_PI*(c1*y+e1));

  return f0;
}
/* -------------------------------------------------------------------------- */

void check_volume(creal vol_n, cint randominput)
{
  double vol_a;

  vol_a = a0 + (b0/(c1*M_PI*c1*M_PI))*(cos(d1*M_PI) - cos((d1+c1)*M_PI))*(cos(e1*M_PI) - cos((e1+c1)*M_PI)); 

  fprintf (stdout,"--------------------------------------------------------------------------------------------------\n");
  fprintf (stdout,"--------------------- C: sinusoidal surface check ------------------------------------------------\n");
  fprintf (stdout," * sinusoidal surface inside the cube in the cube [%.1f,%.1f]x[%.1f,%.1f]x[%.1f,%.1f] in a %dX%dX%d grid  *\n", X0, X0+H, Y0, Y0+H, Z0, Z0+H, NMX, NMY, NMZ);
  fprintf (stdout," * f(x,y,z) = z - a0 - b0*sin(c1*pi*x + pi*d1)*sin(c1*pi*x + pi*e1)                              *\n");
  fprintf (stdout,"--------------------------------------------------------------------------------------------------\n");
  fprintf (stdout,"a0:    %23.16e\n",a0);
  fprintf (stdout,"b0:    %23.16e\n",b0);
  fprintf (stdout,"c1:    %23.16e\n",c1);
  fprintf (stdout,"d1:    %23.16e\n",d1);
  fprintf (stdout,"e1:    %23.16e\n",e1);
  fprintf (stdout,"-----------------------------------------------------------\n");
  fprintf (stdout,"analytical volume: %23.16e\n",vol_a);
  fprintf (stdout,"numerical  volume: %23.16e\n\n",vol_n);
  fprintf (stdout,"absolute error   : %23.16e\n",fabs(vol_a-vol_n));
  fprintf (stdout,"relative error   : %23.16e\n",fabs(vol_a-vol_n)/vol_a); 
  fprintf (stdout,"-----------------------------------------------------------\n");
  if(!randominput) {
    fprintf (stdout,"with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2\n");
    fprintf (stdout,"-----------------------------------------------------------\n");
    fprintf (stdout,"analytical volume:  5.0000000000000000e-01\n");
    fprintf (stdout,"numerical  volume:  5.0000000000000022e-01\n\n");
    fprintf (stdout,"absolute error   :  2.2204460492503131e-16\n");
    fprintf (stdout,"relative error   :  4.4408920985006262e-16\n");
    fprintf (stdout,"------------- C: end sinusoidal surface check -------------\n");
    fprintf (stdout,"-----------------------------------------------------------\n");
  }
  fprintf (stdout,"\n");

  return;
}
