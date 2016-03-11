!****************************************************************************
!* Copyright (C) 2015 by Simone Bnà(a), Sandro Manservisi(a),               *
!* Ruben Scardovelli(a), Philip Yecko(b) and Stephane Zaleski(c,d)          *
!* (a) DIN–Lab. di Montecuccolino, Università di Bologna,                   *
!*     Via dei Colli 16, 40136 Bologna, Italy                               *
!* (b) Physics Department, Cooper Union, New York, NY, USA                  *
!* (c) Sorbonne Universités, UPMC Univ Paris 06, UMR 7190,                  *
!*     Institut Jean Le Rond d’Alembert, F-75005, Paris, France             *
!* (d) CNRS, UMR 7190, Institut Jean Le Rond d’Alembert, F-75005,           *
!*     Paris, France                                                        *
!*                                                                          *
!* You should have received a copy of the CPC license along with Vofi.      *
!* If not, see http://cpc.cs.qub.ac.uk/licence/licence.html.                *
!*                                                                          *
!* e-mail: ruben.scardovelli@unibo.it                                       *
!*                                                                          *
!****************************************************************************

!* -------------------------------------------------------------------------- *
!* DESCRIPTION (reference phase where f(x,y,z) < 0):                          *
!* ellipsoidal cap inside the domain [0,1]x[0,1]x[0,1]                        *
!* f(x,y) = c1*x^2 + c2*y^2 + c3*x*y + c4*x + c5*y - c6                       *
!* f(x,y,z) = f(x,y) + (z-ZC)^2/C2                                            *
!* INPUT PARAMETERS:                                                          *
!* (XC,YC,ZC) center of the ellipsoid; ALPHA: angle between two axes x' and   *
!* x (in the x-y plane); (A1,B1,C1): semiaxis along the three ellipsoid       *
!* (local) axes: x',y',z'                                                     *
!* -------------------------------------------------------------------------- *

REAL(8) FUNCTION IMPL_FUNC (xyz)

  IMPLICIT NONE
  REAL(8), DIMENSION(3), INTENT(IN) :: xyz

  REAL(8), PARAMETER ::  MYPI = 3.141592653589793238462643D0
  REAL(8), PARAMETER ::  ALPHA = MYPI/3.D0
  REAL(8), PARAMETER ::  A1 = 4.D0, B1 = 5.D0, C1 = 6.D0
  REAL(8), PARAMETER ::  XC = 0.50D0, YC = 0.45D0, ZC = -5.97D0

  REAL(8) :: x,y,z,ca,sa,co1,co2,co3,co4,co5,co6,a2,b2,c2

  INTRINSIC DSIN,DCOS

  x = xyz(1)
  y = xyz(2)
  z = xyz(3)

  a2 = A1*A1
  b2 = B1*B1
  c2 = C1*C1
  ca = DCOS(ALPHA)
  sa = DSIN(ALPHA)
  co1 = ca*ca/a2 + sa*sa/b2
  co2 = sa*sa/a2 + ca*ca/b2
  co3 = 2.D0*ca*sa*(b2-a2)/(a2*b2)
  co4 = -(2.D0*co1*XC + co3*YC)
  co5 = -(2.D0*co2*YC + co3*XC)
  co6 = 1.0D0 - (co1*XC*XC + co2*YC*YC + co3*XC*YC)
  
  IMPL_FUNC = co1*x*x + co2*y*y + co3*x*y + co4*x + co5*y - co6 
  IMPL_FUNC = IMPL_FUNC + (z - ZC)*(z - ZC)/c2

  RETURN

END FUNCTION IMPL_FUNC

!* -------------------------------------------------------------------------- *

SUBROUTINE CHECK_VOLUME(volnum)

  IMPLICIT NONE
  REAL(8),INTENT(IN) :: volnum

  REAL(8), PARAMETER ::  MYPI = 3.141592653589793238462643D0
  REAL(8), PARAMETER ::  ALPHA = MYPI/3.D0
  REAL(8), PARAMETER ::  A1 = 4.D0, B1 = 5.D0, C1 = 6.D0
  REAL(8), PARAMETER ::  XC = 0.50D0, YC = 0.45D0, ZC = -5.97D0
 
  REAL(8) :: volana,h0

  INTRINSIC DABS

  h0 = C1 + ZC
  volana = MYPI*A1*B1*h0*h0*(1. - h0/(3.*C1))/C1

  write(*,*) '-----------------------------------------------------------'
  write(*,*) '------------------ F: cap check (1 cell) ------------------'
  write(*,100) volana
  write(*,101) volnum
  write(*,*) ' '
  write(*,102) DABS(volnum-volana)
  write(*,103) DABS(volnum-volana)/volana
  write(*,*) '-----------------------------------------------------------'
  write(*,*) 'with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2'
  write(*,*) '-----------------------------------------------------------'
  write(*,*) 'analytical volume:  9.4090699975015856e-03'
  write(*,*) 'numerical  volume:  9.4090641039290823E-03'
  write(*,*) ' '
  write(*,*) 'absolute error   :  5.8935725032877029E-09'
  write(*,*) 'relative error   :  6.2637141660681008E-07'
  write(*,*) '---------------- F: end cap check (1 cell) ----------------'
  write(*,*) '-----------------------------------------------------------'
  write(*,*) ' '
  100 FORMAT(' analytical volume: ', ES23.16)
  101 FORMAT(' numerical  volume: ', ES23.16)
  102 FORMAT(' absolute error   : ', ES23.16)
  103 FORMAT(' relative error   : ', ES23.16)

END SUBROUTINE CHECK_VOLUME
