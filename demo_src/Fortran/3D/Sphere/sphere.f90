!* -------------------------------------------------------------------------- *
!* DESCRIPTION (reference phase where f(x,y,z) < 0):                          *
!* ellipsoid/sphere inside the cube [0,1]x[0,1]x[0,1]                         *
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
  REAL(8), PARAMETER ::  ALPHA = 0.D0
  REAL(8), PARAMETER ::  A1 = 1.D0, B1 = 1.D0, C1 = 1.D0
  REAL(8), PARAMETER ::  XC = 0.D0, YC = 0.D0, ZC = 0.D0

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
  REAL(8), PARAMETER ::  ALPHA = 0.D0
  REAL(8), PARAMETER ::  A1 = 1.D0, B1 = 1.D0, C1 = 1.D0
  REAL(8), PARAMETER ::  XC = 0.D0, YC = 0.D0, ZC = 0.D0
 
  REAL(8) :: volana,invfrac

  INTRINSIC DABS

  invfrac = 8.D0
  volana = 4.D0*MYPI*A1*B1*C1/(3.D0*invfrac);

  write(*,*) '-----------------------------------------------------------'
  write(*,*) '------------------- F: 1/8 sphere check -------------------'
  write(*,100) volana
  write(*,101) volnum
  write(*,*) ' '
  write(*,102) DABS(volnum-volana)
  write(*,103) DABS(volnum-volana)/volana
  write(*,*) '-----------------------------------------------------------'
  write(*,*) 'with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2'
  write(*,*) '-----------------------------------------------------------'
  write(*,*) 'analytical volume:  5.2359877559829882E-01'
  write(*,*) 'numerical  volume:  5.2359877559829937E-01'
  write(*,*) ' '
  write(*,*) 'absolute error   :  5.5511151231257827E-16'
  write(*,*) 'relative error   :  1.0601848938211723E-15'
  write(*,*) '----------------- F: end 1/8 sphere check -----------------'
  write(*,*) '-----------------------------------------------------------'
  write(*,*) ' '
  100 FORMAT(' analytical volume: ', ES23.16)
  101 FORMAT(' numerical  volume: ', ES23.16)
  102 FORMAT(' absolute error   : ', ES23.16)
  103 FORMAT(' relative error   : ', ES23.16)

END SUBROUTINE CHECK_VOLUME
