!* -------------------------------------------------------------------------- *
!* DESCRIPTION (reference phase where f(x,y) < 0):                            *
!* ellipse inside the square [0,1]x[0,1]                                      *
!* f(x,y) = c1*x^2 + c2*y^2 + c3*x*y + c4*x + c5*y - c6                       *
!* PARAMETERS:                                                                *
!* (XC,YC): center of the ellipse; ALPHA: angle between two axes x' and x;    *
!* (A1,B1): semiaxis along the two ellipse axes (local) x' and y'             *
!* -------------------------------------------------------------------------- *

REAL(8) FUNCTION IMPL_FUNC(xyz)

  IMPLICIT NONE

  REAL(8), DIMENSION(*), INTENT(IN) :: xyz

  REAL(8), PARAMETER :: MYPI = 3.141592653589793238462643D0
  REAL(8), PARAMETER :: ALPHA = 0.48D0
  REAL(8), PARAMETER :: A1 = 0.17D0, B1=0.21D0
  REAL(8), PARAMETER :: XC = 0.523D0, YC = 0.475D0

  REAL(8) :: x,y,a2,b2,ca,sa,c1,c2,c3,c4,c5,c6

  INTRINSIC DCOS,DSIN

  x = xyz(1)
  y = xyz(2)
  
  a2 = A1*A1
  b2 = B1*B1
  ca = DCOS(ALPHA)
  sa = DSIN(ALPHA)
  c1 = ca*ca/a2 + sa*sa/b2
  c2 = sa*sa/a2 + ca*ca/b2
  c3 = 2.D0*ca*sa*(b2-a2)/(a2*b2)
  c4 = -(2.D0*c1*XC + c3*YC)
  c5 = -(2.D0*c2*YC + c3*XC)
  c6 = 1.0D0 - (c1*XC*XC + c2*YC*YC + c3*XC*YC)
  
  IMPL_FUNC = c1*x*x + c2*y*y + c3*x*y + c4*x + c5*y - c6

  RETURN

END FUNCTION IMPL_FUNC

!* -------------------------------------------------------------------------- *

SUBROUTINE CHECK_AREA(areanum)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: areanum

  REAL(8), PARAMETER :: MYPI = 3.141592653589793238462643D0
  REAL(8), PARAMETER :: ALPHA = 0.48D0
  REAL(8), PARAMETER :: A1 = 0.170D0, B1 = 0.210D0
  REAL(8), PARAMETER :: XC = 0.523D0, YC = 0.475D0

  REAL(8) :: areana

  INTRINSIC DABS

  areana = MYPI*A1*B1

  write(*,*) '-----------------------------------------------------------'
  write(*,*) '-------------------- F: ellipse check ---------------------'
  write(*,100) areana
  write(*,101) areanum
  write(*,*) ' '
  write(*,102) DABS(areanum-areana)
  write(*,103) DABS(areanum-areana)/areana
  write(*,*) '-----------------------------------------------------------'
  write(*,*) 'with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2'
  write(*,*) '-----------------------------------------------------------'
  write(*,*) 'analytical area :  1.1215485773315563E-01'
  write(*,*) 'numerical  area :  1.1215485773315677E-01'
  write(*,*) ' '
  write(*,*) 'absolute error  :  1.1379786002407855E-15'
  write(*,*) 'relative error  :  1.0146494081855289E-14'
  write(*,*) '------------------ F: end ellipse check -------------------'
  write(*,*) '-----------------------------------------------------------'
  write(*,*) ' '
  100 FORMAT(' analytical area : ', ES23.16)
  101 FORMAT(' numerical  area : ', ES23.16)
  102 FORMAT(' absolute error  : ', ES23.16)
  103 FORMAT(' relative error  : ', ES23.16)

END SUBROUTINE CHECK_AREA
