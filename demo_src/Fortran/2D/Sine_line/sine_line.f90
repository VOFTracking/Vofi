!* -------------------------------------------------------------------------- *
!* DESCRIPTION (reference phase where f(x,y) < 0):                            *
!* sinusoidal line in the square [0,1]x[0,1]                                  *
!* f(x,y) = y - B0*sin(C0 pi x+ pi/D0) - A0                                   *
!* -------------------------------------------------------------------------- *

REAL(8) FUNCTION IMPL_FUNC(xyz)

  IMPLICIT NONE

  REAL(8), DIMENSION(*), INTENT(IN) :: xyz

  REAL(8), PARAMETER :: MYPI = 3.141592653589793238462643D0
  REAL(8), PARAMETER :: A0 = 0.5D0, B0 = 0.25D0
  REAL(8), PARAMETER :: C0 = 4.0D0, D0 = 14.0D0

  REAL(8) :: x,y

  INTRINSIC DSIN

  x = xyz(1)
  y = xyz(2)
    
  IMPL_FUNC = y - A0 - B0*DSIN(C0*MYPI*x + MYPI/D0)

  RETURN

END FUNCTION IMPL_FUNC

!* -------------------------------------------------------------------------- *

SUBROUTINE CHECK_AREA(areanum)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: areanum

  REAL(8), PARAMETER :: MYPI = 3.141592653589793238462643D0
  REAL(8), PARAMETER :: A0 = 0.5D0, B0 = 0.25D0
  REAL(8), PARAMETER :: C0 = 4.0D0, D0 = 14.0D0

  REAL(8) :: areana

  INTRINSIC DCOS,DABS

  areana = A0 + B0*(-DCOS((C0 + 1./D0)*MYPI) + DCOS(MYPI/D0))/(C0*MYPI)

  write(*,*) '-----------------------------------------------------------'
  write(*,*) '-------------------- F: sine line check -------------------'
  write(*,100) areana
  write(*,101) areanum
  write(*,*) ' '
  write(*,102) DABS(areanum-areana)
  write(*,103) DABS(areanum-areana)/areana
  write(*,*) '-----------------------------------------------------------'
  write(*,*) 'with Intel i7 3.4 GHz + Linux openSUSE 12.3 + gcc 4.7.2 -O3'
  write(*,*) '-----------------------------------------------------------'
  write(*,*) 'analytical area :  5.0000000000000000E-01'
  write(*,*) 'numerical  area :  4.9999999999993749E-01'
  write(*,*) ' '
  write(*,*) 'absolute error  :  6.2505556286396313E-14'
  write(*,*) 'relative error  :  1.2501111257279263E-13'
  write(*,*) '----------------- F: end sine line check ------------------'
  write(*,*) '-----------------------------------------------------------'
  write(*,*) ' '
  100 FORMAT(' analytical area : ', ES23.16)
  101 FORMAT(' numerical area  : ', ES23.16)
  102 FORMAT(' absolute error  : ', ES23.16)
  103 FORMAT(' relative error  : ', ES23.16)

END SUBROUTINE CHECK_AREA

