!* -------------------------------------------------------------------------- *
!* DESCRIPTION (reference phase where f(x,y) < 0):                            *
!* gaussian line in the square [0,1]x[0,1]                                    *
!* f(x,y) = y - YY0 - A0 exp[-GA (x - XX0)^2]                                 *
!* -------------------------------------------------------------------------- *

REAL(8) FUNCTION IMPL_FUNC(xyz)

  IMPLICIT NONE

  REAL(8), DIMENSION(*), INTENT(IN) :: xyz

  REAL(8), PARAMETER :: YY0 = 0.220D0, A0 = 0.51D0
  REAL(8), PARAMETER :: XX0 = 0.541D0, GA = 60.3D0

  REAL(8) :: x,y

  INTRINSIC DEXP

  x = xyz(1)
  y = xyz(2)
    
  IMPL_FUNC = y - YY0 - A0*DEXP(-GA*(x-XX0)*(x-XX0))

  RETURN 

END FUNCTION IMPL_FUNC

!* -------------------------------------------------------------------------- *

SUBROUTINE CHECK_AREA(areanum)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: areanum

  REAL(8), PARAMETER :: YY0 = 0.220D0, A0 = 0.51D0
  REAL(8), PARAMETER :: XX0 = 0.541D0, GA = 60.3D0

  REAL(8) :: areana

  areana = 0.3364089454607542483401167D0

  write(*,*) '-----------------------------------------------------------'
  write(*,*) '-------------------- F: gaussian check --------------------'
  write(*,100) areana
  write(*,101) areanum
  write(*,*) ' '
  write(*,102) DABS(areanum-areana)
  write(*,103) DABS(areanum-areana)/areana
   write(*,*) '----------------------------------------------------------'
  write(*,*) 'with Intel i7 3.4 GHz + Linux openSUSE 12.3 + gcc 4.7.2 -O3 '
  write(*,*) '-----------------------------------------------------------'
  write(*,*) 'analytical area :  3.3640894546075423E-01'
  write(*,*) 'numerical  area :  3.3640894546075717E-01'
  write(*,*) ' '
  write(*,*) 'absolute error  :  2.9420910152566648E-15'
  write(*,*) 'relative error  :  8.7455790190926772E-15'
  write(*,*) '------------------ F: end gaussian check ------------------'
  write(*,*) '-----------------------------------------------------------'
  write(*,*) ' '
  100 FORMAT(' analytical area : ', ES23.16)
  101 FORMAT(' numerical area  : ', ES23.16)
  102 FORMAT(' absolute error  : ', ES23.16)
  103 FORMAT(' relative error  : ', ES23.16)

END SUBROUTINE CHECK_AREA

