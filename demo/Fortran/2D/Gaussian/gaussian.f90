!* -------------------------------------------------------------------------- *
!* DESCRIPTION (reference phase where f(x,y) < 0):                            *
!* gaussian line in the square [0,1]x[0,1]                                    *
!* f(x,y) = y - YY0 - A0 exp[-GA (x - XX0)^2]                                 *
!* -------------------------------------------------------------------------- *

SUBROUTINE CHECK_AREA(areanum)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: areanum

  REAL(8), PARAMETER :: YY0 = 0.220D0, A0 = 0.51D0
  REAL(8), PARAMETER :: XX0 = 0.541D0, GA = 60.3D0

  REAL(8) :: areana

  areana = 0.3364089454607542483401167D0

   write(*,*) '-----------------------------------------------------'
   write(*,*) '----------------- F: gaussian check -----------------'
   write(*,*) ' '
   write(*,100) areana
   write(*,101) areanum
   write(*,102) DABS(areanum-areana)
   write(*,103) DABS(areanum-areana)/areana
   write(*,*) '--------------- F: end gaussian check ---------------'
   write(*,*) '-----------------------------------------------------'
   write(*,*) ' '
   100 FORMAT('analytical area :', ES23.16)
   101 FORMAT('numerical area  :', ES23.16)
   102 FORMAT('absolute error  :', ES23.16)
   103 FORMAT('relative error  :', ES23.16)

  RETURN

END SUBROUTINE CHECK_AREA

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

END FUNCTION IMPL_FUNC

