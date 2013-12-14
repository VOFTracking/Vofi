!* -------------------------------------------------------------------------- *
!* DESCRIPTION (reference phase where f(x,y) < 0):                            *
!* rectangle inside the square [0,1]x[0,1]                                    *
!* PARAMETERS:                                                                *
!* (XC,YC): center of the rectangle; ALPHA: angle between two axes x' and x;  *
!* (A1,B1): half sides along the local x' and y' axes                         *
!* -------------------------------------------------------------------------- *

SUBROUTINE CHECK_AREA(areanum)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: areanum

  REAL(8), PARAMETER :: ALPHA = 0.0D0
  REAL(8), PARAMETER :: A1 = 0.2D0, XC = 0.52D0
  REAL(8), PARAMETER :: B1 = 0.3D0, YC = 0.44D0

  REAL(8) :: areana

  areana = 4.D0*A1*B1

   write(*,*) '-----------------------------------------------------'
   write(*,*) '----------------- F: rectangle check ----------------'
   write(*,*) ' '
   write(*,100) areana
   write(*,101) areanum
   write(*,102) DABS(areanum-areana)
   write(*,103) DABS(areanum-areana)/areana
   write(*,*) '-------------- F: end rectangle check ---------------'
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

  REAL(8), PARAMETER :: ALPHA = 0.0D0
  REAL(8), PARAMETER :: A1 = 0.2D0, XC = 0.52D0
  REAL(8), PARAMETER :: B1 = 0.3D0, YC = 0.44D0

  REAL(8) :: x,y,f0,f1,ca,sa

  INTRINSIC DCOS,DSIN,MAX

  x = xyz(1)
  y = xyz(2)

  ca = DCOS(ALPHA);
  sa = DSIN(ALPHA);
  f1 = - A1 - ((x-XC)*ca+(y-YC)*sa);
  f0 = f1;
  f1 = ((x-XC)*ca+(y-YC)*sa) - A1;
  f0 = MAX(f0,f1);
  f1 = -B1 - ((y-YC)*ca - (x-XC)*sa);
  f0 = MAX(f0,f1);
  f1 = ((y-YC)*ca - (x-XC)*sa) - B1;
  f0 = MAX(f0,f1);
    
  IMPL_FUNC = f0

END FUNCTION IMPL_FUNC

