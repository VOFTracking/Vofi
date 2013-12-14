!* -------------------------------------------------------------------------- *
!* DESCRIPTION (reference phase where f(x,y,z) < 0):                          *
!* sinusoidal surface inside the cube [0,1]x[0,1]x[0,1]                       *
!* f(x,y,z) = z - A0 - B0*sin(C1*pi*x + pi*D1)*sin(C1*pi*x + pi*E1)           *
!* -------------------------------------------------------------------------- *

SUBROUTINE CHECK_VOLUME(volnum)

  IMPLICIT NONE
  REAL(8),INTENT(IN) :: volnum

  REAL(8), PARAMETER ::  MYPI = 3.141592653589793238462643D0
  REAL(8), PARAMETER ::  A0 = 0.5D0, B0 = (1.D0/6.D0), C1 = 1.6D0
  REAL(8), PARAMETER ::  D1 = (1.D0/7.D0), E1 = (1.6D0/5.D0)
 
  REAL(8) :: volana

  INTRINSIC DABS

  volana = 0.5D0;

   write(*,*) ' '
   write(*,*) '-----------------------------------------------------'
   write(*,*) '------------ F: sinusoidal surface check ------------'
   write(*,*) ' '
   write(*,100) volana
   write(*,101) volnum
   write(*,102) DABS(volnum-volana)
   write(*,103) DABS(volnum-volana)/volana
   write(*,*) '---------- F: end sinusoidal surface check ----------'
   write(*,*) '-----------------------------------------------------'
   write(*,*) ' '

   100 FORMAT('analytical volume:', ES23.16)
   101 FORMAT('numerical  volume:', ES23.16)
   102 FORMAT('absolute error   :', ES23.16)
   103 FORMAT('relative error   :', ES23.16)

  RETURN

END SUBROUTINE CHECK_VOLUME

!* -------------------------------------------------------------------------- *

REAL(8) FUNCTION IMPL_FUNC (xyz)

  IMPLICIT NONE
  REAL(8), DIMENSION(3), INTENT(IN) :: xyz

  REAL(8), PARAMETER ::  MYPI = 3.141592653589793238462643D0
  REAL(8), PARAMETER ::  A0 = 0.5D0, B0 = (1.D0/6.D0), C1 = 1.6D0
  REAL(8), PARAMETER ::  D1 = (1.D0/7.D0), E1 = (1.6D0/5.D0)

  REAL(8) :: x,y,z

  INTRINSIC DSIN

  x = xyz(1)
  y = xyz(2)
  z = xyz(3)
  
  IMPL_FUNC = z - A0 - B0*DSIN(MYPI*(C1*x+D1))*DSIN(MYPI*(C1*y+E1))

END FUNCTION IMPL_FUNC
