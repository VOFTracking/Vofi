!* -------------------------------------------------------------------------- *
!* DESCRIPTION (reference phase where f(x,y,z) < 0):                          *
!* sinusoidal surface inside the cube [0,1]x[0,1]x[0,1]                       *
!* f(x,y,z) = z - A0 - B0*sin(C1*pi*x + pi*D1)*sin(C1*pi*x + pi*E1)           *
!* -------------------------------------------------------------------------- *

MODULE SINESURFACE_MOD

  IMPLICIT NONE

  REAL(8), PARAMETER ::  MYPI = 3.141592653589793238462643D0
  REAL(8) ::  A0 = 0.5D0, B0 = (1.D0/6.D0), C1 = 1.6D0
  REAL(8) ::  D1 = (1.D0/7.D0), E1 = (1.D0/5.D0)

  CONTAINS

!  initialize a random seed from the system clock at every run 
SUBROUTINE INIT_RANDOM_SEED()

      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
      INTRINSIC RANDOM_SEED, SYSTEM_CLOCK 

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
      
END SUBROUTINE INIT_RANDOM_SEED


SUBROUTINE INIT()  

  IMPLICIT NONE
  
  REAL(8) :: seed
  REAL(8), PARAMETER :: scalingfactor  = 0.05D0
    
  INTRINSIC RANDOM_NUMBER  
  
  CALL INIT_RANDOM_SEED()
  
  !B0 --> from 1.15 to 1.20 
  CALL RANDOM_NUMBER(seed) 
  B0 = 1.15D0 + seed*scalingfactor
  
  !C1 --> from 1.6 to 1.65
  CALL RANDOM_NUMBER(seed) 
  C1 = 1.6D0 + seed*scalingfactor
   
  !D1 --> from 0.12 to 0.17
  CALL RANDOM_NUMBER(seed) 
  D1 = 0.12D0 + seed*scalingfactor
 
  !E1 --> from 0.15 to 0.25
  CALL RANDOM_NUMBER(seed) 
  E1 = 0.15D0 + seed*scalingfactor
  
END SUBROUTINE INIT

!* -------------------------------------------------------------------------- *
  
REAL(8) FUNCTION IMPL_FUNC (xyz)

  IMPLICIT NONE
  REAL(8), DIMENSION(3), INTENT(IN) :: xyz
  REAL(8) :: x,y,z

  INTRINSIC DSIN

  x = xyz(1)
  y = xyz(2)
  z = xyz(3)
  
  IMPL_FUNC = z - A0 - B0*DSIN(MYPI*(C1*x+D1))*DSIN(MYPI*(C1*y+E1))

  RETURN

END FUNCTION IMPL_FUNC

!* -------------------------------------------------------------------------- *

SUBROUTINE CHECK_VOLUME(volnum)

  IMPLICIT NONE
  REAL(8),INTENT(IN) :: volnum
  REAL(8) :: volana

  INTRINSIC DABS, DCOS

  !volana = 0.5D0;
  volana = A0 + (B0/(C1*MYPI*C1*MYPI))*(DCOS(D1*MYPI) - DCOS((D1+C1)*MYPI))*(DCOS(E1*MYPI) - DCOS((E1+C1)*MYPI)); 

  write(*,*) '-----------------------------------------------------------'
  write(*,*) '--------------- F: sinusoidal surface check ---------------'
  write(*,100) volana
  write(*,101) volnum
  write(*,*) ' '
  write(*,102) DABS(volnum-volana)
  write(*,103) DABS(volnum-volana)/volana
  write(*,*) '-----------------------------------------------------------'
!   write(*,*) 'with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2'
!   write(*,*) '-----------------------------------------------------------'
!   write(*,*) 'analytical volume:  5.0000000000000000E-01'
!   write(*,*) 'numerical  volume:  4.9999999999999994E-01'
!   write(*,*) ' '
!   write(*,*) 'absolute error   :  5.5511151231257827E-17'
!   write(*,*) 'relative error   :  1.1102230246251565E-16'
!   write(*,*) '------------- F: end sinusoidal surface check -------------'
!   write(*,*) '-----------------------------------------------------------'
  write(*,*) ' '
   100 FORMAT(' analytical volume: ', ES23.16)
   101 FORMAT(' numerical  volume: ', ES23.16)
   102 FORMAT(' absolute error   : ', ES23.16)
   103 FORMAT(' relative error   : ', ES23.16)

END SUBROUTINE CHECK_VOLUME

END MODULE SINESURFACE_MOD