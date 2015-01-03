!* -------------------------------------------------------------------------- *
!* DESCRIPTION (reference phase where f(x,y) < 0):                            *
!* sinusoidal line in the square [0,1]x[0,1]                                  *
!* f(x,y) = y - b0*sin(c0 pi x+ pi/d0) - a0                                   *
!* -------------------------------------------------------------------------- *

MODULE SINELINE_MOD

  IMPLICIT NONE

  REAL(8), PARAMETER :: mypi = 3.141592653589793238462643D0
  REAL(8) :: a0 = 0.5D0, b0 = 0.25D0
  REAL(8) :: c0 = 4.0D0, d0 = 14.0D0

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
  REAL(8), PARAMETER :: scalingfactor  = 0.1D0

    
  INTRINSIC RANDOM_NUMBER  
  
  CALL INIT_RANDOM_SEED()

  ! a0 --> from 0.45 to 0.55
  CALL RANDOM_NUMBER(seed) 
  a0 = 0.45D0 + seed*scalingfactor
 
  ! b0 --> from 0.20 to 0.30
  CALL RANDOM_NUMBER(seed)
  b0 = 0.20D0 + seed*scalingfactor

  ! c0 --> from 0.35 to 0.45 
  CALL RANDOM_NUMBER(seed) 
  c0 = 0.35D0 + seed*scalingfactor;

  ! d0 --> from 13.95 to 14.05
  CALL RANDOM_NUMBER(seed)  
  d0 = 13.95D0 + seed*scalingfactor;

  
END SUBROUTINE INIT


REAL(8) FUNCTION IMPL_FUNC(xyz)

  IMPLICIT NONE

  REAL(8), DIMENSION(*), INTENT(IN) :: xyz
  REAL(8) :: x,y

  INTRINSIC DSIN

  x = xyz(1)
  y = xyz(2)
    
  IMPL_FUNC = y - a0 - b0*DSIN(c0*mypi*x + mypi/d0)

  RETURN

END FUNCTION IMPL_FUNC

!* -------------------------------------------------------------------------- *

SUBROUTINE CHECK_AREA(areanum)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: areanum
  REAL(8) :: areana

  INTRINSIC DCOS,DABS

  areana = a0 + b0*(-DCOS((c0 + 1./d0)*mypi) + DCOS(mypi/d0))/(c0*mypi)

  write(*,*) '-----------------------------------------------------------'
  write(*,*) '-------------------- F: sine line check -------------------'
  write(*,100) areana
  write(*,101) areanum
  write(*,*) ' '
  write(*,102) DABS(areanum-areana)
  write(*,103) DABS(areanum-areana)/areana
  write(*,*) '-----------------------------------------------------------'
!   write(*,*) 'with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2'
!   write(*,*) '-----------------------------------------------------------'
!   write(*,*) 'analytical area :  5.0000000000000000E-01'
!   write(*,*) 'numerical  area :  4.9999999999993749E-01'
!   write(*,*) ' '
!   write(*,*) 'absolute error  :  6.2505556286396313E-14'
!   write(*,*) 'relative error  :  1.2501111257279263E-13'
!   write(*,*) '----------------- F: end sine line check ------------------'
!   write(*,*) '-----------------------------------------------------------'
  write(*,*) ' '
  100 FORMAT(' analytical area : ', ES23.16)
  101 FORMAT(' numerical  area : ', ES23.16)
  102 FORMAT(' absolute error  : ', ES23.16)
  103 FORMAT(' relative error  : ', ES23.16)

END SUBROUTINE CHECK_AREA

END MODULE SINELINE_MOD