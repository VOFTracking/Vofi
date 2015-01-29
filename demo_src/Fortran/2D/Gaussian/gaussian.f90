!* -------------------------------------------------------------------------- *
!* DESCRIPTION (reference phase where f(x,y) < 0):                            *
!* gaussian line in the square [0,1]x[0,1]                                    *
!* f(x,y) = y - yy0 - a0 exp[-ga (x - xx0)^2]                                 *
!* -------------------------------------------------------------------------- *

MODULE GAUSSIAN_MOD

  IMPLICIT NONE

  REAL(8), PARAMETER :: mypi = 3.141592653589793238462643D0
  REAL(8) :: yy0 = 0.220D0, a0 = 0.51D0
  REAL(8) :: xx0 = 0.541D0, ga = 60.3D0

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


SUBROUTINE INIT(randominput)  

  IMPLICIT NONE
  
  LOGICAL, INTENT(IN) :: randominput
  REAL(8) :: seed
  REAL(8), PARAMETER :: scalingfactor  = 0.06D0
    
  INTRINSIC RANDOM_NUMBER  

  IF(randominput .EQV. .TRUE.) THEN
  
    CALL INIT_RANDOM_SEED()

    ! yy0 --> from 0.19 to 0.25
    CALL RANDOM_NUMBER(seed) 
    yy0 = 0.19D0 + seed*scalingfactor
 
    ! xx0 --> from 0.511 to 0.561 
    CALL RANDOM_NUMBER(seed)
    xx0 = 0.511D0 + seed*scalingfactor

    ! a0 --> from 0.48 to 0.54
    CALL RANDOM_NUMBER(seed) 
    a0 = 0.48D0 + seed*scalingfactor;

    ! ga --> from 60.00 to 0.60.6 
    CALL RANDOM_NUMBER(seed)  
    ga = 60.D0 + seed*scalingfactor;

  END IF
  
END SUBROUTINE INIT


REAL(8) FUNCTION IMPL_FUNC(xyz)

  IMPLICIT NONE

  REAL(8), DIMENSION(*), INTENT(IN) :: xyz
  REAL(8) :: x,y

  INTRINSIC DEXP

  x = xyz(1)
  y = xyz(2)
    
  IMPL_FUNC = y - yy0 - a0*DEXP(-ga*(x-xx0)*(x-xx0))

  RETURN 

END FUNCTION IMPL_FUNC

!* -------------------------------------------------------------------------- *

SUBROUTINE CHECK_AREA(areanum, randominput)

  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: randominput
  REAL(8), INTENT(IN) :: areanum
  REAL(8) :: areana

  ! analytical integration with x in [0,1] 
  areana = yy0 + 0.5*a0*SQRT(mypi/ga)*(ERF(SQRT(ga)*(1.-xx0) )-ERF(-SQRT(ga)*xx0));

  write(*,*) '-----------------------------------------------------------'
  write(*,*) '-------------------- F: gaussian check --------------------'
  write(*,100) areana
  write(*,101) areanum
  write(*,*) ' '
  write(*,102) DABS(areanum-areana)
  write(*,103) DABS(areanum-areana)/areana
  write(*,*) '-----------------------------------------------------------'
  IF(randominput .EQV. .FALSE.) THEN
    write(*,*) 'with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2'
    write(*,*) '-----------------------------------------------------------'
    write(*,*) 'analytical area :  3.3640894546075423E-01'
    write(*,*) 'numerical  area :  3.3640894546075717E-01'
    write(*,*) ' '
    write(*,*) 'absolute error  :  2.9420910152566648E-15'
    write(*,*) 'relative error  :  8.7455790190926772E-15'
    write(*,*) '------------------ F: end gaussian check ------------------'
    write(*,*) '-----------------------------------------------------------'
  END IF
  write(*,*) ' '
  100 FORMAT(' analytical area : ', ES23.16)
  101 FORMAT(' numerical  area : ', ES23.16)
  102 FORMAT(' absolute error  : ', ES23.16)
  103 FORMAT(' relative error  : ', ES23.16)

END SUBROUTINE CHECK_AREA

END MODULE GAUSSIAN_MOD

