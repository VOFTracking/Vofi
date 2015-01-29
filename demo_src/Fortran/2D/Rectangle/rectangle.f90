!* -------------------------------------------------------------------------- *
!* DESCRIPTION (reference phase where f(x,y) < 0):                            *
!* rectangle inside the square [0,1]x[0,1]                                    *
!* PARAMETERS:                                                                *
!* (XC,YC): center of the rectangle; ALPHA: angle between two axes x' and x;  *
!* (A1,B1): half sides along the local x' and y' axes                         *
!* -------------------------------------------------------------------------- *

MODULE RECTANGLE_MOD

  IMPLICIT NONE

  REAL(8), PARAMETER :: mypi = 3.141592653589793238462643D0
  REAL(8) :: alpha = 0.D0
  REAL(8) :: a1 = 0.2D0, b1=0.3D0
  REAL(8) :: xc = 0.52D0, yc = 0.44D0

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
  REAL(8), PARAMETER :: axesscalingfactor  = 0.06D0
  REAL(8), PARAMETER :: alphascalingfactor  = mypi*0.5D0
  REAL(8), PARAMETER :: centerscalingfactor  = 0.1D0
  REAL(8), PARAMETER :: x0  = 0.45D0
    
  INTRINSIC RANDOM_NUMBER  
  
  IF(randominput .EQV. .TRUE.) THEN

    CALL INIT_RANDOM_SEED()
  
    ! a1 --> from 0.15 to 0.19
    CALL RANDOM_NUMBER(seed) 
    a1 = 0.17D0 + seed*axesscalingfactor
  
    ! b1 --> from 0.19 to 0.23
    CALL RANDOM_NUMBER(seed)
    b1 = 0.27D0 + seed*axesscalingfactor

    ! alpha --> from zero to pi/2
    CALL RANDOM_NUMBER(seed) 
    alpha = seed*alphascalingfactor;

    ! xc & yc --> from 0.4 to 0.6 
    CALL RANDOM_NUMBER(seed)  
    xc = x0 + seed*centerscalingfactor;

    CALL RANDOM_NUMBER(seed)
    yc = x0 + seed*centerscalingfactor;
  
  END IF

END SUBROUTINE INIT


REAL(8) FUNCTION IMPL_FUNC(xyz)

  IMPLICIT NONE

  REAL(8), DIMENSION(*), INTENT(IN) :: xyz
  REAL(8) :: x,y,f0,f1,ca,sa

  INTRINSIC DCOS,DSIN,MAX

  x = xyz(1)
  y = xyz(2)

  ca = DCOS(alpha);
  sa = DSIN(alpha);
  f1 = - a1 - ((x-xc)*ca+(y-yc)*sa);
  f0 = f1;
  f1 = ((x-xc)*ca+(y-yc)*sa) - a1;
  f0 = MAX(f0,f1);
  f1 = -b1 - ((y-yc)*ca - (x-xc)*sa);
  f0 = MAX(f0,f1);
  f1 = ((y-yc)*ca - (x-xc)*sa) - b1;
  f0 = MAX(f0,f1);
    
  IMPL_FUNC = f0

  RETURN

END FUNCTION IMPL_FUNC

!* -------------------------------------------------------------------------- *

SUBROUTINE CHECK_AREA(areanum, randominput)

  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: randominput
  REAL(8), INTENT(IN) :: areanum
  REAL(8) :: areana

  areana = 4.D0*a1*b1

  write(*,*) '-----------------------------------------------------------'
  write(*,*) '-------------------- F: rectangle check -------------------'
  write(*,100) areana
  write(*,101) areanum
  write(*,*) ' '
  write(*,102) DABS(areanum-areana)
  write(*,103) DABS(areanum-areana)/areana
  write(*,*) '-----------------------------------------------------------'
  IF(randominput .EQV. .FALSE.) THEN
    write(*,*) 'with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2'
    write(*,*) '-----------------------------------------------------------'
    write(*,*) 'analytical area :  2.3999999999999999E-01'
    write(*,*) 'numerical  area :  2.3999999999999991E-01'
    write(*,*) ' '
    write(*,*) 'absolute error  :  8.3266726846886741E-17'
    write(*,*) 'relative error  :  3.4694469519536142E-16'
    write(*,*) '----------------- F: end rectangle check ------------------'
    write(*,*) '-----------------------------------------------------------'
  END IF
  write(*,*) ' '
  100 FORMAT(' analytical area : ', ES23.16)
  101 FORMAT(' numerical  area : ', ES23.16)
  102 FORMAT(' absolute error  : ', ES23.16)
  103 FORMAT(' relative error  : ', ES23.16)

END SUBROUTINE CHECK_AREA

END MODULE RECTANGLE_MOD
