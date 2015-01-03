!* -------------------------------------------------------------------------- *
!* DESCRIPTION (reference phase where f(x,y) < 0):                            *
!* ellipse inside the square [0,1]x[0,1]                                      *
!* f(x,y) = c1*x^2 + c2*y^2 + c3*x*y + c4*x + c5*y - c6                       *
!* PARAMETERS:                                                                *
!* (xC,yC): center of the ellipse; alpha: angle between two axes x' and x;    *
!* (a1,b1): semiaxis along the two ellipse axes (local) x' and y'             *
!* -------------------------------------------------------------------------- *

MODULE ELLIPSE_MOD

  IMPLICIT NONE

  REAL(8), PARAMETER :: mypi = 3.141592653589793238462643D0
  REAL(8) :: alpha = 0.48D0
  REAL(8) :: a1 = 0.17D0, b1=0.21D0
  REAL(8) :: xc = 0.523D0, yc = 0.475D0

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
  REAL(8), PARAMETER :: axesscalingfactor  = 0.04D0
  REAL(8), PARAMETER :: alphascalingfactor  = mypi*0.5D0
  REAL(8), PARAMETER :: centerscalingfactor  = 0.1D0
  REAL(8), PARAMETER :: x0  = 0.45D0
    
  INTRINSIC RANDOM_NUMBER  
  
  CALL INIT_RANDOM_SEED()
  
  ! a1 --> from 0.15 to 0.19
  CALL RANDOM_NUMBER(seed) 
  a1 = 0.15D0 + seed*axesscalingfactor
  
  ! b1 --> from 0.19 to 0.23
  CALL RANDOM_NUMBER(seed)
  b1 = 0.19D0 + seed*axesscalingfactor

  ! alpha --> from zero to pi/2
  CALL RANDOM_NUMBER(seed) 
  alpha = seed*alphascalingfactor;

  ! xc & yc --> from 0.4 to 0.6 
  CALL RANDOM_NUMBER(seed)  
  xc = x0 + seed*centerscalingfactor;

  CALL RANDOM_NUMBER(seed)
  yc = x0 + seed*centerscalingfactor;
  
END SUBROUTINE INIT


REAL(8) FUNCTION IMPL_FUNC(xyz)

  IMPLICIT NONE

  REAL(8), DIMENSION(*), INTENT(IN) :: xyz
  REAL(8) :: x,y,a2,b2,ca,sa,c1,c2,c3,c4,c5,c6

  INTRINSIC DCOS,DSIN

  x = xyz(1)
  y = xyz(2)
  
  a2 = a1*a1
  b2 = b1*b1
  ca = DCOS(alpha)
  sa = DSIN(alpha)
  c1 = ca*ca/a2 + sa*sa/b2
  c2 = sa*sa/a2 + ca*ca/b2
  c3 = 2.D0*ca*sa*(b2-a2)/(a2*b2)
  c4 = -(2.D0*c1*xc + c3*yc)
  c5 = -(2.D0*c2*yc + c3*xc)
  c6 = 1.0D0 - (c1*xc*xc + c2*yc*yc + c3*xc*yc)
  
  IMPL_FUNC = c1*x*x + c2*y*y + c3*x*y + c4*x + c5*y - c6

  RETURN

END FUNCTION IMPL_FUNC

!* -------------------------------------------------------------------------- *

SUBROUTINE CHECK_AREA(areanum)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: areanum
  REAL(8) :: areana

  INTRINSIC DABS

  areana = mypi*a1*b1

  write(*,*) '-----------------------------------------------------------'
  write(*,*) '-------------------- F: ellipse check ---------------------'
  write(*,100) areana
  write(*,101) areanum
  write(*,*) ' '
  write(*,102) DABS(areanum-areana)
  write(*,103) DABS(areanum-areana)/areana
  write(*,*) '-----------------------------------------------------------'
!   write(*,*) 'with Intel i7 3.4 GHz + Linux openSUSE 13.1 + gcc 4.8.1 -O2'
!   write(*,*) '-----------------------------------------------------------'
!   write(*,*) 'analytical area :  1.1215485773315563E-01'
!   write(*,*) 'numerical  area :  1.1215485773315677E-01'
!   write(*,*) ' '
!   write(*,*) 'absolute error  :  1.1379786002407855E-15'
!   write(*,*) 'relative error  :  1.0146494081855289E-14'
!   write(*,*) '------------------ F: end ellipse check -------------------'
!   write(*,*) '-----------------------------------------------------------'
  write(*,*) ' '
  100 FORMAT(' analytical area : ', ES23.16)
  101 FORMAT(' numerical  area : ', ES23.16)
  102 FORMAT(' absolute error  : ', ES23.16)
  103 FORMAT(' relative error  : ', ES23.16)

END SUBROUTINE CHECK_AREA

END MODULE ELLIPSE_MOD