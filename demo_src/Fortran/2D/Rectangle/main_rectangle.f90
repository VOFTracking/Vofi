PROGRAM RECTANGLE

  USE RECTANGLE_MOD
  IMPLICIT NONE
  INTEGER, PARAMETER :: NMX = 10, NMY = 10, NDIM = 3, N2D = 2
  REAL(8), PARAMETER :: X0 = 0.D0, Y0 = 0.D0, H = 1.D0    
  INTEGER :: nx=NMX, ny=NMY, ndim0=N2D
  INTEGER :: i,j,itrue
  REAL(8), DIMENSION(NMX,NMY) :: cc
  REAL(8), DIMENSION(3) :: xv,xloc

  REAL(8) :: h0,fh,area_n
  REAl(8), EXTERNAL :: VOFI_GET_CC,VOFI_GET_FH

  ! *********************************************************************
  ! PROGRAM TO INITIALIZE THE COLOR FUNCTION SCALAR FIELD 
  ! *********************************************************************
 
  ! grid spacing 
  h0 = H/nx
  itrue = 1

  CALL INIT()

  ! starting point to get fh
  xv(1) = 0.5D0; xv(2) = 0.5D0; xv(3) = 0.0D0;
  fh = VOFI_GET_FH(IMPL_FUNC,xv,h0,ndim0,itrue)

  ! put now starting point in (X0,Y0,0.) to initialize the color function 
  xv(1) = X0; xv(2) = Y0;

  ! xloc: minor vertex of each cell of the grid 
  xloc(3) = 0.0D0
  DO j=1,ny
    DO i=1,nx
      xloc(1) = xv(1) + (i-1.D0)*h0
      xloc(2) = xv(2) + (j-1.D0)*h0

      cc(i,j) = VOFI_GET_CC(IMPL_FUNC,xloc,h0,fh,ndim0)
    END DO
  END DO 

  ! final global check 
   area_n = 0.0D0

  DO j=1,ny
    DO i=1,nx
      area_n = area_n + cc(i,j)
    END DO
  END DO 
  
  area_n = area_n*h0*h0

  CALL CHECK_AREA(area_n)

END PROGRAM RECTANGLE
