!****************************************************************************
!* Copyright (C) 2015 by Simone Bnà(a), Sandro Manservisi(a),               *
!* Ruben Scardovelli(a), Philip Yecko(b) and Stephane Zaleski(c,d)          *
!* (a) DIN–Lab. di Montecuccolino, Università di Bologna,                   *
!*     Via dei Colli 16, 40136 Bologna, Italy                               *
!* (b) Physics Department, Cooper Union, New York, NY, USA                  *
!* (c) Sorbonne Universités, UPMC Univ Paris 06, UMR 7190,                  *
!*     Institut Jean Le Rond d’Alembert, F-75005, Paris, France             *
!* (d) CNRS, UMR 7190, Institut Jean Le Rond d’Alembert, F-75005,           *
!*     Paris, France                                                        *
!*                                                                          *
!* You should have received a copy of the CPC license along with Vofi.      *
!* If not, see http://cpc.cs.qub.ac.uk/licence/licence.html.                *
!*                                                                          *
!* e-mail: ruben.scardovelli@unibo.it                                       *
!*                                                                          *
!****************************************************************************

PROGRAM GAUSSIAN

  USE GAUSSIAN_MOD
  IMPLICIT NONE
  INTEGER, PARAMETER :: NMX = 10, NMY = 10, NDIM = 3, N2D = 2
  REAL(8), PARAMETER :: X0 = 0.D0, Y0 = 0.D0, H = 1.D0    
  INTEGER :: nx=NMX, ny=NMY, ndim0=N2D
  INTEGER :: i,j,itrue
  REAL(8), DIMENSION(NMX,NMY) :: cc
  REAL(8), DIMENSION(3) :: xv,xloc
  CHARACTER(len=32) :: arg
  LOGICAL :: randominput = .FALSE. 

  REAL(8) :: h0,fh,area_n
  REAl(8), EXTERNAL :: VOFI_GET_CC,VOFI_GET_FH

  ! *********************************************************************
  ! PROGRAM TO INITIALIZE THE COLOR FUNCTION SCALAR FIELD 
  ! *********************************************************************
  DO i = 1, iargc()
    CALL getarg(i, arg)
    IF(arg == '-r' .OR. arg == '--randominput') THEN
      randominput = .TRUE.
    END IF
  END DO

  ! grid spacing 
  h0 = H/nx
  itrue = 1

  CALL INIT(randominput)

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

  CALL CHECK_AREA(area_n, randominput)

END PROGRAM GAUSSIAN
