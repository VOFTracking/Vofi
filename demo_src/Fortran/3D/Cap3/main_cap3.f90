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

PROGRAM CAP3

  IMPLICIT NONE
  INTEGER, PARAMETER :: NMX = 2, NMY = 2, NMZ = 1 
  INTEGER, PARAMETER :: NDIM = 3, N2D = 2, N3D = 3
  REAL(8), PARAMETER :: X0 = -1.D0, Y0 = -1.D0, Z0 = 0.D0, H = 2.D0    
  INTEGER :: nx=NMX, ny=NMY, nz=NMZ, ndim0=N3D
  INTEGER :: i,j,k,itrue
  REAL(8), DIMENSION(NMX,NMY,NMZ) :: cc
  REAL(8), DIMENSION(3) :: xv,xloc
  REAL(8) :: h0,fh,vol_n
  REAl(8), EXTERNAL :: IMPL_FUNC,VOFI_GET_CC,VOFI_GET_FH

  ! *********************************************************************
  ! PROGRAM TO INITIALIZE THE COLOR FUNCTION SCALAR FIELD 
  ! *********************************************************************
 
  h0 = H/nx
  itrue = 1

  ! starting point to get fh
  xv(1) = 0.5D0; xv(2) = 0.5D0; xv(3) = 0.5D0;
  fh = VOFI_GET_FH(IMPL_FUNC,xv,h0,ndim0,itrue)

  ! put now starting point in (X0,Y0,Z0) to initialize the color function 
  xv(1) = X0; xv(2) = Y0; xv(3) = Z0;

  ! xloc: minor vertex of each cell of the grid 
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        xloc(1) = xv(1) + (i-1.D0)*h0
        xloc(2) = xv(2) + (j-1.D0)*h0
        xloc(3) = xv(3) + (k-1.D0)*h0

        cc(i,j,k) = VOFI_GET_CC(IMPL_FUNC,xloc,h0,fh,ndim0)
      END DO
    END DO
  END DO 

  ! final global check 
  vol_n = 0.0D0

  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        vol_n = vol_n + cc(i,j,k)
      END DO
    END DO
  END DO 
  
  vol_n = vol_n*h0*h0*h0

  CALL check_volume(vol_n)

END PROGRAM CAP3
