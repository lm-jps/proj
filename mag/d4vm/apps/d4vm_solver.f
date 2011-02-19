      SUBROUTINE d4vm_solver(S_MTX, U0, V0, W0, UX, VX, WX, UY, VY, WY,
     &                       xsize, ysize, threshold)

C.... Dave4vm matrix solver
C....
C.... Written by Jacob Hageman, jacob.hageman@nasa.gov, (301)286-1803
C....
C.... Solves the A matrix using SVD with a divide-and-conquer algorithm
C....
C.... INPUTS:
C....   S_MTX - Matrix to solve
C....   XSIZE - X dimension size
C....   YSIZE - Y dimension size
C....   THRESHOLD - Threshold to resolve aperture problem
C....
C.... OUTPUTS:
C....   U0, V0, W0, UX, VX, WX, UY, VY, WY - Velocities
C....

      IMPLICIT NONE

C.... INPUTS & PARAMETERS
      INTEGER*8 xsize, ysize
      REAL*4  S_MTX(10, 10, xsize*ysize)
      REAL*8  threshold

C.... OUTPUTS
      REAL*4  U0(xsize*ysize), V0(xsize*ysize), W0(xsize*ysize)
      REAL*4  UX(xsize*ysize), VX(xsize*ysize), WX(xsize*ysize)
      REAL*4  UY(xsize*ysize), VY(xsize*ysize), WY(xsize*ysize)

C.... WORKING VARIABLES
      INTEGER i, j, k, l
      INTEGER pix_count
      REAL*4 total

      INTEGER LWORK, LIWORK, RANK, INFO
      PARAMETER (LWORK = 2000, LIWORK = 500)
      REAL*8 A(9,9), B(9), S(9), WORK(LWORK), RCOND, IWORK(LIWORK)

C.... Initialize
      RCOND = -1.d0
      U0 = 0
      V0 = 0
      W0 = 0
      UX = 0
      VX = 0
      WX = 0
      UY = 0
      VY = 0
      WY = 0

C.... Loop over pixel elements
      pix_count = 0
      DO i = 1, xsize*ysize

        total = S_MTX(1,1,i) + S_MTX(2,2,i) + S_MTX(3,3,i) + 
     &          S_MTX(4,4,i) + S_MTX(5,5,i) + S_MTX(6,6,i) +
     &          S_MTX(7,7,i) + S_MTX(8,8,i) + S_MTX(9,9,i)

        IF (total .GT. threshold) THEN

C.... Count number of good pixels
          pix_count = pix_count + 1

          A = DBLE(S_MTX(1:9, 1:9, i))
          B = DBLE(-S_MTX(1:9, 10, i))

          CALL DGELSD(9, 9, 1, A, 9, B, 9, S, RCOND, RANK, WORK,
     &                LWORK, IWORK, INFO)

          IF (INFO .EQ. 0) THEN
            U0(i) = B(1)
            V0(i) = B(2)
            UX(i) = B(3)
            VY(i) = B(4)
            UY(i) = B(5)
            VX(i) = B(6)
            W0(i) = B(7)
            WX(i) = B(8)
            WY(i) = B(9)
          ENDIF
        ENDIF
      ENDDO

      END
