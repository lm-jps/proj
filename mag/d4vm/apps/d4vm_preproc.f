      SUBROUTINE d4vm_preproc(bx0, bx1, by0, by1, bz0, bz1,
     .                        dx, dy, dt, bx, bxx, bxy,
     .                        by, byx, byy, bz, bzx, bzy, bzt,
     .                        xsize, ysize)

C.... Preprocessing for dave4vm
C....
C.... Written by Jacob Hageman, jacob.hageman@nasa.gov, (301)286-1803
C....
C.... Calculates mag data at midpoint
C....
C.... INPUTS:
C....   BX0 - Frame 0 of X mag data
C....   BX1 - Frame 1 of X mag data
C....   BY0 - Frame 0 of Y mag data
C....   BY1 - Frame 1 of Y mag data
C....   BZ0 - Frame 0 of Z mag data
C....   BZ1 - Frame 1 of Z mag data
C....   DX - X scale factor
C....   DY - Y scale factor
C....   DT - Time between frames (sec)
C....   XSIZE - Image dimension in x
C....   YSIZE - Image dimension in y
C....
C.... OUTPUTS:
C....   Bx - X component of mag data at midpoint
C....   Bxx - X derivative of Bx
C....   Bxy - Y derivative of Bx
C....   By - Y component of mag data at midpoint
C....   Byx - X derivative of By
C....   Byy - Y derivative of By
C....   Bz - Z component of mag data at midpoint
C....   Bzx - X derivative of Bz
C....   Bzy - Y derivative of Bz
C....   Bzt - Time derivative of Z mag data at midpoint
C....

      IMPLICIT NONE

      INTEGER*8 xsize, ysize
      REAL*4  bx0(xsize*ysize), bx1(xsize*ysize)
      REAL*4  by0(xsize*ysize), by1(xsize*ysize)
      REAL*4  bz0(xsize*ysize), bz1(xsize*ysize)
      REAL*8  dx, dy, dt
      REAL*4  bx(xsize*ysize), bxx(xsize*ysize), bxy(xsize*ysize)
      REAL*4  by(xsize*ysize), byx(xsize*ysize), byy(xsize*ysize)
      REAL*4  bz(xsize*ysize), bzx(xsize*ysize), bzy(xsize*ysize)
      REAL*4  bzt(xsize*ysize)

      INTEGER i

      DO i = 1, xsize*ysize
        bx(i) = (bx0(i) + bx1(i))/2.0
        by(i) = (by0(i) + by1(i))/2.0
        bz(i) = (bz0(i) + bz1(i))/2.0
        bzt(i) = (bz1(i) - bz0(i))/dt
      END DO

      CALL d4vm_derivs(Bx, Bxx, Bxy, xsize, ysize)
      CALL d4vm_derivs(By, Byx, Byy, xsize, ysize)
      CALL d4vm_derivs(Bz, Bzx, Bzy, xsize, ysize)

C.... Scale the derivatives
      Bxx = Bxx/dx
      Bxy = Bxy/dy
      Byx = Byx/dx
      Byy = Byy/dy
      Bzx = Bzx/dx
      Bzy = Bzy/dy

      END
