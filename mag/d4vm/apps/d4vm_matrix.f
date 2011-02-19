      SUBROUTINE d4vm_matrix(A, Bx, Bxx, Bxy, By, Byx, Byy,
     &                       Bz, Bzx, Bzy, Bzt,
     &                       k_th, k_x, k_y, k_xx, k_yy, k_xy,
     &                       xsize, ysize, kern_size)

C.... Dave4vm processing routine
C....
C.... Written by Jacob Hageman, jacob.hageman@nasa.gov, (301)286-1803
C....
C.... Calculates plasma velocites from mdag data
C....
C.... INPUTS:
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
C....   K_th - Convolution kernel
C....   K_x - Convolution kernel times x
C....   K_y - Convolution kernel times y
C....   K_xx - Convolution kernel times x^2
C....   K_yy - Convolution kernel times y^2
C....   K_xy - Convolution kernel times x*y
C....   XSIZE - X dimension size
C....   YSIZE - Y dimension size
C....   KERN_SIZE - Kernel size
C....
C.... OUTPUTS:
C....   A - Calculated matrix
C....

      IMPLICIT NONE

C.... INPUTS & PARAMETERS
      INTEGER*8 xsize, ysize, kern_size(2)
      REAL*4  Bx(xsize*ysize), Bxx(xsize*ysize), Bxy(xsize*ysize)
      REAL*4  By(xsize*ysize), Byx(xsize*ysize), Byy(xsize*ysize)
      REAL*4  Bz(xsize*ysize), Bzx(xsize*ysize), Bzy(xsize*ysize)
      REAL*4  Bzt(xsize*ysize)

      REAL*4  k_th(kern_size(1), kern_size(2))
      REAL*4  k_x(kern_size(1), kern_size(2))
      REAL*4  k_y(kern_size(1), kern_size(2))
      REAL*4  k_xx(kern_size(1), kern_size(2))
      REAL*4  k_yy(kern_size(1), kern_size(2))
      REAL*4  k_xy(kern_size(1), kern_size(2))

C.... OUTPUTS
      REAL*4  A(10,10,xsize*ysize)

C.... WORKING VARIABLES
      INTEGER i, j
      REAL*4  temp(xsize*ysize)

C.... DAVE Elements
      REAL*4  GGx(xsize*ysize), GGy(xsize*ysize), GGt(xsize*ysize)
      REAL*4  GGxx(xsize*ysize), GGyy(xsize*ysize), GGxy(xsize*ysize)
      REAL*4  GGtx(xsize*ysize), GGty(xsize*ysize), GGtt(xsize*ysize)

      REAL*4  G(xsize*ysize)
      REAL*4  Gx(xsize*ysize), xGx(xsize*ysize), yGx(xsize*ysize)
      REAL*4  Gy(xsize*ysize), xGy(xsize*ysize), yGy(xsize*ysize)
      REAL*4  Ht(xsize*ysize)
      REAL*4  Gxx(xsize*ysize), Gyy(xsize*ysize), Gxy(xsize*ysize)
      REAL*4  Gtx(xsize*ysize), Gty(xsize*ysize)

      REAL*4  xGxx(xsize*ysize), xGyy(xsize*ysize), xGxy(xsize*ysize)
      REAL*4  xGtx(xsize*ysize), xGty(xsize*ysize)

      REAL*4  yGxx(xsize*ysize), yGyy(xsize*ysize), yGxy(xsize*ysize)
      REAL*4  yGtx(xsize*ysize), yGty(xsize*ysize)

      REAL*4  xxGxx(xsize*ysize), xxGxy(xsize*ysize), xxGyy(xsize*ysize)
      REAL*4  xyGxx(xsize*ysize), xyGxy(xsize*ysize), xyGyy(xsize*ysize)
      REAL*4  yyGxx(xsize*ysize), yyGxy(xsize*ysize), yyGyy(xsize*ysize)

      REAL*4  Gtt(xsize*ysize)

C.... Magnetogram Terms
      REAL*4  BxBx(xsize*ysize), ByBy(xsize*ysize), BxBy(xsize*ysize)
      REAL*4  BzBx(xsize*ysize), BzBy(xsize*ysize)

      REAL*4  BxBxx(xsize*ysize), BxByy(xsize*ysize)
      REAL*4  BxxBxx(xsize*ysize), ByyByy(xsize*ysize)
      REAL*4  BxxByy(xsize*ysize)
      REAL*4  ByBxx(xsize*ysize), ByByy(xsize*ysize)

      REAL*4  BzBxx(xsize*ysize), BzByy(xsize*ysize)
 
      REAL*4  BztBxx(xsize*ysize), BztByy(xsize*ysize)

      REAL*4  BzxBx(xsize*ysize), BzxBy(xsize*ysize)
      REAL*4  BzxBxx(xsize*ysize), BzxByy(xsize*ysize)

      REAL*4  BzyBx(xsize*ysize), BzyBy(xsize*ysize)
      REAL*4  BzyBxx(xsize*ysize), BzyByy(xsize*ysize)

      REAL*4  BztBx(xsize*ysize), BztBy(xsize*ysize)

      REAL*4  xBzxBx(xsize*ysize), xBzxBy(xsize*ysize)
      REAL*4  xBzyBx(xsize*ysize), xBzyBy(xsize*ysize)

      REAL*4  yBzxBx(xsize*ysize), yBzxBy(xsize*ysize)
      REAL*4  yBzyBx(xsize*ysize), yBzyBy(xsize*ysize)

      REAL*4  yBxBxx(xsize*ysize), yBxByy(xsize*ysize)
      REAL*4  yByBxx(xsize*ysize), yByByy(xsize*ysize)
      REAL*4  xByBxx(xsize*ysize), xByByy(xsize*ysize)

      REAL*4  xBzxBxx(xsize*ysize), xBzxByy(xsize*ysize)
      REAL*4  yBzxBxx(xsize*ysize), yBzxByy(xsize*ysize)

      REAL*4  xBxxBxx(xsize*ysize), xBxxByy(xsize*ysize)
      REAL*4  xByyByy(xsize*ysize), yBxxBxx(xsize*ysize)
      REAL*4  yBxxByy(xsize*ysize), yByyByy(xsize*ysize)

      REAL*4  xBxBxx(xsize*ysize), xBxByy(xsize*ysize)
      REAL*4  xBzBxx(xsize*ysize), xBzByy(xsize*ysize)

      REAL*4  xBztBxx(xsize*ysize), xBztByy(xsize*ysize)
      REAL*4  yBztBxx(xsize*ysize), yBztByy(xsize*ysize)

      REAL*4  xyBxxBxx(xsize*ysize), xyBxxByy(xsize*ysize)
      REAL*4  xyByyByy(xsize*ysize)

      REAL*4  xyBzxBxx(xsize*ysize), xyBzxByy(xsize*ysize)
      REAL*4  xyBzyBxx(xsize*ysize), xyBzyByy(xsize*ysize)

      REAL*4  yBzBxx(xsize*ysize), yBzByy(xsize*ysize)

      REAL*4  xBzyBxx(xsize*ysize), xBzyByy(xsize*ysize)
      REAL*4  yBzyBxx(xsize*ysize), yBzyByy(xsize*ysize)

      REAL*4  xxBxxBxx(xsize*ysize), xxBxxByy(xsize*ysize)
      REAL*4  xxByyByy(xsize*ysize)

      REAL*4  xxBzxBxx(xsize*ysize), xxBzxByy(xsize*ysize)
      REAL*4  xxBzyBxx(xsize*ysize), xxBzyByy(xsize*ysize)

      REAL*4  yyBxxBxx(xsize*ysize), yyBxxByy(xsize*ysize)
      REAL*4  yyByyByy(xsize*ysize)

      REAL*4  yyBzxBxx(xsize*ysize), yyBzxByy(xsize*ysize)
      REAL*4  yyBzyBxx(xsize*ysize), yyBzyByy(xsize*ysize)

C.... 2D convolutions
      CALL d4vm_2dconv(Bz*Bz, k_th, G, xsize, ysize, kern_size)

      GGx = Bz*Bzx
      CALL d4vm_2dconv(GGx, k_th, Gx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGx, k_x, xGx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGx, k_y, yGx, xsize, ysize, kern_size)

      GGy = Bz*Bzy
      CALL d4vm_2dconv(GGy, k_th, Gy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGy, k_x, xGy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGy, k_y, yGy, xsize, ysize, kern_size)

      GGt = Bzt*Bz
      CALL d4vm_2dconv(GGt, k_th, Ht, xsize, ysize, kern_size)

      GGxx = Bzx*Bzx
      CALL d4vm_2dconv(GGxx, k_th, Gxx, xsize, ysize, kern_size)

      GGyy = Bzy*Bzy
      CALL d4vm_2dconv(GGyy, k_th, Gyy, xsize, ysize, kern_size)

      GGxy = Bzx*Bzy
      CALL d4vm_2dconv(GGxy, k_th, Gxy, xsize, ysize, kern_size)

      GGtx = Bzt*Bzx
      CALL d4vm_2dconv(GGtx, k_th, Gtx, xsize, ysize, kern_size)

      GGty = Bzt*Bzy
      CALL d4vm_2dconv(GGty, k_th, Gty, xsize, ysize, kern_size)

C.....................

      CALL d4vm_2dconv(GGxx, k_x, xGxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGyy, k_x, xGyy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGxy, k_x, xGxy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGtx, k_x, xGtx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGty, k_x, xGty, xsize, ysize, kern_size)

      CALL d4vm_2dconv(GGxx, k_y, yGxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGyy, k_y, yGyy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGxy, k_y, yGxy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGtx, k_y, yGtx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGty, k_y, yGty, xsize, ysize, kern_size)

      CALL d4vm_2dconv(GGxx, k_xx, xxGxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGxy, k_xx, xxGxy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGyy, k_xx, xxGyy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(GGxx, k_xy, xyGxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGyy, k_xy, xyGyy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGxy, k_xy, xyGxy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(GGxx, k_yy, yyGxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGxy, k_yy, yyGxy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(GGyy, k_yy, yyGyy, xsize, ysize, kern_size)

      GGtt=Bzt*Bzt
      CALL d4vm_2dconv(GGtt, k_th, Gtt, xsize, ysize, kern_size)

C.... Extra vector magnetogram terms

      CALL d4vm_2dconv(Bx*Bx, k_th, BxBx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(By*By, k_th, ByBy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bx*By, k_th, BxBy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bz*Bx, k_th, BzBx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bz*By, k_th, BzBy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bx*Bxx, k_th, BxBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bx*Byy, k_th, BxByy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bxx*Bxx, k_th, BxxBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Byy*Byy, k_th, ByyByy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bxx*Byy, k_th, BxxByy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(By*Bxx, k_th, ByBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(By*Byy, k_th, ByByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bz*Bxx, k_th, BzBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bz*Byy, k_th, BzByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bzt*Bxx, k_th, BztBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzt*Byy, k_th, BztByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bzx*Bx, k_th, BzxBx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzx*By, k_th, BzxBy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzx*Bxx, k_th, BzxBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzx*Byy, k_th, BzxByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bzy*Bx, k_th, BzyBx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzy*By, k_th, BzyBy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzy*Bxx, k_th, BzyBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzy*Byy, k_th, BzyByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bzt*Bx, k_th, BztBx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzt*By, k_th, BztBy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bzx*Bx, k_x, xBzxBx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzx*By, k_x, xBzxBy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzy*Bx, k_x, xBzyBx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzy*By, k_x, xBzyBy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bzy*Bx, k_y, yBzyBx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzy*By, k_y, yBzyBy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bzx*Bx, k_y, yBzxBx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzx*By, k_y, yBzxBy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bx*Bxx, k_y, yBxBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bx*Byy, k_y, yBxByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(By*Bxx, k_y, yByBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(By*Byy, k_y, yByByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(By*Bxx, k_x, xByBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(By*Byy, k_x, xByByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bzx*Bxx, k_x, xBzxBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzx*Byy, k_x, xBzxByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bzx*Bxx, k_y, yBzxBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzx*Byy, k_y, yBzxByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bxx*Bxx, k_x, xBxxBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bxx*Byy, k_x, xBxxByy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Byy*Byy, k_x, xByyByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bxx*Bxx, k_y, yBxxBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bxx*Byy, k_y, yBxxByy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Byy*Byy, k_y, yByyByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bx*Bxx, k_x, xBxBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bx*Byy, k_x, xBxByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bz*Bxx, k_x, xBzBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bz*Byy, k_x, xBzByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bzt*Bxx, k_x, xBztBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzt*Byy, k_x, xBztByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bzt*Bxx, k_y, yBztBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzt*Byy, k_y, yBztByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bxx*Bxx, k_xy, xyBxxBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bxx*Byy, k_xy, xyBxxByy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Byy*Byy, k_xy, xyByyByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bzx*Bxx, k_xy, xyBzxBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzx*Byy, k_xy, xyBzxByy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzy*Bxx, k_xy, xyBzyBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzy*Byy, k_xy, xyBzyByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bz*Bxx, k_y, yBzBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bz*Byy, k_y, yBzByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bzy*Bxx, k_x, xBzyBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzy*Byy, k_x, xBzyByy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzy*Bxx, k_y, yBzyBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzy*Byy, k_y, yBzyByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bxx*Bxx, k_xx, xxBxxBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bxx*Byy, k_xx, xxBxxByy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Byy*Byy, k_xx, xxByyByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bzx*Bxx, k_xx, xxBzxBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzy*Bxx, k_xx, xxBzyBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzx*Byy, k_xx, xxBzxByy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzy*Byy, k_xx, xxBzyByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bxx*Bxx, k_yy, yyBxxBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bxx*Byy, k_yy, yyBxxByy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Byy*Byy, k_yy, yyByyByy, xsize, ysize, kern_size)

      CALL d4vm_2dconv(Bzy*Bxx, k_yy, yyBzyBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzy*Byy, k_yy, yyBzyByy, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzx*Bxx, k_yy, yyBzxBxx, xsize, ysize, kern_size)
      CALL d4vm_2dconv(Bzx*Byy, k_yy, yyBzxByy, xsize, ysize, kern_size)

C.... Form Matrix

      A(1,1,:) = Gxx
      A(1,2,:) = Gxy
      A(1,3,:) = Gx + xGxx
      A(1,4,:) = Gx + yGxy
      A(1,5,:) = yGxx
      A(1,6,:) = xGxy
      A(1,7,:) = -BzxBxx - BzxByy
      A(1,8,:) = -BzxBx - xBzxBxx - xBzxByy
      A(1,9,:) = -BzxBy - yBzxBxx - yBzxByy
      A(1,10,:) = Gtx

      A(2,2,:) = Gyy
      A(2,3,:) = Gy + xGxy
      A(2,4,:) = Gy + yGyy
      A(2,5,:) = yGxy
      A(2,6,:) = xGyy
      A(2,7,:) = -BzyBxx - BzyByy
      A(2,8,:) = -BzyBx - xBzyBxx - xBzyByy
      A(2,9,:) = -BzyBy - yBzyBxx - yBzyByy
      A(2,10,:) = Gty

      A(3,3,:) = G + 2*xGx + xxGxx
      A(3,4,:) = G + xGx + xyGxy + yGy
      A(3,5,:) = xyGxx + yGx
      A(3,6,:) = xGy + xxGxy
      A(3,7,:) = -BzBxx - BzByy - xBzxBxx - xBzxByy
      A(3,8,:) = -BzBx - xBzBxx - xBzByy - xBzxBx - xxBzxBxx - xxBzxByy
      A(3,9,:) = -BzBy - xBzxBy - xyBzxBxx - xyBzxByy - yBzBxx - yBzByy
      A(3,10,:) = Ht + xGtx

      A(4,4,:) = G + 2*yGy + yyGyy
      A(4,5,:) = yGx + yyGxy
      A(4,6,:) = xGy + xyGyy
      A(4,7,:) = -BzBxx - BzByy - yBzyBxx - yBzyByy
      A(4,8,:) = -BzBx - xBzBxx - xBzByy - xyBzyBxx - xyBzyByy - yBzyBx
      A(4,9,:) = -BzBy - yBzBxx - yBzByy - yBzyBy - yyBzyBxx - yyBzyByy
      A(4,10,:) = Ht + yGty

      A(5,5,:) = yyGxx
      A(5,6,:) = xyGxy
      A(5,7,:) = -yBzxBxx - yBzxByy
      A(5,8,:) = -xyBzxBxx - xyBzxByy - yBzxBx
      A(5,9,:) = -yBzxBy - yyBzxBxx - yyBzxByy
      A(5,10,:) = yGtx

      A(6,6,:) = xxGyy
      A(6,7,:) = -xBzyBxx - xBzyByy
      A(6,8,:) = -xBzyBx - xxBzyBxx - xxBzyByy
      A(6,9,:) = -xBzyBy - xyBzyBxx - xyBzyByy
      A(6,10,:) = xGty

      A(7,7,:) = BxxBxx + 2*BxxByy + ByyByy
      A(7,8,:) = BxBxx + BxByy + xBxxBxx + 2*xBxxByy + xByyByy
      A(7,9,:) = ByBxx + ByByy + yBxxBxx + 2*yBxxByy + yByyByy
      A(7,10,:) = -BztBxx - BztByy

      A(8,8,:) = BxBx + 2*xBxBxx + 2*xBxByy + xxBxxBxx + 
     &           2*xxBxxByy + xxByyByy
      A(8,9,:) = BxBy + xByBxx + xByByy + xyBxxBxx +
     &           2*xyBxxByy + xyByyByy + yBxBxx + yBxByy
      A(8,10,:) = -BztBx - xBztBxx - xBztByy

      A(9,9,:) = ByBy + 2*yByBxx + 2*yByByy + yyBxxBxx +
     &           2*yyBxxByy + yyByyByy
      A(9,10,:) = -BztBy - yBztBxx - yBztByy
     
      A(10,10,:) = Gtt

      DO i=2, 10
        DO j = 1, i-1
          A(i,j,:) = A(j,i,:)
        ENDDO
      ENDDO

      END
