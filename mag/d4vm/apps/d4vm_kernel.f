      SUBROUTINE d4vm_kernel(dx, dy, kern_size,
     &                       k_th, k_x, k_y, k_xx, k_yy, k_xy)

C.... Generate convolution kernels
C....
C.... Written by Jacob Hageman, jacob.hageman@nasa.gov, (301)286-1803
C....
C.... INPUTS:
C....   dx - image x scaling
C....   dy - image y scaling
C....   kern_size - Desired kernel size
C....
C.... OUTPUTS:
C....   k_th - Base kernel
C....   k_x - Base kernel times x
C....   k_y - Base kernel times y
C....   k_xx - Base kernel times x^2
C....   k_yy - Base kernel times y^2
C....   k_xy - Base kernel times x*y
C....

      IMPLICIT NONE

      INTEGER*8 kern_size(2)
      REAL*8  dx, dy
      REAL*4  k_th(kern_size(1), kern_size(2))
      REAL*4  k_x(kern_size(1), kern_size(2))
      REAL*4  k_y(kern_size(1), kern_size(2))
      REAL*4  k_xx(kern_size(1), kern_size(2))
      REAL*4  k_yy(kern_size(1), kern_size(2))
      REAL*4  k_xy(kern_size(1), kern_size(2))
      REAL*4  x(kern_size(1), kern_size(2))
      REAL*4  y(kern_size(1), kern_size(2))

      INTEGER i, j

C.... Construct kernels
      k_th = 1.0/(kern_size(1)*kern_size(2))

      DO i = 1, kern_size(1)
        DO j = 1, kern_size(2)
          x(i,j) = (i-(kern_size(1)/2)-1)*real(dx)
          y(i,j) = (j-(kern_size(2)/2)-1)*real(dy)
        ENDDO
      ENDDO

      k_x = x * k_th
      k_y = y * k_th
      k_xx = x * k_x
      k_yy = y * k_y
      k_xy = x * k_y

      END
