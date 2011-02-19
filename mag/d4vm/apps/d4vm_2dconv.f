      SUBROUTINE d4vm_2dconv(array, kernel, conv,
     &                       xsize, ysize, kern_size)

C....
C.... 2D convolution routine
C....
C.... Written by Jacob Hageman, jacob.hageman@nasa.gov, (301)286-1803
C....
C.... INPUTS:
C....   ARRAY - Array to perform convolution on
C....   KERN - Convolution kernel
C....   XSIZE - Array x dimension
C....   YSIZE - Array y dimension
C....   KERN_SIZE - Convolution kernel size
C....
C.... OUTPUTS:
C....   CONV - Convolution result
C....
C.... Notes: Convolution kernel defined consistent with IDL

      IMPLICIT NONE

C.... INPUTS
      INTEGER*8 xsize, ysize, kern_size(2)
      REAL*4 array(xsize,ysize)
      REAL*4 kernel(-kern_size(1)/2:kern_size(1)/2,
     &              -kern_size(2)/2:kern_size(2)/2)

C.... OUTPUTS
      REAL*4 conv(xsize,ysize)

C.... WORKING VARIABLES
      INTEGER i, j, m, n
      INTEGER hx
      INTEGER hy

C.... Initialize hx, hy, conv
      hx = kern_size(1)/2
      hy = kern_size(2)/2
      conv = 0

C.... Run through complete array
      DO j = hy+1, ysize-hy
        DO i = hx+1, xsize-hx

          DO m = -hx, hx
            DO n = -hy, hy
 
              conv(i,j) = conv(i,j) + array(i+m,j+n)*kernel(m,n)
 
            END DO
          END DO
        END DO
      END DO

      END
