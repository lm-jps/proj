      SUBROUTINE d4vm_derivs(array, dx, dy, xsize, ysize)

C....
C.... X and Y derivative calculation
C....
C.... Written by Jacob Hageman, jacob.hageman@nasa.gov, (301)286-1803
C....
C.... INPUTS:
C....   ARRAY - Array to perform convolution on
C....   XSIZE - Array x dimension
C....   YSIZE - Array y dimension
C....
C.... OUTPUTS:
C....   DX - Derivative in x direction
C....   DY - Derivative in y direction
C....
C.... Notes: Convolution kernel defined consistent with IDL

      IMPLICIT NONE

C.... INPUTS
      INTEGER*8 xsize, ysize
      REAL*4 array(xsize,ysize)

C.... OUTPUTS
      REAL*4 dx(xsize,ysize), dy(xsize,ysize)

C.... WORKING VARIABLES
      INTEGER i, j
      REAL*4 array_zeros(-2:(xsize+2),-2:(ysize+2))

C.... Copy into array with zero borders to avoid if's
      array_zeros = 0
      array_zeros(1:xsize,1:ysize) = array(1:xsize,1:ysize)

C.... Run through complete array
      DO j = 1, ysize
        DO i = 1, xsize

C....     Calculate dx
          dx(i,j) = -array_zeros(i-1,j) * 0.74038 +
     &              array_zeros(i-2,j) * 0.12019 +
     &              array_zeros(i+1,j) * 0.74038 -
     &              array_zeros(i+2,j) * 0.12019

C....     Calculate dy
          dy(i,j) = -array_zeros(i,j-1) * 0.74038 +
     &              array_zeros(i,j-2) * 0.12019 +
     &              array_zeros(i,j+1) * 0.74038 -
     &              array_zeros(i,j+2) * 0.12019

        END DO
      END DO

      END

