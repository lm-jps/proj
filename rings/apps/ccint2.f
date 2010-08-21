C------------------------------------------------------------------------
c
c
c	interpolation by cubic convolution in 2 dimensions with 16-bit
c	integer data.  follows r.g. keys in ieee trans. acoustics, speech,
c	and signal processing, vol. 29, 1981, pp. 1153-1160.
c
c	input:
c		f(nx,ny) real    	;input data picture
c		nx,ny	 integer	;x,y dimensions
c		x,y	 real		;position of the desired interpolation
c					;point in units of the x,y sample
c					;spacing.  0-indexed.  
c
c	output:
c		ccint2	real		;interpolated value
c
c	note: currently interpolation within one pixel of edge is not
c	implimented.  if x or y is in this range or off the picture, the
c	function value returned is zero.
c	note: to change to a different data type would only require changing
c	the declaration of f.
c
	real*4 function ccint2(f, nx, ny, x, y)
	integer nx, ny
	real*4 f(nx,ny)
	real*4 x, y
c
	real*4 ux(4), uy(4), sum
	integer ix, iy, ix1, iy1, i, j
c
	if ( x.lt.1. .or. x.gt.float(nx-2) .or.
     &	     y.lt.1. .or. y.gt.float(ny-2) ) then
		ccint2=0.
		return
	endif
c
	ix=int(x)
	call ccker( ux, x-float(ix) )
	iy=int(y)
	call ccker( uy, y-float(iy) )
c
	ix1= ix - 1
	iy1= iy - 1
	sum=0.
	do 100 i=1,4
	do 100 j=1,4
100	sum = sum + f( ix1+i, iy1+j ) * ux(i) * uy(j)
	ccint2=sum
	return
	end
c	calculate the interpolation kernel
c
	subroutine ccker(u, s)
	real*4 u(4), s, s2, s3
c       
	s2= s**2
	s3= s2 * s
	u(1) = -0.5 * (s3+s) + s2
	u(2) = 1.5*s3 - 2.5*s2 + 1.
	u(3) = -1.5*s3 + 2.*s2 + 0.5*s
	u(4) = 0.5*(s3-s2)
	return
	end


