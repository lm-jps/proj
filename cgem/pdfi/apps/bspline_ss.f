      SUBROUTINE bspline_ss(fdata, nx, ny, xnew, nxinterp, ynew, 
     1 nyinterp, finterp, degree)
c
c+
c - - Purpose: To interpolate the array fdata to a new grid defined by
c             xnew,ynew locations (in terms of original index 
c             ranges [1:nx],[1:ny]) 
c
c - - Usage:  call bspline_ss(fdata,nx,ny,xnew,nxinterp,ynew,nyinterp,
c             finterp, degree)
c
c - - Input:  fdata(nx,ny) - real*8 array of data to be used for interpolation.
c             [units arbitrary]
c - - Input:  nx,ny - dimensions of input array fdata
c - - Input:  xnew(nxinterp) - real*8 array of desired x values in [1,nx] range
c - - Input:  nxinterp - dimension of xnew
c - - Input:  ynew(nyinterp) - real*8 array of desired y values in [1,ny] range
c - - Input:  nyinterp - dimension of ynew
c - - Output: finterp(nxinterp,nyinterp) - real*8 array of interpolated values  
c             [units arbitrary]
c - - Input:  degree - integer value of the degree of the bspline (3 <= 9).
c             degree must also be an odd integer: allowed values 3,5,7,9.
c-
c   Origins:  This routine was written in Fortran by Dave Bercik.
c   It is based on software in C, published by Philippe Thevenaz, 
c   http://bigwww.epfl.ch/thevenaz/interpolation/  
c   originally described in IEEE TRANSACTIONS ON MEDICAL IMAGING, 
c   VOL. 19, NO. 7, JULY 2000.  It was subsequently edited for compatibility 
c   with the PDFI_SS library by George Fisher. 
c  
c   PDFI_SS Electric Field Inversion Software
c   http://cgem.ssl.berkeley.edu/cgi-bin/cgem/PDFI_SS/index
c   Copyright (C) 2015,2016 University of California
c  
c   This software is based on the concepts described in Kazachenko et al. 
c   (2014, ApJ 795, 17).  It also extends those techniques to 
c   spherical coordinates, and uses a staggered, rather than a centered grid.
c   If you use the software in a scientific publication, 
c   the authors would appreciate a citation to this paper and any future papers 
c   describing updates to the methods.
c  
c   This is free software; you can redistribute it and/or
c   modify it under the terms of the GNU Lesser General Public
c   License as published by the Free Software Foundation;
c   either version 2.1 of the License, or (at your option) any later version.
c  
c   This software is distributed in the hope that it will be useful,
c   but WITHOUT ANY WARRANTY; without even the implied warranty of
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c   See the GNU Lesser General Public License for more details.
c  
c   To view the GNU Lesser General Public License visit
c   http://www.gnu.org/copyleft/lesser.html
c   or write to the Free Software Foundation, Inc.,
c   59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
c
c
      IMPLICIT NONE
c
      INTEGER, PARAMETER :: r8=KIND(1.d0)
      INTEGER, PARAMETER :: r4=KIND(1.0)
      INTEGER, PARAMETER :: i4=KIND(1)

c     Dummy variable declarations:
      INTEGER(KIND=i4) :: nx,ny,nxinterp,nyinterp
      REAL(KIND=r8) :: fdata(nx,ny)
      REAL(KIND=r8) :: xnew(nxinterp)
      REAL(KIND=r8) :: ynew(nyinterp)
      REAL(KIND=r8) :: finterp(nxinterp,nyinterp)
      INTEGER(KIND=i4) :: degree

c     Local variable declarations
      INTEGER(KIND=i4) :: i, j, k, n, npoles, ix, jy
c     INTEGER(KIND=ir) :: istat
      INTEGER(KIND=i4) :: width2, height2
      INTEGER(KIND=i4), DIMENSION(10) :: xindex, yindex
      REAL(KIND=r8) :: lambda, w, w2, w4, t, t0, t1
      REAL(KIND=r8), DIMENSION(10) :: xweight, yweight
      REAL(KIND=r8), ALLOCATABLE :: c(:)
      REAL(KIND=r8), ALLOCATABLE :: poles(:)
      REAL(KIND=r8), ALLOCATABLE :: coeffs (:,:)
c - - Function declarations:
      REAL(KIND=r8) :: causal_init_ss
      REAL(KIND=r8) :: anticausal_init_ss

c     Determine poles for given spline degree
      npoles = degree/2
c - - check for degree being even, and exit if so; very poor results.
      if(npoles*2 .eq. degree) then
        write(6,*) 'bspline_ss:  degree is even, exiting. degree = ',
     1  degree
        stop
      endif
      ALLOCATE(poles(npoles))
      SELECT CASE (degree)
        CASE (2)
          poles(1) = SQRT(8.0_r8) - 3.0_r8 
        CASE (3)
          poles(1) = SQRT(3.0_r8) - 2.0_r8 
        CASE (4)
          poles(1) = SQRT(664.0_r8 - SQRT(438976.0_r8)) + SQRT(304.0_r8)
     1    - 19.0_r8
          poles(2) = SQRT(664.0_r8 + SQRT(438976.0_r8)) - SQRT(304.0_r8)
     1             - 19.0_r8
        CASE (5)
          poles(1) = SQRT(135.0_r8/2.0_r8 - SQRT(17745.0_r8/4.0_r8)) 
     1               + SQRT(105.0_r8/4.0_r8) - 13.0_r8/2.0_r8
          poles(2) = SQRT(135.0_r8/2.0_r8 + SQRT(17745.0_r8/4.0_r8))
     1               - SQRT(105.0_r8/4.0_r8) - 13.0_r8/2.0_r8
        CASE (6)
          poles(1) = 
     1   -0.48829458930304475513011803888378906211227916123938_r8
          poles(2) = 
     1   -0.081679271076237512597937765737059080653379610398148_r8
          poles(3) = 
     1    -0.0014141518083258177510872439765585925278641690553467_r8
        CASE (7)
          poles(1) = 
     1    -0.53528043079643816554240378168164607183392315234269_r8
          poles(2) = 
     1    -0.12255461519232669051527226435935734360548654942730_r8
          poles(3) = 
     1    -0.0091486948096082769285930216516478534156925639545994_r8
        CASE (8)
          poles(1) = 
     1    -0.57468690924876543053013930412874542429066157804125_r8
          poles(2) = 
     1    -0.16303526929728093524055189686073705223476814550830_r8
          poles(3) = 
     1    -0.023632294694844850023403919296361320612665920854629_r8
          poles(4) = 
     1    -0.00015382131064169091173935253018402160762964054070043_r8
        CASE (9)
          poles(1) = 
     1    -0.60799738916862577900772082395428976943963471853991_r8
          poles(2) = 
     1    -0.20175052019315323879606468505597043468089886575747_r8
          poles(3) = 
     1    -0.043222608540481752133321142979429688265852380231497_r8
          poles(4) = 
     1    -0.0021213069031808184203048965578486234220548560988624_r8
        CASE DEFAULT
          write(6,*) 'bspline_ss: exiting, invalid value of degree = ',
     1     degree
          STOP
      END SELECT

c     Find coefficients of bspline expansion.  Treat as separable.
c       - Mirror boundary conditions
      width2 = 2*(nx-1)
      height2 = 2*(ny-1)
      ALLOCATE(coeffs(nx, ny))
      coeffs = fdata(1:nx,1:ny)
c     Calculate overall gain for one dimension
      lambda = 1.0_r8
      DO k=1,npoles
        lambda = (1.0_r8 - poles(k))*(1.0_r8 - 1.0_r8/poles(k))*lambda
      END DO

c     Filter along x for each y
      ALLOCATE(c(nx))
      coeffs = lambda*coeffs
      DO j=1,ny
        c = coeffs(:, j)
        DO k=1,npoles
c         Causal initiation:
          c(1) = causal_init_ss(c, nx, poles(k))
          DO i=2,nx
            c(i) = c(i) + poles(k)*c(i-1)
          END DO
c         Anticausal initiation:
          c(nx) = anticausal_init_ss(c, nx, poles(k))
          DO i=nx-1,1,-1
            c(i) = poles(k)*(c(i+1) - c(i))
          END DO
        END DO
        coeffs(:, j) = c
      END DO
      DEALLOCATE(c)

c     Filter along y for each x
      ALLOCATE(c(ny))
      coeffs = lambda*coeffs
      DO i=1,nx
        c(1:ny) = coeffs(i, 1:ny)
        DO k=1,npoles
c         Causal initiation
          c(1) = causal_init_ss(c, ny, poles(k))
          DO j=2,ny
            c(j) = c(j) + poles(k)*c(j-1)
          END DO
c         Anticausal initiation
          c(ny) = anticausal_init_ss(c, ny, poles(k))
          DO j=ny-1,1,-1
            c(j) = poles(k)*(c(j+1) - c(j))
          END DO
        END DO
        coeffs(i, :) = c
      END DO
      DEALLOCATE(c)

c     Interpolate to new x,y points
      finterp = 0.0_r8
      DO j = 1,nyinterp
        DO i = 1,nxinterp
          xindex = 0
          yindex = 0
          xweight = 0.0_r8
          yweight = 0.0_r8
          IF (MOD(degree, 2) /= 0) THEN
            ix = INT(FLOOR(xnew(i)), KIND=i4) - degree/2
            jy = INT(FLOOR(ynew(j)), KIND=i4) - degree/2
            xindex(1:degree+1) = (/ (ix+k, k=0,degree) /)
            yindex(1:degree+1) = (/ (jy+k, k=0,degree) /)
          ELSE
            ix = INT(FLOOR(xnew(i) + 0.5_r8), KIND=i4) - degree/2
            jy = INT(FLOOR(ynew(j) + 0.5_r8), KIND=i4) - degree/2
            xindex(1:degree+1) = (/ (ix+k, k=0,degree) /)
            yindex(1:degree+1) = (/ (jy+k, k=0,degree) /)
          END IF

c         Calculate interpolation weights:
          SELECT CASE (degree)
            CASE (2)
              w = xnew(i) - REAL(xindex(2), KIND=r8)
              xweight(2) = 3.0_r8/4.0_r8 - w*w
              xweight(3) = 1.0_r8/2.0_r8*(w - xweight(2) + 1.0_r8)
              xweight(1) = 1.0_r8 - xweight(2) - xweight(3) 
              w = ynew(i) - REAL(yindex(2), KIND=r8)
              yweight(2) = 3.0_r8/4.0_r8 - w*w
              yweight(3) = 1.0_r8/2.0_r8*(w - yweight(2) + 1.0_r8)
              yweight(1) = 1.0_r8 - yweight(2) - yweight(3) 
            CASE (3)
              w = xnew(i) - REAL(xindex(2), KIND=r8)
              xweight(4) = 1.0_r8/6.0_r8*w*w*w
              xweight(1) = 1.0_r8/6.0_r8 + 1.0_r8/2.0_r8*w*(w - 1.0_r8)
     1                   - xweight(4)
              xweight(3) = w + xweight(1) - 2.0_r8*xweight(4)
              xweight(2) = 1.0_r8 - xweight(1) - xweight(3) - xweight(4)
              w = ynew(j) - REAL(yindex(2), KIND=r8)
              yweight(4) = 1.0_r8/6.0_r8*w*w*w
              yweight(1) = 1.0_r8/6.0_r8 + 1.0_r8/2.0_r8*w*(w - 1.0_r8) 
     1                     - yweight(4)
              yweight(3) = w + yweight(1) - 2.0_r8*yweight(4)
              yweight(2) = 1.0_r8 - yweight(1) - yweight(3) - yweight(4)
            CASE (4)
              w = xnew(i) - REAL(xindex(3), KIND=r8)       
              w2 = w*w
              t = 1.0_r8/6.0_r8*w2
              xweight(1) = 1.0_r8/2.0_r8 - w
              xweight(1) = xweight(1)*xweight(1)
              xweight(1) = xweight(1)*1.0_r8/24.0_r8*xweight(1)
              t0 = w*(t - 11.0_r8 / 24.0_r8)
              t1 = 19.0_r8/96.0_r8 + w2*(1.0_r8/4.0_r8 - t)
              xweight(2) = t1 + t0
              xweight(4) = t1 - t0
              xweight(5) = xweight(1) + t0 + 1.0_r8/2.0_r8*w
              xweight(3) = 1.0_r8 - xweight(1) - xweight(2) - xweight(4)
     1                   - xweight(5)
              w = ynew(i) - REAL(yindex(3), KIND=r8)       
              w2 = w*w
              t = 1.0_r8/6.0_r8*w2
              yweight(1) = 1.0_r8/2.0_r8 - w
              yweight(1) = yweight(1)*yweight(1)
              yweight(1) = yweight(1)*1.0_r8/24.0_r8*yweight(1)
              t0 = w*(t - 11.0_r8 / 24.0_r8)
              t1 = 19.0_r8/96.0_r8 + w2*(1.0_r8/4.0_r8 - t)
              yweight(2) = t1 + t0
              yweight(4) = t1 - t0
              yweight(5) = yweight(1) + t0 + 1.0_r8/2.0_r8*w
              yweight(3) = 1.0_r8 - yweight(1) - yweight(2) - yweight(4)
     1                   - yweight(5)
            CASE (5)
              w = xnew(i) - REAL(xindex(3), KIND=r8)
              w2 = w*w
              xweight(6) = 1.0_r8/120.0_r8*w*w2*w2
              w2 = w2 - w
              w4 = w2*w2
              w = w - 1.0_r8/2.0_r8
              t = w2*(w2 - 3.0_r8)
              xweight(1) = 1.0_r8/24.0_r8*(1.0_r8/5.0_r8 + w2 + w4)
     1                   - xweight(6)
              t0 = 1.0_r8/24.0_r8*(w2*(w2 - 5.0_r8) + 46.0_r8/5.0_r8)
              t1 = -1.0_r8/12.0_r8*w*(t + 4.0_r8)
              xweight(3) = t0 + t1
              xweight(4) = t0 - t1
              t0 = 1.0_r8/16.0_r8*(9.0_r8/5.0_r8 - t)
              t1 = 1.0_r8/24.0_r8*w*(w4 - w2 - 5.0_r8)
              xweight(2) = t0 + t1
              xweight(5) = t0 - t1
              w = ynew(j) - REAL(yindex(3), KIND=r8)
              w2 = w*w
              yweight(6) = 1.0_r8/120.0_r8*w*w2*w2
              w2 = w2 - w
              w4 = w2*w2
              w = w - 1.0_r8/2.0_r8
              t = w2*(w2 - 3.0_r8)
              yweight(1) = 1.0_r8/24.0_r8*(1.0_r8/5.0_r8 + w2 + w4)
     1                   - yweight(6)
              t0 = 1.0_r8/24.0_r8*(w2*(w2 - 5.0_r8) + 46.0_r8/5.0_r8)
              t1 = -1.0_r8/12.0_r8*w*(t + 4.0_r8)
              yweight(3) = t0 + t1
              yweight(4) = t0 - t1
              t0 = 1.0_r8/16.0_r8*(9.0_r8/5.0_r8 - t)
              t1 = 1.0_r8/24.0_r8*w*(w4 - w2 - 5.0_r8)
              yweight(2) = t0 + t1
              yweight(5) = t0 - t1
            CASE (6)  
              w = xnew(i) - REAL(xindex(4), KIND=r8)
              xweight(1) = 1.0_r8/2.0_r8 - w
              xweight(1) = xweight(1)*xweight(1)*xweight(1)
              xweight(1) = xweight(1)*xweight(1)/720.0_r8 
              xweight(2) = (361.0_r8/192.0_r8 - w*(59.0_r8/8.0_r8 
     1                    + w*(-185.0_r8/16.0_r8 + w*(25.0_r8/3.0_r8 
     2                   + w*(-5.0_r8/ 2.0_r8 + w)*(1.0_r8/2.0_r8 
     3                   + w)))))/120.0_r8 
              xweight(3) = (10543.0_r8/960.0_r8 + w*(-289.0_r8/16.0_r8 
     1                   + w*(79.0_r8/16.0_r8 + w*(43.0_r8/6.0_r8 
     2                   + w*(-17.0_r8/4.0_r8 + w*(-1.0_r8 + w)))))) 
     3                   /48.0_r8
              w2 = w*w
              xweight(4) = (5887.0_r8/320.0_r8 - w2*(231.0_r8/16.0_r8 
     1                     - w2*(21.0_r8/4.0_r8 - w2)))/36.0_r8
              xweight(5) = (10543.0_r8/960.0_r8 + w*(289.0_r8/16.0_r8 
     1                     + w*(79.0_r8/16.0_r8 + w*(-43.0_r8/6.0_r8 
     2                     + w*(-17.0_r8/4.0_r8 + w*(1.0_r8 + w)))))) 
     3                     /48.0_r8
              xweight(7) = 1.0_r8/2.0_r8 + w
              xweight(7) = xweight(7)*xweight(7)*xweight(7)
              xweight(7) = xweight(7)*xweight(7)/720.0_r8
              xweight(6) = 1.0_r8 - xweight(1) - xweight(2) - xweight(3)
     1                   - xweight(4) - xweight(5) - xweight(7)
              w = ynew(i) - REAL(yindex(4), KIND=r8)
              yweight(1) = 1.0_r8/2.0_r8 - w
              yweight(1) = yweight(1)*yweight(1)*yweight(1)
              yweight(1) = yweight(1)*yweight(1)/720.0_r8 
              yweight(2) = (361.0_r8/192.0_r8 - w*(59.0_r8/8.0_r8 
     1                   + w*(-185.0_r8/16.0_r8 + w*(25.0_r8/3.0_r8 
     2                   + w*(-5.0_r8/ 2.0_r8 + w)*(1.0_r8/2.0_r8 
     3                   + w)))))/120.0_r8 
              yweight(3) = (10543.0_r8/960.0_r8 + w*(-289.0_r8/16.0_r8 
     1                     + w*(79.0_r8/16.0_r8 + w*(43.0_r8/6.0_r8 
     2                     + w*(-17.0_r8/4.0_r8 + w*(-1.0_r8 + w)))))) 
     3                     /48.0_r8
              w2 = w*w
                                                                       
              yweight(4) = (5887.0_r8/320.0_r8 - w2*(231.0_r8/16.0_r8 
     1                     - w2*(21.0_r8/4.0_r8 - w2)))/36.0_r8
              yweight(5) = (10543.0_r8/960.0_r8 + w*(289.0_r8/16.0_r8 
     1                     + w*(79.0_r8/16.0_r8 + w*(-43.0_r8/6.0_r8 
     2                     + w*(-17.0_r8/4.0_r8 + w*(1.0_r8 + w)))))) 
     3                     /48.0_r8
              yweight(7) = 1.0_r8/2.0_r8 + w
              yweight(7) = yweight(7)*yweight(7)*yweight(7)
              yweight(7) = yweight(7)*yweight(7)/720.0_r8
              yweight(6) = 1.0_r8 - yweight(1) - yweight(2) - yweight(3)
     1                     - yweight(4) - yweight(5) - yweight(7)
            CASE (7)
              w = xnew(i) - REAL(xindex(4), KIND=r8)
              xweight(1) = 1.0_r8 - w
              xweight(1) = xweight(1)*xweight(1)
              xweight(1) = xweight(1)*xweight(1)*xweight(1)
              xweight(1) = xweight(1)*(1.0_r8 - w)/5040.0_r8
              w2 = w*w
              xweight(2) = (120.0_r8/7.0_r8 + w*(-56.0_r8 + w*(72.0_r8 
     1                     + w*(-40.0_r8 + w2*(12.0_r8 + w*(-6.0_r8 
     2                     + w))))))/720.0
              xweight(3) = (397.0_r8/7.0_r8 - w*(245.0_r8/3.0_r8 
     1                     + w*(-15.0_r8 + w*(-95.0_r8/3.0_r8 
     2                     + w*(15.0_r8 + w*(5.0_r8 + w*(-5.0_r8 
     3                     + w)))))))/240.0_r8
              xweight(4) = (2416.0_r8/35.0_r8 + w2*(-48.0_r8 
     1                     + w2*(16.0_r8 + w2*(-4.0_r8 + w))))/144.0_r8
              xweight(5) = (1191.0_r8/35.0_r8 - w*(-49.0_r8 +w*(-9.0_r8 
     1                     + w*(19.0_r8 + w*(-3.0_r8 + w)*(-3.0_r8 
     2                     + w2)))))/144.0_r8
              xweight(6) = (40.0_r8/7.0_r8 + w*(56.0_r8/3.0_r8 
     1                     + w*(24.0_r8 + w*(40.0_r8/3.0_r8 
     2                     + w2*(-4.0_r8 + w*(-2.0_r8 + w))))))/240.0_r8
              xweight(8) = w2
              xweight(8) = xweight(8)*xweight(8)*xweight(8)
              xweight(8) = xweight(8)*w/5040.0_r8
              xweight(7) = 1.0_r8 - xweight(1) - xweight(2) 
     1                     - xweight(3) - xweight(4) - xweight(5) 
     2                     - xweight(6) - xweight(8)
              w = ynew(j) - REAL(yindex(4), KIND=r8)
              yweight(1) = 1.0_r8 - w
              yweight(1) = yweight(1)*yweight(1)
              yweight(1) = yweight(1)*yweight(1)*yweight(1)
              yweight(1) = yweight(1)*(1.0_r8 - w)/5040.0_r8
              w2 = w*w
              yweight(2) = (120.0_r8/7.0_r8 + w*(-56.0_r8 + w*(72.0_r8 
     1                     + w*(-40.0_r8 + w2*(12.0_r8 + w*(-6.0_r8 
     2                     + w))))))/720.0_r8
              yweight(3) = (397.0_r8/7.0_r8 - w*(245.0_r8/3.0_r8 
     1                     + w*(-15.0_r8 + w*(-95.0_r8/3.0_r8 
     2                     + w*(15.0_r8 + w*(5.0_r8 + w*(-5.0_r8 
     3                     + w)))))))/240.0_r8
              yweight(4) = (2416.0_r8/35.0_r8 + w2*(-48.0_r8 
     1                     + w2*(16.0_r8 + w2*(-4.0_r8 + w))))/144.0_r8
              yweight(5) = (1191.0_r8/35.0_r8 - w*(-49.0_r8 +w*(-9.0_r8 
     1                     + w*(19.0_r8 + w*(-3.0_r8 + w)*(-3.0_r8 
     2                     + w2)))))/144.0_r8
              yweight(6) = (40.0_r8/7.0_r8 + w*(56.0_r8/3.0_r8 
     1                     + w*(24.0_r8 + w*(40.0_r8/3.0_r8 
     2                     + w2*(-4.0_r8 + w*(-2.0_r8 + w))))))/240.0_r8
              yweight(8) = w2
              yweight(8) = yweight(8)*yweight(8)*yweight(8)
              yweight(8) = yweight(8)*w/5040.0_r8
              yweight(7) = 1.0_r8 - yweight(1) - yweight(2) 
     1                     - yweight(3) - yweight(4) - yweight(5) 
     2                     - yweight(6) - yweight(8)
            CASE (8)
              w = xnew(i) - REAL(xindex(5), KIND=r8)
              xweight(1) = 1.0_r8/2.0_r8 - w
              xweight(1) = xweight(1)*xweight(1)
              xweight(1) = xweight(1)*xweight(1)
              xweight(1) = xweight(1)*xweight(1)/40320.0_r8
              w2 = w*w
              xweight(2) = (39.0_r8/16.0_r8 - w*(6.0_r8 
     1                     + w*(-9.0_r8/2.0_r8 + w2)))*(21.0_r8/16.0_r8 
     2                     + w*(-15.0_r8/4.0_r8 + w*(9.0_r8/2.0_r8 
     3                     + w*(-3.0_r8 + w))))/5040.0_r8
              xweight(3) = (82903.0_r8/1792.0_r8 + w*(-4177.0_r8/32.0_r8
     1                     + w*(2275.0_r8/16.0_r8 + w*(-487.0_r8/8.0_r8 
     2                     + w*(-85.0_r8/8.0_r8 + w*(41.0_r8/2.0_r8 
     3                     + w*(-5.0_r8 + w*(-2.0_r8 + w))))))))
     4                     /1440.0_r8
              xweight(4) = (310661.0_r8/1792.0_r8 - w*(14219.0_r8/
     1                     64.0_r8 
     2                     + w*(-199.0_r8/8.0_r8 + w*(-1327.0_r8/16.0_r8
     3                     + w*(245.0_r8/8.0_r8 + w*(53.0_r8/4.0_r8 
     4                     + w*(-8.0_r8 + w*(-1.0_r8 + w))))))))
     5                     /720.0_r8
              xweight(5) = (2337507.0_r8/8960.0_r8 + w2*(-2601.0_r8
     1                     /16.0_r8 
     2                     + w2*(387.0_r8/8.0_r8 + w2*(-9.0_r8 + w2)))) 
     3                     /576.0_r8
              xweight(6) = (310661.0_r8/1792.0_r8 - w*(-14219.0_r8
     1                     /64.0_r8 
     2                     + w*(-199.0_r8/8.0_r8 + w*(1327.0_r8/16.0_r8 
     3                     + w*(245.0_r8/8.0_r8 + w*(-53.0_r8/4.0_r8 
     4                     + w*(-8.0_r8 + w * (1.0_r8 + w))))))))
     5                     / 720.0_r8
              xweight(8) = (39.0_r8/16.0_r8 - w*(-6.0_r8 
     1                     + w*(-9.0_r8/2.0_r8 + w2)))*(21.0_r8/16.0_r8 
     2                     + w*(15.0_r8/4.0_r8 + w*(9.0_r8/2.0_r8 
     3                     + w*(3.0_r8 + w))))/5040.0_r8
              xweight(9) = 1.0_r8/2.0_r8 + w
              xweight(9) = xweight(9)*xweight(9)
              xweight(9) = xweight(9)*xweight(9)
              xweight(9) = xweight(9)*xweight(9)/40320.0_r8
              xweight(7) = 1.0_r8 - xweight(1) - xweight(2) - xweight(3)
     1                     - xweight(4) - xweight(5) - xweight(6) 
     2                     - xweight(8) - xweight(9)                  
              w = ynew(i) - REAL(yindex(5), KIND=r8)
              yweight(1) = 1.0_r8/2.0_r8 - w
              yweight(1) = yweight(1)*yweight(1)
              yweight(1) = yweight(1)*yweight(1)
              yweight(1) = yweight(1)*yweight(1)/40320.0_r8
              w2 = w*w
              yweight(2) = (39.0_r8/16.0_r8 - w*(6.0_r8 
     1                     + w*(-9.0_r8/2.0_r8 + w2)))*(21.0_r8/16.0_r8 
     2                     + w*(-15.0_r8/4.0_r8 + w*(9.0_r8/2.0_r8 
     3                     + w*(-3.0_r8 + w))))/5040.0_r8
              yweight(3) = (82903.0_r8/1792.0_r8 + w*(-4177.0_r8/32.0_r8
     1                     + w*(2275.0_r8/16.0_r8 + w*(-487.0_r8/8.0_r8 
     2                     + w*(-85.0_r8/8.0_r8 + w*(41.0_r8/2.0_r8 
     3                     + w*(-5.0_r8 + w*(-2.0_r8 + w))))))))
     4                     /1440.0_r8
              yweight(4) = (310661.0_r8/1792.0_r8 - w*(14219.0_r8/
     1                     64.0_r8 
     2                     + w*(-199.0_r8/8.0_r8 + w*(-1327.0_r8/
     3                     16.0_r8
     4                     + w*(245.0_r8/8.0_r8 + w*(53.0_r8/4.0_r8 
     5                     + w*(-8.0_r8 + w*(-1.0_r8 + w))))))))
     6                     /720.0_r8
              yweight(5) = (2337507.0_r8/8960.0_r8 + w2*(-2601.0_r8/
     1                     16.0_r8 
     2                     + w2*(387.0_r8/8.0_r8 + w2*(-9.0_r8 + w2)))) 
     3                     /576.0_r8
              yweight(6) = (310661.0_r8/1792.0_r8 - w*(-14219.0_r8/
     1                     64.0_r8
     2                     + w*(-199.0_r8/8.0_r8 + w*(1327.0_r8/16.0_r8 
     3                     + w*(245.0_r8/8.0_r8 + w*(-53.0_r8/4.0_r8 
     4                     + w*(-8.0_r8 + w * (1.0_r8 + w))))))))/ 
     5                     720.0_r8
              yweight(8) = (39.0_r8/16.0_r8 - w*(-6.0_r8 
     1                     + w*(-9.0_r8/2.0_r8 + w2)))*(21.0_r8/16.0_r8 
     2                     + w*(15.0_r8/4.0_r8 + w*(9.0_r8/2.0_r8 
     3                     + w*(3.0_r8 + w))))/5040.0_r8
              yweight(9) = 1.0_r8/2.0_r8 + w
              yweight(9) = yweight(9)*yweight(9)
              yweight(9) = yweight(9)*yweight(9)
              yweight(9) = yweight(9)*yweight(9)/40320.0_r8
              yweight(7) = 1.0_r8 - yweight(1) - yweight(2) - yweight(3)
     1                     - yweight(4) - yweight(5) - yweight(6) 
     2                     - yweight(8) - yweight(9)                  
            CASE (9)
              w = xnew(i) - REAL(xindex(5), KIND=r8)
              xweight(1) = 1.0_r8 - w
              xweight(1) = xweight(1)*xweight(1)
              xweight(1) = xweight(1)*xweight(1)
              xweight(1) = xweight(1)*xweight(1)*(1.0_r8 - w)
     1                     /362880.0_r8
              xweight(2) = (502.0_r8/9.0_r8 + w*(-246.0_r8 + w*(472.0_r8
     1                     + w*(-504.0_r8 + w*(308.0_r8 + w*(-84.0_r8 
     2                     + w*(-56.0_r8/3.0_r8 + w*(24.0_r8 
     3                     + w*(-8.0_r8 + w)))))))))/40320.0_r8
              xweight(3) = (3652.0_r8/9.0_r8 - w*(2023.0_r8/2.0_r8 
     1                     + w*(-952.0_r8 + w*(938.0_r8/3.0_r8 
     2                     + w*(112.0_r8 + w*(-119.0_r8 
     3                     + w*(56.0_r8/3.0_r8 + w*(14.0_r8 + w*(-7.0_r8
     4                     + w)))))))))/10080.0_r8
              xweight(4) = (44117.0_r8/42.0_r8 + w*(-2427.0_r8/2.0_r8 
     1                     + w*(66.0_r8 + w*(434.0_r8 + w*(-129.0_r8 
     2                     + w*(-69.0_r8 + w*(34.0_r8 + w*(6.0_r8 
     3                     + w*(-6.0_r8 + w)))))))))/4320.0_r8
              w2 = w*w
              xweight(5) = (78095.0_r8/63.0_r8 - w2*(700.0_r8 
     1                     + w2*(-190.0_r8 + w2*(100.0_r8/3.0_r8 
     2                     + w2*(-5.0_r8 + w)))))/2880.0_r8
              xweight(6) = (44117.0_r8/63.0_r8 + w*(809.0_r8 
     1                     + w*(44.0_r8 + w*(-868.0_r8/3.0_r8 
     2                     + w*(-86.0_r8 + w*(46.0_r8 
     3                     + w*(68.0_r8/3.0_r8 + w*(-4.0_r8 + w*(-4.0_r8
     4                     + w)))))))))/2880.0_r8
              xweight(7) = (3652.0_r8/21.0_r8 - w*(-867.0_r8/2.0_r8 
     1                     + w*(-408.0_r8 + w*(-134.0_r8 + w*(48.0_r8 
     2                     + w*(51.0_r8 + w*(-4.0_r8 + w)*(-1.0_r8 
     3                     + w)*(2.0_r8 + w)))))))/4320.0_r8
              xweight(8) = (251.0_r8/18.0_r8 + w*(123.0_r8/2.0_r8 
     1                     + w*(118.0_r8 + w*(126.0_r8 + w*(77.0_r8 
     2                     + w*(21.0_r8 + w*(-14.0_r8/3.0_r8 
     3                     + w*(-6.0_r8 + w*(-2.0_r8 + w)))))))))/ 
     4                     10080.0_r8
              xweight(10) = w2*w2
              xweight(10) = xweight(10)*xweight(10)*w/362880.0_r8
              xweight(9) = 1.0_r8 - xweight(1) - xweight(2) 
     1                     - xweight(3) - xweight(4) - xweight(5) 
     2                     - xweight(6) - xweight(7) - xweight(8) 
     3                     - xweight(10)
              w = ynew(j) - REAL(yindex(5), KIND=r8)
              yweight(1) = 1.0_r8 - w
              yweight(1) = yweight(1)*yweight(1)
              yweight(1) = yweight(1)*yweight(1)
              yweight(1) = yweight(1)*yweight(1)*(1.0_r8 - w)
     1                     /362880.0_r8
              yweight(2) = (502.0_r8/9.0_r8 + w*(-246.0_r8 + w*(472.0_r8
     1                     + w*(-504.0_r8 + w*(308.0_r8 + w*(-84.0_r8 
     2                     + w*(-56.0_r8/3.0_r8 + w*(24.0_r8 
     3                     + w*(-8.0_r8 + w)))))))))/40320.0_r8
              yweight(3) = (3652.0_r8/9.0_r8 - w*(2023.0_r8/2.0_r8 
     1                     + w*(-952.0_r8 + w*(938.0_r8/3.0_r8 
     2                     + w*(112.0_r8 + w*(-119.0_r8 
     3                     + w*(56.0_r8/3.0_r8 + w*(14.0_r8 + w*(-7.0_r8
     4                     + w)))))))))/10080.0_r8
              yweight(4) = (44117.0_r8/42.0_r8 + w*(-2427.0_r8/2.0_r8 
     1                     + w*(66.0_r8 + w*(434.0_r8 + w*(-129.0_r8 
     2                     + w*(-69.0_r8 + w*(34.0_r8 + w*(6.0_r8 
     3                     + w*(-6.0_r8 + w)))))))))/4320.0_r8
              w2 = w*w
              yweight(5) = (78095.0_r8/63.0_r8 - w2*(700.0_r8 
     1                     + w2*(-190.0_r8 + w2*(100.0_r8/3.0_r8 
     2                     + w2*(-5.0_r8 + w)))))/2880.0_r8
              yweight(6) = (44117.0_r8/63.0_r8 + w*(809.0_r8 
     1                     + w*(44.0_r8 + w*(-868.0_r8/3.0_r8 
     2                     + w*(-86.0_r8 + w*(46.0_r8 
     3                     + w*(68.0_r8/3.0_r8 + w*(-4.0_r8 + w*(-4.0_r8
     4                     + w)))))))))/2880.0_r8
              yweight(7) = (3652.0_r8/21.0_r8 - w*(-867.0_r8/2.0_r8 
     1                     + w*(-408.0_r8 + w*(-134.0_r8 + w*(48.0_r8 
     2                     + w*(51.0_r8 + w*(-4.0_r8 + w)*(-1.0_r8 
     3                     + w)*(2.0_r8 + w)))))))/4320.0_r8
              yweight(8) = (251.0_r8/18.0_r8 + w*(123.0_r8/2.0_r8 
     1                     + w*(118.0_r8 + w*(126.0_r8 + w*(77.0_r8 
     2                     + w*(21.0_r8 + w*(-14.0_r8/3.0_r8 
     3                     + w*(-6.0_r8 + w*(-2.0_r8 + w)))))))))/ 
     4                     10080.0_r8
              yweight(10) = w2*w2
              yweight(10) = yweight(10)*yweight(10)*w/362880.0_r8
              yweight(9) = 1.0_r8 - yweight(1) - yweight(2) 
     1                     - yweight(3) - yweight(4) - yweight(5) 
     2                     - yweight(6) - yweight(7) - yweight(8) 
     3                     - yweight(10)
            CASE DEFAULT
              write (6,*) 'bspline_ss: Invalid bspline degree'
              stop
          END SELECT
        
c         Mirror boundary conditions
          DO k=1,degree+1
            IF (nx == 1) THEN
              xindex(k) = 1
            ELSE IF (xindex(k) < 1) THEN
              xindex(k) = -xindex(k) - width2*((-xindex(k))/width2) + 2
            ELSE
              xindex(k) = xindex(k) - width2*(xindex(k)/width2)
            END IF
            IF (nx < xindex(k)) THEN
              xindex(k) = width2 - xindex(k) + 2
            END IF
            IF (ny == 1) THEN
              yindex(k) = 1
            ELSE IF (yindex(k) < 1) THEN
              yindex(k) = -yindex(k) - height2*((-yindex(k))/height2)+2
            ELSE
              yindex(k) = yindex(k) - height2*(yindex(k)/height2)
            END IF
            IF (ny < yindex(k)) THEN
              yindex(k) = height2 - yindex(k) + 2 
            END IF
          END DO

c         Interpolate
          DO n=1,degree+1
            w = 0.0_r8
            DO k=1,degree+1
              w = w + xweight(k)*coeffs(xindex(k), yindex(n))
            END DO
            finterp(i, j) = finterp(i, j) + yweight(n)*w
          END DO
        END DO
      END DO
      DEALLOCATE(poles, coeffs)
      END SUBROUTINE bspline_ss
