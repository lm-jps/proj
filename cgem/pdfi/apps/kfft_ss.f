      subroutine kfft_ss(n,k)
c
c+
c - - Purpose:  Compute wavenumbers for use with FFTPACK's fft 
c               routines rfftf and rfftb.
c
c - - Usage:    call kfft_ss(n,k)
c
c - - Input:    n - integer value of length of 1-d vector to be transformed.
c
c - - Output:   k(n) - real*8 array of wavenumbers for each
c               fourier mode
c-
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
      implicit none
c
c - - input variables:
c
      integer :: n
c
c - - output variables:
c
      real*8 :: k(n)
c
c - - local variables:
c
      integer :: l,i
c
c - - code below inferred from FFTPACK documentation contained within
c - - FISHPACK 4.1 documentation
c
      if((n/2)*2 .eq. n) then
        l=n/2
        k(n)=l
      else
        l=(n+1)/2
      endif
      k(1)=0
      do i=2,l
         k(2*i-1)=i-1
         k(2*i-2)=i-1
      end do
      return
      end
