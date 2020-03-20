      subroutine pad_int_gen_ss(m_orig,n_orig,mpad0,npad0,mpadb,mpadt,
     1 npadl,npadr)
c
c+
c - - Purpose: Compute padding on left, right, bottom, top such that
c              n=n_orig+npadl+npadr is divisible by 12, and m=m_orig+mpadb+mpadt
c              is also divisible by 12.  mpad0,npad0 are initial guesses as to
c              the amounts of padding.
c
c - - Usage:   call pad_int_gen_ss(m_orig,n_orig,mpad0,npad0,mpadb,mpadt,
c              npadl,npadr)
c
c - - Input:   m_orig,n_orig: integer values for the number of unpadded cells
c              in the latitude and longitude directions, respectively.
c
c - - Input:   mpad0,npad0: integer values for a "first guess" as to the 
c              amount of desired padding in longitude and latitude directions,
c              respectively.  These values will usually be upper limits to 
c              the chosen values.
c
c - - Output:  mpadb,mpadt:  Integer values of number of padded cells in 
c              latitude on the bottom and top sides of the array, respectively.
c
c - - Output:  npadl, npadr:  Integer values of number of padded cells in 
c              longitude on the left and right sides of the array, respectively.
c-
c - -  PDFI_SS electric field inversion software
c - -  http://cgem.ssl.berkeley.edu/~fisher/public/software/PDFI_SS
c - -  Copyright (C) 2015-2019 Regents of the University of California
c 
c - -  This software is based on the concepts described in Kazachenko et al.
c - -  (2014, ApJ 795, 17).  A detailed description of the software is in
c - -  Fisher et al. (2019, arXiv:1912.08301 ).
c - -  If you use the software in a scientific 
c - -  publication, the authors would appreciate a citation to these papers 
c - -  and any future papers describing updates to the methods.
c
c - -  This is free software; you can redistribute it and/or
c - -  modify it under the terms of the GNU Lesser General Public
c - -  License as published by the Free Software Foundation,
c - -  version 2.1 of the License.
c
c - -  This software is distributed in the hope that it will be useful,
c - -  but WITHOUT ANY WARRANTY; without even the implied warranty of
c - -  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c - -  See the GNU Lesser General Public License for more details.
c
c - -  To view the GNU Lesser General Public License visit
c - -  http://www.gnu.org/copyleft/lesser.html
c - -  or write to the Free Software Foundation, Inc.,
c - -  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
c
      implicit none
c
c - - input variable declarations:
c
      integer :: n_orig,m_orig,npad0,mpad0
c
c - - output variable declarations:
c
      integer :: npadl,npadr,mpadb,mpadt
c
c - - local variable declarations: "cf", common-factor, will be set to 12.
c
      integer :: cf
c
c - - If we ever decide to use something other than 12, we can set it here.
c
      cf=12
c
c - - Following expressions originally from Dave Bercik.
c - - Latitude padding:
c

      mpadb=(((m_orig + 2*mpad0 + mod(m_orig,2))/cf)*cf-m_orig)/2
      mpadt=mpadb+mod(m_orig,2)
c
c - - Longitude padding:
c
      npadl=(((n_orig + 2*npad0 + mod(n_orig,2))/cf)*cf-n_orig)/2
      npadr=npadl+mod(n_orig,2)
c
c - - We're done
c
      return
      end
