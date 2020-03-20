      subroutine get_pils_rad_ss(m,n,brad,bmag,pilmap,thresh_brad, 
     1 thresh_bmag, dilation_param)
c
c+
c - - Purpose:  Given input "brad" array (m-1,n-1) corresponding to radial
c               magnetogram pixels, this returns "pilmap" an (m-1, n-1) 
c               bitmap (of integer type) for all pixels closer to opposite 
c               polarity than dilation_param (integer).  
c               Magnetized regions must satisfy both radial component and 
c               magnitude thresholds to be considered.
c
c - - Usage:    call get_pils_rad_ss(m,n,brad,bmag,pilmap,thresh_brad,
c               thresh_bmag,dilation_param)
c
c - - Input:    m,n - integers equal to the number of cell centers in colat,lon.
c
c - - Input:    brad(m-1,n-1) - real*8 array of radial magnetic field values
c               on CO grid [G]
c
c - - Input:    bmag(m-1,n-1) - real*8 array of magnetic field amplitudes on
c               CO grid [G]
c
c - - Output:   pilmap(m-1,n-1) - integer array of CO PIL locations (1 in PIL
c               regions, 0 outside them)
c
c - - Input:    thresh_brad - real*8 variable - radial field threshold to be
c               considered magnetized (abs. value) [G]
c
c - - Input:    thresh_bmag - real*8 variable - threshold for total field 
c               strength for a pixel to be considered magnetized [G]
c
c - - Input:    dilation_param - integer value defining width of PILs.  Width
c               is 2*dilation_param + 1.
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
c HISTORY: 2016/10/25, BT Welsch: get_pils.f90 started
c          2017/08/01, BTW: copied from get_pils, added bmag input
c          2017/08/03, GHF:  Modified to make compatible with PDFI_SS library
c
c
      implicit none
c
c - - Input calling arguments:
c
c no. cell centers in colat, lon:
      integer :: m,n
c input magnetogram, B_los or B_r:
      real*8 :: brad(m-1,n-1) 
c input mag. vector field strength, |B|:
      real*8 :: bmag(m-1,n-1) 
c max. dist. btwn. +/- to be a PIL:
      integer :: dilation_param 
c threshold in |B_r| below which pixel values are ignored:
      real*8 thresh_brad 
c  threshold in |B| below which pixel values are ignored:
      real*8 thresh_bmag 
c
c - - output arguments:
c
c output map of PIL pixels (integer):
      integer :: pilmap(m-1,n-1) 
c  
c - - local variable declarations:
c
      integer :: posmap(m-1,n-1),negmap(m-1,n-1) 
      integer :: dilpos(m-1,n-1),dilneg(m-1,n-1) 
c
c - - now actually doing stuff:
c
      where ((brad .ge. thresh_brad) .and. (bmag .ge. thresh_bmag))
         posmap = 1
      elsewhere
         posmap = 0
      end where

      where ((brad .le. -thresh_brad) .and. (bmag .ge. thresh_bmag))
         negmap = 1
      elsewhere
         negmap = 0
      end where
c
      call dilate_ss(m,n,posmap, dilpos, dilation_param)
      call dilate_ss(m,n,negmap, dilneg, dilation_param)
c
c print *,dilpos ! for diagnostics of dilation subroutine
c
      pilmap = dilpos*dilneg
c
      return
      end 
