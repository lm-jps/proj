      subroutine dilate_ss(m,n,map,dilmap,dilation_param)
c
c
c+
c     Purpose:  Given a 2-d integer mask array, map, and a single integer,
c               dilation_param, this subroutine computes a new mask, dilmap, 
c               that is dilated from the original mask on all sides by 
c               dilation_param pixels.  The width in each direction is thus 
c               expanded by 2*dilation_param +1 pixels.
c
c               NOTE: - this subroutine is used by subroutines get_pils_ss and 
c               get_pils_rad_ss.
c
c - - Usage:    call dilate_ss(m,n,map,dilmap,dilation_param)
c
c - - Input:    m,n - integer no. of cell centers in colatitude, longitude, 
c               resp.
c
c - - Input:    map(m-1,n-1) - an integer array containing the original mask,
c               assumed to consist of 0s or 1s.  map is assumed a CO array.
c
c - - Output:   dilmap(m-1,n-1) - an integer array on the CO grid
c               containing the dilated mask array, consisting of 0s or 1s
c
c - - Input:    dilation_param - a single integer indicating the amount of
c               desired dilation.
c-
c HISTORY: 2016/10/25, BT Welsch: started
c          2017/08/01 GHF:  Modified for compatibility with pdfi_ss library
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
      implicit none
c
c input variables:
c
      integer :: m,n,map(m-1,n-1) 
      integer :: dilation_param 
c
c output dilated map 
c
      integer :: dilmap(m-1,n-1) 
c
c - - local variables:
c
c - - shifted array
      integer :: shifted(m-1,n-1) 
c - - counting integers:
      integer :: i,j 
c - - array bounds on input map
      integer :: naxis1,naxis2
      integer :: mcoli,mcolf,mrowi,mrowf 
c - - array bounds on shifted map
      integer :: scoli,scolf,srowi,srowf 
c
c - - retain the use of naxis1,naxis2 notation for compatibility with Brian's
c - - original source code, but express these in terms of m and n.
c
      naxis1=m-1
      naxis2=n-1
c
      dilmap(:,:) = 0
c
c Shift pos- & negmap array rows
      do i=-dilation_param,dilation_param   

        shifted(:,:) = 0

        if (i.lt.0) then 
           srowi = -i+1
           srowf = naxis1
           mrowi = 1
           mrowf = naxis1 + i
        else 
           srowi = 1
           srowf = naxis1 - i
           mrowi = i+1
           mrowf = naxis1 
        end if
c       
c  - -  Shift pos- & negmap array columns
        do j=-dilation_param,dilation_param   
c
          if (j.lt.0) then 
             scoli = -j+1
             scolf = naxis2
             mcoli = 1
             mcolf = naxis2 + j
           else 
             scoli = 1
             scolf = naxis2 - j
             mcoli = j+1
             mcolf = naxis2 
           end if
              
c     print *,''
c     print *,i,j,shape(shifted)
c     print *,'S row & col:',srowf-srowi,scolf-scoli
c     print *,'M row & col:',mrowf-mrowi,mcolf-mcoli
c     print *,'========================================'

         shifted(srowi:srowf,scoli:scolf) = map(mrowi:mrowf,mcoli:mcolf)
         dilmap = dilmap + shifted
        
c - - end of j loop
        end do
c - - end of i loop
      end do

      where (dilmap .gt. 0) 
        dilmap = 1
      end where
c
      return
      end 
