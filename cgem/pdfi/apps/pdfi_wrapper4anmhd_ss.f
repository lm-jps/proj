      subroutine pdfi_wrapper4anmhd_ss(m,n,rsun,a,b,c,d,
     1 bloncoe0,blatcoe0,brllcoe0,bloncoe1,blatcoe1,brllcoe1,
     2 vloncoe0,vlatcoe0,vlosllcoe0,vloncoe1,vlatcoe1,vlosllcoe1,
     3 lloncoe0,llatcoe0,lrllcoe0,lloncoe1,llatcoe1,lrllcoe1,
     4 tjul0,tjul1,blon0,blat0,brll0,blon1,blat1,brll1,
     5 elonpdfi,elatpdfi,erllpdfi,srll,srtot,hmll,hmtot,tjulhalf)
c    
c    
c
c - - documentation below can be seen using DOC_LIB,'pdfi_wrapper4jsoc_ss.f'
c - - in an IDL session.
c+
c    Purpose:  This is the subroutine that computes the entire PDFI_SS 
c              solution for the ANMHD test case, given the input data in HMI
c              format for a pair of observation
c              times.  It is identical to pdfi_wrapper4jsoc_ss, except for
c              the value of bthr, which is set to 370G instead of 200G.
c
c    Usage:    call pdfi_wrapper4anmhd_ss(m,n,rsun,a,b,c,d,
c              bloncoe0,blatcoe0,brllcoe0,bloncoe1,blatcoe1,brllcoe1,
c              vloncoe0,vlatcoe0,vlosllcoe0,vloncoe1,vlatcoe1,vlosllcoe1,
c              lloncoe0,llatcoe0,lrllcoe0,lloncoe1,llatcoe1,lrllcoe1,
c              tjul0,tjul1,blon0,blat0,brll0,blon1,blat1,brll1,
c              elonpdfi,elatpdfi,erllpdfi,srll,srtot,hmll,hmtot,tjulhalf)
c
c     Input:   m,n - number of cell centers in the theta (lat), and phi (lon)
c              directions, respectively. 
c
c     Input:   rsun - assumed radius of the Sun [km].  Normally set to 6.96d5.
c
c     Input:   a,b - real*8 Minimum and maximum values of co-latitude [radians]
c              corresponding to the range of the theta (colatitude) edge values.
c
c     Input:   c,d - real*8 Minimum and maximum values of longitude [radians]
c              corresponding to the range of longitude (azimuth) edge values.
c
c     Input:   bloncoe0(n+1,m+1),blatcoe0(n+1,m+1),brllcoe0(n+1,m+1) - real*8
c              arrays of longitudinal, latitudinal and radial components of the 
c              magnetic field evaluated at COE locations at time t0 
c              (corners plus exterior corners on boundary). [G]
c
c     Input:   bloncoe1(n+1,m+1),blatcoe1(n+1,m+1),brllcoe1(n+1,m+1) - 
c              real*8 arrays of the longitudinal, latitudinal and radial 
c              components of the magnetic field evaluated at COE locations at 
c              time t1 (corners plus exterior corners on boundary). [G]
c
c     Input:   vloncoe0(n+1,m+1),vlatcoe0(n+1,m+1),vlosllcoe0(n+1,m+1) - real*8
c              arrays of the longitudinal, latitudinal and LOS components of 
c              the velocity field evaluated at COE locations at time t0 
c              (corners plus exterior corners on boundary).
c              vloncoe0, vlatcoe0 have units of [km/s]; vlosllcoe1 has units 
c              of [m/s] and is directed such that positive values 
c              correspond to redshift. 
c
c     Input:   vloncoe1(n+1,m+1),vlatcoe1(n+1,m+1),vlosllcoe1(n+1,m+1) - real*8
c              arrays of the longitudinal, latitudinal and LOS components of 
c              the velocity field evaluated at COE locations at time t1 
c              (corners plus exterior corners on boundary). 
c              vloncoe1 and vlatcoe1 have units of [km/s]; vlosllcoe1 has 
c              units of [m/s] and is 
c              directed such that positive values correspond to redshift. 
c
c     Input:   lloncoe0(n+1,m+1),llatcoe0(n+1,m+1),lrllcoe0(n+1,m+1) - real*8
c              arrays of the longitudinal, latitudinal and radial components of 
c              the LOS unit vector evaluated at COE locations at time t0 
c              (corners plus exterior corners on boundary).
c
c     Input:   lloncoe1(n+1,m+1),llatcoe1(n+1,m+1),lrllcoe1(n+1,m+1) - real*8
c              arrays of the longitudinal, latitudinal and radial components 
c              of the LOS unit vector evaluated at COE locations at time t1 
c              (corners plus exterior corners on boundary).
c
c     Input:   tjul0,tjul1 - real*8 values of times t0 and t1 [days].
c
c     Output:  elonpdfi(n,m+1),elatpdfi(n+1,m),erllpdfi(n+1,m+1) - real*8 
c              arrays of longitudinal, latidudinal and radial components of 
c              electric field, stored in lon,lat index order [V/cm].
c              If desired, these electric fields can be converted to cE values
c              in [G km/sec] by multiplying by 1d3.
c
c     Output:  blon0(n+1,m),blat0(n,m+1),brll0(n,m),blon1(n+1,m),
c              blat1(n,m+1),brll1(n,m) - real*8 arrays of magnetic field 
c              values at the staggered grid locations at both the t0 and t1 
c              times, in (lon,lat) order. blon0,blon1 on PE grid, blat0,blat1
c              on TE grid, brll0,brll1 on CE grid. [G]
c
c     Output:  srll(n,m) - real*8 array of radial component of Poynting flux 
c              values in lon,lat order on CE grid [erg/(cm^2-sec)]
c
c     Output:  srtot - real*8 value of the area integral of srll [erg/sec]
c
c     Output:  hmll(n,m) - real*8 array of Helicity flux density in lon,lat 
c              order.  [Mx^2/(cm^2-sec)]
c
c     Output:  hmtot - real*8 value of the area integral of hmll. [Mx^2/sec]
c
c     Output:  tjulhalf - real*8 value of time halfway between t0 and t1 [days].
c
c - - MKD, January 2016
c - - GHF, mods OCT 2016 to add ci,di = 0,d-c and add output var tjulhalf
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
c input variables in lon, lat order
c
      integer :: n,m
      real*8 rsun,a,b,c,d 
      real*8 tjul0,tjul1
c
      real*8 :: bloncoe0(n+1,m+1),blatcoe0(n+1,m+1),brllcoe0(n+1,m+1)
      real*8 :: bloncoe1(n+1,m+1),blatcoe1(n+1,m+1),brllcoe1(n+1,m+1)
c
      real*8 :: vloncoe0(n+1,m+1),vlatcoe0(n+1,m+1),vlosllcoe0(n+1,m+1)
      real*8 :: vloncoe1(n+1,m+1),vlatcoe1(n+1,m+1),vlosllcoe1(n+1,m+1)
c
      real*8 :: lloncoe0(n+1,m+1),llatcoe0(n+1,m+1),lrllcoe0(n+1,m+1)
      real*8 :: lloncoe1(n+1,m+1),llatcoe1(n+1,m+1),lrllcoe1(n+1,m+1)
c
c output: electric field at TE and PE edges   
c
      real*8 :: elonpdfi(n,m+1),elatpdfi(n+1,m)
      real*8 :: erllpdfi(n+1,m+1)
c
c output: Poynting flux, helicity flux density arrays at CE grid, lon,lat order
c
      real*8 :: srll(n,m),hmll(n,m)
c
c output: horiz magnetic fields at TE, PE edges and Br at cell centers
c
      real*8 :: blon0(n+1,m),blat0(n,m+1),brll0(n,m)
      real*8 :: blon1(n+1,m),blat1(n,m+1),brll1(n,m)
c
c output: time halfway between magnetogram times for electric fields
c
      real*8 :: tjulhalf
c
c output:  area integrals of Poynting flux, Helicity injection rate density:
c
      real*8 :: srtot,hmtot
c      
c - - local variables 
c 
c b-magnitude and field changes in COE     
      real*8 :: bmagllcoe0(n+1,m+1),bmagllcoe1(n+1,m+1)
      real*8 :: dbloncoe(n+1,m+1),dblatcoe(n+1,m+1),dbrllcoe(n+1,m+1)
c
c variables at half (h) time values (t1-t0)/2  in lon, lat order in COE
      real*8 :: bloncoeh(n+1,m+1),blatcoeh(n+1,m+1),brllcoeh(n+1,m+1)
      real*8 :: vloncoeh(n+1,m+1),vlatcoeh(n+1,m+1),vlosllcoeh(n+1,m+1)
      real*8 :: lloncoeh(n+1,m+1),llatcoeh(n+1,m+1),lrllcoeh(n+1,m+1)
c
c variables at half time in theta, phi order in COE
      real*8 :: btcoeh(m+1,n+1),bpcoeh(m+1,n+1),brtpcoeh(m+1,n+1)
      real*8 :: vtcoeh(m+1,n+1),vpcoeh(m+1,n+1),vlostpcoeh(m+1,n+1)
      real*8 :: ltcoeh(m+1,n+1),lpcoeh(m+1,n+1),lrtpcoeh(m+1,n+1)
c
c bmag in theta, phi order at t0 and t1      
      real*8 :: bmagtpcoe0(m+1,n+1),bmagtpcoe1(m+1,n+1)
c
c change in dB vector in theta, phi order     
      real*8 :: dbltcoe(m+1,n+1),dblpcoe(m+1,n+1),dbrtpcoe(m+1,n+1)
c
c variables to calculate TE and PE at t0 and t1
      real*8 :: btcoe0(m+1,n+1),bpcoe0(m+1,n+1)
      real*8 :: brtpcoe0(m+1,n+1)
      real*8 :: btcoe1(m+1,n+1),bpcoe1(m+1,n+1)
      real*8 :: brtpcoe1(m+1,n+1)
      real*8 :: bt0(m+1,n),bp0(m,n+1),br0(m,n)
      real*8 :: bt1(m+1,n),bp1(m,n+1),br1(m,n)
c      
c additional internal variables      
      real*8 :: ci,di
c  ci and di are values of c,d with c subtracted from both, for use in fishpack
      real*8 :: maskcoe(m+1,n+1)
      real*8 :: bt(m+1,n),bp(m,n+1),br(m,n)
      real*8 :: btt(m+1,n),bpt(m,n+1),brt(m,n)
      real*8 :: brte(m+1,n),brpe(m,n+1)
      real*8 :: bhte(m+1,n),bhpe(m,n+1)
      real*8 :: vt(m+1,n),vp(m,n+1)
      real*8 :: mskte(m+1,n),mskpe(m,n+1),mskr(m,n)
      real*8 :: sinth(m+1),sinth_hlf(m)
      real*8 :: scrb(m+1+1,n+1+1),dscrbdr(m+1+1,n+1+1)
      real*8 :: scrbs(m+1+1,n+1+1),dscrbdrs(m+1+1,n+1+1)
      real*8 :: scrj(m+1,n+1),scrjs(m+1,n+1)
      real*8 :: et(m,n+1),ep(m+1,n),er(m+1,n+1)     
      real*8 :: ezetat(m,n+1),ezetap(m+1,n)
      real*8 :: edt(m,n+1),edp(m+1,n),edr(m+1,n+1)
      real*8 :: eti(m,n+1),epi(m+1,n),eri(m+1,n+1)
      real*8 :: etpdf(m,n+1),eppdf(m+1,n),erpdf(m+1,n+1)
      real*8 :: etpdfi(m,n+1),eppdfi(m+1,n),erpdfi(m+1,n+1)
      real*8 :: sr(m,n),hm(m,n)
c
c keywords: bthr -|B| threshold for masking 
      real*8 :: bthr
c keywords: bthr_relax -|B| threshold for relaxation
      real*8 :: bthr_relax
c number of iterations, verbose (1) or mute (0), show debug (1) or not
      integer :: max_iter,verbose,debug
c maskflag to determine mask values when interpolated value .ne. 0 or 1
c     maskflag should be set to either 0 or 1
      integer :: maskflag
c
c additional constants
      real*8 dphi,dtheta
      real*8 dt
      real*8 pi,dum
c - - debug variables (for sdf) - uncomment if needed
c     integer*8 :: dims(20)
c     
c Function in FISHPACK (FFTPACK) that computes pi
c
      real*8 pimach
c
c---------------------------------------------------------------------------c
c      
      max_iter=25
      verbose=0
      debug=1
c     bthr=200.d0
c - - Set bthr to 370G, to be consistent with KFW2014 threshold
      bthr=370.d0
      bthr_relax=0.d0
      maskflag=0
c      
      pi=pimach(dum)
c - - radius of the Sun [ in km] is now set as an input parameter
c     rsun=6.96d5
      dt=(tjul1-tjul0)*(24.*3600)
      dtheta=(b-a)/(m)
      dphi=(d-c)/(n)
c - - ci and di are set to c,d with c subtracted from both, for use in fishpack
      ci=0.d0
      di=d-c
c      
      bmagllcoe0=sqrt(bloncoe0*bloncoe0+blatcoe0*blatcoe0+
     1 brllcoe0*brllcoe0)
      bmagllcoe1=sqrt(bloncoe1*bloncoe1+blatcoe1*blatcoe1+
     2 brllcoe1*brllcoe1)
c           
      bloncoeh=(bloncoe1+bloncoe0)/2.0d0
      blatcoeh=(blatcoe1+blatcoe0)/2.0d0
      brllcoeh=(brllcoe1+brllcoe0)/2.0d0
      vloncoeh=(vloncoe1+vloncoe0)/2.0d0
      vlatcoeh=(vlatcoe1+vlatcoe0)/2.0d0
c
c - - In next statement we convert LOS velocities from m/sec to km/sec,
c - - and also change sign convention from positive=> redshift to 
c     positive => blueshift.
c
      vlosllcoeh=-(vlosllcoe1+vlosllcoe0)/(2000.0d0)
      lloncoeh=(lloncoe1+lloncoe0)/2.0d0
      llatcoeh=(llatcoe1+llatcoe0)/2.0d0
      lrllcoeh=(lrllcoe1+lrllcoe0)/2.0d0 
c
      dbloncoe=bloncoe1-bloncoe0
      dblatcoe=blatcoe1-blatcoe0
      dbrllcoe=brllcoe1-brllcoe0    
c
c Transpose B_h,v_h,l_h data arrays from lon,lat to theta,phi order and 
c flip sign to get B_theta.
      call bhll2tp_ss(m,n,bloncoeh,blatcoeh,btcoeh,bpcoeh)
      call bhll2tp_ss(m,n,vloncoeh,vlatcoeh,vtcoeh,vpcoeh)
      call bhll2tp_ss(m,n,lloncoeh,llatcoeh,ltcoeh,lpcoeh)
c      
c Transpose B_r,br,lr data arrays from lon,lat to theta,phi order.
      call brll2tp_ss(m,n,brllcoeh,brtpcoeh)
      call brll2tp_ss(m,n,vlosllcoeh,vlostpcoeh)
      call brll2tp_ss(m,n,lrllcoeh,lrtpcoeh)   
c         
c Transpose magnetic field time derivatives from lon,lat order to t,p order
      call bhll2tp_ss(m,n,dbloncoe,dblatcoe,dbltcoe,dblpcoe)
      call brll2tp_ss(m,n,bmagllcoe0,bmagtpcoe0)
      call brll2tp_ss(m,n,bmagllcoe1,bmagtpcoe1)
      call brll2tp_ss(m,n,dbrllcoe,dbrtpcoe)
c 
c Calculate mask using bmag at t0 and t1.  Note normally find_mask_ss uses
c magnetograms at 3 times, but since we've switched to two adjacent times,
c we call it with the magnetic field at the first time repeated.
      call find_mask_ss(m,n,bmagtpcoe0,bmagtpcoe0,
     1 bmagtpcoe1,bthr,maskcoe)
c
c Interpolate db/dt from COE corner/edges to theta and phi edge locations
c Input data: dbltcoe(m+1,n+1),dblpcoe(m+1,n+1),dbrtpcoe(m+1,n+1)
c output data: btt(m+1,n),bpt(m,n+1),brt(m,n)
c output data: brte(m+1,n),brpe(m,n+1), bhte(m+1,n),bhpe(m,n+1)
c output data: vt(m+1,n),vp(m,n+1)
      call interp_data_ss(m,n,dbltcoe/dt,dblpcoe/dt,
     1 dbrtpcoe/dt,vtcoeh,vpcoeh,a,b,btt,bpt,brt,brte,brpe,bhte,bhpe,
     2 vt,vp)
c
c     
      call interp_data_ss(m,n,btcoeh,bpcoeh,brtpcoeh,
     1 vtcoeh,vpcoeh,a,b,bt,bp,br,brte,brpe,bhte,bhpe,
     2 vt,vp)
c
      call interp_var_ss(m,n,maskcoe,maskcoe,maskcoe,
     1 a,b,mskte,mskpe,mskr)
c
c - - fix mask values in cases where interpolated value is between 0 and 1
c
      call fix_mask_ss(m+1,n,maskflag,mskte)
      call fix_mask_ss(m,n+1,maskflag,mskpe)
      call fix_mask_ss(m,n,maskflag,mskr)
c         
c Compute arrays of sin(theta) for the edge locations and at the
c half-cell locations
      call sinthta_ss(a,b,m,sinth,sinth_hlf)
c
c Solve for time derivatives of scrb, dscrbdr, and scrj     
c Note that the input time derivatives are masked
c
c     call ptdsolve_ss(m,n,btt*mskte,bpt*mskpe,brt*mskr,
c    1 rsun,sinth,sinth_hlf,a,b,ci,di,scrb,dscrbdr,scrj)
      call ptdsolve_ss(m,n,btt,bpt,brt,
     1 rsun,sinth,sinth_hlf,a,b,ci,di,scrb,dscrbdr,scrj)
c
c Solve for static values of scrb,dscrbdr,scrj using magnetic field values
c to compute vector potential [if needed]
      call ptdsolve_ss(m,n,bt,bp,br,
     1 rsun,sinth,sinth_hlf,a,b,ci,di,scrbs,dscrbdrs,scrjs)
c
c Calculate PTD contribution into electric field
      call e_ptd_ss(m,n,scrb,scrj,rsun,sinth_hlf,dtheta,dphi,
     1 et,ep,er) 
c
c Calculate FLCT contribution into electric field
c Note the input magnetic fields are masked
c
c     call e_flct_ss(m,n,brte*mskte,brpe*mskpe,
c    1 bhte*mskte,bhpe*mskpe,vt,vp,rsun,sinth,sinth_hlf,
c    2 a,b,ci,di,ezetat,ezetap)
      call e_flct_ss(m,n,brte,brpe,
     1 bhte,bhpe,vt,vp,rsun,sinth,sinth_hlf,
     2 a,b,ci,di,ezetat,ezetap)
c
c Calculate Doppler contribution into electric field
c Note magnetic field arrays are masked
c
      call e_doppler_ss(m,n,btcoeh,
     1 bpcoeh,brtpcoeh,vlostpcoeh,
     2 ltcoeh,lpcoeh,lrtpcoeh,rsun,sinth,
     3 a,b,ci,di,edt,edp,edr)
c
c     call e_doppler_ss(m,n,btcoeh*maskcoe,
c    1 bpcoeh*maskcoe,brtpcoeh*maskcoe,vlostpcoeh*maskcoe,
c    2 ltcoeh,lpcoeh,lrtpcoeh,rsun,sinth,
c    3 a,b,ci,di,edt,edp,edr)
c
c Calculate PDF electric field
      etpdf=et+ezetat+edt
      eppdf=ep+ezetap+edp
      erpdf=er+edr
c
c Calculate non-inductive electric field components    
c Note that magnetic field points (corners) are masked on input.
c
      call e_ideal_ss(m,n,btcoeh,
     1 bpcoeh,brtpcoeh,etpdf,eppdf,erpdf,rsun,
     2 sinth,a,b,ci,di,eti,epi,eri)
c     call e_ideal_ss(m,n,btcoeh*maskcoe,
c    1 bpcoeh*maskcoe,brtpcoeh*maskcoe,etpdf,eppdf,erpdf,rsun,
c    2 sinth,a,b,ci,di,eti,epi,eri)
c
c Calculate PDFI electric field and then convert from cE in G-km/s to E in V/cm:
c
      etpdfi=(etpdf+eti)*1.0d-3
      eppdfi=(eppdf+epi)*1.0d-3
      erpdfi=(erpdf+eri)*1.0d-3
c
c Calculate radial Poynting flux from PDFI electric field and B_h:
      call sr_ss(m,n,etpdfi,eppdfi,bt,bp,sr)
c
c Integrate radial Poynting flux over area:
      call srtot_ss(m,n,rsun,sinth_hlf,dtheta,dphi,sr,srtot)
c
c Calculate Helicity flux density from scrbs and PDFI electric field:
      call hm_ss(m,n,rsun,sinth_hlf,dtheta,dphi,etpdfi,eppdfi,scrbs,hm)
c
c Integrate Helicity flux density over area to get helicity injection rate:
      call hmtot_ss(m,n,rsun,sinth_hlf,dtheta,dphi,hm,hmtot)
c
c Transpose E_h data arrays from theta,phi to lon,lat order 
      call ehyeetp2ll_ss(m,n,etpdfi,eppdfi,elonpdfi,elatpdfi)
      call eryeetp2ll_ss(m,n,erpdfi,erllpdfi)
c
c Calculate Bvec at t0 at TE and PE in lon, lat order
      call bhll2tp_ss(m,n,bloncoe0,blatcoe0,btcoe0,bpcoe0)
      call brll2tp_ss(m,n,brllcoe0,brtpcoe0)
      call interp_var_ss(m,n,btcoe0,bpcoe0,brtpcoe0,
     1 a,b,bt0,bp0,br0)
      call bhyeetp2ll_ss(m,n,bt0,bp0,blon0,blat0)
      call bryeetp2ll_ss(m,n,br0,brll0)
c - - debug statements to create diagnostic .sdf file
c     dims(1)=n
c     dims(2)=m
c     call sdf_rm_f77('debug_wrapper.sdf')
c     call sdf_write_f77('debug_wrapper.sdf','brll0','f',8,2,dims,brll0)
c     dims(1)=n+1
c     dims(2)=m
c     call sdf_write_f77('debug_wrapper.sdf','blon0','f',8,2,dims,blon0)
c     dims(1)=n
c     dims(2)=m+1
c     call sdf_write_f77('debug_wrapper.sdf','blat0','f',8,2,dims,blat0)
c
c
c Calculate Bvec at t1 at TE and PE in lon, lat order      
      call bhll2tp_ss(m,n,bloncoe1,blatcoe1,btcoe1,bpcoe1)
      call brll2tp_ss(m,n,brllcoe1,brtpcoe1)            
      call interp_var_ss(m,n,btcoe1,bpcoe1,brtpcoe1,
     1 a,b,bt1,bp1,br1)
      call bhyeetp2ll_ss(m,n,bt1,bp1,blon1,blat1)
      call bryeetp2ll_ss(m,n,br1,brll1)
c
c - - rotate sr and hm into lon, lat order using brtp2ll_ss subroutine:
      call bryeetp2ll_ss(m,n,sr,srll)
      call bryeetp2ll_ss(m,n,hm,hmll)
c
c Compute time at half-way between vector magnetograms:
c
      tjulhalf=0.5d0*(tjul0+tjul1)
c    
c Debugging/Diagnostics.  If you get tired of this output, can set debug to 0
      if(debug .ne. 0) then     
        write(6,*) 'total mskte,pe,r',sum(mskte),sum(mskpe),sum(mskr)
        write(6,*) 'total bloncoeh, btcoeh',sum(bloncoeh),sum(btcoeh)
        write(6,*) 'total blatcoeh, bpcoeh',sum(blatcoeh),sum(bpcoeh)
        write(6,*) 'total vlosllcoeh, vlostpcoeh',sum(vlosllcoeh),
     1              sum(vlostpcoeh)
        write(6,*) 'total bt, bp',sum(bt),sum(bp)
        write(6,*) 'total br,te,pe',sum(br),sum(brte),sum(brpe)
        write(6,*) 'total vtcoeh,vpcoeh',sum(vtcoeh),sum(vpcoeh)
        write(6,*) 'total vt,vp',sum(vt),sum(vp)
        write(6,*) 'total bhte,bhpe',sum(bhte),sum(bhpe)
        write(6,*) 'total dbltcoe/dt',sum(dbltcoe/dt)
        write(6,*) 'total dblpcoe/dt',sum(dblpcoe/dt)
        write(6,*) 'total dbrtpcoe/dt',sum(dbrtpcoe/dt)
        write(6,*) 'sum btt,bpt,brt',sum(btt),sum(bpt),sum(brt)
        write(6,*) 'a,b,m+1,tot(sinth),tot(sinth_hlf)',a,b,m+1,
     1             sum(sinth),sum(sinth_hlf)
        write(6,*) 'after ptdsolve_ss: sum(srb), dscrbdr,scrj ',
     1             sum(scrb),sum(dscrbdr),sum(scrj)
        write(6,*) 'dtheta,dphi ', dtheta,dphi
        write(6,*) 'sum(et),sum(ep),sum(er) ',sum(et),sum(ep),
     1             sum(er)
        write(6,*) 'total(ezetat),total(ezetap) ',sum(ezetat),
     1             sum(ezetap)
        write(6,*) 'sum(edt),sum(edp),sum(edr) ',
     1             sum(edt),sum(edp),sum(edr)
        write(6,*) 'sum(etpdf),sum(eppdf),sum(erpdf) ',
     1             sum(etpdf),sum(eppdf),sum(erpdf) 
        write(6,*) 'sum(etpdfi)',sum(etpdfi)
        write(6,*) 'sum(eppdfi)',sum(eppdfi)
        write(6,*) 'sum(erpdfi)',sum(erpdfi)   
        write(6,*) 'sum(bpcoeh)',sum(bpcoeh(2:m,2:n))
        write(6,*) 'sum(btcoeh)',sum(btcoeh(2:m,2:n))
        write(6,*) 'integral Poynting flux:',srtot
        write(6,*) 'Helicity injection rate:',hmtot
        write(6,*) 'tjulhalf = ',tjulhalf
        write(6,*) 'Dt=',dt
      endif
c     
      return
      end
