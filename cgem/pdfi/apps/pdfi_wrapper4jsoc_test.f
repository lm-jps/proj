      subroutine pdfi_wrapper4jsoc_test(m,n,a,b,c,d,
     1 bloncoe0,blatcoe0,brllcoe0,bloncoe1,blatcoe1,brllcoe1,
     2 vloncoe0,vlatcoe0,vrllcoe0,vloncoe1,vlatcoe1,vrllcoe1,
     3 lloncoe0,llatcoe0,lrllcoe0,lloncoe1,llatcoe1,lrllcoe1,
     4 tjul0,tjul1,blon0,blat0,brll0,blon1,blat1,brll1,
     5 elonpdfi,elatpdfi,erllpdfi)
c    
c    
c
c - - documentation below can be seen using DOC_LIBRARY,'pdfi_wrapper4jsoc_ss'
c+
c    Usage:  call pdfi_wrapper4jsoc_ss(m,n,a,b,c,d,
c     1 bloncoe0,blatcoe0,brllcoe0,bloncoe1,blatcoe1,brllcoe1,
c     2 vloncoe0,vlatcoe0,vrllcoe0,vloncoe1,vlatcoe1,vrllcoe1,
c     3 lloncoe0,llatcoe0,lrllcoe0,lloncoe1,llatcoe1,lrllcoe1,tjul0,tjul1,
c     4 blon0,blat0,brll0,blon1,blat1,brll1,elonpdfi,elatpdfi,erllpdfi)
c             This subroutine derives PDFI_SS electric fields from vector 
c             magnetogram sequences.
c     Input:  m,n - number of cell centers in the theta (lat), and phi (lon)
c             directions, respectively. 
c     Input:  a,b - Minimum and maximum values of co-latitude in radians
c             corresponding to the range of the theta (colatitude) edge values.
c     Input:  c,d - Minimum and maximum values of longitude in radians 
c             corresponding to the range of longitude (azimuth) edge values.
c
c     Input:  bloncoe0(n+1,m+1),blatcoe0(n+1,m+1),brllcoe0(n+1,m+1) - arrays 
c             of the longitudinal, latitudinal and radial components of the 
c             magnetic field evaluated at COE locations at time t0 
c             (corners plus exterior corners on boundary).
c     Input:  bloncoe1(n+1,m+1),blatcoe1(n+1,m+1),brllcoe1(n+1,m+1) - 
c             arrays of the longitudinal, latitudinal and radial components of 
c             the magnetic field evaluated at COE locations at time t1 
c             (corners plus exterior corners on boundary).
c     Input:  vloncoe0(n+1,m+1),vlatcoe0(n+1,m+1),vrllcoe0(n+1,m+1) - arrays 
c             of the longitudinal, latitudinal and radial components of 
c             the velocity field evaluated at COE locations at time t0 
c             (corners plus exterior corners on boundary).
c     Input:  vloncoe1(n+1,m+1),vlatcoe1(n+1,m+1),vrllcoe1(n+1,m+1) - arrays 
c             of the longitudinal, latitudinal and radial components of 
c             the velocity field evaluated at COE locations at time t1 
c             (corners plus exterior corners on boundary). vloncoe1 and vlatcoe1
c              have units of km/s; vrllcoe1 has units of m/s and is directed such
c              that positive values correspond to redshift. 
c     Input:  lloncoe0(n+1,m+1),llatcoe0(n+1,m+1),lrllcoe0(n+1,m+1) - arrays 
c              of the longitudinal, latitudinal and radial components of 
c             the LOS vector evaluated at COE locations at time t0 (corners 
c             plus exterior corners on boundary).
c     Input:  lloncoe1(n+1,m+1),llatcoe1(n+1,m+1),lrllcoe1(n+1,m+1) - arrays 
c             of the longitudinal, latitudinal and radial components of the 
c             LOS vector evaluated at COE locations at time t1 (corners plus 
c             exterior corners on boundary).
c     Input:  tjul0,tjul1 - times, t0 and t1, in julian days.
c
c     Output: elonpdfi(n,m+1),elatpdfi(n+1,m),erllpdfi(n+1,m+1) - arrays 
c             of lontitudinal, latidudinal and radial components of electric 
c             field, stored in lon,lat index order.
c     Output: blon0(n+1,m),blat0(n,m+1),brll0(n,m),blon1(n+1,m),
c             blat1(n,m+1),brll1(n,m) - magnetic field values at the staggered
c             grid locations at both the t0 and t1 times, in (lon,lat) order
c
c - - MKD, January 2016
c-
      implicit none      
c
c input variables in lon, lat order
c
      integer :: n,m
      real*8 a,b,c,d 
      real*8 tjul0,tjul1
c
      real*8 :: bloncoe0(n+1,m+1),blatcoe0(n+1,m+1),brllcoe0(n+1,m+1)
      real*8 :: bloncoe1(n+1,m+1),blatcoe1(n+1,m+1),brllcoe1(n+1,m+1)
c
      real*8 :: vloncoe0(n+1,m+1),vlatcoe0(n+1,m+1),vrllcoe0(n+1,m+1)
      real*8 :: vloncoe1(n+1,m+1),vlatcoe1(n+1,m+1),vrllcoe1(n+1,m+1)
c
      real*8 :: lloncoe0(n+1,m+1),llatcoe0(n+1,m+1),lrllcoe0(n+1,m+1)
      real*8 :: lloncoe1(n+1,m+1),llatcoe1(n+1,m+1),lrllcoe1(n+1,m+1)
c
c output: electric field at TE and PE edges   
c
      real*8 :: elonpdfi(n,m+1),elatpdfi(n+1,m)
      real*8 :: erllpdfi(n+1,m+1)
c
      real*8 :: blon0(n+1,m),blat0(n,m+1),brll0(n,m)
      real*8 :: blon1(n+1,m),blat1(n,m+1),brll1(n,m)
c
      return
      end
