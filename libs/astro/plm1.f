c
c  plm1.f                            ~soi/(version)/src/functions/plm1.f
c
c  Description:
c     source: Jesper Schou
c     sets plm(i,indx(l))=plm(x(i),l) for l=lmin,lmax and i=1,nx
c     the plm's are normalized such that int from -1 to 1 plm sup 2 dx=1
c     all 0<=m<=l<=1000 tested with nx=5001, max err always at m=0
c     all 0<=m<=10 and m<=l<=2500 tested with nx=12501
c     every 50 l and m with 0<=m<=l<=3000 tested with nx=10001
c     fails at l \approx 1800 and m=700 because of underflow in setting
c     up P_m^m
c     if lmin>m then the Plm's for l=lmin-2 (if applicable) and
c     l=lmin-1 must be set on input
c
c  Responsible:  Kay Leibrand                   KLeibrand@solar.Stanford.EDU
c
c  Bugs:
c     contains writes
c
c  Revision history is at end of file
c


      subroutine setplm(lmin,lmax,m,nx,indx,x,nplm,plm)


      implicit double precision (a-h,o-z)

      parameter (nmax=12501,eps=1.0D-12,lbig=10000)

c     lbig is used not to get indexing errors when debugging

      double precision x(nx),plm(nplm,0:lbig)
      double precision x1(nmax)
      integer indx(0:lbig)

      if (nx.gt.nmax) then
        write (6,*) 'setplm only dimensioned for nmax=',nmax
        write (6,*) 'but called with nx=',nx
        stop
      end if
      if (lmax.gt.1800) then
        write (6,*) 'setplm called with lmax larger than 1800'
        stop
      end if
c     lmin must be >= m
      lmin=max(m,lmin)
      do 10,i=1,nx
        x1(i)=min(1.0D0-eps,max(-1.0D0+eps,x(i)))
 10   continue
      if (lmin.le.(m+1)) then
        im=indx(m)
        c=sqrt((2*m+1)/2.0D0)
        do 15,i=1,m
          c=-c*sqrt(1-0.5D0/i)
 15     continue
      end if
      if (lmin.eq.m) then
        do 20,i=1,nx
          plm(i,im)=c*sqrt(1.0D0-x1(i)*x1(i))**m
 20     continue
      end if
      if ((lmax.gt.m).and.(lmin.le.(m+1))) then
        im1=indx(m+1)
        c=sqrt(2*m+3.0D0)
        do 50,i=1,nx
          plm(i,im1)=x1(i)*c*plm(i,im)
 50     continue
      end if
      m2=m**2
c*$* assert do (serial)
      do 110,l0=max(lmin,m+2),lmax
        l02=l0*l0
        il0=indx(l0)
        il1=indx(l0-1)
        il2=indx(l0-2)
        c1=sqrt((4*l02-1.0D0)/(l02-m2))
c       c2=sqrt(((2*l0+1.0D0)*(l0+m-1)*(l0-m-1))/((2*l0-3)*(l02-m2)))
c changed to double precision to avoid integer overflow
        c2=sqrt(((2*l0+1.0D0)*(l0+m-1)*(l0-m-1))/
     c          ((2*l0-3.0D0)*(l02-m2)))
c*$* assert do prefer (concurrent)
        do 100,i=1,nx
          plm(i,il0)=c1*x1(i)*plm(i,il1)-c2*plm(i,il2)
 100    continue
 110  continue
      end

      subroutine dscopy(n,a,b)

c utility routine to copy double precision array to single precision
 
      double precision a(n)
      real b(n)
 
c$doacross if (n.ge.2000) local(i),share(a,b)
      do i=1,n
        b(i)=a(i)
      end do
 
      end

c $Id: plm1.f,v 1.1 2009/07/15 02:30:55 tplarson Exp $
c $Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/libs/astro/Attic/plm1.f,v $
c $Author: tplarson $
c $Log: plm1.f,v $
c Revision 1.1  2009/07/15 02:30:55  tplarson
c functions needed for jqdotprod
c
c Revision 1.3  1997/05/15  20:08:55  schou
c Added dscopy utility routine.
c
c Revision 1.2  1994/10/26  21:37:08  kay
c added comments
c
