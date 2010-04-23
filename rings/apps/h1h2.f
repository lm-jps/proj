      subroutine h1h2(n, lal, om, dom, err, w, nmode, no, nw, nsol, 
     @    filcsq, filsw, filout, qverb)
c This version determines the frequency and w limits
c from the data. However, the frequency is still forced to 
c be > 2 mHz and < 4.5 mHz.
       implicit real*8(a-h,o-z)
       external fin
       character*100 fildif, filcsq, filsw, filout
       parameter(nk=4600)

c---------------------------------------
       dimension radfw(500), bb(nk)
       dimension fb(10)
       dimension csol(400), esol(400)
       dimension n(nmode),lal(nmode),om(nmode),dom(nmode),err(nmode)
       dimension w(nmode)

C      COMMON BLOCKS
       common/int/aa,xw(100),cw(4,100),b(nk)
       common/into/xo(100),co(4,100)
       common/inp/swi(nk)
       common/svdinp/ai(nk,100),bi(nk),m
       common/svdout/a(nk,100), v(nk,100)
       common/solr/ct,radsol(400)
       common/inpc/r(6000),c(6000),rs,np
       common/fw/sw(500),wf(500),ntabb
       common/finblk/nwi

C      DEFAULT VALUES OF PARAMETERS
c-------------------------------------
c      constants
       data pi,rs/3.1415926535d0,6.9599d10/
c      number of freq. and w knots
C       data no,nw/6,10/
c      number of solution points
C       data nsol/50/
c      reference model file-names
C       data filcsq,filsw/'csq.dat','fw.dat'/
c      input data file name
       data fildif/'freqdif.in'/
c      output file name
C       data filout/'inv.out'/
c      number of errorrealizations for MC
       data nmc/500/

51       format(' spline coef:'/(1p5e12.3))
52       format(1p9e15.7)
57       format(2i6,1p9e15.7)

c ----- END DECLARATIVE STATEMENTS ----------

C       read(*,*)fildif
C       read(*,*)filcsq
C       read(*,*)filsw
C       read(*,*)filout
c       print*,'type number of frequency knots'
C       read(*,*)no
       if(no.lt.2)no=2
c       print*,'type number of w knots'
C       read(*,*)nw
       if(nw.lt.2)nw=2
c      print*,'type number of points where solution needed'
C       read(*,*)nsol
       if(nsol.gt.400)nsol=400
       if(nsol.lt.2)nsol=2

C       open(unit=8,file=fildif)
       open(unit=11,file=filsw)
       open(unit=12, file=filcsq)
       open(unit=19,file=filout)

       nwi=nw

c-------------------------
c	READ INPUT MODEL QUANTITIES
cThis part reads in the two mdoel files, one with the
c sound speed and the other with w, r, and S(w)
       np=0
       do 800 i=1,6000
            read(12,*,end=900) rad,cs
            np=np+1
            r(np)=rad*rs
            c(np)=cs/r(np)
800       continue
900       continue
       ct=c(1)
       ntab=0
       do 1500 i=1,500
            read(11,*,end=1501)wff,rr, smh
            ntab=ntab+1
            sw(ntab)=smh
            radfw(ntab)=rr
            wf(ntab)=wff
1500       continue
1501       continue
       ntabb=ntab
       reps=1.d-7

c-------------------------
c      READ DATA
c Reading in the data. Assumes ferquency, differences
c and errors are in microHerz. 
       omax=0.d0
       omin=10.d0
       wmax=-4.d0
       wmin=4.d0
       do 2000 i=1,nmode
         if(om(i).gt.omax)omax=om(i)
         if(om(i).lt.omin)omin=om(i)
         if(w(i).gt.wmax)wmax=w(i)
         if(w(i).lt.wmin)wmin=w(i)
2000    continue
C       nmode=0
C       do 2000 i=1,6000
C       read(8,*,end=2200)li, ni, fi, dfi, erri
C       omi=fi/1000.d0
C       if(omi.gt.4.5d0.or.omi.lt.2.0d0)go to 1999
C       llll=li
C       dif=dfi/1000.d0
C         if(omi.gt.omax)omax=omi
C         if(omi.lt.omin)omin=omi
C       wi=dlog10(omi/(llll+0.5))
C         if(wi.gt.wmax)wmax=wi
C         if(wi.lt.wmin)wmin=wi
C       omo=omi
C       nmode=nmode+1
C           if(nmode.gt.nk)then
C             print*, 'cannot read all data'
C             print*, 'increase nk to read all data'
C             go to 2200
C           endif
C       n(nmode)=ni
C       lal(nmode)=llll
C       dom(nmode)=dif
C       om(nmode)=omi
C       err(nmode)=erri/1000.d0
C       w(nmode)=wi
C1999       continue
C2000       continue
C2200       continue

c------------------------------
c      SET KNOTS
       oe=omax
       ob=omin
       do=(oe-ob)/(no-1.)
       do 1000 i=1,no
            xo(i)=ob+(i-1)*do
1000       continue
       we=wmax
       wb=wmin
       dw=(we-wb)/(nw-1.)
       do 1200 i=1,nw
            xw(i)=wb+(i-1)*dw
1200       continue
c       print*,'ob,oe ', ob,oe
c       print*,'wb,we ', wb,we

       ic=4
       call bspcof(no,xo,co,ic)
       call bspcof(nw,xw,cw,ic)
       m=nw+no+4
c--------------------------------------

c SET UP SOLUTION LIMITS, and array at which solution will be found
       nuse=4
      call DIVDIF(xw(1),wf,radfw,NUSE,NTABb,FB,REPS,IER,DFB,DDFB)
       rad1=fb(nuse)
       nuse=4
      call DIVDIF(xw(nw),wf,radfw,NUSE,NTABb,FB,REPS,IER,DFB,DDFB)
       rad2=fb(nuse)
       if(rad2.lt.0.975d0)rad2=0.975d0
       if(rad1.gt.0.997d0)rad1=0.997d0
       dr=(rad1-rad2)/(nsol-1.)
       do i=1,nsol
           radsol(i)=rad2+(i-1)*dr
       enddo

c--------------------------------------
c      SET UP COMMON SVD MATRICES
       CALL SETMATRIX(n,lal,om,dom,err,w,nmode,nw,no,nsol)
c      EXACT SOLUTION
       CALL EXACT(n,lal,om,dom,err,csol,nmode,nw,no,nsol)
C      MC exercise for errors
       CALL MC(n,lal,om,dom,err,esol,nmc,nmode,nw,no,nsol)
c-----------------------------------
c 	WRITE OUT SOLUTION AND ERRORS
       do i=1,nsol
       write(19,52)radsol(i),csol(i),esol(i)
C       print*,i,radsol(i),filout
       enddo
       close(19)

       close(11)
       close(12)

       RETURN
       END

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       SUBROUTINE SETMATRIX(n,lal,om,dom,err,w,nmode,nw,no,nsol)
c Sets the matrix for the fits. Called just once by the
c program
       implicit real*8(a-h,o-z)
       external fin
       parameter(nk=4600)
       dimension fb(10)
       dimension n(nmode),lal(nmode),om(nmode),dom(nmode),err(nmode)
       dimension w(nmode)
       common/int/aa,xw(100),cw(4,100),b(nk)
       common/into/xo(100),co(4,100)
       common/inp/swi(nk)
       common/svdinp/ai(nk,100),bi(nk),m
       common/svdout/a(nk,100), v(nk,100)
       common/solr/ct,radsol(400)
       common/inpc/r(6000),c(6000),rs,np
       common/fw/sw(500),wf(500),ntab

       reps=1.d-7

              m=nw+no+4
       do 3000 i=1,nmode
              nuse=4
             call DIVDIF(w(i),wf,sw,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
              swi(i)=fb(nuse)
              erri=err(i)*swi(i)/om(i)
              bi(i)=fb(nuse)/(erri*om(i))
              do 2400 j=1,nw+2
                     ai(i,j)=bspev(nw,xw,cw,ic,w(i),j)/erri
2400              continue
              do 2600 j=1,no+2
                     ai(i,j+nw+2)=bspev(no,xo,co,ic,om(i),j)/erri
2600              continue
3000       continue

       do 3050 i=nmode+1,nmode+4
               bi(i)=0.0
              do 3050 j=1,m
3050              ai(i,j)=0.0

       wi=xw(nw)
       wj=xw(1)
       do 3080 j=1,nw+2
              ai(nmode+1,j)=bspev(nw,xw,cw,ic,wi,j)
              ai(nmode+4,j)=bspd2(nw,xw,cw,ic,wj,j)
3080   continue
       oi=xo(no)
       oj=xo(1)
       do 3090 j=1,no+2
              ai(nmode+2,j+nw+2)=bspd2(no,xo,co,ic,oi,j)
              ai(nmode+3,j+nw+2)=bspd2(no,xo,co,ic,oj,j)
3090       continue

       return
       end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       SUBROUTINE EXACT(n,lal,om,dom,err,csol,nmode,nw,no,nsol)
       implicit real*8(a-h,o-z)
        parameter(nk=4600)
       dimension csol(nsol), bb(nk)
       dimension n(nmode),lal(nmode),om(nmode),dom(nmode),err(nmode)
c      ---------------
       common/int/aa,xw(100),cw(4,100),b(nk)
       common/into/xo(100),co(4,100)
       common/inp/swi(nk)
       common/svdinp/ai(nk,100),bi(nk),m
       common/svdout/a(nk,100), v(nk,100)
       common/solr/ct,radsol(400)
       common/inpc/r(6000),c(6000),rs,np
c       ---------------
       data pi/3.1415926535d0/

       do 200 i=1,nmode
       bb(i)=bi(i)*dom(i)
200    continue
       do 300 i=nmode+1,nmode+4
        bb(i)=0.0
300    continue

       CALL DECOMP(bb,nmode,nw,no,nsol)
       CALL SOLUTION (csol,iera,nmode,nsol)

            if(iera.gt.0)then
            print*, 'error in exact solution'
            stop
            endif

       return
       end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


       SUBROUTINE MC(n,lal,om,dom,err,esol,nmc,nmode,nw,no,nsol)
       implicit real*8(a-h,o-z)
       external fin
       parameter(nk=4600)
       dimension bb(nk)
       dimension n(nmode),lal(nmode),om(nmode),dom(nmode),err(nmode)
       dimension sumc(nsol),sumcc(nsol),esol(nsol), dsol(nsol)
c      ---------------
       common/int/aa,xw(100),cw(4,100),b(nk)
       common/into/xo(100),co(4,100)
       common/inp/swi(nk)
       common/svdinp/ai(nk,100),bi(nk),m
       common/svdout/a(nk,100), v(nk,100)
       common/solr/ct,radsol(400)
       common/inpc/r(6000),c(6000),rs,np
c      ---------------
       data pi/3.1415926535d0/

       do j=1,nsol
       sumc(j)=0.
       sumcc(j)=0.
       enddo

c      random number key
       idum=-20
c      loop over number of realizations
       numm=0
       do 100 kk=1,nmc

c     first get the new "data"
       do 200 i=1,nmode
       dd=dom(i)+gasdev(idum)*err(i)
       bb(i)=bi(i)*dd
c       write(33,*)n(i),lal(i),dom(i),dd
200    continue
       do 300 i=nmode+1,nmode+4
        bb(i)=0.0
300    continue
       CALL DECOMP(bb,nmode,nw,no,nsol)
       CALL SOLUTION (dsol,iera,nmode,nsol)

c      now populate the sun amd sum of squares
       if(iera.eq.0)then
       numm=numm+1
       do 400 i=1,nsol
       ds=dsol(i)
       sumc(i)=sumc(i)+ds
       sumcc(i)=sumcc(i)+ds*ds
400    continue
       endif

100    continue
c      now find average and sigma and populate returned array

       do 500 i=1,nsol
       av=sumc(i)/numm
       sq=sumcc(i)/numm
       sig=sq-av*av
       sig=sqrt(sig)
       esol(i)=sig*1.5d0
500    continue
53     format(1p3e13.5)

       return
       end



c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       SUBROUTINE DECOMP(bb,nmode,nw,no,nsol)
c       decomposes data into H_1 and H_2
       implicit real*8(a-h,o-z)
       external fin
       parameter(nk=4600)
       dimension sigma(100),e(100),wk(nk)
       dimension bb(nk)
       common/int/aa,xw(100),cw(4,100),b(nk)
       common/into/xo(100),co(4,100)
       common/inp/swi(nk)
       common/svdinp/ai(nk,100),bi(nk),m
       common/svdout/a(nk,100), v(nk,100)



       do 3000 i=1,nmode
       b(i)=bb(i)

       do 2400 j=1,nw+2
       a(i,j)=ai(i,j)
2400       continue
       do 2600 j=1,no+2
       a(i,j+nw+2)=ai(i,j+nw+2)
2600       continue

3000       continue

       do 3050 i=nmode+1,nmode+4
        b(i)=bb(i)
       do 3050 j=1,m
3050       a(i,j)=ai(i,j)

       do 3080 j=1,nw+2
       a(nmode+1,j)=ai(nmode+1,j)
c       a(nmode+4,j)=ai(nmode+2,j)
3080       continue
       
c       do 3090 j=1,no+2
c       a(nmode+2,j+nw+2)=ai(nmode+3,j+nw+2)
c       a(nmode+3,j+nw+2)=ai(nmode+4,j+nw+2)
c3090       continue


       neq=nmode+1
       la=nk
       lv=nk       
       reps=1.e-14
c       print *,' call svd'
      call SVD(m,neq,A,V,SIGMA,LA,LV,E,REPS,IER)
       if(ier.gt.0) print *,'err in svd, ier=',ier
c       print *,(sigma(i),i=1,m)
      call SVDEVL(m,nmode,a,V,SIGMA,La,LV,B,WK,REPS)


       RETURN
       END


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE SOLUTION(sound,iera,nmode,nsol)
c This routine does the real inversion.
        implicit real*8(a-h,o-z)
       external fin
       parameter(nk=4600)
       dimension sound(nsol), fb(10)
       common/int/aa,xw(100),cw(4,100),b(nk)
       common/into/xo(100),co(4,100)
       common/inp/swi(nk)
       common/svdinp/ai(nk,100),bi(nk),m
       common/svdout/a(nk,100), v(nk,100)
       common/solr/ct,radsol(400)
       common/inpc/r(6000),c(6000),rs,np
       data pi/3.1415926535d0/

c       sound speed inversion
       reps=1.e-7
       aeps=1.e-10
       ntab=np

c       print *,' surface a',ct

       iera=0
       nmax=0
       do 6000 i=1,nsol
       rr=radsol(i)*rs
       nuse=4
      call DIVDIF(rr,r,c,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
       aa=fb(nuse)
       apr=dfb
       xl=ct
       xu=aa
      call ADPINT(RIN,XL,XU,REPS,AEPS,DIF,Fin,IER,NPT,NMAX)
       if(ier.gt.0)then
       iera=ier
       endif
       delc=-2.*rr*apr*rin/pi
       sound(i)=delc
6000       continue

       end


C       -------------------------------------------------

       function fin(w)
c  The function being integrated.
       implicit real*8(a-h,o-z)
       parameter(nk=4600)
       common/int/aa,xw(100),cw(4,100),b(nk)
       common/finblk/nwi
       data pi/3.1415926535d0/

       nw=nwi
       ic=4
       x0=dlog10(500.*w/pi)
       h1= bspd1(nw,xw,cw,ic,b,x0)/2.3025851d0
       if(aa.gt.w) then
       fin=h1/sqrt(aa**2-w**2)
       else
       fin=0.0
       endif
       end

C       -------------------------------------------------



       function bspev(n,x,c,ic,x0,i)
       implicit real*8(a-h,o-z)
       dimension x(n+3),c(ic,n+2)

       if(i.gt.3) then
       a0=x(i-3)
       else
       a0=x(1)-(4-i)*(x(2)-x(1))
       endif
       if(i.gt.2) then
       a1=x(i-2)
       else
       a1=x(1)-(3-i)*(x(2)-x(1))
       endif
       if(i.gt.1) then
       a2=x(i-1)
       else
       a2=2.*x(1)-x(2)
       endif
       a3=x(i)
       a4=x(i+1)
       if(x0.le.a0) then
       bspev=0.0
       else if(x0.le.a1) then
       bspev=c(1,i)*(x0-a0)**3
       else if(x0.le.a2) then
       bspev=c(3,i)*(x0-a1)**3+c(1,i)*(a1-a0)*(3.*(x0-a1)**2
     1   +(a1-a0)*(3.*(x0-a1)+(a1-a0)))
       else if(x0.le.a3) then
       bspev=c(4,i)*(x0-a3)**3+c(2,i)*(a3-a4)*(3.*(x0-a3)**2
     1   +(a3-a4)*(3.*(x0-a3)+(a3-a4)))
       else if(x0.le.a4) then
       bspev=c(2,i)*(x0-a4)**3
       else
       bspev=0.0
       endif
       end

C       -------------------------------------------------

       function bspde1(n,x,c,ic,x0,i)
       implicit real*8(a-h,o-z)
       dimension x(n+3),c(ic,n+2)

       if(i.gt.3) then
       a0=x(i-3)
       else
       a0=x(1)-(4-i)*(x(2)-x(1))
       endif
       if(i.gt.2) then
       a1=x(i-2)
       else
       a1=x(1)-(3-i)*(x(2)-x(1))
       endif
       if(i.gt.1) then
       a2=x(i-1)
       else
       a2=2.*x(1)-x(2)
       endif
       a3=x(i)
       a4=x(i+1)
       if(x0.le.a0) then
       bspde1=0.0
       else if(x0.le.a1) then
       bspde1=3.*c(1,i)*(x0-a0)**2
       else if(x0.le.a2) then
       bspde1=3.*c(3,i)*(x0-a1)**2+c(1,i)*(a1-a0)*(6.*(x0-a1)
     1   +3.*(a1-a0))
       else if(x0.le.a3) then
       bspde1=3.*c(4,i)*(x0-a3)**2+c(2,i)*(a3-a4)*(6.*(x0-a3)
     1   +3.*(a3-a4))
       else if(x0.le.a4) then
       bspde1=3.*c(2,i)*(x0-a4)**2
       else
       bspde1=0.0
       endif
       end

C       -------------------------------------------------

       function bspd2(n,x,c,ic,x0,i)
       implicit real*8(a-h,o-z)
       dimension x(n+3),c(ic,n+2)

       if(i.gt.3) then
       a0=x(i-3)
       else
       a0=x(1)-(4-i)*(x(2)-x(1))
       endif
       if(i.gt.2) then
       a1=x(i-2)
       else
       a1=x(1)-(3-i)*(x(2)-x(1))
       endif
       if(i.gt.1) then
       a2=x(i-1)
       else
       a2=2.*x(1)-x(2)
       endif
       a3=x(i)
       a4=x(i+1)
       if(x0.le.a0) then
       bspd2=0.0
       else if(x0.le.a1) then
       bspd2=6.*c(1,i)*(x0-a0)
       else if(x0.le.a2) then
       bspd2=6.*c(3,i)*(x0-a1)+6.*c(1,i)*(a1-a0)
       else if(x0.le.a3) then
       bspd2=6.*c(4,i)*(x0-a3)+6.*c(2,i)*(a3-a4)
       else if(x0.le.a4) then
       bspd2=6.*c(2,i)*(x0-a4)
       else
       bspd2=0.0
       endif
       end

C       -------------------------------------------------

       subroutine bspcof(n,x,c,ic)
       implicit real*8(a-h,o-z)
       dimension c(ic,n+2),a(4,4),int(4),x(n+3)

       x(n+1)=2.*x(n)-x(n-1)
       x(n+2)=3.*x(n)-2.*x(n-1)
       x(n+3)=4.*x(n)-3.*x(n-1)
       do 1000 i=1,n+2
       if(i.gt.3) then
       a10=x(i-2)-x(i-3)
       else
       a10=x(2)-x(1)
       endif
       a34=x(i)-x(i+1)
       if(i.gt.2) then
       a21=x(i-1)-x(i-2)
       else
       a21=x(2)-x(1)
       endif
       if(i.gt.1) then
       a23=x(i-1)-x(i)
       else
       a23=x(1)-x(2)
       endif
       a(1,1)=a10*(a10**2+3.*a10*a21+3.*a21**2)
       a(1,2)=0.0
       a(1,3)=a21**3
       a(1,4)=0.0
       c(1,i)=1.0
       a(2,1)=0.0
       a(2,2)=a34*(a34**2+3.*a34*a23+3.*a23**2)
       a(2,3)=0.0
       a(2,4)=a23**3
       c(2,i)=1.0
       a(3,1)=a10*(a10+2.*a21)
       a(3,2)=-a34*(a34+2.*a23)
       a(3,3)=a21**2
       a(3,4)=-a23**2
       c(3,i)=0.0
       a(4,1)=a10
       a(4,2)=-a34
       a(4,3)=a21
       a(4,4)=-a23
       c(4,i)=0.0

       m=4
       num=1
       lj=4
       iflg=0
      call GAUELM(m,NUM,A,c(1,i),DET,INT,LJ,IER,IFLG)
       if(ier.gt.0) print *,'error in gauelm, ier=',ier,x(i)
1000       continue
       end

C       -------------------------------------------------

       function bspevl(n,x,c,ic,wt,x0)
       implicit real*8(a-h,o-z)
       dimension x(n),c(ic,n),wt(n)

       bspevl=0.0
       do 1000 i=1,n+2
       if(i.gt.3) then
       a0=x(i-3)
       else
       a0=x(1)-(4-i)*(x(2)-x(1))
       endif
       if(i.gt.2) then
       a1=x(i-2)
       else
       a1=x(1)-(3-i)*(x(2)-x(1))
       endif
       if(i.gt.1) then
       a2=x(i-1)
       else
       a2=2.*x(1)-x(2)
       endif
       a3=x(i)
       a4=x(i+1)
       if(x0.le.a0) then
       bsp=0.0
       else if(x0.le.a1) then
       bsp=c(1,i)*(x0-a0)**3
       else if(x0.le.a2) then
       bsp=c(3,i)*(x0-a1)**3+c(1,i)*(a1-a0)*
     1     (3.*(x0-a1)**2+(a1-a0)*(3.*(x0-a1)+(a1-a0)))
       else if(x0.le.a3) then
       bsp=c(4,i)*(x0-a3)**3+c(2,i)*(a3-a4)*
     1      (3.*(x0-a3)**2+(a3-a4)*(3.*(x0-a3)+(a3-a4)))
       else if(x0.le.a4) then
       bsp=c(2,i)*(x0-a4)**3
       else
       bsp=0.0
       endif
       bspevl=bspevl+wt(i)*bsp
1000       continue

       end

C       -------------------------------------------------

       function bspd1(n,x,c,ic,wt,x0)
       implicit real*8(a-h,o-z)
       dimension x(n),c(ic,n),wt(n)

       bspd1=0.0
       do 1000 i=1,n+2
       if(i.gt.3) then
       a0=x(i-3)
       else
       a0=x(1)-(4-i)*(x(2)-x(1))
       endif
       if(i.gt.2) then
       a1=x(i-2)
       else
       a1=x(1)-(3-i)*(x(2)-x(1))
       endif
       if(i.gt.1) then
       a2=x(i-1)
       else
       a2=2.*x(1)-x(2)
       endif
       a3=x(i)
       a4=x(i+1)
       if(x0.le.a0) then
       bsp=0.0
       else if(x0.le.a1) then
       bsp=3.*c(1,i)*(x0-a0)**2
       else if(x0.le.a2) then
       bsp=3.*c(3,i)*(x0-a1)**2+c(1,i)*(a1-a0)*
     1     (6.*(x0-a1)+3.*(a1-a0))
       else if(x0.le.a3) then
       bsp=3.*c(4,i)*(x0-a3)**2+c(2,i)*(a3-a4)*
     1      (6.*(x0-a3)+3.*(a3-a4))
       else if(x0.le.a4) then
       bsp=3.*c(2,i)*(x0-a4)**2
       else
       bsp=0.0
       endif
       bspd1=bspd1+wt(i)*bsp
1000       continue

       end

C       -------------------------------------------------

      SUBROUTINE GAUELM(N,NUM,A,X,DET,INT,LJ,IER,IFLG)
       implicit real*8(a-h,o-z)
      DIMENSION A(LJ,N),INT(N),X(LJ,NUM)

      IF(N.LE.0.OR.N.GT.LJ) THEN
        IER=111
        RETURN
      ENDIF

      IER=122
      IF(IFLG.LE.1) THEN
        DET=1.0
        DO 2600 K=1,N-1
          R1=0.0
          DO 2200 L=K,N
            IF(ABS(A(L,K)).GT.R1) THEN
              R1=ABS(A(L,K))
              KM=L
            ENDIF
2200      CONTINUE

          INT(K)=KM
          IF(KM.NE.K) THEN
            DO 2300 L=K,N
              T1=A(K,L)
              A(K,L)=A(KM,L)
2300        A(KM,L)=T1
            DET=-DET
          ENDIF

          DET=DET*A(K,K)
          IF(A(K,K).EQ.0.0) RETURN
C         IF(ABS(A(K,K)).LT.REPS) RETURN
          DO 2500 L=K+1,N
            A(L,K)=A(L,K)/A(K,K)
            DO 2500 L1=K+1,N
2500      A(L,L1)=A(L,L1)-A(L,K)*A(K,L1)
2600    CONTINUE
        DET=DET*A(N,N)
        INT(N)=N
        IF(A(N,N).EQ.0.0) RETURN
C         IF(ABS(A(N,N)).LT.REPS) RETURN

        IER=0
        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
        IFLG=2
      ENDIF

      IER=0
      DO 5000 J=1,NUM
        DO 3000 K=1,N-1
          IF(K.NE.INT(K)) THEN
            T1=X(K,J)
            X(K,J)=X(INT(K),J)
            X(INT(K),J)=T1
          ENDIF
          DO 3000 L=K+1,N
3000    X(L,J)=X(L,J)-A(L,K)*X(K,J)

        X(N,J)=X(N,J)/A(N,N)
        DO 3300 K=N-1,1,-1
          DO 3200 L=N,K+1,-1
3200      X(K,J)=X(K,J)-X(L,J)*A(K,L)
3300    X(K,J)=X(K,J)/A(K,K)
5000  CONTINUE
      END

C       --------------------------------------

      SUBROUTINE SVD(N,M,A,V,SIGMA,LA,LV,E,REPS,IER)
       implicit real*8(a-h,o-z)
      PARAMETER(ITMAX=30)
      DIMENSION A(LA,N),V(LV,N),SIGMA(N),E(N)

      IF(N.GT.M.OR.N.LE.0.OR.M.LE.0.OR.M.GT.LA.OR.N.GT.LV) THEN
        IER=111
        RETURN
      ENDIF

      IER=0
      G=0
      RMAX=0

      DO 3000 I=1,N
        E(I)=G
        S=0
        DO 1200 J=I,M
1200    S=S+A(J,I)**2

        IF(S.LE.0.0) THEN
          G=0
        ELSE
          F=A(I,I)
          G=SQRT(S)
          IF(F.GE.0.0) G=-G
          H=F*G-S
          A(I,I)=F-G

          DO 1800 J=I+1,N
            S=0
            DO 1400 K=I,M
1400        S=S+A(K,I)*A(K,J)
            F=S/H
            DO 1600 K=I,M
1600        A(K,J)=A(K,J)+F*A(K,I)
1800      CONTINUE
        ENDIF

        SIGMA(I)=G
        S=0
        DO 2000 J=I+1,N
2000    S=S+A(I,J)**2

        IF(S.LE.0.0) THEN
          G=0
        ELSE
          F=A(I,I+1)
          G=SQRT(S)
          IF(F.GE.0.0) G=-G
          H=F*G-S
          A(I,I+1)=F-G
          DO 2200 J=I+1,N
2200      E(J)=A(I,J)/H

          DO 2800 J=I+1,M
            S=0
            DO 2400 K=I+1,N
2400        S=S+A(J,K)*A(I,K)
            DO 2600 K=I+1,N
2600        A(J,K)=A(J,K)+S*E(K)
2800      CONTINUE
        ENDIF
        R1=ABS(SIGMA(I))+ABS(E(I))
        IF(R1.GT.RMAX) RMAX=R1
3000  CONTINUE

c       print *,'3000 over'
      DO 4000 I=N,1,-1
        IF(G.NE.0) THEN
          H=A(I,I+1)*G
          DO 3200 J=I+1,N
3200      V(J,I)=A(I,J)/H

          DO 3800 J=I+1,N
            S=0
            DO 3400 K=I+1,N
3400        S=S+A(I,K)*V(K,J)
            DO 3600 K=I+1,N
3600        V(K,J)=V(K,J)+S*V(K,I)
3800      CONTINUE
        ENDIF

        DO 3900 J=I+1,N
          V(I,J)=0.0
          V(J,I)=0.0
3900    CONTINUE
        V(I,I)=1
        G=E(I)
4000  CONTINUE

c       print *,'4000 over'
      DO 5000 I=N,1,-1
        G=SIGMA(I)
        DO 4200 J=I+1,N
4200    A(I,J)=0
        IF(G.NE.0.0) THEN
          H=A(I,I)*G

          DO 4700 J=I+1,N
            S=0
            DO 4400 K=I+1,M
4400        S=S+A(K,I)*A(K,J)
            F=S/H
            DO 4600 K=I,M
4600        A(K,J)=A(K,J)+F*A(K,I)
4700      CONTINUE

          DO 4800 J=I,M
4800      A(J,I)=A(J,I)/G
        ELSE
          DO 4900 J=I,M
4900      A(J,I)=0.0
        ENDIF
        A(I,I)=A(I,I)+1
5000  CONTINUE

c       print *,'5000 over'
      AEPS=REPS*RMAX
      DO 8000 K=N,1,-1
c       print *,k
        DO 7500 ITR=1,ITMAX
c       print *,itr

          DO 5200 Ll=K,1,-1
       l=ll
            IF(ABS(E(L)).LT.AEPS) GO TO 6000
            IF(ABS(SIGMA(L-1)).LT.AEPS) GO TO 5400
5200      CONTINUE

5400      C=0.0
          S=1.0
c       print *,l,k
          DO 5800 I=L,K
c       print *,i
            F=S*E(I)
            E(I)=C*E(I)
            IF(ABS(F).LT.AEPS) GO TO 6000
            G=SIGMA(I)
            SIGMA(I)=SQRT(F*F+G*G)
            C=G/SIGMA(I)
            S=-F/SIGMA(I)

            DO 5600 J=1,M
              R1=A(J,L-1)
              R2=A(J,I)
              A(J,L-1)=R1*C+R2*S
              A(J,I)=C*R2-S*R1
5600        CONTINUE
5800      CONTINUE

6000      Z=SIGMA(K)
c       print *,k,z
          IF(L.EQ.K) THEN
            IF(Z.LT.0.0) THEN
              SIGMA(K)=-Z
              DO 6200 J=1,N
6200          V(J,K)=-V(J,K)
            ENDIF
            GO TO 8000
          ENDIF

          IF(ITR.EQ.ITMAX) THEN
            IER=11
            GO TO 7500
          ENDIF

          X=SIGMA(L)
          Y=SIGMA(K-1)
          G=E(K-1)
          H=E(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.*H*Y)
          G=SQRT(1.+F*F)
          IF(F.LT.0.0) G=-G
          F=((X-Z)*(X+Z)+H*(Y/(F+G)-H))/X

          C=1.0
          S=1.0
          DO 7000 I=L+1,K
c       print *,i,l,k
            G=E(I)
            Y=SIGMA(I)
            H=S*G
            G=C*G
            E(I-1)=SQRT(F*F+H*H)
            C=F/E(I-1)
            S=H/E(I-1)
            F=C*X+S*G
            G=C*G-S*X
            H=S*Y
            Y=C*Y

            DO 6400 J=1,N
              X=V(J,I-1)
              Z=V(J,I)
              V(J,I-1)=C*X+S*Z
              V(J,I)=C*Z-S*X
6400        CONTINUE

            SIGMA(I-1)=SQRT(F*F+H*H)
            IF(SIGMA(I-1).NE.0.0) THEN
              C=F/SIGMA(I-1)
              S=H/SIGMA(I-1)
            ENDIF
            F=C*G+S*Y
            X=C*Y-S*G
            DO 6600 J=1,M
              Y=A(J,I-1)
              Z=A(J,I)
              A(J,I-1)=C*Y+S*Z
              A(J,I)=C*Z-S*Y
6600        CONTINUE
7000      CONTINUE

          E(L)=0
          E(K)=F
          SIGMA(K)=X
7500    CONTINUE
8000  CONTINUE
      END

C       --------------------------------------

      SUBROUTINE SVDEVL(N,M,U,V,SIGMA,LU,LV,B,WK,REPS)
       implicit real*8(a-h,o-z)
      DIMENSION U(LU,N),V(LV,N),SIGMA(N),B(*),WK(N)

      SMAX=0.0
      DO 2000 I=1,N
        IF(SIGMA(I).GT.SMAX) SMAX=SIGMA(I)
2000  CONTINUE

      AEPS=SMAX*REPS
      DO 3000 I=1,N
        S=0.0
        IF(SIGMA(I).GT.AEPS) THEN
          DO 2400 J=1,M
2400      S=S+U(J,I)*B(J)
          S=S/SIGMA(I)
        ENDIF
        WK(I)=S
3000  CONTINUE

      DO 4000 I=1,N
        S=0.0
        DO 3400 J=1,N
3400    S=S+V(I,J)*WK(J)
        B(I)=S
4000  CONTINUE
      END

C       ---------------------------------------------

      SUBROUTINE DIVDIF(XB,X,F,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
       implicit real*8(a-h,o-z)
      PARAMETER(NMAX=10)
      DIMENSION X(NTAB),F(NTAB),FB(*),XN(NMAX),XD(NMAX)

      NEXT=NEARST(XB,X,NTAB)
      FB(1)=F(NEXT)
      XD(1)=F(NEXT)
      XN(1)=X(NEXT)
      IER=0
      PX=1.0

      DFB=0.0
      DDFB=0.0
      DPX=0.0
      DDPX=0.0

      IP=NEXT
      IN=NEXT

      NIT=MIN0(NMAX,NUSE,NTAB)
      IF(NUSE.GT.NMAX.OR.NUSE.GT.NTAB) IER=12
      IF(NUSE.LT.1) THEN
        IER=11
        NIT=MIN0(6,NTAB,NMAX)
      ENDIF
      NUSE=1

      DO 5000 J=2,NIT

        IF(IN.LE.1) GO TO 2200
        IF(IP.GE.NTAB) GO TO 2000
        IF(ABS(XB-X(IP+1)).LT.ABS(XB-X(IN-1))) GO TO 2200
2000    IN=IN-1
        NEXT=IN
        GO TO 2800
2200    IP=IP+1
        NEXT=IP

2800    XD(J)=F(NEXT)
        XN(J)=X(NEXT)
        DO 3000 K=J-1,1,-1
3000    XD(K)=(XD(K+1)-XD(K))/(XN(J)-XN(K))

        DDPX=DDPX*(XB-XN(J-1))+2.*DPX
        DPX=DPX*(XB-XN(J-1))+PX
        DFB=DFB+DPX*XD(1)
        DDFB=DDFB+DDPX*XD(1)

        PX=PX*(XB-XN(J-1))
        ERR=XD(1)*PX
        FB(J)=FB(J-1)+ERR
        NUSE=J

        IF(ABS(ERR).LT.REPS) RETURN
5000  CONTINUE

      IER=24
      END

c       -------------------------------------------

      FUNCTION NEARST(XB,X,NTAB)
       implicit real*8(a-h,o-z)
      DIMENSION X(NTAB)

      LOW=1
      IGH=NTAB
      IF(.NOT.(XB.LT.X(LOW).EQV.XB.LT.X(IGH))) THEN


1500    IF(IGH-LOW.GT.1) THEN
          MID=(LOW+IGH)/2
          IF(XB.LT.X(MID).EQV.XB.LT.X(LOW)) THEN
            LOW=MID
          ELSE
            IGH=MID
          ENDIF
          GO TO 1500
        ENDIF
      ENDIF

      IF(ABS(XB-X(LOW)).LT.ABS(XB-X(IGH))) THEN
        NEARST=LOW
      ELSE
        NEARST=IGH
      ENDIF
      END

c       ---------------------------------------------------

      SUBROUTINE ADPINT(RINT,XL,XU,REPS,AEPS,DIF,F,IER,NPT,NMAX)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      PARAMETER(IPMAX=100,IFMAX=5,MAXPT=100000)
      EXTERNAL F
      DIMENSION XU1(IPMAX)
      SAVE

      IFAIL=0
      RINT=0.0
      DIF=0.0
      IF(XL.EQ.XU) RETURN
      IF(NMAX.LE.0) NMAX=MAXPT
      AEPSL=AEPS
      IER=0
      NPT=0
      RL=XL
      RU=XU
      IU=0

1000  CALL KRONRD(FINT,RL,RU,DIF0,NP,F)
      NPT=NPT+NP
      RM=0.5*(RL+RU)
      Q=IU.GE.IPMAX.OR.RM.EQ.RL.OR.RM.EQ.RU
      IF(DIF0.LT.DMAX1(ABS(FINT)*REPS,AEPSL).OR.Q) THEN
        RINT=RINT+FINT
        DIF=DIF+DIF0
        IF(Q.AND.DIF0.GT.DMAX1(ABS(RINT)*REPS,AEPSL)) THEN
          IER=11
          IFAIL=IFAIL+1
          IF(IFAIL.GT.IFMAX) THEN
            IER=22
            AEPSL=DIF*0.5
          ENDIF
        ENDIF
        IF(ABS(RU-XU).LT.REPS.OR.IU.LE.0) RETURN
        RL=RU
        RU=XU1(IU)
        IU=IU-1
      ELSE
        IU=IU+1
        XU1(IU)=RU
        RU=RM
      ENDIF

      IF(NPT.LT.NMAX) GO TO 1000
      IER=13
      RU=XU
      CALL KRONRD(FINT,RL,RU,DIF0,NP,F)
      NPT=NPT+NP
      RINT=RINT+FINT
      DIF=DIF+DIF0
      END

c       ---------------------------------------------------

      SUBROUTINE KRONRD(RI,A,B,DIF,N,F)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  W7(4),A7(4),WK7(4),WK15(4),AK15(4)
      SAVE


      DATA W7  /0.12948496616886969327D0, 0.27970539148927666790D0,
     *          0.38183005050511894495D0, 0.41795918367346938775D0/
      DATA A7  /0.94910791234275852452D0, 0.74153118559939443986D0,
     *          0.40584515137739716690D0, 0.0/
      DATA WK7 /0.06309209262997855329D0, 0.14065325971552591874D0,
     *          0.19035057806478540991D0, 0.20948214108472782801D0/
      DATA WK15/0.02293532201052922496D0, 0.10479001032225018383D0,
     *          0.16900472663926790282D0, 0.20443294007529889241D0/
      DATA AK15/0.99145537112081263920D0, 0.86486442335976907278D0,
     *          0.58608723546769113029D0, 0.20778495500789846760D0/

      AT=(B-A)/2.
      BT=(B+A)/2.
      FBT=F(BT)
      R1=W7(4)*FBT
      RI=WK7(4)*FBT
      DO 2000 K=1,3
        F1=F(AT*A7(K)+BT)
        F2=F(BT-AT*A7(K))
        R1=R1+W7(K)*(F1+F2)
        RI=RI+WK7(K)*(F1+F2)
2000  CONTINUE

      DO 2500 K=1,4
2500  RI=RI+WK15(K)*(F(AT*AK15(K)+BT)+F(BT-AT*AK15(K)))

      RI=RI*AT
      R1=R1*AT
      DIF=ABS(RI-R1)
      N=15
      END
c       --------------
      real*8 FUNCTION GASDEV(IDUM)
       implicit real*8(a-h,o-z)
       save
      DATA ISET/0/
      IF (ISET.EQ.0) THEN
1       V1=2.*RAN1(IDUM)-1.
        V2=2.*RAN1(IDUM)-1.
        R=V1**2+V2**2
        IF(R.GE.1.)GO TO 1
        FAC=SQRT(-2.*LOG(R)/R)
        GSET=V1*FAC
        GASDEV=V2*FAC
        ISET=1
      ELSE
        GASDEV=GSET
        ISET=0
      ENDIF
      RETURN
      END

c       ---------------------------


      real*8 FUNCTION RAN1(IDUM)
       implicit real*8(a-h,o-z)
       save
      DIMENSION R(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
      DATA IFF /0/
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11      CONTINUE
        IDUM=1
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
      IF(J.GT.97.OR.J.LT.1)PAUSE
      RAN1=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
      END
