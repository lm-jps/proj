        subroutine ola(dl_inp, n_inp, f_inp, ux_inp, euxinp,
     1    uy_inp, euyinp, nmodes, filek1, x0, rcgx, rcgy,
     1    quartx, quarty, solnx, solny, errx, erry, filker, filcoef,
     1    qave, qcoef, qverb, ob, oe, num, beg, endd, amu, retval)
	IMPLICIT NONE
        character*100 filek1, filobs, filker, filcoef, filsol
	character*100 filkerx, filkery, filcoefx, filcoefy
	integer*4 retval
	REAL*8 GG, amu, anorml, beg, cc, dif1, dif2, dthcx, dthcy, dx
	REAL*8 endd, er1x, er1y, erri1, erri2
	REAL*8 f, ob, oe, omax, omin, ot, pi
	REAL*8 s1x, s1y, sa, sb, sc, sum1, sum1x, sum1y, sum2, sum3
	REAL*8 sumc1x, sumc1y, sumc2x
	REAL*8 sumc2y, sumccx, sumccy, sumcx, sumcy, sume1, sume2
	REAL*8 term
	INTEGER*4 i, ierr, ii, inum, itnum, ipiv, j, ji, jj
	INTEGER*4 k, kk, l, lll
	INTEGER*4 ndim, ni, nmd, nmode, nmodes, nnn, np, npt, nrhs, num
	INTEGER*4 ntab, numker, nx0
	LOGICAL qverb, QAVE, QCOEF
        parameter(npt=4000,nmd=3100,ndim=300,pi=3.1415926535d0)
        parameter(nx0=100)
        parameter(GG=6.67232d-8)
c       kernels are read in real*4
        real*4 ak(npt), rd(npt), ooo,rs,rms
c       radius, etc. in read*8
        REAL*8 rad(npt), fakc(npt,nmd)
c       arrays relating to data
	REAL*8 fdif(50,2000), w(2,nmd), error(2,nmd), rl(nmd)
	INTEGER*4 n(nmd)
	REAL*8 err(2,50,2000), freq(2,50,2000)
        REAL*8 om(nmd), dom(2,nmd)
c       integration weights
        REAL*8 weight(npt)
c       inversion related arrays
        REAL*8 aa1x(nmd,nmd), aa1y(nmd,nmd), suml(nmd,nmd,3)
	REAL*8 ov1(nmd,nmd), vv1x(nmd,nx0), cov(2,nmd,nmd)
	REAL*8 x0(num),rcgx(num),rcgy(num),quartx(3,num),quarty(3,num)
	REAL*8 solnx(num),solny(num),errx(num),erry(num)
        REAL*8  vv1y(nmd,nx0)
        dimension  ipiv(nmd)
c       solution related arrays
        REAL*8 delta(nx0), ctil1x(nmd,nx0), ctil2x(nmd,nx0), cencx(nx0)
        REAL*8 ctil1y(nmd,nx0), ctil2y(nmd,nx0), cency(nx0)
        REAL*8 avcx(npt), avccx(nx0,npt)
        dimension dthcx(3), dthcy(3)
        REAL*8 avcy(npt), avccy(nx0,npt)
        REAL*8 ckintx(nx0), ckinty(nx0)
c       interpolation related array
        REAL*8  fb(10)
        REAL*8 dl_inp(nmodes), f_inp(nmodes)
        INTEGER*4 n_inp(nmodes)
        REAL*8 ux_inp(nmodes), uy_inp(nmodes)
        REAL*8 euxinp(nmodes), euyinp(nmodes)

51      format('# radius ', 101e13.5)
53      format('#', 101e13.5)
52      format(101e14.6)
C 54      format(f7.5,e13.5, 5(1x,f6.4),1p10e13.5)
55      format(1x, 2i6,1x,9e14.6)

C							checks added 2012.07.16
C        PRINT*,'begin ola'
	retval = 0
	if (num.gt.nx0 .or. nmodes.gt.90000) then
	  retval = 1
	  return
	endif

C        if(qverb)then
C         print*, ob,oe
C	  print*, amu
C        endif

        open(unit=21, file=filek1, status='old')
C FILE WITH AVERAGING KERNELS
        if(QAVE)then
        filkerx='Ux_'//filker
        filkery='Uy_'//filker
        open(19,file=filkerx)
        open(29,file=filkery)
        endif
c FILE WITH INVERSION COEFFICIENTS
       if(QCOEF)then
        filcoefx='Ux_'//filcoef
        filcoefy='Uy_'//filcoef
        open(27,file=filcoefx)
        open(37,file=filcoefy)
        endif

c       Constants to normalise various quantities to dimensionless units
        anorml=4.*pi/3.

c%%%%%%%%%%%%%%%%%INITIALIZE MATRICES %%%%%%

        do 1400 i=1,50
        do 1400 j=1,2000
        fdif(i,j)=0.0
        freq(1,i,j)=0.0
        err(1,i,j)=0.0
        freq(2,i,j)=0.0
        err(2,i,j)=0.0
1400    continue

        do 9092 i=1,nmd

                do 9091 j=1,nmd
                ov1(j,i)=0.
                cov(1,j,i)=0.
                cov(2,j,i)=0.
                  do 9083 kk=1,3
                  suml(j,i,kk)=0.d0
9083    continue
9091    continue

	do 9090 j=1,nx0
	  vv1x(i,j)=0.
	  vv1y(i,j)=0.
C			       additional initializations added added 2012.07.16
	  ctil1x(i,j)=0.0
	  ctil2x(i,j)=0.0
	  ctil1y(i,j)=0.0
	  ctil2y(i,j)=0.0
9090    continue

9092    continue

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C       nsp=0
c       READ IN OBSERVED errors
        do 1600 i=1,nmodes
        ni=n_inp(i)
	l=dl_inp(i)
        f=f_inp(i)/1000.
        IF(F.GT.OB.AND.F.LT.OE.and.l.lt.2000.and.ni.lt.50)THEN
                err(1,ni+1,l+1)=euxinp(i)
                freq(1,ni+1,l+1)=ux_inp(i)
                err(2,ni+1,l+1)=euyinp(i)
                freq(2,ni+1,l+1)=uy_inp(i)
                fdif(ni+1, l+1)=f_inp(i)
	
C                nsp=nsp+1
        ENDIF
1600    continue
C        if(qverb)then
C          print*,'nsp ',nsp
C          print*,'i frequencies and velocities', i
C          print*, 'MESSAGE: READ IN DATA'
C        endif

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

c       READ IN THE KERNELS

c READ IN RADIUS AND MASS
CC - CSB - added a '*' format char
        read(21, *)rs,rms
c READ IN MESH
        read(21, *)np,(rd(ji),ji=1,np)
C        if(qverb)then
C          print*, 'R, M ', rs, rms
C	endif
        nmode=0
        omax=0.
        omin=20000
        do 8100 i=1,20000
        read(21, *, end=999)lll,nnn, ooo,(ak(ii),ii=1,np)   
        if(lll.ge.2000) go to 8110
        if(nnn.ge.50) go to 8110

        dif1=freq(1,nnn+1,lll+1)
        dif2=freq(2,nnn+1,lll+1)
        erri1=err(1,nnn+1,lll+1)
        erri2=err(2,nnn+1,lll+1)
        

        if(erri1.le.0.) go to 8110
        if(abs(erri1).le.1e-11) go to 8110
        if(erri2.le.0.) go to 8110
        if(abs(erri2).le.1e-11) go to 8110


81      nmode=nmode+1
        ot=fdif(nnn+1,lll+1)

                inum=0
c       TURNING REAL*4 kernels to real*8 ONES
                do 8111 jj=1,np-1
                inum=inum+1
                fakc(inum,nmode)=ak(jj)
                rad(inum)=rd(jj)
8111            continue

        n(nmode)=nnn
        rl(nmode)=lll
        om(nmode)=ot/1000.d0
        error(1,nmode)=erri1
        error(2,nmode)=erri2
        dom(1,nmode)=dif1
        dom(2,nmode)=dif2
        go to 8100

8110    continue
8100    continue
999     continue

        np=inum
C        if(qverb)then
C          print*,'np ',np
C        endif

        sume1=0.d0
        sume2=0.d0
        do 8101 i=1,nmode
        sume1=sume1+1./(error(1,i)**2)
        sume2=sume2+1./(error(2,i)**2)
8101    continue
        
C        if(qverb)then
C          print*,'read kernels'
C          print*,'read kernels',nmode
C        endif
c %%%%%%%%%%%%%%%%%%%%%%



c       set covariance matrix
        sumcx=0.
        sumcy=0.
        do 840 j=1,nmode
        sumcx=sumcx+error(1,j)*error(1,j)
        cov(1,j,j)=error(1,j)*error(1,j)
        sumcy=sumcy+error(2,j)*error(2,j)
        cov(2,j,j)=error(2,j)*error(2,j)
840     continue

        sumcx=sumcx/nmode
        sumcy=sumcy/nmode
C        if(qverb)then
C          print*, 'sumc', sumcx, sumcy
C        endif

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c        do 9010 i=1,nmode
c        w(1,i)=dom(1,i)
c9010    continue

c %%%%%%%%%%%%%%%%%%%%
C        if(qverb)then
C          print*,'calc. integration weights'
C	endif
c       get integration weights

        weight(1)=0.5d0*(rad(1)-rad(2)) 
        do 100 j=1,np-2
        weight(j+1)=0.5d0*(rad(j)-rad(j+2))
100     continue
        weight(np)=0.5*(rad(np-1)-rad(np))


c%%%%%%%%%%%%%%%%%%%%%%
c NOW THE OLA RELATED STUFF
C        if(qverb)then
C          print*,'calc overlap matrix'
C	endif

c       set the overlap matrix
c       ----------------------------------------

        do 200 i=1,nmode
c       print*,'overlap ', i
           do 200 j=i,nmode

                sum2=0.d0
                sum3=0.d0
                sa=0.
                sb=0.
                sc=0.
                
                do 202 k=1,np
        term=weight(k)*fakc(k,i)*fakc(k,j)
        sa=sa+term*(rad(k)**2)
        sb=sb+term*rad(k)
        sc=sc+term
202     continue

        suml(i,j,1)=sa
        suml(j,i,1)=sa
        suml(i,j,2)=sb
        suml(j,i,2)=sb
        suml(i,j,3)=sc
        suml(j,i,3)=sc
200     continue

c       The lagrange multiplier
        do 203 i=1,nmode
                sum1=0.
                do 204 k=1,np
                sum1=sum1+weight(k)*fakc(k,i)
204     continue
        ov1(i,nmode+1)=0.5*sum1
        ov1(nmode+1,i)=ov1(i,nmode+1)
        suml(i,nmode+1,1)=0.5*sum1
        suml(nmode+1,i,1)=ov1(i,nmode+1)
203     continue

        numker=nmode+1

C        if(qverb)then
C          print*,'numker ',numker
C        endif


c       set target radii
c       ------------------
C        if(qverb)then
C          print*,'calc. targets'
C	endif
        ntab=np
        dx=(endd-beg)/(num-1.)
        do 360 i=1,num
        x0(i)=beg+(i-1)*dx
360     continue

c       NOW FOR THE SOLUTIONS
c       =====================

C NEED LOOP OVER X AND Y HERE
	do 9000 itnum=1,num

c       set the overlap matrix
c       -----------

C          if(qverb)then
C            print*,'calc. ov1'
C          endif
          do 250 i=1,nmode
            do 250 j=i,nmode
              sum1=0
              term=(rad(k)-x0(itnum))**2
              sum1=suml(i,j,1)-2*x0(itnum)*suml(i,j,2)
     1          +(x0(itnum)**2)*suml(i,j,3)
              ov1(i,j)=sum1
              ov1(j,i)=sum1
250       continue

c       set inhomogeneous part of rhs
c       -----------------------------
        do 400 j=1,num
        vv1x(nmode+1,j)=0.5
        vv1y(nmode+1,j)=0.5
400     continue

c       set matrix for lhs
c       ------------------
C        if(qverb)then
C          print*,'calc. lhs'
C        endif

        do 710 j=1,nmode
                do 720 i=1,nmode
                cc=cov(1,i,j)/sumcx
        aa1x(i,j)=ov1(i,j)+amu*cc
                cc=cov(2,i,j)/sumcy
        aa1y(i,j)=ov1(i,j)+amu*cc
720     continue
710     continue

C        if(qverb)then
C          print*,'numker ', numker, nmode
C	endif
        do 711 j=1,numker
        do 712 i=nmode+1, numker
        aa1x(i,j)=ov1(i,j)
        aa1x(j,i)=ov1(j,i)
        aa1y(i,j)=ov1(i,j)
        aa1y(j,i)=ov1(j,i)
712     continue
                do 730 k=1,num
                ctil1x(j,k)=vv1x(j,k)
                ctil1y(j,k)=vv1y(j,k)
730     continue
711     continue

c       solving the equations
c       --------------------

C        if(qverb)then
C          print*,'solving eqn'
C	endif
        nrhs=1
        call solve(nmd,nmd,numker,nrhs,aa1x,ctil1x,ierr,ipiv,qverb)
C	if(qverb)then
C          print*,'done Ux'
C	endif

        nrhs=1
        call solve(nmd,nmd,numker,nrhs,aa1y,ctil1y,ierr,ipiv,qverb)
C	if(qverb)then
C          print*,'done Uy'
C	endif

        jj=1
        do 7181 i=1,numker
        ctil2x(i,itnum)=ctil1x(i,jj)
        ctil2y(i,itnum)=ctil1y(i,jj)
7181    continue

c       Calculating averaging kernel
c       -------------------------------

        do 1112 j=1,np

                jj=1
                sumccx=0.
                sumccy=0.
                do 1114 i=1,nmode
                sumccx=sumccx+ctil1x(i,jj)*fakc(j,i)
                sumccy=sumccy+ctil1y(i,jj)*fakc(j,i)
1114    continue
        avcx(j)=sumccx
        avccx(itnum,j)=sumccx
        avcy(j)=sumccy
        avccy(itnum,j)=sumccy
1113    continue


1112    continue

c       Finding cg of averaging kernel
c       ------------------------------

                sumc1x=0.
                sumc2x=0.
                sumc1y=0.
                sumc2y=0.
                do 1116 j=1,np
                sumc1x=sumc1x+weight(j)*rad(j)*avcx(j)
                sumc2x=sumc2x+weight(j)*avcx(j)
                sumc1y=sumc1y+weight(j)*rad(j)*avcy(j)
                sumc2y=sumc2y+weight(j)*avcy(j)
1116            continue
C        if(qverb)then
C          print*, 'int. ', sumc1x
C	endif
        cencx(jj)=sumc1x/sumc2x
        ckintx(jj)=sumc2x
        cency(jj)=sumc1y/sumc2y
        ckinty(jj)=sumc2y
1115    continue

C        if(qverb)then
C          print*,'calling width for Ux'
C	endif
        call fwidth(npt,np,num,rad,avcx,dthcx,qverb)
C	if(qverb)then
C          print*,'width done for Ux', dthcx
C          print*,'calling width for Uy'
C        endif
        call fwidth(npt,np,num,rad,avcy,dthcy,qverb)
C	if(qverb)then
C          print*,'width done for Uy', dthcy
C	endif

c       the solution
c       ------------


        do 415 j=1,1
        sum1x=0.
        sum1y=0.
        er1x=0.
        er1y=0.

        do 420 i=1,nmode
        sum1x=sum1x+ctil1x(i,j)*dom(1,i)
        er1x=er1x+ctil1x(i,j)**2*cov(1,i,i)
        sum1y=sum1y+ctil1y(i,j)*dom(2,i)
        er1y=er1y+ctil1y(i,j)**2*cov(2,i,i)
420     continue

                s1x=abs(dthcx(1)-dthcx(3))
                s1y=abs(dthcy(1)-dthcy(3))

        er1x=sqrt(er1x)
        er1y=sqrt(er1y)


C	if(qverb)then
C          print*, x0(itnum), sum1x, er1x
C	endif
	rcgx(itnum) = cencx(j)
	rcgy(itnum) = cency(j)
	do 421 i = 1,3
	  quartx(i,itnum) = dthcx(i)
	  quarty(i,itnum) = dthcy(i)
421	continue
	solnx(itnum) = sum1x
	solny(itnum) = sum1y
	errx(itnum) = er1x
	erry(itnum) = er1y

415     continue


700     continue
9000    continue

        if(QAVE)then
c       WRITE OUT THE AVERAGING KERNELS
        
        write(19,51)  (x0(i),i=1,num)
        write(29,51)  (x0(i),i=1,num)
        do 8010 j=1,np,2
        write(19,52)rad(j),(avccx(k,j),k=1,num)
        write(29,52)rad(j),(avccy(k,j),k=1,num)
8010    continue
        endif
        if(QCOEF)then
c       WRITE OUT INVERSION COEFFICIENTS
        write(27,53)  (x0(i),i=1,num)
        write(37,53)  (x0(i),i=1,num)
        do 8012 j=1,numker-1
        write(27,52) (ctil2x(j,i),i=1,num)
        write(37,52) (ctil2y(j,i),i=1,num)
8012    continue
        endif

C	close(10)
C	close(20)
	close(21)
	if(QAVE)then
	  close(19)
	  close(29)
	endif
	if(QCOEF)then
	  close(27)
	  close(37)
	endif

        end
