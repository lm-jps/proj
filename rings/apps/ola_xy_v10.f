        subroutine ola(dl_inp, n_inp, f_inp, ux_inp, euxinp, 
     1    uy_inp, euyinp, nmodes, filek1, filsolx, filsoly, 
     1    filker, filcoef, 
     1    qave, qcoef, qverb, ob, oe, num, beg, endd, amu)
        implicit real*8(a-h,o-p,r-z)
        implicit logical(Q)
        character*100 filek1, filobs, filker, filcoef, filsol
       character*100 filkerx,filkery, filcoefx, filcoefy
        character*100 filsolx,filsoly
        parameter(npt=4000,nmd=3100,ndim=300,pi=3.1415926535d0)
        parameter(nx0=100)
        parameter(GG=6.67232d-8)
c       kernels are read in real*4
        real*4 ak(npt), rd(npt), ooo,rs,rms
c       radius, etc. in read*8
        dimension rad(npt),  fakc(npt,nmd)
c       arrays relating to data
        dimension  fdif(50,2000), w(2,nmd), error(2,nmd), n(nmd),rl(nmd)
        dimension om(nmd),dom(2,nmd),err(2,50,2000), freq(2, 50,2000)
c       integration weights
        dimension weight(npt)
c       inversion related arrays
        dimension  aa1x(nmd,nmd), aa1y(nmd,nmd), suml(nmd,nmd,3)
        dimension  ov1(nmd,nmd),vv1x(nmd,nx0),x0(nx0), cov(2,nmd,nmd)
        dimension  vv1y(nmd,nx0)
        dimension  ipiv(nmd)
c       solution related arrays
        dimension delta(nx0), ctil1x(nmd,nx0),ctil2x(nmd,nx0),cencx(nx0)
        dimension  ctil1y(nmd,nx0),ctil2y(nmd, nx0),  cency(nx0)
        dimension dthcx(3),  avcx(npt), avccx(nx0,npt)
        dimension dthcy(3),  avcy(npt), avccy(nx0,npt)
        dimension ckintx(nx0), ckinty(nx0)
c       interpolation related array
        dimension  fb(10)
        dimension dl_inp(nmodes), n_inp(nmodes), f_inp(nmodes)
        dimension ux_inp(nmodes), uy_inp(nmodes)
        dimension euxinp(nmodes), euyinp(nmodes)
	dimension errt(2,50,2000),freqt(2,50,2000),fdift(50,2000)
c
c
C        data ob,oe/1.0,3.01/
C        data num/76/

49      format('#target rad, err-suppr. param, CG of av. kernel, ')
50      format('# 1st,2nd,3rd quart. pt, (3rd-1st)quart.pt., soln, err')
51      format('# radius ', 101e13.5)
53      format('#', 101e13.5)
52      format(101e14.6)
54      format(f7.5,e13.5, 5(1x,f6.4),1p10e13.5)
55      format(1x, 2i6,1x,9e14.6)

C        print*, 'Pl. type in name of file with kernels'
C        read(*,*)filek1
C        print*, 'Pl. type in name of file with observations'
C        read(*,*)filobs
C        print*, 'Pl. type file-name for solutions'
C        read(*,*)filsol
C        print*,'Do you want averaging kernels?'
C        print*,'Type .T. for yes, .F. for no'
C        read(*,*)QAVE
C        if(QAVE)then
C        print*, 'Pl. type file-name for averaging kernels'
C        read(*,*)filker
C        endif
C        print*,'Do you want Inversion coefficients?'
C        print*,'Type .T. for yes, .F. for no'
C        read(*,*)QCOEF
C        if(QCOEF)then
C        print*, 'Pl. type file-name for inversion coefficients'
C        read(*,*)filcoef
C        endif
C        print*, 'Pl. type in frequency limits (in mHz)'
C        read(*,*)ob,oe
C        print*, ob,oe
C        print*, 'Pl. type in number target points'
C        read (*,*)num
C                if(num.gt.100)then
C                print*,' Too many target points, upper limit is 100'
C                stop
C                endif
C        print*, 'Pl. type first and last target points'
C        read (*,*)beg, endd
C        print *,num, beg, endd
c       print *, 'type 1 for Ux inversion  or 2 for Uy inversion'
c       read *, icase
c       print*, icase
C        print*,' Pl. type in error-suppression parameter'
C        read(*,*)amu
        if(qverb)then
          print*, ob,oe
	  print*, amu
        endif

c FILE WITH KERNELS, SHOULD THIS BE HARWIRED?
CC - CSB - changed to formatted
CC        open(unit=21, file=filek1,form='unformatted')
        open(unit=21, file=filek1, status='old')
c FILES WITH SOLUTION
c        filsolx='Ux_'//filsol
c        filsoly='Uy_'//filsol
c        filsolx=filsol//'_ux'
c        filsoly=filsol//'_uy'
        open(10,file=filsolx)
        open(20,file=filsoly)

        
C FILE WITT AVERAGING KERNELS
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

c       %np%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

        do 9090 j=1,num
        vv1x(i,j)=0.
        vv1y(i,j)=0.
9090    continue

9092    continue

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nsp=0
c       READ IN OBSERVED errors
        do 1600 i=1,90000
        ni=n_inp(i)
	lltemp=l
	l=dl_inp(i)
        f=f_inp(i)/1000.
        IF(F.GT.OB.AND.F.LT.OE)THEN
                err(1,ni+1,l+1)=euxinp(i)
                freq(1,ni+1,l+1)=ux_inp(i)
                err(2,ni+1,l+1)=euyinp(i)
                freq(2,ni+1,l+1)=uy_inp(i)
                fdif(ni+1, l+1)=f_inp(i)


		errt(1,ni+1,l+1)=ex1
                freqt(1,ni+1,l+1)=ux1
                errt(2,ni+1,l+1)=ey1
                freqt(2,ni+1,l+1)=uy1
                fdift(ni+1, l+1)=fre

                nsp=nsp+1
        ENDIF
	if(i.gt.nmodes) goto 1700
1600    continue
1700    continue
        if(qverb)then
          print*,'nsp ',nsp
          print*,'i frequencies and velocities', i
          print*, 'MESSAGE: READ IN DATA'
        endif



c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


c       READ IN THE KERNELS

c READ IN RADIUS AND MASS
CC - CSB - added a '*' format char
        read(21, *)rs,rms
c READ IN MESH
        read(21, *)np,(rd(ji),ji=1,np)
        if(qverb)then
          print*, 'R, M ', rs, rms
	endif
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
C        write(23,55)lll,nnn,ooo, dif1, erri1, dif2, erri2
        go to 8100

8110    continue
8100    continue
999     continue

        np=inum
        if(qverb)then
          print*,'np ',np
        endif

        sume1=0.d0
        sume2=0.d0
        do 8101 i=1,nmode
        sume1=sume1+1./(error(1,i)**2)
        sume2=sume2+1./(error(2,i)**2)
8101    continue
        
        if(qverb)then
          print*,'read kernels'
          print*,'read kernels',nmode
        endif
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
        if(qverb)then
          print*, 'sumc', sumcx, sumcy
        endif

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c        do 9010 i=1,nmode
c        w(1,i)=dom(1,i)
c9010    continue



c %%%%%%%%%%%%%%%%%%%%
        if(qverb)then
          print*,'calc. integration weights'
	endif
c       get integration weights

        weight(1)=0.5d0*(rad(1)-rad(2)) 
        do 100 j=1,np-2
        weight(j+1)=0.5d0*(rad(j)-rad(j+2))
100     continue
        weight(np)=0.5*(rad(np-1)-rad(np))


c%%%%%%%%%%%%%%%%%%%%%%
c NOW THE OLA RELATED STUFF
        if(qverb)then
          print*,'calc overlap matrix'
	endif

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

        if(qverb)then
          print*,'numker ',numker
        endif


c       set target radii
c       ------------------
        if(qverb)then
          print*,'calc. targets'
	endif
        ntab=np
        dx=(endd-beg)/(num-1.)
        do 360 i=1,num
        x0(i)=beg+(i-1)*dx
360     continue


        write(10,49)
        write(10,50)
        write(20,49)
        write(20,50)
c       write(10,*)

        

c       NOW FOR THE SOLUTIONS
c       =====================

C NEED LOOP OVER X AND Y HERE
                do 9000 itnum=1,num

c       set the overlap matrix
c       -----------

        if(qverb)then
          print*,'calc. ov1'
        endif
        do 250 i=1,nmode
           do 250 j=i,nmode
                sum1=0
                term=(rad(k)-x0(itnum))**2
        sum1=suml(i,j,1)-2*x0(itnum)*suml(i,j,2)
     1    +(x0(itnum)**2)*suml(i,j,3)
        ov1(i,j)=sum1
        ov1(j,i)=sum1
250     continue

c       set inhomogeneous part of rhs
c       -----------------------------
        do 400 j=1,num
        vv1x(nmode+1,j)=0.5
        vv1y(nmode+1,j)=0.5
400     continue

c       set matrix for lhs
c       ------------------
        if(qverb)then
          print*,'calc. lhs'
        endif

        do 710 j=1,nmode
                do 720 i=1,nmode
                cc=cov(1,i,j)/sumcx
        aa1x(i,j)=ov1(i,j)+amu*cc
                cc=cov(2,i,j)/sumcy
        aa1y(i,j)=ov1(i,j)+amu*cc
720     continue
710     continue

        if(qverb)then
          print*,'numker ', numker, nmode
	endif
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

        if(qverb)then
          print*,'solving eqn'
	endif
        nrhs=1
        call solve(nmd,nmd,numker,nrhs,aa1x,ctil1x,ierr,ipiv,qverb)
	if(qverb)then
          print*,'done Ux'
	endif

        nrhs=1
        call solve(nmd,nmd,numker,nrhs,aa1y,ctil1y,ierr,ipiv,qverb)
	if(qverb)then
          print*,'done Uy'
	endif

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
        if(qverb)then
          print*, 'int. ', sumc1x
	endif
        cencx(jj)=sumc1x/sumc2x
        ckintx(jj)=sumc2x
        cency(jj)=sumc1y/sumc2y
        ckinty(jj)=sumc2y
1115    continue

        if(qverb)then
          print*,'calling width for Ux'
	endif
        call fwidth(npt,np,num,rad,avcx,dthcx,qverb)
	if(qverb)then
          print*,'width done for Ux', dthcx

          print*,'calling width for Uy'
        endif
        call fwidth(npt,np,num,rad,avcy,dthcy,qverb)
	if(qverb)then
          print*,'width done for Uy', dthcy
	endif

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


	if(qverb)then
          print*, x0(itnum), sum1x, er1x
	endif
        write(10,54)x0(itnum), amu, cencx(j), (dthcx(i),i=3,1,-1),
     1  s1x, sum1x, er1x
        write(20,54)x0(itnum), amu, cency(j), (dthcy(i),i=3,1,-1),
     1  s1y, sum1y, er1y

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

	close(10)
	close(20)
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
