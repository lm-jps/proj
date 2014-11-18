c------------------------------------------------------------------------------------------------
c	#define CODE_NAME 		"limbfit"
c	#define CODE_VERSION 	"V5.2r0" 
c	#define CODE_DATE 		"Tue Mar 26 16:30:08 HST 2013" 
c------------------------------------------------------------------------------------------------
c Revision 1.0  2009/01/08  17:20:00  Marcelo Emilio
c changed assumed-size array declaration from "real xxx(1)" to "real xxx(3)"
c jreg=1024, jang=129,jprf=160, jpt=4000000, ahi=30000.0
c variable funcs changed to funcsb
c Revision 2.0 2010/18/12 11:00:00 Marcelo Emilio
c 2010/01/05 JK:Here's a version of limb_5.f that shouldn't return the center guess if clmt is
c               input too large...lets try this at clmt of e-6, e-7, and e-8... -- J
c modified with input shift function to iteratively generate LDF

c Subroutine Limb
c Input: 
c   1) Position and intensity values for limb annulus pixels in the real array
c      "anls(3,*)" x,y coordinates: anls(1,*),anls(2,*); intensity anls(3,*).
c      Expect the number of annulus pixels, integer "ndat", to be <= 40000
c   2) Intial guessimate of the center position, real "cmx,cmy"; try (511.5, 
c      511.5)
c   3) Number of angular bins in which to find average profiles, integer "nang" 
c      to be <=16 (16 gives good resolution for a 7 pixel annulus) 
c   4) Number of radial bins, integer "nprf" <=40 (if this parameter is set to
c      much less than 40, the order of the Chebyschev polynomial "jc" should be
c      changed)
c   5) Number of angular bins in which to compute alpha & beta, the limb scale
c      factor and offset, integer "nreg" <= 512
c Output:
c   1) New improved center position
c   2) Mean radius of all pixels, real "rmean"
c   3) Center-finder convergence diagnostic: number iterations, integer "nitr"
c   4) LDF profile-fitting diagnostic: number of points cut, integer "ncut"
c   5) Average profile for "nang" angular bins:  real array "lprf(nprf,nang+1)"
c      where the last (nang+1) profile is the overall mean limb profile to which
c      the Chebyschev polynomial is fit.  The radial information is stored
c      in real array "rprf(nprf)".
c   6) LDF scale factor and offset per angular bin, real arrays "alph(nreg)" 
c      and "beta(nreg)" 	
c   7) A general success/failure flag: integer "ifail".  For successful 
c      operation ifail=0; any other value corresponds to a failure at a specific
c      numbered site in the code which can be determined by searching for 
c      "ifail=#".  With one exception, a failure (ifail >< 0) is coupled with
c      a "return" to the limb wrapper calling procedure.
c   8) b0 is the "preshift" add this to the output beta to get the net
c	*) centyp to specify if cmx/cmy are guess center(=0) or to be used as center(=1)
	Subroutine limb(anls,npts,cmx,cmy,rguess,nitr,ncut,rprf,lprf, 
     &                   rsi, rso, dx, dy, 
     &                   alph,beta,ifail,b0, centyp, ahi)
	
c      parameter(jpt=40000,jreg=4096,jb=50,jc=50,jang=129,jprf=160)
c	values for HMI rolls
	implicit none
	integer jpt,jb,jc,jreg,jang,jprf,nang,nprf,nreg
	parameter(jpt=8000000,jb=200,jc=200,jreg=1024,jang=181,jprf=64)
	parameter(nang=180,nprf=64,nreg=1024)
	integer i,j,itr,itr2,ilt,ilt2,lmt,ir,lista(3),setn,nremv,centyp
	integer n,nb,ind,nc,spflag,nbad,wbad,cut,jnd
	integer npts,nitr,ncut,ifail,nrem
	real rguess,ain,aout,chebev
	real anls(3,jpt),an(3,jpt),cmx,cmy,alph(jreg),beta(jreg),rsi,rso
	real pi,tpi,alo,ahi,flag,asi,aso,rni,rno,clmt,rem,dreg,dx,dy,fbrem
	real rmax,rmin,rminp,rmaxp,dr,pix(jb+2),rad(jb+2),bin(jb+2)
	real a,x,y,r,savr(jpt),savi(jpt),rmean,imean,rout,rin,theta
	real inte(jpt),thet(jpt),sig(jpt),scmx,scmy,scrit,stpscl
	real cf(3),cov(3,3),chi,cfo(3),cxo,cyo 
	real c(jc),cp(jc),cpp(jc),lsq,ln,dev,err,vr
	real d0,d1,d2,expn(6,jreg),da,db,dc,dd,de,df,dbb,daa
	real crit,svx,svy,xscal,yscal
	real lprf(jprf,jang),prf(jprf,jang),rprf(jprf),dprf,dang
	real b0(jreg)	
	common /ldfcom/ spflag,nb,rad,bin
	common /model/ nc,c,cp,cpp

	external funcs
c	print*,("in fortran")	
c	print*,("cc:"),centyp,ahi,rsi,rso,dx,dy
c	print*,("--:"),ifail,ncut,nitr,cmx,cmy,rguess
c	print*,("--:"),jprf,jang,jreg,jpt,npts
c	print*,(">>"),rprf(1),rprf(jprf),lprf(1,1),lprf(jprf,jang)
c	print*,(">>"),anls(1,1),anls(2,1),anls(3,1),anls(1,npts),anls(2,npts),anls(3,npts)
c	print*,(">>:"),alph(1),alph(jreg),beta(1),beta(jreg),b0(1),b0(jreg)
	pi=3.141592654
	tpi=2.0*pi
	dreg=tpi/float(nreg)	! angular bin size for alpha & beta 
	dang=tpi/float(nang)	! angular bin size for profiles
	alo=1.0	! lower limit for acceptable data
c	ahi=70000.0	! upper limit ...
	flag=-2147483648.0	! flag bad pixels
c	rsi=18.0		! rmean-rsi lower limit to use in center-finder
c	rso=18.0		! rmean+rso upper limit ...
	asi=4.0		! Imean*asi upper limit to use in center-finder
	aso=0.5		! Imean*aso lower limit ...
	fbrem=2.0	! number of std.dev. to cut data in sin/cos fit 
c	dx=-4.0		! x increment to step by in center-finder
c	dy=-4.0		! y increment ...
	stpscl=0.0001	! limit on how small to scale the steps
	ilt=10		! "itr" loop cutoff for stepping dx,dy to center
	ilt2=1		! "itr2" cutoff for finding <r> w/ new center
c	clmt=0.0000001	! coefficient convergence criteria for center-finder
	clmt=1E-7  ! coefficient convergence criteria for center-finder
	rni=rsi	! as "rsi", but for remaining calculations  
	rno=rso		! as "rso", ...  
	nb=nprf		! number of radial bins (jb=50)
	nc=16		! order of the Chebyschev polynomial (jc=50)
	rem=20.0		! number of std.dev. to cut data in the Chebyschev fit
	lmt=20		! limit on # of cuts of rem*sigma; loop counter "cut"
c Initialize stuff; find the min & max radii
	ifail=0
	svx=cmx
	svy=cmy
	do i=1,3
	   lista(i)=i
	enddo
	do i=1,jpt
	   inte(i)=0.0
	   thet(i)=0.0
	   savr(i)=0.0
	   savi(i)=0.0
	enddo
	rmax=0.0
	rmin=100000.0
	do 2 i=1,npts 
	   a=anls(3,i) 
	   if((a.le.alo).or.(a.gt.ahi)) goto 2
	   r=((anls(1,i)-cmx)**2+(anls(2,i)-cmy)**2)**0.5
	   if(r.gt.rmax) rmax=r
	   if(r.lt.rmin) rmin=r
	   ir=ir+1
2       enddo 
	if(centyp.eq.1) goto 40 
c	PRINT*, ("center determination"),cmx,cmy
c Begin center-finding routine.  Starting from the rough center input, find
c the mean radius of the annulus.  Define a sub-annulus; find intensity(angle).
c Iteratively, find a new center that minimizes the difference in the measured
c intensity and the function cf(1)*sin(theta)+cf(2)*cos(theta)+cf(3).
c The convergence criteria is cf(1)**2+cf(2)**2 <<< cf(3)**2.  The code loops
c until this criteria is met or the counter "itr" reaches limit "ilt".
c After convergence, a new annulus and I(theta) is found with the new center.
c For robustness, this new array (with new mean radius) is tested for 
c convergence until counter "itr2" reaches limit "ilt2"

	itr2=0
	scrit=1.0
	nitr=0
10      continue	! loop ilt2 times s.t. test center w/ new radius
	itr2=itr2+1

c Compute radial and intensity means
	ir=0	
	do 20 i=1,npts
	   a=anls(3,i)
	   if((a.le.alo).or.(a.gt.ahi)) goto 20
	   r=((anls(1,i)-cmx)**2+(anls(2,i)-cmy)**2)**0.5
	   ir=ir+1
	   savr(ir)=r
	   savi(ir)=a
20     enddo 
c	print*,ir
	call Hpsort(ir,savr)
	call Hpsort(ir,savi)
	i=ir/2
	rmean=savr(i)
	imean=savi(i)
c	PRINT*, ("rmean rguess"), rmean, rguess

c Define annulus radial range
	rin=rguess-rsi
	rout=rguess+rso
c Define annulus with intensity range
	ain=imean*asi
	aout=imean*aso

	itr=0
30      continue	! loop ilt times s.t. converge on a center
	itr=itr+1

c       PRINT*,"rmean, imean, ain, aout",rmean, imean, ain, aout
	do i=1,jpt
	   sig(i)=1.0
	enddo
c Put sub-annulus data into intensity and angle arrays 
	call darr(anls,npts,ain,aout,cmx,cmy,rin,rout,thet,inte,setn)
c	call gplot(setn, thet, inte, 1, hc, it)	

c Fit sin/cos function (found in "funcs") to intensity(angle)
	call lfit(thet,inte,sig,setn,cf,3,lista,3,cov,3,chi,funcs,ifail)
	if(ifail.gt.0)then
c add by IS Oct5	
		ifail=20
		return
	endif
c Remove outliers from fit 
	call fitb(thet,inte,setn,fbrem,chi,sig,cf,nremv)
c Fit again
	call lfit(thet,inte,sig,setn,cf,3,lista,3,cov,3,chi,funcs,ifail)
	if(ifail.gt.0)then
c add by IS Oct5	
		ifail=21
		return
	endif 
	crit=(cf(1)*cf(1)+cf(2)*cf(2))/cf(3)/cf(3)
c Save the center with the best convergence
c	PRINT*,"scrit,criti, cmx, cmy:",scrit,crit, cmx, cmy	
    	if(crit.lt.scrit)then
	   scrit=crit
	   scmx=cmx
	   scmy=cmy
	endif
c Test convergence; Save number of interations
	if(crit.lt.clmt)then
	  nitr=nitr+itr
	  if(itr2.le.ilt2) then
c save file to debug
c	OPEN(3,FILE='output.txt')
c 	WRITE(3,*)setn
c	WRITE(3,'(I8, F18.8, F18.8)'),(j, thet(J), inte(j), J=1,setn)	
c	WRITE(3,*),(j, thet(J), inte(j),char(13)//char(10), J=1,setn)
c	WRITE(3,*),(j, thet(J), inte(j), J=1,setn)
    	
	goto 10	! loop again
	  endif 
c-1/2010 JRK added if condition to force one center calculation with big clmt
	  if(nitr.gt.2)goto 40	! continue on only if one center correction
	endif
	
c Recompute the coefficients for a center position offset by dx,dy
	cfo(1)=cf(1)
	cfo(2)=cf(2)
	cxo=cmx
	cmx=cmx+dx
	cyo=cmy
	cmy=cmy+dy
	do i=1,jpt
	   sig(i)=1.0
	enddo
	call darr(anls,npts,ain,aout,cmx,cmy,rin,rout,thet,inte,setn)
	call lfit(thet,inte,sig,setn,cf,3,lista,3,cov,3,chi,funcs,ifail)
	if(ifail.gt.0)then
c add by IS Oct5	
		ifail=22
		return
	endif 
	call fitb(thet,inte,setn,fbrem,chi,sig,cf,nremv)
	call lfit(thet,inte,sig,setn,cf,3,lista,3,cov,3,chi,funcs,ifail)
c	PRINT*, "ITR, CHI, NREMV", itr, chi, nremv	
	if(ifail.gt.0)then
c add by IS Oct5	
		ifail=23
		return
	endif
	
c Compute the factor with two sets of coefficients to scale the step size
	xscal=cfo(1)/(cf(1)-cfo(1))
	yscal=cfo(2)/(cf(2)-cfo(2))
	xscal=dx*xscal
	yscal=dy*yscal
	if ((abs(xscal).lt.stpscl).and.(abs(yscal).lt.stpscl)) goto 35
	cmx=cxo-xscal
	cmy=cyo-yscal
c	PRINT*, "cmx, cmy, xscal, yscal, nitr", cmx, cmy, xscal, yscal, nitr
	if(itr.le.ilt) goto 30	! try new center
35	continue		! no convergence
	if(itr2.le.ilt2) then
	  nitr=nitr+itr
	  goto 10		! try new radius
	endif
c No convergence; Reset center to that with the best criteria
c This is flagged in the output file by the parameter "nitr" being its maximum
c value and by "ifail"=1.  Good data converges with under 5 iterations.
	cmx=scmx
	cmy=scmy	
	nitr=nitr+itr
c Do not return here because convergence criteria could be too strict unless
c the center does not change or it changes alot
	cxo=abs(cmx-svx)
	cyo=abs(cmy-svy)
	if(((cxo.eq.0.0).and.(cyo.eq.0.0)).or.
     &    ((cxo.gt.10.0).and.(cyo.gt.10.0)))then
	  ifail=1
 	  return
	endif

c Begin limb figure determination.  Define a sub-annulus about the mean radius
c with the new center coordinates.  Radially bin intensity data.  Fit a Cheby-
c schev polynomial of order "nc".  Iteratively, refit polynomial to data after
c outliers (> "rem"*sigma) have been removed.  The limit to iteration counter 
c "cut" is "lmt".  Note: if cut > lmt, the code continues without action.	

40	continue
c Find the min and max radii
	rmax=0.0
	rmin=100000.0
	do 42 i=1,npts 
	   a=anls(3,i)
	   if((a.le.alo).or.(a.gt.ahi)) goto 42
	   r=((anls(1,i)-cmx)**2+(anls(2,i)-cmy)**2)**0.5
	   if(r.gt.rmax) rmax=r
	   if(r.lt.rmin) rmin=r
42      enddo 

c Find a sub-annulus about the mean radius.  This may not be necessary with the 
c SOI data, but was for data where larger annuli (or full disk was available).
	n=0
	rin=rguess-rni
	rout=rguess+rno
	do 50 i=1,npts
	   a=anls(3,i) 
	   if((a.le.alo).or.(a.gt.ahi)) goto 50
	   x=anls(1,i)
	   y=anls(2,i)
	   r=((x-cmx)**2+(y-cmy)**2)**0.5
	   if((r.lt.rin).or.(r.gt.rout)) goto 50
	   n=n+1
	   savr(n)=r
	   savi(n)=a
	   an(1,n)=x
	   an(2,n)=y
	   an(3,n)=a
50      enddo 
c Loop to cut points deviating from the Chebyschev fit.
	nbad=0		! current number of bad points
	cut=0		! loop counter
	ncut=0		! total points cut
60      continue

c Find min and max radii
	rmax=0.0
	rmin=100000.0
	do 70 i=1,n  
	   a=an(3,i) 
	   if((a.le.alo).or.(a.gt.ahi)) goto 70
	   r=((an(1,i)-cmx)**2+(an(2,i)-cmy)**2)**0.5
	   if(r.gt.rmax) rmax=r
	   if(r.lt.rmin) rmin=r
70	enddo	
c Put data into "nb" radial bins
	dr=(rmax-rmin)/float(nb)
	rad(1)=rmin-dr/2.0
	do i=2,nb+2
	   rad(i)=rad(i-1)+dr
c	(IS: 1st bin -> rmin-dr , 2nd -> rmin)
	enddo
c	print*,rad
	do i=1,nb+2
	   pix(i)=0.0
	   bin(i)=0.0
	enddo
	nrem=0
c-	print*,"-",alo,ahi,tpi,dreg,cmx,cmy,nrem
	do 80 i=1,n
	   a=an(3,i)
	   if((a.le.alo).or.(a.gt.ahi)) goto 80
	   nrem=nrem+1
	   r=((an(1,i)-cmx)**2+(an(2,i)-cmy)**2)**0.5
c  this added to preshift limb distortion to get cleaner mean LDF (IS: maximum distorsion = 2 pixels)
	   theta=atan2((an(2,i)-cmy),(an(1,i)-cmx))
	   if(theta.lt.0.0) theta=theta+tpi
	   ind=int(theta/dreg+1.0)
c-	   print*,"->",i,a,an(1,i),an(2,i),r,rmin,dr,ind,theta,b0(ind)
	   r = r - b0(ind)
c  end of radius correction
	   ind=int((r-rmin)/dr+2.0)
c-		print*," ",r,rmin,dr,ind
c	   if(ind.lt.1) goto 80
c changed w/Jeff Oct 6,2011
	   if(ind.lt.1.or.ind.gt.nb) goto 80
	   bin(ind)=bin(ind)+a
	   pix(ind)=pix(ind)+1.0
c	   print*,"=>",i,ind,bin(ind),pix(ind)
80	enddo 
	do i=2,nb+1
	   if (pix(i).ne.0.0) then
	   	bin(i)=bin(i)/pix(i)
c	   print*,"Z=",i,pix(i),bin(i)
		endif
	enddo
	bin(1)=2.0*bin(2)-bin(3)
c nb+2 bin from above loop only contains intensity values for r=rmax
	bin(nb+2)=2.0*bin(nb+1)-bin(nb)
c	print*,"X>",nb,bin(nb+2)

c Make the Chebyschev fit
	rminp=rmin-dr
	rmaxp=rmax+dr
	spflag=0

c	print*,"3->",bin

c	print*,"before: ",spflag,nb,rad,bin,cp,cpp
c	print*,rminp," ",rmaxp," ",dr,c,nc
	call chebft(rminp,rmaxp,c,nc,ifail)

c	print*,"after: ",spflag,nb,rad,bin,cp,cpp
	if(ifail.gt.0)then
c add by IS Oct5	
		ifail=24
		return
	endif
	call chder(rminp,rmaxp,c,cp,nc)
	call chder(rminp,rmaxp,cp,cpp,nc)

c Do a least-squares calculation
	lsq=0.0
	ln=0.0
	do 90 i=1,n       
	   a=an(3,i)
	   if((a.le.alo).or.(a.gt.ahi)) goto 90
	   ln=ln+1.0
	   r=((an(1,i)-cmx)**2+(an(2,i)-cmy)**2)**0.5
	   lsq=lsq+(chebev(rminp,rmaxp,c,nc,r,ifail)-a)**2
90      enddo
	if(ifail.gt.0) then
c add by IS Oct5	
		ifail=25
		return
	endif
	if(ln.eq.0.0)then
c This should not happen with normal data.
	  ifail=2
	  return
	endif
	dev=(lsq/ln)**0.5

c Cut points from data if they are "rem" sigma away from the fit
	wbad=nbad	! save the number cut to compare to the next pass
	nbad=0
	cut=cut+1
c If "lmt" is exceeded, the parameter passed to output file "nbad" should be 
c obviously large and "ifail"=3. 
	if(cut.gt.lmt)then
	  ifail=3
	  return
	endif
	
	vr=(dev*rem)**2
	do 100 i=1,n
	   a=an(3,i)
	   if((a.le.alo).or.(a.gt.ahi)) goto 100
	   r=((an(1,i)-cmx)**2+(an(2,i)-cmy)**2)**0.5
	   err=(a-chebev(rminp,rmaxp,c,nc,r,ifail))**2
	   if(err.ge.vr)then
	     nbad=nbad+1
	     an(3,i)=flag
	   endif
100     enddo
	if(ifail.gt.0) then
c add by IS Oct5	
		ifail=26
		return 
	endif
	ncut=ncut+nbad
c This criteria works ok for reasonable data.
	if((wbad.eq.0).and.(nbad.eq.0)) goto 110 	! continue on
	goto 60 	! take another cut

110	continue
c Compute the radial profile for "nang" angular bins
	do j=1,jang
	   do i=1,jprf
	      prf(i,j)=0.0
	   enddo
	enddo
	dprf=dr  ! =(rmax-rmin)/float(nb)
	do 118 i=1,n
 	   a=an(3,i)
	   if((a.le.alo).or.(a.gt.ahi)) goto 118
	   x=an(1,i)
	   y=an(2,i)
	   r=((x-cmx)**2+(y-cmy)**2)**0.5
c  this added to preshift limb distortion to get cleaner mean LDF
	   theta=atan2((y-cmy),(x-cmx))
	   if(theta.lt.0.0) theta=theta+tpi
	   ind=int(theta/dreg+1.0)
	   r = r - b0(ind)
c  end of radius correction
	   jnd=int((r-rmin)/dprf+1.0)
	   if(jnd.lt.1)goto 118
	   if(jnd.gt.nprf)then
c Just in case the common radial range doesn't include all profiles
	     if(jnd.gt.nprf+2)then
c Need to consider making a new radial range...or there is problems in the data
		ifail=4
	        return
	     endif
	     jnd=nprf
	   endif
	   theta=atan2((y-cmy),(x-cmx))
	   if(theta.lt.0.0) theta=theta+tpi
	   ind=int(theta/dang+1.0)
	   if(ind.gt.nang)then
c This should not happen with normal data.
c	     ifail=5
		 ind=nang
c	     return
	   endif
	   lprf(jnd,ind)=lprf(jnd,ind)+a
	   prf(jnd,ind)=prf(jnd,ind)+1.0
118	enddo 
	do i=1,nang
	   do j=1,nprf
	      if(prf(j,i).gt.0.1) lprf(j,i)=lprf(j,i)/prf(j,i)
	   enddo
	end do
c The last profile in the array is the average.  "bin" & "rad" were offset 
c one array element to insure the Chebyschev fit do a good job near the edges
c	print*,"jang=",jang
	do i=1,nprf
	   lprf(i,jang)=bin(i+1)
	   rprf(i)=rad(i+1)
c		print*,"avg #:",i,lprf(i,jang),rprf(i)
	enddo
c Determine quantities per angular bin for alpha, beta calculation 
	do j=1,jreg
	   do i=1,6
	      expn(i,j)=0.0
	   enddo
	enddo
	do 120 i=1,n
	   a=an(3,i)
	   if((a.le.alo).or.(a.gt.ahi)) goto 120
	   x=an(1,i)
	   y=an(2,i)
	   r=((x-cmx)**2+(y-cmy)**2)**0.5
	   theta=atan2((y-cmy),(x-cmx))
	   if(theta.lt.0.0) theta=theta+tpi
	   ind=int(theta/dreg+1.0)
c		IS Tue Sep 21 14:54:47 PDT 2010: switched the 2 following paragraphs
	   if(ind.gt.nreg)then
		 ind=nang
c	     ifail=6
c	     return
	   endif
c	correct radius bin
	   r = r - b0(ind)
c	end correction
	   D0=chebev(rminp,rmaxp,c,nc,R,ifail)
	   if(ifail.eq.11) then
c	    print*,"ifail=11!"
	   	ifail=0
	    goto 120
	   endif
	   D1=chebev(rminp,rmaxp,cp,nc,R,ifail)
	   D2=chebev(rminp,rmaxp,cpp,nc,R,ifail)
	   expn(1,ind)=expn(1,ind)+D0*a 
	   expn(2,ind)=expn(2,ind)+D1*a
	   expn(3,ind)=expn(3,ind)+D2*a
	   expn(4,ind)=expn(4,ind)+D0*D1
	   expn(5,ind)=expn(5,ind)+D0*D2+D1**2
	   expn(6,ind)=expn(6,ind)+D0**2
120	enddo
	if(ifail.gt.0)then
c add by IS Oct5	
		ifail=27
		return
	endif
c Determine alpha (the LDF scale factor) and beta (the offset)
	do i=1,nreg
	   thet(i)=float(i)*360.0/float(nreg)
	   DA=expn(1,i)
	   DB=expn(2,i)
	   DC=expn(3,i)
	   DD=expn(4,i)
	   DE=expn(5,i)
	   DF=expn(6,i)
	   DBB=DA*DE-DC*DF-DB*DD
	   if(DBB.eq.0.0)then
	     beta(i)=-10.0
	     alph(i)=-10.0
	     ifail=7
	     goto 130 
c If too many points are cut, DBB (and DAA) can be zero.  The flag this with
c ifail=7 but continue on.  Alpha,betas with value -10 in a angular particular
c bin will show up in the output profiles.  Any further software will have check
c for this unnatural phenomenon.
	   endif
	   beta(i)=(DA*DD-DB*DF)/DBB
	   DBB=beta(i)
	   DAA=DD-DBB*DE
	   if(DAA.eq.0.0)then
	     alph(i)=-10.0
	     ifail=7
	     goto 130
	   endif
	   alph(i)=(DB-DBB*DC)/DAA-1.0
130	enddo
c	print*,"end... ",cmx,cmy	
	return
	end	
	

c Convert data into angle vs. radius, intensity, or I*r using trial center, 
c cuts in radius, and intensity cuts
	Subroutine Darr(anls,npts,ain,aout,cmx,cmy,rin,rout,thet,inte,n)
	implicit none
	integer jpt
	parameter(jpt=8000000)

	integer i,n,npts
	real anls(3,jpt),aout,ain,rout,rin,thet(1),inte(1)
	real a,x,y,theta,r,tpi,cmx,cmy

	tpi=6.283185307
	n=0

	do 10 i=1,npts
	   a=anls(3,i)
	   if((a.le.aout).or.(a.ge.ain)) goto 10 
	   x=anls(1,i)
	   y=anls(2,i)
	   r=((x-cmx)**2+(y-cmy)**2)**0.5
	   if((r.lt.rin).or.(r.gt.rout)) goto 10
	   theta=atan2(cmx-x,cmy-y)
	   if(theta.lt.0.0) theta=theta+tpi
	   n=n+1
	   thet(n)=theta
c	   inte(n)=a*r
c	   inte(n)=a
	   inte(n)=r
10      enddo 

	return
	end


c Flags points to cut for a better fit
	Subroutine Fitb(thet,inte,n,rem,chi,sig,cf,npts)
	implicit none
	integer i,n,npts
	real thet(1),inte(1),rem,chi,sig(1),cf(3),dev,err,tp

	npts=0
	dev=rem*(chi/float(n))**0.5
	do i=1,n
	   tp=cf(1)*sin(thet(i))+cf(2)*cos(thet(i))+cf(3)
	   err=abs(tp-inte(i))
	   if(err.ge.dev)then
	     npts=npts+1
	     sig(i)=99999.0	! sig(i) is normally 1.0
	   endif
	enddo

	return
	end	


c Computes the limb darkening function at any radius with a cubic spline	
	Real Function Ldf(r,ifail)
	implicit none
	integer jb
	parameter(jb=200)

	integer nb,nb2,spflag,ifail
	real r,bin(jb+2),rad(jb+2),snd(jb)

	common /ldfcom/ spflag,nb,rad,bin

	nb2=nb+2 
	ldf=0.0
	if(spflag.eq.0)then
	  call spline(rad,bin,nb2,2.0e30,2.0e30,snd)
	  spflag=1
	endif
	call splint(rad,bin,snd,nb2,r,ldf,ifail)
	
	return
	end


c Performs a Heap Sort (see Numerical Recipes) to find the mean of array "ra"
	Subroutine Hpsort(n,ra)
	implicit none
	integer i,ir,j,l,n
	real ra(1),rra

	if(n.lt.2) return

	l=n/2+1
	ir=n
10	continue
 
	if(l.gt.1)then
	  l=l-1
	  rra=ra(l)
	else
	  rra=ra(ir)
	  ra(ir)=ra(l)
	  ir=ir-1
	  if(ir.eq.1)then
	    ra(1)=rra
	    return
  	  endif 
	endif
	i=l
	j=l+1
20	if(j.le.ir)then
	  if(j.lt.ir)then
	    if(ra(j).lt.ra(j+1)) j=j+1
	  endif
	  if(rra.lt.ra(j))then
	    ra(i)=ra(j)
	    i=j
	    j=j+j
	  else
	    j=ir+1
	  endif
	  goto 20
	endif
	ra(i)=rra
	goto 10

	end
   

	SUBROUTINE COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
	  implicit none
	  integer ncvmp,ma,mfit
      parameter(ncvmp=3)
      real COVAR(ncvmp,ncvmp)
      real swap
      integer LISTA(3),j,i,ncvm
      
      DO 12 J=1,MA-1
        DO 11 I=J+1,MA
          COVAR(I,J)=0.
11      CONTINUE
12    CONTINUE

      DO 14 I=1,MFIT-1
        DO 13 J=I+1,MFIT
          IF(LISTA(J).GT.LISTA(I)) THEN
            COVAR(LISTA(J),LISTA(I))=COVAR(I,J)
          ELSE
            COVAR(LISTA(I),LISTA(J))=COVAR(I,J)
          ENDIF
13      CONTINUE
14    CONTINUE

      SWAP=COVAR(1,1)
      DO 15 J=1,MA
        COVAR(1,J)=COVAR(J,J)
        COVAR(J,J)=0.
15    CONTINUE
      COVAR(LISTA(1),LISTA(1))=SWAP
      DO 16 J=2,MFIT
        COVAR(LISTA(J),LISTA(J))=COVAR(1,J)
16    CONTINUE
      DO 18 J=2,MA
        DO 17 I=1,J-1
          COVAR(I,J)=COVAR(J,I)
17      CONTINUE
18    CONTINUE

      RETURN
      END


	SUBROUTINE GAUSSJ(A,N,NP,B,M,MP,ifail)
	implicit none
	integer nmax,npp,mpp,n,mp,ifail,m,np
      PARAMETER(NMAX=200,npp=3,mpp=1)
      real B(npp,mpp),A(npp,npp),big,dum,pivinv
      integer IPIV(NMAX),INDXR(NMAX),INDXC(NMAX),j,k,i,irow,icol,l,ll

      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE

      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG) THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
	           ifail=9
		   return
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1

        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF

        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) THEN
	   ifail=9
	   return
	endif
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END


C FUNCTION FOR LFIT
	SUBROUTINE FUNCS(X,AFUNC)
	implicit none
	REAL AFUNC(3) ,x

	AFUNC(1)=SIN(X)
	AFUNC(2)=COS(X)
	AFUNC(3)=1.0

	RETURN
	END


	SUBROUTINE LFIT(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,COVAR,NCVM,CHISQ,
     & FUNCS,ifail)
	implicit none
	integer mmax,ncvmp,jpt
	PARAMETER (MMAX=3,ncvmp=3,jpt=8000000)
	integer lista(3),kk,j,k,ihit,ifail,i,ncvm,ndata,mfit,ma
	real X(jpt),Y(jpt),SIG(jpt),A(3),COVAR(ncvmp,ncvmp)
	real BETA(MMAX),AFUNC(MMAX),CHISQ,SUM,WT,SIG2I,YM
	external funcs

      KK=MFIT+1
      DO 12 J=1,MA
        IHIT=0
        DO 11 K=1,MFIT
          IF (LISTA(K).EQ.J) IHIT=IHIT+1
11      CONTINUE
        IF (IHIT.EQ.0) THEN
          LISTA(KK)=J
          KK=KK+1
        ELSE IF (IHIT.GT.1) THEN
	     ifail=8
	     return
        ENDIF
12    CONTINUE
      IF (KK.NE.(MA+1)) THEN
         ifail=8
         return
      ENDIF

      DO 14 J=1,MFIT
        DO 13 K=1,MFIT
          COVAR(J,K)=0.
13      CONTINUE
        BETA(J)=0.
14    CONTINUE

      DO 18 I=1,NDATA
        CALL FUNCS(X(I),AFUNC)
        YM=Y(I)
        IF(MFIT.LT.MA) THEN
          DO 15 J=MFIT+1,MA
            YM=YM-A(LISTA(J))*AFUNC(LISTA(J))
15        CONTINUE
        ENDIF
        SIG2I=1./SIG(I)**2
        DO 17 J=1,MFIT
          WT=AFUNC(LISTA(J))*SIG2I
          DO 16 K=1,J
            COVAR(J,K)=COVAR(J,K)+WT*AFUNC(LISTA(K))
16        CONTINUE
          BETA(J)=BETA(J)+YM*WT
17      CONTINUE
18    CONTINUE

      IF (MFIT.GT.1) THEN
        DO 21 J=2,MFIT
          DO 19 K=1,J-1
            COVAR(K,J)=COVAR(J,K)
19        CONTINUE
21      CONTINUE
      ENDIF

      CALL GAUSSJ(COVAR,MFIT,NCVM,BETA,1,1,ifail)
      if(ifail.gt.0)  then
c add by IS Oct5	
		ifail=30
		return
	  endif
      DO 22 J=1,MFIT
        A(LISTA(J))=BETA(J)
22    CONTINUE
      CHISQ=0.
      DO 24 I=1,NDATA
        CALL FUNCS(X(I),AFUNC)
        SUM=0.
        DO 23 J=1,MA
          SUM=SUM+A(J)*AFUNC(J)
23      CONTINUE
        CHISQ=CHISQ+((Y(I)-SUM)/SIG(I))**2
24    CONTINUE

      CALL COVSRT(COVAR,NCVM,MA,LISTA,MFIT)

	RETURN
	END


	SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
	implicit none
	integer nmax,jb2
	PARAMETER (NMAX=200,jb2=202)
	real X(jb2),Y(jb2),Y2(jb2),U(NMAX),yp1,YPN,p,qn,un,k,sig
	integer n,i

      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF

      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE

      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE

	RETURN
	END


	SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y,ifail)
	implicit none
	integer jb2
	PARAMETER (jb2=202)
	
	real XA(jb2),YA(jb2),Y2A(jb2),x,y,a,b
	integer k,klo,khi,h,ifail,n

      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF(H.EQ.0.)then
		ifail=10
c		print*,"ifail=10"
		return
      endif
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.

	RETURN
	END


c Chebyschev Approximation to function expressed in func on interval
c [a,b].  Numerical Recipes pp 184-190.


c Compute coefficients c(1:n)
	Subroutine CHEBFT(a,b,c,n,ifail)
	implicit none
	integer nmax
	real pi
	parameter(nmax=200,pi=3.141592653589793D0)	
	integer n,j,k,ifail
	real a,b,c(nmax),bma,bpa,fac,y,f(nmax),ldf
	double precision sum
c	external func
	
	bma=0.5*(b-a)
	bpa=0.5*(b+a)
	do k=1,n
	   y=cos(pi*(k-0.5)/n)
	   f(k)=ldf(y*bma+bpa,ifail)
	enddo
	if(ifail.gt.0) then
c add by IS Oct5	
		ifail=31
c		print*,"ifail=31"
		return
	endif
	fac=2.0/n
	do j=1,n
	   sum=0.D0
	   do k=1,n
	      sum=sum+f(k)*cos((pi*(j-1))*((k-0.5D0)/n))
	   enddo	
	   c(j)=fac*sum
	enddo

	return
	end


c Evaluate the approximation to f(x)
	Function CHEBEV(a,b,c,m,x,ifail)
	implicit none

	integer m,j,ifail,jb
	parameter(jb=200)
	real chebev,a,b,x,c(jb),d,dd,sv,y,y2

	if((x-a)*(x-b).gt.0.0) then
c CHEBEV out of range
	  ifail=11
	  return
	endif
	d=0.0
	dd=0.0
	y=(2.0*x-a-b)/(b-a)
	y2=2.0*y
	do j=m,2,-1
	   sv=d
	   d=y2*d-dd+c(j)
	   dd=sv
	enddo
	chebev=y*d-dd+0.5*c(1)

	return
	end


c Compute the coefficients of the derivative of the function whose
c coefficients are c(n)
	Subroutine CHDER(a,b,c,cder,n)
	implicit none
	integer jc
	parameter(jc=200)
	integer n,j
	real a,b,c(jc),cder(jc),con

	cder(n)=0.0
	cder(n-1)=2.0*(n-1)*c(n)
	if(n.ge.3)then
	  do j=n-2,1,-1
	     cder(j)=cder(j+2)+2.0*j*c(j+1)
	  enddo
	endif
	con=2.0/(b-a)
	do j=1,n
	   cder(j)=cder(j)*con
	enddo

	return
	end

