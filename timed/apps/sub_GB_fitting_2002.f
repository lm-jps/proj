C SUBROUTINE TO COMPUTE THE TRAVEL TIMES
c BASED ON THE DEFINITION
c BY GIZON & BIRCH (2002) EQUATION A1
c VERSION 1.3, 07/26/2010
c AUTHOR SEBASTIEN COUVIDAT

c After computing the travel-time maps (positive and negative times)
c the code computes the mean value and rms variation of these times.
c Times that are larger than a threshold (in absolute values)
c are deemed dubious and are re-computed based on a neighbor average of the
c cross-covariance. The threshold is equal to the rms variation time
c a certain amount. This amount can be changed in the code

c ISSUES TO SOLVE
c the format of the cross-covariance file is assumed to be negative 
c times followed by positive times
c MAKE SURE THAT kmmin-thresh >= 1, kpmin-thresh>=1, kmmax+thresh<=nt, 
c kmmin+thresh<=nt
c COMPILE ON n00 WITH exe77

        SUBROUTINE gbtimes02(correl2,correl3,nx,ny,nt,iii,tau)

        IMPLICIT REAL(a-h,o-z)
        PARAMETER (pi=3.141592653589793116)
        PARAMETER (naxis=3,nmax=4096)
        PARAMETER (dt=45.,nt2=1791,dt2=5.,nt3=63,nthresh=31,nt4=5)
        REAL correl2(nt,ny,nx),correl3(nt,ny,nx)
        REAL correl(nx,ny,nt),correltemp(nt),correlint(nt2),amaxi
        REAL taup(nx,ny),tp(nt3),tm(nt3),sinc(nt2,nt)
        REAL taum(nx,ny),correlref(nt),correlrefint(nt2),tau(nx,ny,2)
        REAL c0,c1,c2,alags(nt3),tpc(nt4),tmc(nt4),alagscp(nt4)
        REAL aminimump,aminimumm,ameanp,ameanm,rmsp,rmsm,threshold
	REAL time(nt2),time0(nt),temp,alagscm(7)
        CHARACTER*100 filename,filewrite

        DO 5 i=1,nx
          DO 5 j=1,ny
           DO 5 k=1,nt
             correl(i,j,k)=correl3(k,j,i)
 5      CONTINUE

        threshold=5.0
        nedge=8

c       compute the time window to select the first-bounce skip
        CALL window(kpmin,kpmax,kmmin,kmmax,nt,nt2,iii,dt,dt2)

c       compute the reference cross-covariance
c       and avoid the edges (NB: can change the parameters)
        DO k=1,nt
           correlref(k) = 0.
           DO i=nedge,nx-nedge
              DO j=nedge,ny-nedge
                 correlref(k)=correlref(k)+correl(i,j,k)
              ENDDO
           ENDDO
           correlref(k)=correlref(k)/(nx-2*nedge+1)/(ny-2*nedge+1) 
        ENDDO

c       normalize the reference cross-covariance
        amaxi=MAXVAL(correlref)

        DO k=1,nt
           correlref(k)=correlref(k)/amaxi
        ENDDO

        DO 10 i=1,nx
          DO 10 j=1,ny
           DO 10 k=1,nt
             correl(i,j,k)=correl2(k,j,i)
 10      CONTINUE

c       initialize some interpolation parameters
        DO i=1,nt2
           time(i)=(i-1.)*dt2
        ENDDO
        DO i=1,nt
           time0(i)=(i-1.)*dt
        ENDDO

        DO i=1,nt2 
           DO j=1,nt
              temp=pi*(time(i)-time0(j))/dt
              IF (temp .EQ. 0.) THEN 
                 sinc(i,j)=1. 
              ELSE 
                 sinc(i,j)=SIN(temp)/temp
              ENDIF
           ENDDO
        ENDDO

c       interpolate in time the reference cross-covariance
        CALL interpolation(correlref,nt,correlrefint,sinc,nt2,kpmin,
     +       kpmax,kmmin,kmmax,nthresh)

        DO i=1,nt3
           alags(i)=(i-1-nthresh)*dt2
        ENDDO

c       compute the travel times
        DO i=1,nx
           DO j=1,nx

c             normalize the cross-covariance
              DO k=1,nt
                 correltemp(k)=correl(i,j,k)/amaxi
              ENDDO

c             interpolate in time the cross-covariance
              CALL interpolation(correltemp,nt,correlint,sinc,nt2,
     +             kpmin,kpmax,kmmin,kmmax,nthresh)

              DO lag=1,nt3
                 tp(lag)=0.
                 tm(lag)=0.
              ENDDO

              DO lag=-nthresh,nthresh
                 DO k=kpmin,kpmax
                    tp(lag+nthresh+1) = tp(lag+nthresh+1)+(correlint(k)
     +                 -correlrefint(k-lag))**2
                 ENDDO
                 DO k=kmmin,kmmax
                    tm(lag+nthresh+1) = tm(lag+nthresh+1)+(correlint(k)
     +                 -correlrefint(k-lag))**2
                 ENDDO
              ENDDO

              aminimump=10000.
              minipl=0
              DO k=1,nt3
                 IF (tp(k) .LT. aminimump) THEN 
                    aminimump=tp(k)
                    minipl=k
                 ENDIF
              ENDDO

              aminimumm=10000.
              miniml=0
              DO k=1,nt3
                 IF (tm(k) .LT. aminimumm) THEN 
                    aminimumm=tm(k)
                    miniml=k
                 ENDIF
              ENDDO

c             to avoid issues at the edges:
              IF(minipl .EQ. 1 .OR. minipl .EQ. 2) THEN 
                 minipl=(nt4+1)/2
              ENDIF
              IF(miniml .EQ. 1 .OR. miniml .EQ. 2) THEN 
                 miniml=(nt4+1)/2
              ENDIF
              IF(minipl .EQ. nt3 .OR. minipl .EQ. nt3-1) THEN 
                 minipl=nt3-2
              ENDIF
              IF(miniml .EQ. nt3 .OR. miniml .EQ. nt3-1) THEN 
                 miniml=nt3-2
              ENDIF


              DO k=1,nt4
                 tpc(k)=tp(minipl-(nt4+1)/2+k)
                 tmc(k)=tm(miniml-(nt4+1)/2+k)
                 alagscp(k)=alags(minipl-(nt4+1)/2+k)
                 alagscm(k)=alags(miniml-(nt4+1)/2+k)
              ENDDO

c             fit the parabola
              CALL parabola(alagscp,tpc,nt4,c0,c1,c2)
              IF(c2 .EQ.0) c2=1.e-5
              taup(i,j)  = -c1/2./c2
              CALL parabola(alagscm,tmc,nt4,c0,c1,c2)
              IF(c2 .EQ.0) c2=1.e-5
              taum(i,j)  =  c1/2./c2
              tau(i,j,1) = (taup(i,j)+taum(i,j))/2
              tau(i,j,2) =  taup(i,j)-taum(i,j)

              IF((ABS(tau(i,j,1)).GT.2*dt).AND.(i.GT.1).AND.(j.GT.1))
     +           tau(i,j,1)=0.5*(tau(i-1,j,1)+tau(i,j-1,1))
              IF((ABS(tau(i,j,2)).GT.2*dt).AND.(i.GT.1).AND.(j.GT.1)) 
     +           tau(i,j,2)=0.5*(tau(i-1,j,1)+tau(i,j-1,1))
              IF(ABS(tau(i,j,1)).GT.2*dt) tau(i,j,1)=0.   
              IF(ABS(tau(i,j,2)).GT.2*dt) tau(i,j,2)=0.   
           ENDDO

        ENDDO


c       compute the mean and square root of variance of the travel times
        ameanp=0.
        ameanm=0.
        DO i=1,nx
           DO j=1,ny
              ameanp=ameanp+taup(i,j)
              ameanm=ameanm+taum(i,j)
           ENDDO
        ENDDO
        ameanp=ameanp/nx/ny
        ameanm=ameanm/nx/ny

        rmsp=0.
        rmsm=0.
        DO i=1,nx
           DO j=1,ny
              rmsp=rmsp+(taup(i,j)-ameanp)**2
              rmsm=rmsm+(taum(i,j)-ameanm)**2
           ENDDO
        ENDDO
        rmsp=SQRT(rmsp/nx/ny)
        rmsm=SQRT(rmsm/nx/ny)

        threshold=threshold*MAX(rmsp,rmsm)
C        PRINT *,'mean values',ameanp,ameanm
C        PRINT *,'rms variations and threshold:',rmsp,rmsm,threshold

c       re-compute the outlier travel times except at the map edges
        DO i=2,nx-1
           DO j=2,nx-1
              IF( (ABS(taup(i,j)-ameanp) .GT. threshold) .OR. 
     +            (ABS(taum(i,j)-ameanm) .GT. threshold))THEN

                 DO k=1,nt
                    correltemp(k)=(correl(i,j,k)+correl(i-1,j,k)
     +           +correl(i+1,j,k)+correl(i,j-1,k)+correl(i-1,j-1,k)
     +           +correl(i+1,j-1,k)+correl(i,j+1,k)+correl(i-1,j+1,k)
     +           +correl(i+1,j+1,k))/9./amaxi
                 ENDDO
                 
                 CALL interpolation(correltemp,nt,correlint,sinc,nt2,
     +                kpmin,kpmax,kmmin,kmmax,nthresh)
                 
                 DO lag=1,nt3
                    tp(lag)=0.
                    tm(lag)=0.
                 ENDDO
                 
                 DO lag=-nthresh,nthresh
                    DO k=kpmin,kpmax
                       tp(lag+nthresh+1) = tp(lag+nthresh+1)+
     +                   (correlint(k)-correlrefint(k-lag))**2
                    ENDDO
                    DO k=kmmin,kmmax
                       tm(lag+nthresh+1) = tm(lag+nthresh+1)+
     +                   (correlint(k)-correlrefint(k-lag))**2
                    ENDDO
                 ENDDO

                 aminimump=10000.
                 minipl=0
                 DO k=1,nt3
                    IF (tp(k) .LT. aminimump) THEN 
                       aminimump=tp(k)
                       minipl=k
                    ENDIF
                 ENDDO

                 aminimumm=10000.
                 miniml=0
                 DO k=1,nt3
                    IF (tm(k) .LT. aminimumm) THEN 
                       aminimumm=tm(k)
                       miniml=k
                    ENDIF
                 ENDDO
                 
c                to avoid issues at the edges:
                 IF(minipl .EQ. 1 .OR. minipl .EQ. 2) THEN 
                    minipl=(nt4+1)/2
                 ENDIF
                 IF(miniml .EQ. 1 .OR. miniml .EQ. 2) THEN 
                    miniml=(nt4+1)/2
                 ENDIF
                 IF(minipl .EQ. nt3 .OR. minipl .EQ. nt3-1) THEN 
                    minipl=nt3-2
                 ENDIF
                 IF(miniml .EQ. nt3 .OR. miniml .EQ. nt3-1) THEN 
                    miniml=nt3-2
                 ENDIF


                 DO k=1,nt4
                    tpc(k)=tp(minipl-(nt4+1)/2+k)
                    tmc(k)=tm(miniml-(nt4+1)/2+k)
                    alagscp(k)=alags(minipl-(nt4+1)/2+k)
                    alagscm(k)=alags(miniml-(nt4+1)/2+k)
                 ENDDO

c                fit the parabola
                 CALL parabola(alagscp,tpc,nt4,c0,c1,c2)
                 IF(c2 .EQ. 0) c2=1.e-5
                 taup(i,j)  = -c1/2./c2
                 CALL parabola(alagscm,tmc,nt4,c0,c1,c2)
                 IF(c2 .EQ. 0) c2=1.e-5
                 taum(i,j)  =  c1/2./c2
                 tau(i,j,1) = (taup(i,j)+taum(i,j))/2.
                 tau(i,j,2) =  taup(i,j)-taum(i,j)

              IF((ABS(tau(i,j,1)).GT.2*dt).AND.(i.GT.1).AND.(j.GT.1))
     +           tau(i,j,1)=0.5*(tau(i-1,j,1)+tau(i,j-1,1))
              IF((ABS(tau(i,j,2)).GT.2*dt).AND.(i.GT.1).AND.(j.GT.1)) 
     +           tau(i,j,2)=0.5*(tau(i-1,j,1)+tau(i,j-1,1))
              IF(ABS(tau(i,j,1)).GT.2*dt) tau(i,j,1)=0.   
              IF(ABS(tau(i,j,2)).GT.2*dt) tau(i,j,2)=0.   
              ENDIF     
           ENDDO
        ENDDO

        RETURN

        END

c -------------------------------------------------------------------------------
c       uses the Whittaker-Shannon interpolation formula:
c       optimal interpolation for 1) band-limited signals and 2) when the Nyquist frequency
c       is larger than the bandwidth

        SUBROUTINE interpolation(f,nt0,fi,sinc,nt,kpmin,kpmax,kmmin,
     +             kmmax,nthresh)

        IMPLICIT REAL(a-h,o-z)
        INTEGER nt0,nt,i,j
        INTEGER kpmin,kpmax,kmmin,kmmax,nthresh
        REAL f(nt0),fi(nt)
        REAL sinc(nt,nt0)
        
c       we only care about the interpolated value around the time window
        DO i=kpmin-nthresh,kpmax+nthresh 
           fi(i)=0.
           DO j=1,nt0
              fi(i)=fi(i)+f(j)*sinc(i,j)
           ENDDO
        ENDDO
        DO i=kmmin-nthresh,kmmax+nthresh 
           fi(i)=0.
           DO j=1,nt0
              fi(i)=fi(i)+f(j)*sinc(i,j)
           ENDDO
        ENDDO

        RETURN

        END



c -------------------------------------------------------------------------------
c routine to fit for a parabola

        SUBROUTINE parabola(x,y,n,c0,c1,c2)

        IMPLICIT REAL(a-h,o-z)
        INTEGER n,i
        REAL x(n),y(n),a11,a12,a13,a21,a22,a23,a31,a32,a33
        REAL determinant,am11,am12,am13,am21,am22,am23,am31,am32,am33
        REAL c0,c1,c2,sy,sxy,sxxy

        a11=n*1.
        a12=0.
        a13=0.
        a23=0.
        a33=0.
        sy =0.
        sxy=0.
        sxxy=0.
        do i=1,n 
           a12=a12+x(i)
           a13=a13+x(i)**2
           a23=a23+x(i)**3
           a33=a33+x(i)**4
           sy =sy+y(i)
           sxy=sxy+x(i)*y(i)
           sxxy=sxxy+x(i)*x(i)*y(i)
        enddo
        a21=a12
        a22=a13
        a31=a13
        a32=a23

c       Compute the minors/cofactors/adjoint matrix M
        am11=a22*a33-a32*a23
        am21=(a21*a33-a31*a23)*(-1.)
        am31=a21*a32-a31*a22
        am12=(a12*a33-a32*a13)*(-1.)
        am22=a11*a33-a31*a13
        am32=(a11*a32-a31*a12)*(-1.)
        am13=a12*a23-a22*a13
        am23=(a11*a23-a21*a13)*(-1.)
        am33=a11*a22-a21*a12

        determinant=a11*a22*a33-a11*a32*a23-a12*a21*a33+a12*a31*a23
     +              +a13*a21*a32-a13*a31*a22

        IF(determinant .EQ. 0) determinant=1.e-5
        c0 = (am11*sy+am12*sxy+am13*sxxy)/determinant
        c1 = (am21*sy+am22*sxy+am23*sxxy)/determinant
        c2 = (am31*sy+am32*sxy+am33*sxxy)/determinant

        RETURN
        
        END

c ----------------------------------------------------------------------------

        SUBROUTINE window(kpmin,kpmax,kmmin,kmmax,nt,nt2,iii,dt,dt2)

        IMPLICIT REAL(a-h,o-z)
        INTEGER iii,k,nt,nt2,kpmin,kpmax,kmmin,kmmax
        REAL centrage,dt,dt2,time(nt2),width

        width=700.0

        DO k=1,nt2 
           time(k) = ((k-1)-(nt-1)/2*nt2/nt)*dt2
        ENDDO

        IF (iii .EQ. 1) THEN 
           centrage=810.
        ELSE IF (iii .EQ. 2) THEN
           centrage=1130.
        ELSE IF (iii .EQ. 3) THEN 
           centrage=1480.
        ELSE IF (iii .EQ. 4) THEN
           centrage=1750.
        ELSE IF (iii .EQ. 5) THEN
           centrage=1900.
        ELSE IF (iii .EQ. 6) THEN
           centrage=2050.
        ELSE IF (iii .EQ. 7) THEN
           centrage=2300.
        ELSE IF (iii .EQ. 8) THEN
           centrage=2520.
        ELSE IF (iii .EQ. 9) THEN
           centrage=2740.
        ELSE IF (iii .EQ. 10) THEN
           centrage=2990.
        ELSE IF (iii .EQ. 11) THEN
           centrage=3200.
        ENDIF
        
C        PRINT *,'CENTER OF TIME WINDOW',centrage

        kpmin=10000
        kmmin=10000
        kpmax=0
        kmmax=0

        DO k=1,nt2 
           IF ((ABS(time(k)-centrage) .LE. width) .AND. 
     +        (k .GE.(nt-1.)/2.*nt2/nt) .AND. (k .GT. kpmax)) THEN
              kpmax = k   
           ENDIF
           IF ((ABS(time(k)-centrage) .LE. width) .AND. 
     +        (k .GE.(nt-1.)/2.*nt2/nt) .AND. (k .LT. kpmin)) THEN
              kpmin = k   
           ENDIF
           IF ((ABS(time(k)+centrage) .LE. width) .AND. 
     +        (k .LT.(nt-1.)/2.*nt2/nt) .AND. (k .GT. kmmax)) THEN
              kmmax = k
           ENDIF
           IF ((ABS(time(k)+centrage) .LE. width) .AND. 
     +        (k .LT.(nt-1.)/2.*nt2/nt) .AND. (k .LT. kmmin)) THEN
              kmmin = k
           ENDIF
        ENDDO

        RETURN
        END

