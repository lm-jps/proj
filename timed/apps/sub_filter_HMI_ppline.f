      SUBROUTINE FILTER(dat,n1,n2,n3,velo,wid,spat_res)

      IMPLICIT REAL*4(a-h,o-z)
C      INCLUDE '/usr/local/include/fftw3.f'
C      INCLUDE '/usr/local/include.old/fftw3.f'
      INTEGER FFTW_R2HC
      PARAMETER (FFTW_R2HC=0)
      INTEGER FFTW_HC2R
      PARAMETER (FFTW_HC2R=1)
      INTEGER FFTW_DHT
      PARAMETER (FFTW_DHT=2)
      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)
      INTEGER FFTW_RODFT00
      PARAMETER (FFTW_RODFT00=7)
      INTEGER FFTW_RODFT01
      PARAMETER (FFTW_RODFT01=8)
      INTEGER FFTW_RODFT10
      PARAMETER (FFTW_RODFT10=9)
      INTEGER FFTW_RODFT11
      PARAMETER (FFTW_RODFT11=10)
      INTEGER FFTW_FORWARD
      PARAMETER (FFTW_FORWARD=-1)
      INTEGER FFTW_BACKWARD
      PARAMETER (FFTW_BACKWARD=+1)
      INTEGER FFTW_MEASURE
      PARAMETER (FFTW_MEASURE=0)
      INTEGER FFTW_DESTROY_INPUT
      PARAMETER (FFTW_DESTROY_INPUT=1)
      INTEGER FFTW_UNALIGNED
      PARAMETER (FFTW_UNALIGNED=2)
      INTEGER FFTW_CONSERVE_MEMORY
      PARAMETER (FFTW_CONSERVE_MEMORY=4)
      INTEGER FFTW_EXHAUSTIVE
      PARAMETER (FFTW_EXHAUSTIVE=8)
      INTEGER FFTW_PRESERVE_INPUT
      PARAMETER (FFTW_PRESERVE_INPUT=16)
      INTEGER FFTW_PATIENT
      PARAMETER (FFTW_PATIENT=32)
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)
      INTEGER FFTW_ESTIMATE_PATIENT
      PARAMETER (FFTW_ESTIMATE_PATIENT=128)
      INTEGER FFTW_BELIEVE_PCOST
      PARAMETER (FFTW_BELIEVE_PCOST=256)
      INTEGER FFTW_NO_DFT_R2HC
      PARAMETER (FFTW_NO_DFT_R2HC=512)
      INTEGER FFTW_NO_NONTHREADED
      PARAMETER (FFTW_NO_NONTHREADED=1024)
      INTEGER FFTW_NO_BUFFERING
      PARAMETER (FFTW_NO_BUFFERING=2048)
      INTEGER FFTW_NO_INDIRECT_OP
      PARAMETER (FFTW_NO_INDIRECT_OP=4096)
      INTEGER FFTW_ALLOW_LARGE_GENERIC
      PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
      INTEGER FFTW_NO_RANK_SPLITS
      PARAMETER (FFTW_NO_RANK_SPLITS=16384)
      INTEGER FFTW_NO_VRANK_SPLITS
      PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
      INTEGER FFTW_NO_VRECURSE
      PARAMETER (FFTW_NO_VRECURSE=65536)
      INTEGER FFTW_NO_SIMD
      PARAMETER (FFTW_NO_SIMD=131072)
      INTEGER FFTW_NO_SLOW
      PARAMETER (FFTW_NO_SLOW=262144)
      INTEGER FFTW_NO_FIXED_RADIX_LARGE_N
      PARAMETER (FFTW_NO_FIXED_RADIX_LARGE_N=524288)
      INTEGER FFTW_ALLOW_PRUNING
      PARAMETER (FFTW_ALLOW_PRUNING=1048576)
      INTEGER FFTW_WISDOM_ONLY
      PARAMETER (FFTW_WISDOM_ONLY=2097152)

      PARAMETER(PI=3.1415926)
      DIMENSION dat(n1,n2,n3)
      INTEGER n1,n2,n3
      REAL nu

      REAL ps_filt(n3)
      COMPLEX datc(n1/2+1,n2,n3)
      INTEGER*8 plan,pinver

      DO 5 i=1,20 
        DO 5 j=1,n2
          DO 5 k=1,n3
            dat(i,j,k)=dat(i,j,k)*SIN(i/20.*PI/2)
            dat(n1+1-i,j,k)=dat(n1+1-i,j,k)*SIN(i/20.*PI/2)
 5    CONTINUE
      DO 10 i=1,n1
        DO 10 j=1,20
          DO 10 k=1,n3
            dat(i,j,k)=dat(i,j,k)*SIN(j/20.*PI/2)
            dat(i,n2+1-j,k)=dat(i,n2+1-j,k)*SIN(j/20.*PI/2)
 10   CONTINUE
      DO 15 i=1,n1
        DO 15 j=1,n2
          DO 15 k=1,30
            dat(i,j,k)=dat(i,j,k)*SIN(k/30.*PI/2)
            dat(i,j,n3+1-k)=dat(i,j,n3+1-k)*SIN(k/30.*PI/2)
 15   CONTINUE

      wid=SQRT(wid/2.35482)
C      WRITE(*,*) 'Forward FFT begins......'
      CALL SFFTW_PLAN_DFT_R2C_3D(plan,n1,n2,n3,dat,datc,FFTW_ESTIMATE)
      CALL SFFTW_EXECUTE(plan)
      CALL SFFTW_DESTROY_PLAN(plan)
      
C      WRITE(*,*) 'Phaspeed Filtering begins......'
      DO 50 i=1,n1/2+1
        DO 40 j=1,n2
          hi=i-1.
          IF(j.LE.n2/2+1) THEN
            hj=j-1.
          ELSE
            hj=n2+1.-j
          ENDIF
          hk=SQRT((hi/n1)**2+(hj/n2)**2)/(spat_res*12.15)*697*2*PI
C 0.12 degr/pixel is spatial resolution. need to change for difference case
          DO 30 l=1,n3
            IF(l.LE.n3/2+1) THEN
              nu=(l-1.)/n3/45.*1000*1000
            ELSE
              nu=(n3+1.-l)/n3/45.*1000*1000
            ENDIF
            IF(nu.LE.1200) THEN
              ps_filt(l)=0.
            ELSE
              IF(hk.EQ.0.) THEN
                ps_filt(l)=0.
              ELSE
                ps_filt(l)=EXP(-(nu/hk-velo)**2/(2.*wid**2))
              ENDIF
              q0=SQRT(hk/100.)*1000.
              q1=SQRT(hk/80.)*1000.+350
              IF(nu.LT.(q0+q1)/2.) THEN
                ps_filt(l)=EXP(-(nu-(q0+q1)/2.)**2*8.*LOG(10.)/(q1
     *             -q0)**2)*EXP(-(nu/hk-velo)**2/(2.*wid**2))
              ENDIF
            ENDIF
            IF((nu.GT.1200).AND.(nu.LE.1800)) 
     *        ps_filt(l)=ps_filt(l)*SIN((nu-1200.)/600.*PI/2)
            IF((nu.GT.6000).AND.(nu.LE.6500))
     *        ps_filt(l)=ps_filt(l)*SIN((nu-6000)/500.*PI/2)
            IF(nu.GT.6500)
     *        ps_filt(l)=0.
C	write(*,*) i,j,l,hk,nu,ps_filt(l)
            datc(i,j,l)=datc(i,j,l)*ps_filt(l)
 30       CONTINUE
 40     CONTINUE
 50   CONTINUE

C      WRITE(*,*) 'Inverse FFT begins......'
      CALL SFFTW_PLAN_DFT_C2R_3D(pinv,n1,n2,n3,datc,dat,FFTW_ESTIMATE)
      CALL SFFTW_EXECUTE(pinv)
      CALL SFFTW_DESTROY_PLAN(pinv)

c	open(3,file='filtered_data2.dat',status='new',form='unformatted')
c	write(3) dat
c	close(3)

      RETURN
      END
