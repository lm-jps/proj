CC This is the program to do inversions for sound-speed perturbation
CC   by use of the MultiChannel Deconvolution, using only t_oi 
CC   observation                                      ......July 24, 2002
C
C  Modified 2011.03.18 by R Bogart to accept kernel input from calling
C    routine
C
      SUBROUTINE INVER_CS_SUB(aa_oi, ac, cs)
      IMPLICIT REAL*4(a-h,o-z)
      PARAMETER(nrays=11,nx1=157,ny1=157,nz1=11,nnx=256,nny=256,nnz=11,
     *          lm=nnx,ln=nny,lwork=400)
      DIMENSION ac(nx1,ny1,nz1,nrays),dtau(lm,ln)
      DIMENSION ac_exp(lm,ln,nz1,nrays),actt(lm,ln)
      DIMENSION cs_exp(lm,ln),vv(nrays),v(nz1,nrays)
      INTEGER ipiv(nrays),i,j,k,ii,jj,kk
      DIMENSION aa_oi(nnx,nny,nrays),cs(nnx,nny,nnz)
      COMPLEX stau(lm/2+1,ln),s_oi(lm/2+1,ln,nrays)
      COMPLEX sactt(lm/2+1,ln),s_ac(lm/2+1,ln,nz1,nrays)
      COMPLEX aa(nz1,nrays),yy(nrays),scs(lm/2+1,ln,nz1)
      COMPLEX bb(nz1,nrays)
      COMPLEX work(lwork),temp(nz1,nrays),h(nz1,nrays),res(nz1)
      COMPLEX ONE, ZERO
      DATA vv/8.25, 13.10,18.25,23.82,31.28,38.79,45.41,52.82,60.51,
     +        67.50,73.84/
      ONE=CMPLX(1.,0.)
      ZERO=CMPLX(0.,0.)

      scale=SQRT(1./lm/ln)

      DO 20 k=1,nrays
        DO 18 i=1,lm
        DO 18 j=1,ln
          dtau(i,j)=aa_oi(i,j,k)
 18     CONTINUE
        CALL fft_forw(dtau,stau,lm,ln)
        DO 24 i=1,lm/2+1
        DO 24 j=1,ln
          s_oi(i,j,k)=stau(i,j)*scale
 24     CONTINUE
 20   CONTINUE

CC Expanding ac into lm and ln from nx1 and ny1 
      CALL expand_ac(ac,ac_exp,nx1,ny1,nz1,nrays,lm,ln)
      DO 30 k=1,nrays
      DO 30 kk=1,nz1
        DO 32 i=1,lm
        DO 32 j=1,ln
          actt(i,j)=ac_exp(i,j,kk,k)
 32     CONTINUE
        CALL fft_forw(actt,sactt,lm,ln)
        DO 34 i=1,lm/2+1
        DO 34 j=1,ln
          s_ac(i,j,kk,k)=sactt(i,j)*scale
 34     CONTINUE 
 30   CONTINUE

      DO 50 i=1,lm/2+1
      DO 50 j=1,ln

        IF(j.LE.ln/2) THEN
          ppak=((1.*i/lm)**2+(1.*j/ln)**2)**0.5
        ELSE
          ppak=((1.*i/lm)**2+(1-1.*j/ln)**2)**0.5
        ENDIF
        epson_v=0.02
        epson_h=0.3
        DO 40 ii=1,nz1
          DO 40 jj=1,nz1
            v(ii,jj)=0.
 40     CONTINUE
        DO 41 ii=1,nz1
C          v(i,i)=epson_v**2*(vv(i)/vv(1))**3
          v(ii,ii)=epson_v**2*(vv(ii)/vv(1))**0.5+(epson_h*ppak)**2
 41     CONTINUE

        DO 52 jj=1,nrays
          yy(jj)=s_oi(i,j,jj)
          DO 54 ii=1,nz1
            aa(ii,jj)=s_ac(i,j,jj,ii)
 54       CONTINUE
 52     CONTINUE
CC  Calculate temp(*,*)=CONJ(TRANSPOSE(aa))*aa
        CALL cgemm('C','N',nz1,nrays,nz1,ONE,aa,nz1,aa,nz1,ZERO,
     *               temp,nz1)
        DO 55 ii=1,nrays
        DO 55 jj=1,nz1
          temp(ii,jj)=temp(ii,jj)+v(ii,jj)
 55     CONTINUE
CC Calculate the inverse of temp
        CALL cgetrf(nrays,nz1,temp,nrays,ipiv,info)
        CALL cgetri(nz1,temp,nrays,ipiv,work,lwork,info)
CC Calculate the h=(temp)^(-1)*aa^h
        CALL cgemm('N','C',nz1,nrays,nz1,ONE,temp,nz1,aa,nz1,ZERO,h,nz1)
CC Calculate the res=h*yy
        CALL cgemv('N',nz1,nrays,ONE,h,nz1,yy,1,ZERO,res,1)
        DO 56 ii=1,nz1
          scs(i,j,ii)=res(ii)
 56     CONTINUE
 50   CONTINUE 

      DO 60 k=1,nnz
        DO 62 i=1,lm/2+1
        DO 62 j=1,ln
          sactt(i,j)=scs(i,j,k)
 62     CONTINUE
        DO 63 i=1,nnx
        DO 63 j=1,nny
          cs_exp(i,j)=0.
 63     CONTINUE
        CALL fft_inver(sactt,cs_exp,lm,ln)
        DO 64 i=1,nnx
        DO 64 j=1,nny
          cs(i,j,k)=cs_exp(i,j)*scale*scale
 64     CONTINUE
 60   CONTINUE

      RETURN
      END
        

CCCCCCC       Subroutine for the forward FFT       CCCCCCCCCCCCCC
      SUBROUTINE fft_forw(a,c,m,n)
      IMPLICIT REAL*4(a-h,o-z), INTEGER(i-n)
C      INCLUDE 'fftw3.f'
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)

      DIMENSION a(m,n)
      COMPLEX c(m/2+1,n)
      INTEGER m,n
      INTEGER*8 PLAN

      CALL SFFTW_PLAN_DFT_R2C_2D(plan,m,n,a,c,FFTW_ESTIMATE)
      CALL SFFTW_EXECUTE_DFT_R2C(plan,a,c)
      CALL SFFTW_DESTROY_PLAN(plan)

      RETURN 
      END

CCCCCCC       Subroutine for the inverse FFT       CCCCCCCCCCCCCC
      SUBROUTINE fft_inver(c,a,m,n)
      IMPLICIT REAL*4(a-h,o-z), INTEGER(i-n)
C      INCLUDE 'fftw3.f'
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)

      DIMENSION a(m,n)
      COMPLEX c(m/2+1,n)
      INTEGER m,n
      INTEGER*8 PLAN

      CALL SFFTW_PLAN_DFT_C2R_2D(plan,m,n,c,a,FFTW_ESTIMATE)
      CALL SFFTW_EXECUTE_DFT_C2R(plan,c,a)
      CALL SFFTW_DESTROY_PLAN(plan)

      RETURN
      END


CCCCCC       subroutine to expand the ac into lm x ln    CCCCCCCCCCCCC
      SUBROUTINE expand_ac(temp,b,mm,nn,nz,nrays,lm,ln)
      IMPLICIT REAL(a-h,o-z)
      DIMENSION temp(mm,nn,nz,nrays),b(lm,ln,nz,nrays)

      DO 10 k=1,nrays
      DO 10 kk=1,nz
        DO 11 i=1,mm/2+1
        DO 11 j=1,nn/2+1
          b(i,j,kk,k)=temp(mm/2+2-i,nn/2+2-j,kk,k)
 11     CONTINUE
        DO 12 i=1,mm/2+1
        DO 12 j=1,nn/2
          b(i,ln+1-j,kk,k)=temp(mm/2+2-i,nn/2+1+j,kk,k)
 12     CONTINUE
        DO 13 i=1,mm/2
        DO 13 j=1,nn/2+1
          b(lm+1-i,j,kk,k)=temp(mm/2+1+i,nn/2+2-j,kk,k)
 13     CONTINUE
        DO 14 i=1,mm/2
        DO 14 j=1,nn/2
          b(lm+1-i,ln+1-j,kk,k)=temp(mm/2+1+i,nn/2+1+j,kk,k)
 14     CONTINUE
 10   CONTINUE
      
      RETURN
      END
