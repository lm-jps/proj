CC This is te program to do inversions for 3D velocities by use of the 
CC MultiChannel Deconvolution                    ......July 24, 2002

      SUBROUTINE INVER_V_SUB(aa_oi,aa_we,aa_ns,alon0,alat0,aavx,aavy,
     +           aavz,vx,vy,vz,kn)
C  kn = kernel number. kn=1 is for ray kernel and kn=2 is for Born kernel.
      IMPLICIT REAL*4(a-h,o-z)
      PARAMETER(nrays=11,nx1=157,ny1=157,nz1=11,nnx=256,nny=256,nnz=11,
     *      lm=nnx,ln=nny,lwork=400,nset=3,mm=nset*nz1,nn=nset*nrays)
      DIMENSION aavx(nx1,ny1,nz1,nrays,nset),aavy(nx1,ny1,nz1,
     *      nrays,nset),aavz(nx1,ny1,nz1,nrays,nset)
      DIMENSION avxt(lm,ln),avyt(lm,ln),avzt(lm,ln)
      DIMENSION doi(lm,ln),dwe(lm,ln),dns(lm,ln),delv(nnx,nny,nnz)
      DIMENSION avx_exp(lm,ln,nz1,nrays,nset)
      DIMENSION avy_exp(lm,ln,nz1,nrays,nset)
      DIMENSION avz_exp(lm,ln,nz1,nrays,nset)
      DIMENSION vx_exp(lm,ln),vy_exp(lm,ln),vz_exp(lm,ln)
      DIMENSION vv(nrays),v(mm,nn),velo_param(nnz),trans_born(nnz)
      DIMENSION velo_param2(nnz)
      DIMENSION ang1(nnx,nny),ang2(nnx,nny)
      INTEGER ipiv(mm)
      REAL*4 aa_oi(nnx,nny,nrays),aa_we(nnx,nny,nrays)
      REAL*4 aa_ns(nnx,nny,nrays),vx(nnx,nny,nnz)
      REAL*4 vy(nnx,nny,nnz),vz(nnx,nny,nnz)
      COMPLEX soi(lm/2+1,ln),swe(lm/2+1,ln),sns(lm/2+1,ln)
      COMPLEX s_oi(lm/2+1,ln,nrays)
      COMPLEX s_we(lm/2+1,ln,nrays),s_ns(lm/2+1,ln,nrays)
      COMPLEX savxt(lm/2+1,ln),s_avx(lm/2+1,ln,nz1,nrays,nset)
      COMPLEX savyt(lm/2+1,ln),s_avy(lm/2+1,ln,nz1,nrays,nset)
      COMPLEX savzt(lm/2+1,ln),s_avz(lm/2+1,ln,nz1,nrays,nset)
      COMPLEX aa(mm,nn),aat(mm,nn),yy(nn),svx(lm/2+1,ln,nz1)
      COMPLEX bb(nz1,nrays),svy(lm/2+1,ln,nz1),svz(lm/2+1,ln,nz1)
      COMPLEX work(lwork),temp(mm,nn),h(mm,nn),res(mm)
      COMPLEX ONE, ZERO
C      DATA vv/6.25,8.83,10.9,14.2,18.3,20.8,22.7,24.6,27.,30.,34./
      DATA vv/8.25, 13.10,18.25,23.82,31.28,38.79,45.41,52.82,60.51,
     +        67.50,73.84/
      DATA velo_param/15675.0,11266.0,20440.0,32157.0,29716.0,48487.5,
     +        1.,1.,1.,1.,1./
      DATA trans_born/1.5,1.15,0.8,1.,1.,1.,1.,1.,1.,1.,1./
      DATA velo_param2/23512.5,12955.9,16352.0,32157.0,29716.0,48487.5,
     +        1.,1.,1.,1.,1./
      ONE=CMPLX(1.,0.)
      ZERO=CMPLX(0.,0.)

      scale=SQRT(1./lm/ln)

      DO 20 k=1,nrays
        DO 18 i=1,lm
        DO 18 j=1,ln
          doi(i,j)=aa_oi(i,j,k)
          dwe(i,j)=aa_we(i,j,k)
          dns(i,j)=aa_ns(i,j,k)
 18     CONTINUE
        CALL fft_forw(doi,soi,lm,ln)
        CALL fft_forw(dwe,swe,lm,ln)
        CALL fft_forw(dns,sns,lm,ln)
        DO 24 i=1,lm/2+1
        DO 24 j=1,ln
          s_oi(i,j,k)=soi(i,j)*scale
          s_we(i,j,k)=swe(i,j)*scale
          s_ns(i,j,k)=sns(i,j)*scale
 24     CONTINUE
 20   CONTINUE

CC Expanding ac into lm and ln from nx1 and ny1 
      CALL expand_av(aavx,avx_exp,nx1,ny1,nz1,nrays,nset,lm,ln)
      CALL expand_av(aavy,avy_exp,nx1,ny1,nz1,nrays,nset,lm,ln)
      CALL expand_av(aavz,avz_exp,nx1,ny1,nz1,nrays,nset,lm,ln)
      DO 30 iset=1,nset
      DO 30 k=1,nrays
      DO 30 kk=1,nz1
        DO 32 i=1,lm
        DO 32 j=1,ln
          avxt(i,j)=avx_exp(i,j,kk,k,iset)
          avyt(i,j)=avy_exp(i,j,kk,k,iset)
          avzt(i,j)=avz_exp(i,j,kk,k,iset)
 32     CONTINUE
        CALL fft_forw(avxt,savxt,lm,ln)
        CALL fft_forw(avyt,savyt,lm,ln)
        CALL fft_forw(avzt,savzt,lm,ln)
        DO 34 i=1,lm/2+1
        DO 34 j=1,ln
          s_avx(i,j,kk,k,iset)=savxt(i,j)*scale
          s_avy(i,j,kk,k,iset)=savyt(i,j)*scale
          s_avz(i,j,kk,k,iset)=savzt(i,j)*scale
 34     CONTINUE 
 30   CONTINUE

      DO 50 i=1,lm/2+1
      DO 50 j=1,ln

        IF(j.LE.ln/2) THEN
          ppak=((1.*i/lm)**2+(1.*j/ln)**2)**0.5
        ELSE
          ppak=((1.*i/lm)**2+(1-1.*j/ln)**2)**0.5
        ENDIF
        epson1_v=0.06
        epson2_v=0.06
        epson1_h=0.3
        epson2_h=0.3
        DO 40 ii=1,mm
          DO 40 jj=1,nn
            v(ii,jj)=0.
 40     CONTINUE
        DO 45 ii=1,nz1
C          v(ii,ii)=epson1_v**2*(vv(ii)/vv(1))**0.5+(epson1_h*ppak)**2
C          v(ii+nz1,ii+nz1)=v(ii,ii)
C          v(ii+2*nz1,ii+2*nz1)=epson2_v**2*(vv(ii)/vv(1))**0.5+
C     +                         (epson2_h*ppak)**2
          v(ii,ii)=epson1_v**2*(vv(ii)/vv(1))**1+(epson1_h*ppak)**2
          v(ii+nz1,ii+nz1)=v(ii,ii)
          v(ii+2*nz1,ii+2*nz1)=epson2_v**2*(vv(ii)/vv(1))**1+
     +                         (epson2_h*ppak)**2
 45     CONTINUE

        DO 52 jj=1,nrays
          yy(jj)=s_we(i,j,jj)
          yy(jj+nrays)=s_ns(i,j,jj)
          yy(jj+2*nrays)=s_oi(i,j,jj)
          DO 53 ii=1,nz1
            aat(ii,jj)=s_avx(i,j,jj,ii,1)
            aat(ii+nz1,jj)=s_avx(i,j,jj,ii,2)
            aat(ii+2*nz1,jj)=s_avx(i,j,jj,ii,3)
            aat(ii,jj+nrays)=s_avy(i,j,jj,ii,1)
            aat(ii+nz1,jj+nrays)=s_avy(i,j,jj,ii,2)
            aat(ii+2*nz1,jj+nrays)=s_avy(i,j,jj,ii,3)
            aat(ii,jj+2*nrays)=s_avz(i,j,jj,ii,1)
            aat(ii+nz1,jj+2*nrays)=s_avz(i,j,jj,ii,2)
            aat(ii+2*nz1,jj+2*nrays)=s_avz(i,j,jj,ii,3)
 53       CONTINUE
 52     CONTINUE
        DO 54 ii=1,nn
        DO 54 jj=1,mm
          aa(ii,jj)=aat(ii,jj)
 54     CONTINUE
CC  Calculate temp(*,*)=CONJ(TRANSPOSE(aa))*aa
        CALL cgemm('C','N',mm,nn,mm,ONE,aa,mm,aa,mm,ZERO,temp,mm)
        DO 55 ii=1,nn
        DO 55 jj=1,mm
          temp(ii,jj)=temp(ii,jj)+v(ii,jj)
 55     CONTINUE
CC Calculate the inverse of temp
        CALL cgetrf(nn,mm,temp,nn,ipiv,info)
        CALL cgetri(mm,temp,nn,ipiv,work,lwork,info)
CC Calculate the h=(temp)^(-1)*aa^h
        CALL cgemm('N','C',mm,nn,mm,ONE,temp,mm,aa,mm,ZERO,h,mm)
CC Calculate the res=h*yy
        CALL cgemv('N',mm,nn,ONE,h,mm,yy,1,ZERO,res,1)
        DO 56 ii=1,nz1
          svx(i,j,ii)=res(ii)
          svy(i,j,ii)=res(ii+nz1)
          svz(i,j,ii)=res(ii+2*nz1)
 56     CONTINUE
 50   CONTINUE 

      CALL correct_direc(nnx,nny,alon0,alat0,ang1,ang2)
      DO 60 k=1,nnz
        DO 62 i=1,lm/2+1
        DO 62 j=1,ln
          savxt(i,j)=svx(i,j,k)
          savyt(i,j)=svy(i,j,k)
          savzt(i,j)=svz(i,j,k)
 62     CONTINUE
        CALL fft_inver(savxt,vx_exp,lm,ln)
        CALL fft_inver(savyt,vy_exp,lm,ln)
        CALL fft_inver(savzt,vz_exp,lm,ln)
        DO 64 i=1,nnx
        DO 64 j=1,nny
        IF(kn .EQ. 1) THEN
          vx(i,j,k)=(vx_exp(i,j)*COS(ang1(i,j))+
     +           vy_exp(i,j)*SIN(ang1(i,j)))*scale*scale*velo_param(k)
          vy(i,j,k)=(vy_exp(i,j)*COS(-ang2(i,j))+
     +           vx_exp(i,j)*SIN(-ang2(i,j)))*scale*scale*velo_param(k)
          vz(i,j,k)=vz_exp(i,j)*scale*scale*velo_param(k)
        ELSE 
          vx(i,j,k)=(vx_exp(i,j)*COS(ang1(i,j))+
     +           vy_exp(i,j)*SIN(ang1(i,j)))*scale*scale*velo_param2(k)
          vy(i,j,k)=(vy_exp(i,j)*COS(-ang2(i,j))+
     +           vx_exp(i,j)*SIN(-ang2(i,j)))*scale*scale*velo_param2(k)
          vz(i,j,k)=vz_exp(i,j)*scale*scale*velo_param2(k)
        ENDIF
 64     CONTINUE
 60   CONTINUE

      END
        

CCCCCCC       Subroutine for the forward FFT       CCCCCCCCCCCCCC
C      SUBROUTINE fft_forw(a,c,m,n)
C      IMPLICIT REAL*4(a-h,o-z)
C      INCLUDE '/usr/local/include.old/fftw3.f'
C
C      INTEGER*8 PLAN
C      REAL a(m,n)
C      COMPLEX c(m/2+1,n)
C
C      CALL SFFTW_PLAN_DFT_R2C_2D(plan,m,n,a,c,FFTW_ESTIMATE)
C      CALL SFFTW_EXECUTE(plan)
C      CALL SFFTW_DESTROY_PLAN(plan)
C
C      RETURN 
C      END

CCCCCCC       Subroutine for the inverse FFT       CCCCCCCCCCCCCC
C      SUBROUTINE fft_inver(c,a,m,n)
C      IMPLICIT REAL*4 (a-h,o-z)
C      INCLUDE '/usr/local/include.old/fftw3.f'
C
C      INTEGER*8 PLAN
C      REAL a(m,n)
C      COMPLEX c(m/2+1,n)
C
C      CALL SFFTW_PLAN_DFT_C2R_2D(plan,m,n,c,a,FFTW_ESTIMATE)
C      CALL SFFTW_EXECUTE(plan)
C      CALL SFFTW_DESTROY_PLAN(plan)
C
C      RETURN
C      END


CCCCCC       subroutine to expand the ac into lm x ln    CCCCCCCCCCCCC
      SUBROUTINE expand_av(temp,b,mm,nn,nz,nrays,nset,lm,ln)
      IMPLICIT REAL*4(a-h,o-z)
      DIMENSION temp(mm,nn,nz,nrays,nset),b(lm,ln,nz,nrays,nset)

      DO 10 iset=1,nset
      DO 10 k=1,nrays
      DO 10 kk=1,nz
        DO 11 i=1,mm/2+1
        DO 11 j=1,nn/2+1
          b(i,j,kk,k,iset)=temp(mm/2+2-i,nn/2+2-j,kk,k,iset)
 11     CONTINUE
        DO 12 i=1,mm/2+1
        DO 12 j=1,nn/2
          b(i,ln+1-j,kk,k,iset)=temp(mm/2+2-i,nn/2+1+j,kk,k,iset)
 12     CONTINUE
        DO 13 i=1,mm/2
        DO 13 j=1,nn/2+1
          b(lm+1-i,j,kk,k,iset)=temp(mm/2+1+i,nn/2+2-j,kk,k,iset)
 13     CONTINUE
        DO 14 i=1,mm/2
        DO 14 j=1,nn/2
          b(lm+1-i,ln+1-j,kk,k,iset)=temp(mm/2+1+i,nn/2+1+j,kk,k,iset)
 14     CONTINUE
 10   CONTINUE
      
      RETURN
      END

CCCCC	subroutine to correct directions after inversions   CCCCCCC
      SUBROUTINE correct_direc(nx,ny,alon0,alat0,ang1,ang2)
      IMPLICIT REAL*4(a-h,o-z)
      PARAMETER(PI=3.1415926)
      DIMENSION alon(nx,ny),alat(nx,ny),ang1(nx,ny),ang2(nx,ny)
CC ang1 is corresponding to angle relative to horizontal, and ang2 is
CC relative to the vertical

CC R_sun is the solar radius in a unit of 0.12 heliographic degree.
      Rsun=478. 
      alon1=alon0*PI/180.
      alat1=alat0*PI/180.
      DO 10 i=1,nx
        DO 10 j=1,ny
          xx=(i-128.5)/Rsun 
          yy=(j-128.5)/Rsun
          cc=SQRT(xx**2+yy**2)
          IF(cc.LE.PI/2) THEN
            tt=ASIN(COS(cc)*SIN(alat1)+yy*SIN(cc)*COS(alat1)/cc)
            pp=alon1+ATAN(xx*SIN(cc)/(cc*COS(alat1)*COS(cc)-
     *              yy*SIN(alat1)*SIN(cc)))
          ELSE
            tt=0.
            pp=0.
          ENDIF
          alat(i,j)=tt 
          alon(i,j)=pp
 10   CONTINUE

      DO 20 i=2,nx-1
        DO 20 j=2,ny-1
          dx=(alat(i+1,j)-alat(i-1,j))*180/PI
          dy=(alon(i,j+1)-alon(i,j-1))*180/PI
          ang1(i,j)=ATAN(-dx/0.12/2)
          ang2(i,j)=ATAN(-dy/0.12/2)
  20  CONTINUE
      DO 30 j=1,ny
        ang1(1,j)=ang1(2,j)
        ang1(nx,j)=ang1(nx-1,j)
        ang2(1,j)=ang2(2,j)
        ang2(nx,j)=ang2(nx-1,j)
  30  CONTINUE
      DO 35 i=1,nx
        ang1(i,1)=ang1(i,2)
        ang1(i,ny)=ang1(i,ny-1)
        ang2(i,1)=ang2(i,2)
        ang2(i,ny)=ang2(i,ny-1)
  35  CONTINUE

      RETURN
      END 
