      SUBROUTINE DO_FITTING(rr,n1,n2,num,n_cor,n_start,coef_temp,n_coe,
     *                      ia,apos,aneg)

      IMPLICIT REAL*4(a-h,o-z)
      PARAMETER (nca=8,tol=1.E-4)
      DIMENSION rr(num,n2,n1),apos(n1,n2),aneg(n1,n2),ia(n_coe)
      DIMENSION x(n_cor),y(n_cor),coef_temp(n_coe),coef(n_coe)

      DO 5 ll=1,n_cor
        x(ll)=n_start-1+ll
 5    CONTINUE

      DO 10 i=1,n1
        DO 10 j=1,n2 

          CALL SCOPY(n_cor,rr(100+n_start,j,i),1,y,1)
          amax=GET_MAX(y,n_cor)
          CALL SSCAL(n_cor,1./amax,y,1)
          CALL SCOPY(n_coe,coef_temp,1,coef,1)
          CALL LMFIT(x,y,n_cor,coef,ia,n_coe,nca,tol)
          k=1
          DO 15 WHILE((coef(5).GT.coef_temp(5)+1.7) .OR. (coef(5).LT.
     +                 coef_temp(5)-1.7))
            CALL SSCAL(n_cor,0.,y,1)
            DO 20 ii=i-k,i+k
            DO 20 jj=j-k,j+k
              IF((ii.GE.1).AND.(ii.LE.n1).AND.(jj.GE.1).AND.(jj.LE.n2))
     +          CALL SAXPY(n_cor,1.,rr(100+n_start,jj,ii),1,y,1)
 20         CONTINUE
            amax=GET_MAX(y,n_cor)
            CALL SSCAL(n_cor,1./amax,y,1)
            CALL SCOPY(n_coe,coef_temp,1,coef,1)
            CALL LMFIT(x,y,n_cor,coef,ia,n_coe,nca,tol)
            k=k+1
 15       CONTINUE
          apos(i,j)=coef(5)

          CALL SCOPY(n_cor,rr(100-n_start-n_cor+1,j,i),-1,y,1)
          amax=GET_MAX(y,n_cor)
          CALL SSCAL(n_cor,1./amax,y,1)
          CALL SCOPY(n_coe,coef_temp,1,coef,1)
          CALL LMFIT(x,y,n_cor,coef,ia,n_coe,nca,tol)
          k=1
          DO 25 WHILE((coef(5).GT.coef_temp(5)+1.7) .OR. (coef(5).LT.
     +                 coef_temp(5)-1.7))
            CALL SSCAL(n_cor,0.,y,1)
            DO 30 ii=i-k,i+k
            DO 30 jj=j-k,j+k
              IF((ii.GE.1).AND.(ii.LE.n1).AND.(jj.GE.1).AND.(jj.LE.n2))
     +         CALL SAXPY(n_cor,1.,rr(100-n_start-n_cor+1,jj,ii),-1,y,1)
 30         CONTINUE
            amax=GET_MAX(y,n_cor)
            CALL SSCAL(n_cor,1./amax,y,1)
            CALL SCOPY(n_coe,coef_temp,1,coef,1)
            CALL LMFIT(x,y,n_cor,coef,ia,n_coe,nca,tol)
            k=k+1
 25       CONTINUE
          aneg(i,j)=coef(5)

 10   CONTINUE

      RETURN
      END
