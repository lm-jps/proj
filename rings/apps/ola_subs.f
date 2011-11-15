c       =======================================================
c              SUBROUTINES AND FUNCTIONS
c       =======================================================

c       -----------------------------
        subroutine fwidth(npt,np,num,rad,avc,widthc,bverb)
cc        call fwidth(npt,np,num,rad,avc,idthc)
        implicit real*8(a,c-h,o-z)
CC      IMP. log. now 'b', not 'q' ('q' was already taken)
        implicit logical(b)
        parameter(nn=4000)
        parameter(nx0=100)
        dimension rad(npt), avc(npt), widthc(3)
        dimension dx(nn), sumker(nn), qq(3)
        dimension fb(10),y(nn)
        
        q1=0.25d0
        q2=0.5d0
        q3=0.75d0


        do 100 j=1,np-1
                dx(j)=rad(j)-rad(j+1)
100     continue

        k=1
        sum=0.d0
        sumker(1)=0.d0
                do 300 j=1,np-1
                sker=0.5d0*(avc(j)+avc(j+1))
                sum=sum+dx(j)*sker
                sumker(j+1)=sum
cc      write(72,*)rad(j), avc(j),sumker(j+1)
300     continue
200     continue

        j=1
        anor=sumker(np)
        if(bverb)then
          print*, 'anor ', anor
        endif
        do 400 i=1,np
        sumker(i)=sumker(i)/anor
cc      write(82,*)rad(i), sumker(i)
400     continue
        if(bverb) print*, sumker(i),sumker(np)
410     continue

        ntab=np
        reps=1.d-5

        j=1
                do 1010 i=1,np
                y(i)=sumker(i)
1010    continue
        nuse=4
      call DIVDIF(q1,y,rad,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
        widthc(1)=fb(nuse)
        nuse=4
      call DIVDIF(q2,y,rad,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
        widthc(2)=fb(nuse)
        nuse=4
      call DIVDIF(q3,y,rad,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
        widthc(3)=fb(nuse)
1000    continue
        if(bverb) print*,'wid ', widthc

        return
        end
c       -------------------------
c       SUBROUTINE WHICH CALLS LAPACK ROUTINES
c       =====================================

        subroutine solve(la,lb,norder,nrhs,a,b,ierr,ipiv,qverb)
        implicit real*8(a-h,o-p,r-z)
        implicit logical(q)
        character*1 trans, uplo

        dimension a(la,norder), b(lb,nrhs), ipiv(norder)
        dimension work (10000)
        
        lwork=10000
        uplo='u'
        call dsytrf(uplo, norder, a ,la, ipiv, work, lwork, ier)
        if(ier.ne.0)then
        print*,'error in transformation', ier
        endif

        trans='N'
        uplo='u'
        call dsytrs(uplo,norder,nrhs,a,la,ipiv,b,lb,ierr)
        if(qverb) print*, 'solution ier', ierr
        return
        end


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

c       --------------------------------------------------

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

