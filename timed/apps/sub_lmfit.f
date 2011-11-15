      FUNCTION  get_max(a,n)
CC     to get the maximum value of the array a(n)
      REAL get_max, a(n)
      b=0.
      DO 3 i=1, n
         IF(a(i) .GT. b) b=a(i)
 3    ENDDO
      get_max=b
      RETURN
      END
      
      SUBROUTINE lmfit(x,y,ndata,a,ia,na,nca,tol)
      REAL x(ndata), y(ndata), a(na)
      DIMENSION sig(ndata), ia(na),covar(nca,nca), alpha(nca,nca)
      alamda=-0.1
      DO 8 i=1, nca
	DO 8 j=1, nca
	  covar(i,j)=0.
	  alpha(i,j)=0.
 8    CONTINUE
      DO 9 i=1, ndata
	sig(i)=1.
 9    CONTINUE
      itmax=0
      chisq=1.
      DO 10 WHILE(SQRT(ABS(chisq)).GT.tol)
      CALL mrqmin(x,y,sig,ndata,a,ia,na,covar,alpha,nca,
     *     chisq,alamda)
      itmax=itmax+1
      IF(itmax.GE.35) GOTO 13
 10   CONTINUE
 13   alamda=0.
      CALL mrqmin(x,y,sig,ndata,a,ia,na,covar,alpha,nca,
     *     chisq,alamda)
      RETURN
      END

      SUBROUTINE funcs(x,a,y,dyda,na)
      REAL x,y,a(na),dyda(na)
      bx=EXP(-a(2)**2/4.*((x-a(3))**2))
      cx=COS(a(4)*(x-a(5)))
      y=a(1)*bx*cx
      dyda(1)=bx*cx
      dyda(2)=-0.5*a(1)*a(2)*((x-a(3))**2)*bx*cx
      dyda(3)=a(1)*bx*cx*(x-a(3))*(a(2)**2)/2.
      dyda(4)=a(1)*bx*(x-a(5))*(-SIN(a(4)*(x-a(5))))
      dyda(5)=a(1)*bx*a(4)*SIN(a(4)*(x-a(5)))
      RETURN
      END
      

C      SUBROUTINE funcs(x,a,y,dyda,na)
C      REAL x,y,a(na),dyda(na)
C      y=a(1)*EXP(a(2)*x)+a(3)+a(4)*SIN(x)
C      dyda(1)=EXP(a(2)*x)
C      dyda(2)=a(1)*EXP(a(2)*x)*x
C      dyda(3)=1.
C      dyda(4)=SIN(x)
C      RETURN
C      END

      SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,
     *     chisq,alamda)
      INTEGER ma,nca,ndata,ia(ma),MMAX
      REAL alamda,chisq,a(ma),alpha(nca,nca),covar(nca,nca),
     * sig(ndata),x(ndata),y(ndata)
      PARAMETER (MMAX=20) 
      INTEGER j,k,l,mfit
      REAL ochisq,atry(MMAX),beta(MMAX),da(MMAX)
      SAVE ochisq,atry,beta,da,mfit
     
      if(alamda.lt.0.)then 
         mfit=0
         do 11 j=1,ma
            if (ia(j).ne.0) mfit=mfit+1
 11      enddo 
         alamda=0.001
         call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq)

c	write(*,*) 'alpha=', alpha

         ochisq=chisq
         do 12 j=1,ma
            atry(j)=a(j)
 12      enddo 
      endif
      do 14 j=1,mfit 
         do 13 k=1,mfit
            covar(j,k)=alpha(j,k)
 13      enddo 
         covar(j,j)=alpha(j,j)*(1.+alamda)
         da(j)=beta(j)
 14   enddo 
c	write(*,*) 'covar=',covar
      call gaussj(covar,mfit,nca,da,1,1)
      if(alamda.eq.0.)then 
         call covsrt(covar,nca,ma,ia,mfit)
         call covsrt(alpha,nca,ma,ia,mfit) 
         return
      endif
      j=0
      do 15 l=1,ma 
         if(ia(l).ne.0) then
            j=j+1
            atry(l)=a(l)+da(j)
         endif
 15   enddo 
      call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq)
      if(chisq.lt.ochisq)then 
         alamda=0.1*alamda
         ochisq=chisq
         do 17 j=1,mfit
            do 16 k=1,mfit
               alpha(j,k)=covar(j,k)
 16         enddo 
            beta(j)=da(j)
 17      enddo 
         do 18 l=1,ma
            a(l)=atry(l)
 18      enddo 
      else 
         alamda=10.*alamda
         chisq=ochisq
      endif
      return
      END

      SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,
     * chisq)
      INTEGER ma,nalp,ndata,ia(ma),MMAX
      REAL chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata),
     * y(ndata)
      EXTERNAL funcs
      PARAMETER (MMAX=20)

      INTEGER mfit,i,j,k,l,m
      REAL dy,sig2i,wt,ymod,dyda(MMAX)
	
c	write(*,*) 'In sub 1, alpha=',alpha

      mfit=0
      do 11 j=1,ma
         if (ia(j).ne.0) mfit=mfit+1
 11   enddo 
      do 13 j=1,mfit 
         do 12 k=1,j
            alpha(j,k)=0.
 12      enddo 
         beta(j)=0.
 13   enddo 
      chisq=0.
      do 16 i=1,ndata 
         call funcs(x(i),a,ymod,dyda,ma)
         sig2i=1./(sig(i)*sig(i))
         dy=y(i)-ymod
         j=0
         do 15 l=1,ma
            if(ia(l).ne.0) then
               j=j+1
               wt=dyda(l)*sig2i
               k=0
               do 14 m=1,l
                  if(ia(m).ne.0) then
                     k=k+1
                     alpha(j,k)=alpha(j,k)+wt*dyda(m)
                  endif
 14            enddo 
               beta(j)=beta(j)+dy*wt
            endif
 15      enddo 
         chisq=chisq+dy*dy*sig2i
 16   enddo 
      do 18 j=2,mfit 
         do 17 k=1,j-1
            alpha(k,j)=alpha(j,k)
 17      enddo 
 18   enddo 
c	write(*,*) 'In sub 2, alpha=',alpha
      return
      END

      SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
      INTEGER ma,mfit,npc,ia(ma)
      REAL covar(npc,npc)
      INTEGER i,j,k
      REAL swap
      do 12 i=mfit+1,ma
         do 11 j=1,i
            covar(i,j)=0.
            covar(j,i)=0.
 11      enddo 
 12   enddo 
      k=mfit
      do 15 j=ma,1,-1
         if(ia(j).ne.0)then
            do 13 i=1,ma
               swap=covar(i,k)
               covar(i,k)=covar(i,j)
               covar(i,j)=swap
 13         enddo 
            do 14 i=1,ma
               swap=covar(k,i)
               covar(k,i)=covar(j,i)
               covar(j,i)=swap
 14         enddo 
            k=k-1
         endif
 15   enddo 
      return
      END


      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),
     * ipiv(NMAX) 
      REAL big,dum,pivinv
      do 11 j=1,n
         ipiv(j)=0
 11   enddo 
      do 22 i=1,n 
         big=0.
         do 13 j=1,n 
            if(ipiv(j).ne.1)then
               do 12 k=1,n
                  if (ipiv(k).eq.0) then
                     if (abs(a(j,k)).ge.big)then
                        big=abs(a(j,k))
                        irow=j
                        icol=k
                     endif
                  else if (ipiv(k).gt.1) then
		     return
c                     pause '1,singular matrix in gaussj'
                  endif
 12            enddo 
            endif
 13      enddo 
         ipiv(icol)=ipiv(icol)+1
         if (irow.ne.icol) then
            do 14 l=1,n
               dum=a(irow,l)
               a(irow,l)=a(icol,l)
               a(icol,l)=dum
 14         enddo 
            do 15 l=1,m
               dum=b(irow,l)
               b(irow,l)=b(icol,l)
               b(icol,l)=dum
 15         enddo 
         endif
         indxr(i)=irow 
         indxc(i)=icol
         if (a(icol,icol).eq.0.) return
c pause '2,singular matrix in gaussj'
         pivinv=1./a(icol,icol)
         a(icol,icol)=1.
         do 16 l=1,n
            a(icol,l)=a(icol,l)*pivinv
 16      enddo 
         do 17 l=1,m
            b(icol,l)=b(icol,l)*pivinv
 17      enddo 
         do 21 ll=1,n 
            if(ll.ne.icol)then 
               dum=a(ll,icol)
               a(ll,icol)=0.
               do 18 l=1,n
                  a(ll,l)=a(ll,l)-a(icol,l)*dum
 18            enddo 
               do 19 l=1,m
                  b(ll,l)=b(ll,l)-b(icol,l)*dum
 19            enddo 
            endif
 21     enddo 
 22   enddo 
      do 24 l=n,1,-1 
         if(indxr(l).ne.indxc(l))then
            do 23 k=1,n
               dum=a(k,indxr(l))
               a(k,indxr(l))=a(k,indxc(l))
               a(k,indxc(l))=dum
 23         enddo 
         endif
 24   enddo 
      return 
      END
