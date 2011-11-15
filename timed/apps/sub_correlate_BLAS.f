CC This is a new subroutine to compute correlation by calling the BLAS
CC library to compute the dot-product. However, the computation of 
CC averages seems stupid. But no BLAS FUNCTION does that.  ... Mar 2, 2005
      SUBROUTINE c_correlate(x,y,n,lag,m,ans)
CC    x and y are data series with dimension of n, and
CC    lag's dimension is m. ans(m) is the results
      DIMENSION x(n), y(n), ans(m),lag(m)

CC    first is to get the average of both x and y
      total_x=0.
      total_y=0.
      DO 10 i=1, n
         total_x=total_x+x(i)
         total_y=total_y+y(i)
 10   CONTINUE
      ave_x=total_x/n
      ave_y=total_y/n
      DO 20 i=1,n
         x(i)=x(i)-ave_x
         y(i)=y(i)-ave_y
 20   CONTINUE

CC    The following is to begin the calculation of
CC    cross-covariance function
      DO 100 l=1, m
         ilag=lag(l)
         IF(ilag.LT.0) THEN
            tot=0.
            tot=SDOT(n+ilag,x(1-ilag),1,y,1)
            ans(l)=tot/n
         ELSE
            tot=0.
            tot=SDOT(n-ilag,x,1,y(1+ilag),1)
            ans(l)=tot/n
         ENDIF
 100  ENDDO

      RETURN
      END

