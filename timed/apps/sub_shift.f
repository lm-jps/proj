      SUBROUTINE shift(x,n,delta)
CC     x(n) is the original data array for shift.
CC     delta is the value to be shifted, and the 
CC     routine returns x(n) after shifting.
      DIMENSION x(n), y1((n+1)/2), y2((n+1)/2)
      DIMENSION y1_temp((n+1)/2),y2_temp((n+1)/2)
      REAL delta
CC     first is to divide x(n) into 2 sub-arrays
      DO 5 i=1, (n+1)/2
         y1(i)=x((n+1)/2+i-1)
         y2(i)=x((n+1)/2-i+1)
 5    CONTINUE
      
      IF(delta.GE.0.) THEN
 	 m=INT(delta)
	 del=delta-m
         DO 15 i=m+2, (n+1)/2
            y1_temp(i)=y1(i-m-1)+(y1(i-m)-y1(i-m-1))*(1-del)
            y2_temp(i)=y2(i-m-1)+(y2(i-m)-y2(i-m-1))*(1-del)
 15      ENDDO
	 DO 16 i=1, m+1
           y1_temp(i)=y1(1)-(y1(2)-y1(1))*(i+del-1)
           y2_temp(i)=y2(1)-(y2(2)-y2(1))*(i+del-1)
 16 	 ENDDO
         
      ELSE
	 m=INT(-delta)
	 del=-delta-m
         DO 20 i=1, (n+1)/2-(m+1)
            y1_temp(i)=y1(i+m)+(y1(i+m+1)-y1(i+m))*(del)
            y2_temp(i)=y2(i+m)+(y2(i+m+1)-y2(i+m))*(del)
 20      ENDDO
      ENDIF

      DO 25 i=1,(n+1)/2
         x(i)=y2_temp((n+1)/2-i+1)
         x(i+(n+1)/2-1)=y1_temp(i)
 25   ENDDO
      x((n+1)/2)=y1(1)
      
      RETURN
      END
