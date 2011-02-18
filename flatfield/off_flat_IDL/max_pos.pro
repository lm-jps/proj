
;	max_pos.pro
;
;
;	This function, given an array with a single well
;	defined maximum, finds the maximum position 
;	of that array to decimal accuracy using an interpolation
;	technique to find the maximum of an array from, Niblack, W., 
;	"An Introduction to Digital Image Processing", p 139.
;
;
;           ------------INPUT-------------			
;
;	f	-	An array.
;			type: array,any type
;
;           ------------OUTPUT------------			
;
;	max	-	An array containing the position of the maximum of an array.
;			type: vector,floating point,fltarr(2)
;
;	max(0)	-	The decimal position of the x-maximum.
;			type: scalar,floating point
;
;	max(1)	-	The decimal position of the y-maximum.
;			type: scalar,floating point

function max_pos,f,loc=loc

  ssz=size(f)
  fmax=max(f,loc)
  xmax=loc mod ssz(1)
  ymax=loc/ssz(1)

  ;A more complicated interpolation technique to find the maximum of an array.
  ;from, Niblack, W., "An Introduction to Digital Image Processing", p 139.
  	if (xmax*ymax gt 0) and (xmax lt (ssz(1)-1)) and (ymax lt (ssz(2)-1)) then begin
        denom = fmax*2 - f(xmax-1,ymax) - f(xmax+1,ymax)
        xfra = (xmax-.5) + (fmax-f(xmax-1,ymax))/denom
        denom = fmax*2 - f(xmax,ymax-1) - f(xmax,ymax+1)
        yfra = (ymax-.5) + (fmax-f(xmax,ymax-1))/denom
        xmax=xfra
        ymax=yfra
      endif

  rmax=[xmax,ymax]

return,rmax

end
