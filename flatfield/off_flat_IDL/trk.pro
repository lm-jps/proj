
;+
; NAME:		TRK
; PURPOSE:
;	This function, given an 2 images, will return the shifts in
;   	x and y direction that will align the second
;	image b with respect to a. 
; CATEGORY:
;	image processing
; CALLING SEQUENCE:
;	shift = trk(img,ref,mask=mask)
; INPUTS:
;
;	img	-	The image that you want to track with 
;			respect to a.
;			type: array,any type,arr(nx,ny)
;	ref	-	The reference image.
;			type: array,any type,arr(nx,ny)
;
; OUTPUTS:
;
;	shifts = vector containing  shifts
;		TYPE: float(2)
;
; KEYWORDS
;	mask = fourier filter is applied to image and reference
;	       currently available,
;	       mask = 1 -> highpass
;	       mask = 2 -> butterwoth highpass
;	
; SIDE EFFECTS:	None.
; RESTRICTIONS: None.
; PROCEDURE:
; MODIFICATION HISTORY:
;	Dec. 91 TR
;       Feb. 07 RW (avoids sticking to flat field)
;-
;

function trk,img,ref,mask=mask,sub=sub,sf=sf, ccf=ccf, flatfield=flatfield

  ;print,'track is running'

  ;Compute sizes of a and b and nx and ny
  sza=size(img) & szb=size(ref)
  nx=sza(1) & ny=sza(2)
  shifts = fltarr(2)
  ma = 1

  ;Array size : ERROR CHECK.
  if (sza(0) ne szb(0))or(sza(1) ne szb(1))or(sza(2) ne szb(2)) then begin 
    print,'Array sizes not equal' 
    return,0
  endif

	if keyword_set(sub) then begin
	  subz = sub
	  if keyword_set(mask) then begin
		case mask  of
		  1: ma = hipass(subz,subz)
		  2: ma = bwth(subz,subz)
		endcase
	  endif
	  lox = nx/2-subz/2
	  hix = lox+subz-1
	  loy = ny/2-subz/2
	  hiy = loy+subz-1

	endif else begin
	  lox=0
	  hix=nx-1
	  loy=0
	  hiy=ny-1
	  if keyword_set(mask) then begin
		case mask  of
		  1: ma = hipass(nx,ny)
		  2: ma = bwth(nx,ny)
		endcase
	  endif
	endelse
	  
	sref=ref(lox:hix,loy:hiy)	 
	simg=img(lox:hix,loy:hiy)	 
  if keyword_set(sf) then begin
	sref=sref-sfit(sref,1)
	simg=simg-sfit(simg,1)
  endif

  ;Compute fft's
  fref=fft(sref,-1) & fimg=fft(simg,-1)

  ;Compute the correlation between the two images (ccf)
  ccf=fft(fimg*conj(fref)*ma,1)


  ;Compute nx/2 and ny/2
  sz = size(sref)
  sx=sz(1)/2 & sy=sz(2)/2

  ;Convert ccf to float and shift.
  ccf=float(ccf) & ccf=shift(ccf,sx,sy)

  ;interpolate ccf at center point to avoid sticking to the flat field


  if keyword_set(flatfield) then begin
  c4x=(spline([findgen(7), 8+findgen(7)], [ccf[sx-7:sx-1, sy], ccf[sx+1:sx+7, sy]], findgen(15)))[7]
  c4y=(spline([findgen(7), 8+findgen(7)], [reform(ccf[sx, sy-7:sy-1]), reform(ccf[sx, sy+1:sy+7])], findgen(15)))[7]
  ccf[sx,sy]=0.5*(c4x+c4y)
  endif

  ;Find position of max of ccf.
  ps=max_pos(ccf)

  ;Determine shift and shift function.
  sx=ps(0)-sx
  sy=ps(1)-sy
  shifts(0)=-sx
  shifts(1)=-sy

  ;print,'Shifts in x and y',-sx,-sy



return, shifts

end
