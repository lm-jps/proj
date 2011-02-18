pro dark_rem_inorbit, imx, headers, time=time, imr, headim, filename=filename, xlow=xlow, xhigh=xhigh, ylow=ylow, yhigh=yhigh, drk=drk



ss=size(imx)

nbin=ss[1]

ns=ss[3]

if not keyword_set(xlow) then xlow=0
if not keyword_set(xhigh) then xhigh=nbin-1

if not keyword_set(ylow) then ylow=0
if not keyword_set(yhigh) then yhigh=nbin-1

nbinx=xhigh-xlow+1
nbiny=yhigh-ylow+1

; Identify dark images
exposure=getpar(headers, 'HSHIEXP')


; give exposure by hand


idx_image=where(exposure ne 0)
idx_dark=where(exposure eq 0, n_dark)





headim=headers[idx_image]


;;;;
n_image=n_elements(idx_image)


imr=fltarr(nbinx, nbiny, n_image)
dark=fltarr(nbinx,nbiny)
if keyword_set(drk) then dark=drk

dark_e=fltarr(nbin/16,nbin/16, n_image)



time=getshs(headers)
;;;







  for i=xlow, xhigh do begin
     
 
   for j=ylow, yhigh do begin

    if not keyword_set(drk) then dark[i-xlow,j-ylow]=median(imx[i,j,idx_dark])

  endfor

endfor





  for i=0, n_image-1 do begin

   imr[*,*,i]=imx[xlow:xhigh,ylow:yhigh,idx_image[i]] - dark

  endfor





end
