pro fstop_center, imr, xx, yy, xxn, yyn

nk=4
kerl=fltarr(2*nk+1,2*nk+1)-1
kerl(nk,nk)=(2*nk+1)^2-1

imk=imr



nleg=(size(imr))[3]
nx=(size(imr))[1]

samp=4096/nx


fac=502./1024. ; changed from 540 !!

for i=0, nleg-1 do  imk[*,*,i]=convol(imr[*,*,i], kerl)

dist=shift(dist(nx), nx/2, nx/2)
idx=where(dist gt fac*nx)


im=fltarr(2*nx,2*nx,nleg)

xy=fltarr(nleg, 2)

for i=0, nleg-1 do begin


dd=imk[*,*,i]
dd[idx]=0.0


im[nx/2:3*nx/2-1, nx/2:3*nx/2-1,i]=dd
im[*,*,i]=shift(im[*,*,i], -fix(xx[i]/samp), -fix(yy[i]/samp))

endfor

;tvim, avg(im,2)

xy0=(trk(avg(im,2), rotate(avg(im,2),2))+1.)/2.



for i=0, nleg-1 do begin

xy[i,*]=(trk(avg(im,2), im[*,*,i])-xy0)*samp + [fix(xx[i]), fix(yy[i])]


;spcoord_qs, im[nx/2:3*nx/2-1, nx/2:3*nx/2-1,i], rphi, xc=nx/2+xy[i,0]-xy0[0], yc=nx/2+xy[i,0]-xy0[1], frad=1.2 
; tvim, rphi
;oplot, [460,460], [0, 512]



endfor


xxn=xy[*,0]
yyn=xy[*,1]

end

