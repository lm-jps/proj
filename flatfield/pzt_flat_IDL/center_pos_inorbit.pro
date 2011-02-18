pro center_pos_inorbit, im, x0, y0, write=write




nx=(size(im))[1]
ny=(size(im))[2]
nl=(size(im))[3]

rad=nx/2

apod=apod_circ(nx, rad, 0)
bin=float(4096/nx)
imr=fltarr(nx,ny,nl)

for i=0, nl-1 do imr[*,*,i]=im[*,*,i]*apod


dx=imr[0:nx-3,1:nx-2,*]-imr[2:nx-1,1:nx-2,*] ;centred first difference
dy=imr[1:nx-2,0:nx-3,*]-imr[1:nx-2,2:nx-1,*]

dv=fltarr(nx,ny,nl)
dv[1:nx-2,1:ny-2,*]=sqrt(dx^2+dy^2)
xy=fltarr(nl,2)



for i=0, nl-1 do begin
print, i
xy[i,*]=-(trk(dv[*,*,i], rotate(dv[*,*,i],2))+1.0)/2.0
endfor

x0=xy[*,0]*bin
y0=xy[*,1]*bin

if keyword_set(write) then begin
openw,1,'/scr20/richard/hmi/legpos2'
printf, 1, x0
printf, 1, y0
close,1
endif


end
