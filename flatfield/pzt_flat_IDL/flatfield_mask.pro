pro flatfield_mask, imr, offset, level_spot, level_plage, radius,  masks, flat=flat

nx=(size(imr))[1]
ny=(size(imr))[2]
nl=(size(imr))[3]

lb=fltarr(nx,ny,nl)
apod=lb
masks=intarr(nx,ny,nl)

sigma=9.0*(nx/1024)
bin=nx/1024

x=rebin(findgen(nx),nx,ny)
y=transpose(rebin(findgen(ny),ny,nx))
;a=[0.340000,      1.37000,     -2.04000,      2.70000,     -1.94000,     0.559000]
;a=[-1.36813,  30.9044, -232.513, 999.353, -2613.42, 4211.68,
;-3921.77, 1512.53, 591.547, -802.067, 226.122]
a=[0.229947, 1.99002,-4.01991, 5.96688, -4.53731, 1.37037]

for k=0, nl-1 do begin
lb[*,*,k]=limb_darkening(nx, radius, offset[k,*], a)  ; limb darkening function: current limb darkening function inaccurate at the limb\
apod[*,*,k]=apod_circ(nx, radius*0.99, 0, offx=offset[k,0], offy=offset[k,1])
endfor
apodc=apod_circ(nx, 1.09*nx/2.0, 0)
;if not keyword_set(flat) then begin
;c1=readfits('/home/soi/CM/tables/calib/flat/hr/vers_0/1.fits') ;get approximate (large scale) flatfield -- this one is definitely not erfect
;flat=1./c1
;endif



if not keyword_set(flat) then flat=1.0+fltarr(nx,ny)

for k=0L, nl-1 do begin
imsm=smooth(imr[*,*,k]/flat/lb[*,*,k],bin*2)



idx=where(apod[*,*,k] eq 1.0)

print, k, avg(imsm[idx])
threshold=avg(imsm[idx])*level_spot
threshold_plage=avg(imsm[idx])*level_plage

ixx=where((imsm lt threshold or imsm gt threshold_plage) and apod[*,*,k] eq 1.0 and apodc eq 1.0, count)  ; masks out dark regions and bight regions (plages).
mm=replicate(1.0, nx, ny)
nn=n_elements(ixx)
print, nn


if count gt 0 then begin
for i=0L, nn-1 do begin 
idd_x=indgen(21)-10+x[ixx[i]]
idd_y=indgen(21)-10+y[ixx[i]]
idd_x=idd_x[where(idd_x ge 0 and idd_x lt nx)]
idd_y=idd_y[where(idd_y ge 0 and idd_y lt ny)]

ibb=transpose(rebin(idd_y,n_elements(idd_y),n_elements(idd_x)))*nx+rebin(idd_x,n_elements(idd_x),n_elements(idd_y))


mm[ibb]=mm[ibb]*(1.0-exp(-((x[ibb]-x[ixx[i]])^2+(y[ibb]-y[ixx[i]])^2)/2.0/sigma))



endfor
endif


; obtain neigborhood for missing pixels
idx=where(mm lt 0.2, ct)
msk=intarr(nx,ny)
;i0x=where(apod[*,*,k] ne 1.0 or apodc ne 1.0)
;msk[i0x]=-1
if ct gt 0 then msk[idx]=1
masks[*,*,k]=msk
print, 'mask number', k

endfor



end
