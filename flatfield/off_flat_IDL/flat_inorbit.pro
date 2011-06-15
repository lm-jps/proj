;drk: 4096 x 4096 x 2



pro flat_inorbit, year, month, day, cam, focus,ftsid,flatfield,time_stamp=time_stamp,drk=drk, write=write
; arguments
;mandatory input:
;year, month, day: day of offpoint
;cam: camera ; 1: side, 2:front
;focus: focus position
;ftsid

;output:
;flatfield


;optional input
;time_stamp: time stamp for input, e.g. "2010.10.15_17:58:59_TAI"
;drk: array of dark image

;flag:
;write: if set, flatfield is written to hmi.offpoint_flatfield

print, systime()

if not keyword_set(drk) then begin
make_dark, year, month, day, cam, ftsid, time_stamp_dark, drk, write=write

endif


ftsstring=strcompress(string(ftsid), /rem)

;constants
minoffpos=10
std=1500.
n_focpos=10
flat_radius=1.08
camerastr=["side", "front"]
focstr=string(focus, format="(i2.2)")

spawn, 'show_info -q ds="hmi.lev1_nrt['+string(year, format='(i4.4)')+'.'+string(month, format='(i2.2)')+'.'+string(day, format='(i2.2)')+'/1d][?HFTSACID = '+ftsstring+'?][?CAMERA='+string(cam, format='(i1.1)')+'?]"', result
nf=long(strmid(result, 0, 4))

if nf eq 0 then begin
print, 'use lev 1 series'
spawn, 'show_info -q ds="hmi.lev1['+string(year, format='(i4.4)')+'.'+string(month, format='(i2.2)')+'.'+string(day, format='(i2.2)')+'/1d][?HFTSACID = '+ftsstring+'?][?CAMERA='+string(cam, format='(i1.1)')+'?]"', result
nf=long(strmid(result, 0, 4))
series_string="hmi.lev1"
endif else begin
series_string="hmi.lev1_nrt"
endelse





spawn, 'show_info -q ds="'+series_string+'['+string(year, format='(i4.4)')+'.'+string(month, format='(i2.2)')+'.'+string(day, format='(i2.2)')+'/1d][?HFTSACID = '+ftsstring+'?][?CAMERA='+string(cam, format='(i1.1)')+'?]" key=FSN,SAT_Y0,SAT_Z0,HCFTID,RSUN_LF,HCAMID > test.txt'

 
ll=dblarr(6, nf)
openr,1,'test.txt'
readf,1,ll
close,1

n_pos=nf/n_focpos
ll=reform(ll[*,0:n_pos*n_focpos-1], 6, n_focpos, n_pos)

fsn=reform(long(ll[0,*,*]))
rads=reform(ll[4,*,*])
sat_y0=reform(ll[1,*,*])
sat_z0=reform(ll[2,*,*])
foc=reform(ll[3,*,0])
hcam=reform(fix(ll[5,*,0]))

ar=sqrt(sat_y0^2+sat_z0^2)

selectnames_fsn, nm, 'hmi.lev0a', min(long(fsn)), max(long(fsn)), nfiles
nm=nm[fsn-fsn[0]]

if (nfiles-1) ne (max(long(fsn))-min(long(fsn))) then begin & print, 'missing image' & stop & endif



nm=reform(nm, n_focpos, n_pos)

print, n_focpos*n_pos, " images"

focind=where(foc eq focus and hcam gt 1)
nfoc=n_elements(focind)

ipp=ptrarr(nfoc, /all)
fsn_list=lonarr(1)
flatfield=fltarr(4096,4096,nfoc)

print, 'get reference flat'
spawn, 'show_info -q ds="hmi.offpoint_flatfield['+string(cam, format='(i1.1)')+'][2010.03.01-'+string(year, format='(i4.4)')+'.'+string(month, format='(i2.2)')+'.'+string(day, format='(i2.2)')+']" -p', result
ff=fitsio_read_image(result[n_elements(result)-1]+'/offpoint_flatfield.fits')




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

for kk=0, nfoc-1 do begin ; begin loop over identical focus positions // 


hhh=intarr(nf)
idd=where(abs(ar[focind[kk], *]-ar[focind[kk], 1:*]) lt 3.0)
hhh[idd]=1
wo, hhh, sek, sv=1
nsek=(size(sek))[2]
ipx=ptrarr(nsek, /all)
count=0
for i=0, nsek-1 do begin
if (sek[1,i]-sek[0,i] gt 2) then begin
*ipx[count]=sek[0,i]+1+indgen(sek[1,i]-sek[0,i]-1)
count=count+1
endif
endfor
ipx=ipx[0:count-1]



nl=n_elements(ipx)


stdy=fltarr(count)
stdz=fltarr(count)
meany=fltarr(count)
meanz=fltarr(count)

for i=0, count-1 do begin
stdy[i]=max(sat_y0[focind[kk], *ipx[i]])-min(sat_y0[focind[kk], *ipx[i]])
stdz[i]=max(sat_z0[focind[kk], *ipx[i]])-min(sat_z0[focind[kk], *ipx[i]])
meany[i]=avg(sat_y0[focind[kk], *ipx[i]])
meanz[i]=avg(sat_z0[focind[kk], *ipx[i]])
print, i, meany[i], meanz[i], stdy[i], stdz[i], n_elements(*ipx[i])
endfor

iks=where(stdy lt 2.0 and stdz lt 2.0)
if n_elements(iks) lt minoffpos then begin & print, "not enough good offpoint positions" & stop & endif
meany=meany[iks]
meanz=meanz[iks]
ipx=ipx[iks]
*ipp[kk]=ipx

nl=n_elements(ipx)





minim=min(sqrt(meany^2+meanz^2), mi)
if minim gt 20.0 then begin & print, "no good center position" & stop & endif
ind0=mi


print, "centered:", mi


;print, "number of focus sweeps"
;read, nl


;print, "focus"
;read, focus






radius=median(rads[focind[kk],*])
rad=radius*0.96



im=fltarr(4096,4096,nl)

off=fltarr(nl, 2)
imm=[0]
for i=0, nl-1 do imm=[imm, (*ipx[i])[1]]
imm=imm[1:*]
readimages, nm[focind[kk], imm], imx, head, nbin=1024, /noshow
apod=apod_circ(1024, 440, 72)

print, 'track images'
for i=0, nl-1 do off[i,*]=trk((imx[*,*,ind0]-mean(imx[*,*,ind0]))*apod, (imx[*,*,i]-mean(imx[*,*,i]))*apod)*4.0



print, 'determine center positions'
fstop_center_cm, imx, off[*,0], off[*,1], xx0, yy0, nl-1
fstop_center, imx, xx0, yy0, x0, y0


for k=0, nl-1 do begin 


index=*ipx[k]
ixx=index
nk=n_elements(index)
imr=fltarr(4096,4096,nk)

print, 'position', k, x0[k], y0[k]



readimages, nm[focind[kk],index], imx, head, /noshow
;readimages, nm[focind[kk],index], imx1, head
;readimages, nm[18+cam,index], imx2, head
;imr[*,*,indgen(n_elements(index))*2]=imx1
;imr[*,*,indgen(n_elements(index))*2+1]=imx2



for l=0, nk-1 do imr[*,*,l]=imx[*,*,l]-drk  ;dark subtraction



fsna=fsn[focind[kk],ixx]
fsna=fsna[sort(fsna)]


fsn_list=[fsn_list, fsna]


count=0L
for i=0, 4095 do for j=0, 4095 do begin & idx=where((imr[i,j,*]-median(imr[i,j,*])) gt std, cnt) & if sqrt((i-2048-x0[k])^2+(j-2048-y0[k])^2) lt rad and cnt gt 0 then begin & imr[i,j,idx]=!values.f_nan & count=count+1 & endif & endfor
print, "cosmic rays", count


dd=avg(imr, 2, /nan)
indexx=where(finite(dd) eq 0, count)
if count gt 0 then dd[indexx]=-1.0
im[*,*,k]=dd





endfor






;creating masks

;lower_limit=8000 ;front camera
;upper_limit=10000 ; front camera

lower_limit=0.85 ;side camera
upper_limit=1.15 ; side camera

offset=[[x0],[y0]]/4.0

print, 'mask active regions'

flatfield_mask, rebin(im,1024,1024,nl), offset, lower_limit,upper_limit, radius/4., masks, flat=rebin(ff, 1024,1024)

for i=0, 1023 do for j=0, 1023 do if n_elements(where(abs(masks[i,j,*]) eq 1 and sqrt((i-511.5)^2+(j-511.5)^2) lt 512.0*flat_radius)) eq nl then begin & print, "change mask limit" & stop & endif


for i=0, nl-1 do begin & idx=where(rebin(masks[*,*,i],4096,4096) eq 1, ct) & dd=im[*,*,i] & if ct gt 0 then dd[idx]=-1.0 & im[*,*,i]=dd  & endfor
;;;

openw,1,'legpos2'
printf, 1, x0
printf, 1, y0
close,1


openw,1,'imr_hmi0.bin'
writeu,1,im
close,1

print, 'calculate flatfield'
print, "./flatfield_iter_int 4096 "+string(nl, format='(i2.2)')+" "+string(cam, format='(i1.1)')+" 110 "+string(fix(radius/2048.0*100.0), format='(i2.2)')

spawn, "./flatfield_iter_int 4096 "+string(nl, format='(i2.2)')+" "+string(cam, format='(i1.1)')+" 110 "+string(fix(radius/2048.0*100.0), format='(i2.2)')

print, 'read flatfield'

ffield=fltarr(4096,4096)
openr,1,'flatfield_out0.bin'
readu,1,ffield
close,1
flatfield[*,*,kk]=ffield



endfor




fsn_list=fsn_list[1:*]
if nfoc gt 1 then begin
print, stddev(flatfield[*,*,0]-flatfield[*,*,1])
flatfield=avg(flatfield, 2)
endif

print, 'write out flatfield'
filename_flatfield='flat_'+camerastr[cam-1]+'_'+focstr+'_'+string(fsn_list[0],format='(i8.8)')+'.bin'
openw, 1, filename_flatfield
writeu,1,flatfield
close,1



if not keyword_set(time_stamp) then begin
spawn, 'show_info -q ds="'+series_string+'[]['+strcompress(string(min(fsn_list)),/rem)+']" key=T_OBS', result
time_stamp=result[0]
endif



ingest_string='write_offpoint instrument="HMI" file_offpoint="'+filename_flatfield+'" series_offpoint="hmi.offpoint_flatfield" camera='+string(cam, format='(i1.1)')+' pztflag=0 focus='+string(focus, format='(i2.2)')+' t_obs="'+time_stamp+'" fsn_list_offpoint='+strcompress(string(min(fsn_list)),/rem)+','+strcompress(string(max(fsn_list)),/rem)

print, ingest_string

openw,1,'ingest_command_'+camerastr[cam-1]+'_'+focstr+'_'+string(fsn_list[0],format='(i8.8)')+'.csh'
printf,1,ingest_string
close,1



print, min(fsn_list), max(fsn_list)



if keyword_set(write) then begin
print, "writing flatfield into database"
spawn, ingest_string

endif


print, "FLATFIELD COMPLETED"

end
