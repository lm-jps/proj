pro pzt_flatfield, year, month, day, flatfield_front, flatfield_side, ff_front, ff_side, x0, y0, offpoint_front, offpoint_side, msk_front, msk_side, nowrite=nowrite, plots=plots


ftsstring=['3001', '3021']

print, systime()


;jul=julday(month, day, year)
;caldat, jul-1, month1, day1, year1

;daystring1=string(year1, format='(i4.4)')+'.'+string(month1, format='(i2.2)')+'.'+string(day1, format='(i2.2)')
daystring=string(year, format='(i4.4)')+'.'+string(month, format='(i2.2)')+'.'+string(day, format='(i2.2)')

spawn, "show_info ds='hmi.offpoint_flatfield[][1977.01.01_TAI-"+daystring+"][?PZTFLAG=1?]' key=T_OBS; echo $status", result

if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure1" & return & endif
t_obs_last_update=strmid(result[n_elements(result)-2], 0, 10)

k=0
repeat begin
spawn, "show_info -q ds='hmi.lev0a[][?T_OBS >= $("+t_obs_last_update+"_23:59:59.99_TAI) and T_OBS <= $("+daystring+"_23:59:59.99_TAI)?][?HFTSACID="+ftsstring[k]+"?][?HCAMID=1?]' key=T_OBS; echo $status", result_time
if result_time[n_elements(result_time)-1] ne '0' then begin & print, "access to database failure2" & return & endif


if n_elements(result_time) lt 3 then begin
print, "not enough data with FTSID ", ftsstring[k]
endif

spawn, "show_info -q ds='hmi.lev0a[][?T_OBS >= $("+t_obs_last_update+"_23:59:59.99_TAI) and T_OBS <= $("+daystring+"_23:59:59.99_TAI)?][?HFTSACID="+ftsstring[k]+"?][?HCAMID=1?]' key=FSN; echo $status", result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure3" & return & endif

spawn, "show_info -q ds='hmi.lev0a[][?T_OBS >= $("+t_obs_last_update+"_23:59:59.99_TAI) and T_OBS <= $("+daystring+"_23:59:59.99_TAI)?][?HFTSACID="+ftsstring[k]+"?][?HCAMID=1?]' key=HCFTID; echo $status", result_focus
if result_focus[n_elements(result_focus)-1] ne '0' then begin & print, "access to database failure" & return & endif

k=k+1
endrep until (n_elements(result_time) ge 3 or k eq n_elements(ftsstring)) 


if (n_elements(result_time) lt 3) then begin & print, 'not enough data' & return & endif

focus=long(result_focus[0])
focstr='['+strcompress(string(focus),/remove_all)+']'
print, 'focus: ', focstr

minute=fltarr(n_elements(result_time)-1)
for i=0, n_elements(minute)-1 do begin
minute[i]=date_obs2min(result_time[i])
print, minute[i]-minute[0]
endfor

fsn=[0L,0L]
for i=0, n_elements(minute)-1 do begin
for j=i+1, n_elements(minute)-1 do begin
if minute[j]-minute[i] lt 36.0*60.0 then begin
print, (minute[j]-minute[i]), 'minutes'
fsna=[long(result[i]), long(result[j])]
if (fsna[1]-fsna[0]) gt (fsn[1]-fsn[0]) then fsn=fsna
endif
endfor
endfor





print, fsn[0], fsn[1]




if fsn[0] eq 0 then begin
print, "not enough data with correct FTSID"
return
endif


selectnames_fsn,nm1,'hmi.lev0a', fsn[0], fsn[0]+23
selectnames_fsn,nm2,'hmi.lev0a', fsn[1], fsn[1]+23

nm=[[nm1],[nm2]]




size=4096
nl=(size(nm))[2]
normrad=1970.0
dst=shift(dist(4096),2048,2048)
flatfield_front=fltarr(size,size,nl)
flatfield_side=fltarr(size,size,nl)

keychange_front=strarr(2)
keychange_side=strarr(2)


msk_front=fltarr(4096,4096,nl)
msk_side=fltarr(4096,4096,nl)

k=nl-1
;front cam
spawn, 'show_info ds="hmi.lev1_nrt[]['+string(fsn[k]+1,format="(i8.8)")+'-'+string(fsn[k]+22,format="(i8.8)")+'][?CAMERA=2?]" key=RSUN_LF,T_OBS,HCFTID; echo $status', resultp
if resultp[n_elements(resultp)-1] ne '0' then begin & print, "access to database failure4" & return & endif

if n_elements(resultp)-1 ne 12 then begin & print, "lev1_nrt missing" & return & endif



solrad_front=fix(float(strmid(resultp[n_elements(resultp)-2],0, 10))/2048.*100.) 
radius=float(strmid(resultp[n_elements(resultp)-2],0, 10))

t_start_front=strmid(strmid(resultp[n_elements(resultp)-2], 12, 26),0,19)+'_UTC'
t_start_front_search=strmid(strmid(resultp[n_elements(resultp)-3], 12, 26),0,19)+'_UTC'
t_0_front='1976.12.31_23:59:45_UTC'
focstr_front=strmid(resultp[n_elements(resultp)-2], 39,2)


spawn, 'show_info ds="hmi.dark[2]['+t_0_front+'-'+t_start_front_search+']" key=T_OBS; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure5" & return & endif

t_obs_dark_front=result[n_elements(result)-2]


print, 'show_info ds="hmi.dark[2]['+t_0_front+'-'+t_start_front_search+']" -P'
spawn, 'show_info ds="hmi.dark[2]['+t_0_front+'-'+t_start_front_search+']" -P; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure6" & return & endif

dark_path_front=result[n_elements(result)-2]
dark_front=fitsio_read_image(dark_path_front+'/dark.fits')

spawn, 'show_info ds="hmi.bad_pixel_list[2]['+t_0_front+'-'+t_start_front_search+']" key=T_OBS; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure7" & return & endif
t_obs_badpix_front=result[n_elements(result)-2]

spawn, 'show_info ds="hmi.offpoint_flatfield[2]['+t_0_front+'-'+t_start_front_search+']'+focstr+'[?PZTFLAG=1?]" key=T_OBS; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure" & return & endif
t_obs_offpoint_front=result[n_elements(result)-2]

spawn, 'show_info ds="hmi.offpoint_flatfield[2]['+t_0_front+'-'+t_start_front_search+']'+focstr+'[?PZTFLAG=1?]" -q -P; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure8" & return & endif
last_flat_front=fitsio_read_image(result[n_elements(result)-2]+'/offpoint_flatfield.fits')

spawn, 'show_info ds="hmi.flatfield[2]['+t_0_front+'-'+t_start_front_search+']" key=T_START; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure9" & return & endif
t_start_last_front=result[n_elements(result)-2]

spawn, 'show_info ds="hmi.flatfield[2]['+t_0_front+'-'+t_start_front+']" key=T_STOP; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure10" & return & endif
t_stop_last_front=result[n_elements(result)-2]


spawn, 'show_info ds="hmi.flatfield[2]['+t_0_front+'-'+t_start_front_search+']" key=FLATFIELD_VERSION; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure" & return & endif
version_last_front=result[n_elements(result)-2]


keychange_front[0]='set_keys ds="hmi.flatfield[2]['+t_start_last_front+']" FLATFIELD_VERSION=1'
keychange_front[1]='set_keys ds="hmi.flatfield[2]['+t_start_last_front+']" T_STOP="'+t_start_front+'"'

print, 'show_info ds="hmi.offpoint_flatfield[2]['+t_0_front+'-'+t_start_front_search+']'+focstr+'[?PZTFLAG=0?]" -P'
spawn, 'show_info ds="hmi.offpoint_flatfield[2]['+t_0_front+'-'+t_start_front_search+']'+focstr+'[?PZTFLAG=0?]" -P; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure" & return & endif

offpoint_path_front=result[n_elements(result)-2]
offpoint_front=fitsio_read_image(offpoint_path_front+'/offpoint_flatfield.fits')

spawn, 'show_info ds="hmi.offpoint_flatfield[2]['+t_0_front+'-'+t_start_front_search+']'+focstr+'[?PZTFLAG=0?]" key="FSN_INPUT"; echo $status', result_fsn_input
if result_fsn_input[n_elements(result_fsn_input)-1] ne '0' then begin & print, "access to database failure" & return & endif
string_offpoint_front=strcompress(result_fsn_input[n_elements(result_fsn_input)-2], /remove_all)
apos=strpos(string_offpoint_front, 'PZT_FSN')
if apos ne -1 then string_offpoint_front=strmid(string_offpoint_front, 0, apos-1)

print, string_offpoint_front

print, 'solrad front' ,solrad_front




; side cam

spawn, 'show_info ds="hmi.lev1_nrt[]['+string(fsn[k]+1,format="(i8.8)")+'-'+string(fsn[k]+22,format="(i8.8)")+'][?CAMERA=1?]" key=RSUN_LF,T_OBS,HCFTID; echo $status', resultp
if resultp[n_elements(resultp)-1] ne '0' then begin & print, "access to database failure" & return & endif
if n_elements(resultp)-1 ne 12 then begin & print, "lev1_nrt missing" & return & endif

solrad_side=fix(float(strmid(resultp[n_elements(resultp)-2],0, 10))/2048.*100.)

t_start_side=strmid(strmid(resultp[n_elements(resultp)-2], 12, 26),0,19)+'_UTC'
t_start_side_search=strmid(strmid(resultp[n_elements(resultp)-3], 12, 26),0,19)+'_UTC'
t_0_side='1976.12.31_23:59:45_UTC'
focstr_side=strmid(resultp[n_elements(resultp)-2], 39,2)

spawn, 'show_info ds="hmi.dark[1]['+t_0_side+'-'+t_start_side_search+']" key=T_OBS; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure" & return & endif
t_obs_dark_side=result[n_elements(result)-2]

print, 'show_info ds="hmi.dark[1]['+t_0_side+'-'+t_start_side_search+']" -P'
spawn, 'show_info ds="hmi.dark[1]['+t_0_side+'-'+t_start_side_search+']" -P; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure" & return & endif
dark_path_side=result[n_elements(result)-2]
dark_side=fitsio_read_image(dark_path_side+'/dark.fits')


spawn, 'show_info ds="hmi.bad_pixel_list[1]['+t_0_side+'-'+t_start_side_search+']" key=T_OBS; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure" & return & endif
t_obs_badpix_side=result[n_elements(result)-2]

spawn, 'show_info ds="hmi.offpoint_flatfield[1]['+t_0_side+'-'+t_start_side_search+']'+focstr+'[?PZTFLAG=1?]" key=T_OBS; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure" & return & endif
t_obs_offpoint_side=result[n_elements(result)-2]

spawn, 'show_info ds="hmi.offpoint_flatfield[1]['+t_0_side+'-'+t_start_side_search+']'+focstr+'[?PZTFLAG=1?]" -q -P; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure" & return & endif
last_flat_side=fitsio_read_image(result[n_elements(result)-2]+'/offpoint_flatfield.fits')

spawn, 'show_info ds="hmi.flatfield[1]['+t_0_side+'-'+t_start_side_search+']" key=T_START; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure" & return & endif
t_start_last_side=result[n_elements(result)-2]


spawn, 'show_info ds="hmi.flatfield[1]['+t_0_side+'-'+t_start_side+']" key=T_STOP; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure" & return & endif
t_stop_last_side=result[n_elements(result)-2]

spawn, 'show_info ds="hmi.flatfield[1]['+t_0_side+'-'+t_start_side_search+']" key=FLATFIELD_VERSION; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure" & return & endif
version_last_side=result[n_elements(result)-2]


keychange_side[0]='set_keys ds="hmi.flatfield[1]['+t_start_last_side+']" FLATFIELD_VERSION=1'
keychange_side[1]='set_keys ds="hmi.flatfield[1]['+t_start_last_side+']" T_STOP="'+t_start_side+'"'

print, 'show_info ds="hmi.offpoint_flatfield[1]['+t_0_side+'-'+t_start_side_search+']'+focstr+'[?PZTFLAG=0?]" -P'
spawn, 'show_info ds="hmi.offpoint_flatfield[1]['+t_0_side+'-'+t_start_side_search+']'+focstr+'[?PZTFLAG=0?]" -P; echo $status', result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure" & return & endif
offpoint_path_side=result[n_elements(result)-2]
offpoint_side=fitsio_read_image(offpoint_path_side+'/offpoint_flatfield.fits')

spawn, 'show_info ds="hmi.offpoint_flatfield[1]['+t_0_side+'-'+t_start_side_search+']'+focstr+'[?PZTFLAG=0?]" key="FSN_INPUT"; echo $status', result_fsn_input
if result_fsn_input[n_elements(result_fsn_input)-1] ne '0' then begin & print, "access to database failure" & return & endif
string_offpoint_side=strcompress(result_fsn_input[n_elements(result_fsn_input)-2], /remove_all)
apos=strpos(string_offpoint_side, 'PZT_FSN')
if apos ne -1 then string_offpoint_side=strmid(string_offpoint_side, 0, apos-1)

print, string_offpoint_side

print, 'solrad side', solrad_side


print, t_start_front
print, t_stop_last_front







nim=12
strat_front=''
strat_side=''

; loop over PZT flatfields
for k=0, nl-1 do begin

readimages, nm[indgen(nim)*2, k], imx, head, nbin=size, /noshow
dark_rem_inorbit, imx, head, imr, headim, drk=dark_front

fsn_front=getpar(headim, 'FSN')

for i=0, n_elements(fsn_front)-1 do strat_front=strat_front+string(fsn_front[i], format='(i8.7)')+","
strat_front=strcompress(strat_front, /remove_all)





for i=0, nim-2 do imr[*,*,i]=imr[*,*,i]/offpoint_front

center_pos_inorbit, rebin(imr,1024,1024,nim-1), x0, y0

x00=round(x0[0])
y00=round(y0[0])



flatfield_mask, rebin(imr, 1024,1024,nim-1), [[x0],[y0]]/4.0, 0.9, 1.1, radius/4., masks
dd=total(masks,3)
dm=bytarr(1024,1024)
idx=where(dd gt 0)
if idx[0] ne -1 then dm[idx]=1
dm=rebin(dm, 4096,4096)


msk_front[*,*,k]=dm


std=1500
;cosmic ray finder
cnt=0L
for i=0, 4095 do for j=0, 4095 do begin & idx=where(abs(imr[i,j,*]-median(imr[i,j,*])) gt std, count) & if sqrt((i-2048-x00)^2+(j-2048-y00)^2) lt float(solrad_front)/100.*2048.0*0.99 and count gt 0 then begin & cnt=cnt+1 & imr[i,j,idx]=-1.0 & endif & endfor


print, "number of cosmic rays", cnt

if keyword_set(plots) then plot, x0, y0, psym=2






openw,1,'legpos2'
printf, 1, x0
printf, 1, y0
close,1

if keyword_set(plots) then  plot, x0, y0,psym=-2
openw,1,'imr_hmi0.bin'
writeu,1,imr
close,1


 spawn, "./flatfield_iter_next 4096 "+string(nim-1,format="(i2.2)")+" 2 "+string(solrad_front+1,format="(i2.2)")+" "+string(solrad_front,format="(i2.2)")+"; echo $status", result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure" & return & endif


ff=fltarr(size,size)
openr,1,'flatfield_out0.bin'
readu,1,ff
close,1

apod=apod_circ(4096, 2048.0*solrad_front/100*0.99, 2048.*solrad_front/100*0.02)
ff=(ff-1.0)*apod+1.0

;idx=where(dst gt solrad/101.*2048.)
;ff[idx]=1.0

;idx=where(dm ne 0)
;if idx[0] ne -1 then ff[idx]=1.0




flatfield_front[*,*,k]=ff







readimages, nm[indgen(nim)*2+1,k], imx, head, nbin=size, /noshow
dark_rem_inorbit, imx, head, imr, headim, drk=dark_side

fsn_side=getpar(headim, 'FSN')
for i=0, n_elements(fsn_side)-1 do strat_side=strat_side+string(fsn_side[i], format='(i8.7)')+","
strat_side=strcompress(strat_side, /remove_all)


for i=0, nim-2 do imr[*,*,i]=imr[*,*,i]/offpoint_side

center_pos_inorbit, rebin(imr,1024,1024,nim-1), x0, y0
x00=round(x0[0])
y00=round(y0[0])

flatfield_mask, rebin(imr, 1024,1024,nim-1), [[x0],[y0]]/4.0, 0.9, 1.1, radius/4., masks
dd=total(masks,3)
dm=bytarr(1024,1024)
idx=where(dd ne 0)
if idx[0] ne -1 then dm[idx]=1
dm=rebin(dm, 4096,4096)

msk_side[*,*,k]=dm

if keyword_set(plots) then plot, x0, y0, psym=2

cnt=0
for i=0, 4095 do for j=0, 4095 do begin & idx=where(abs(imr[i,j,*]-median(imr[i,j,*])) gt std, count) & if sqrt((i-2048-x00)^2+(j-2048-y00)^2) lt float(solrad_side)/100.*0.99*2048.0 and count gt 0 then begin & cnt=cnt+1 & imr[i,j,idx]=-1.0 & endif & endfor
print, "number of cosmic rays", cnt

if keyword_set(plots) then plot, x0, y0, psym=-2


openw,1,'legpos2'
printf, 1, x0
printf, 1, y0
close,1


openw,1,'imr_hmi0.bin'
writeu,1,imr
close,1



 spawn, "./flatfield_iter_next 4096 "+string(nim-1,format="(i2.2)")+" 1 "+string(solrad_side+1,format="(i2.2)")+" "+string(solrad_side,format="(i2.2)")+"; echo $status", result
if result[n_elements(result)-1] ne '0' then begin & print, "access to database failure" & return & endif

ff=fltarr(size,size)
openr,1,'flatfield_out0.bin'
readu,1,ff
close,1

apod=apod_circ(4096, 2048.0*solrad_side/100.0*.99, 2048.*solrad_side/100*0.02)
ff=(ff-1.0)*apod+1.0

;idx=where(dst gt solrad_side/101.*2048.)
;ff[idx]=1.0

;idx=where(dm ne 0)
;if idx[0] ne -1 then ff[idx]=1.0


flatfield_side[*,*,k]=ff

endfor                 ; end loop over PZT flatfields



for j=0, 4095 do begin
for i=0, 4095 do begin
idm=where(msk_front[i,j,*] eq 0)
if idm[0] ne -1 then ff[i,j]=avg(flatfield_front[i,j,idm],2) else ff[i,j]=1.0
endfor
endfor

ff_front=ff*offpoint_front
ff_front[where(offpoint_front ne 1.0)]=ff_front[where(offpoint_front ne 1.0)]/avg(ff_front[where(dst lt normrad)])



for j=0, 4095 do begin
for i=0, 4095 do begin
idm=where(msk_side[i,j,*] eq 0)
if idm[0] ne -1 then ff[i,j]=avg(flatfield_side[i,j,idm],2) else ff[i,j]=1.0
endfor
endfor





ff_side=ff*offpoint_side
ff_side[where(offpoint_side ne 1.0)]=ff_side[where(offpoint_side ne 1.0)]/avg(ff_side[where(dst le normrad)])







string_front='write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="'+t_start_front+'" camera=2 focus='+focstr_front+' t_obs_badpix="'+t_obs_badpix_front+'" t_obs_dark="'+t_obs_dark_front+'" t_obs_offpoint="'+t_start_front+'"  file_flatfield="pzt_flat_front_'+string(fsn[nl-1],format='(i8.8)')+'.bin"'

if version_last_front eq 1 then string_front=string_front+' flatfield_version=1 t_stop="'+t_stop_last_front+'"'


string_side='write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="'+t_start_side+'" camera=1 focus='+focstr_side+' t_obs_badpix="'+t_obs_badpix_side+'" t_obs_dark="'+t_obs_dark_side+'" t_obs_offpoint="'+t_start_side+'"  file_flatfield="pzt_flat_side_'+string(fsn[nl-1],format='(i8.8)')+'.bin"'

if version_last_side eq 1 then string_side=string_side+' flatfield_version=1 t_stop="'+t_stop_last_side+'"'

string_off_front='write_offpoint instrument="HMI" file_offpoint="pzt_flat_front_'+string(fsn[nl-1],format='(i8.8)')+'.bin"' + ' series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="'+t_start_front+'" fsn_list_offpoint="'+string_offpoint_front+'" fsn_list_pzt="'+strmid(strat_front,0,strlen(strat_front)-1)+'" nx=4096 ny=4096 focus='+focstr_front


string_off_side='write_offpoint instrument="HMI" file_offpoint="pzt_flat_side_'+string(fsn[nl-1],format='(i8.8)')+'.bin"' + ' series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="'+t_start_side+'" fsn_list_offpoint="'+string_offpoint_side+'" fsn_list_pzt="'+strmid(strat_side,0,strlen(strat_side)-1)+'" nx=4096 ny=4096 focus='+focstr_side


err_front=stddev(flatfield_front[2048-256:2048+255,2048-256:2048+255,0]-flatfield_front[2048-256:2048+255,2048-256:2048+255,1])/2.0
print, 'error front:', err_front

min_front=min(flatfield_front)
std_front=min([stddev(flatfield_front[*,*,0]),stddev(flatfield_front[*,*,1])])

print, "min front std front", min_front, std_front

err_side=stddev(flatfield_side[2048-256:2048+255,2048-256:2048+255,0]-flatfield_side[2048-256:2048+255,2048-256:2048+255,1])/2.0
print, 'error side:', err_side

min_side=min(flatfield_side)
std_side=min([stddev(flatfield_side[*,*,0]),stddev(flatfield_side[*,*,1])])

print, "min side, std side", min_side, std_side

openw,1,'pzt_flat_side_'+string(fsn[nl-1],format='(i8.8)')+'.bin'
writeu,1, ff_side
close,1

openw,1,'pzt_flat_front_'+string(fsn[nl-1],format='(i8.8)')+'.bin'
writeu,1, ff_front
close,1

!p.multi=[0,2,1]
if keyword_set(plots) then  tvim, rebin(ff_front/offpoint_front,256,256)
if keyword_set(plots) then tvim, rebin(ff_side/offpoint_side,256,256)

script_filename='ingest_commands_'+string(fsn[nl-1],format='(i8.8)')+'.csh'

print, "output file"
print, script_filename

openw,1,script_filename
printf,1,'#!/bin/csh'
printf,1,''

;if version_last_side eq 0 then begin
;printf,1,keychange_side[0]
;printf,1,keychange_side[1]
;endif

;if version_last_front eq 0 then begin
;printf,1,keychange_front[0]
;printf,1,keychange_front[1]
;endif


printf,1,''
printf,1,string_off_side
printf,1,''
printf,1,string_off_front
;printf,1,''
;printf,1,string_side
;printf,1,''
;printf,1,string_front
close,1

if not keyword_set(nowrite) then begin
if err_front lt 0.002 and err_side lt 0.002 and min_front gt 0.0 and min_side gt 0.0 and std_front lt 0.05 and std_side lt 0.05 then begin
spawn, 'chmod 744 '+script_filename
spawn, './'+script_filename
endif else begin
print, "errors too large: Please check!"
endelse
endif


end

