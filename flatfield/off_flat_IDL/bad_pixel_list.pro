pro bad_pixel_list, year, month, day, cam, focus, write=write


;**************************************
;List of pixels with unstable gain / hard coded - they can not be
;                                                 detected from the
;                                                 flat field

hard_bad_list=ptrarr(2, /allocate)
;*hard_bad_list[0]=[ ... ] ; side camera pixel with unstable gain ;for now  - none
*hard_bad_list[1]=[9856346, 9860442, 8389978] ; front camera pixel with unstable gain

;************************************************************


camstr=["side","front"]





spawn, 'show_info -q ds="hmi.offpoint_flatfield['+string(cam,format='(i1.1)')+']['+string(year, format='(i4.4)')+'.'+string(month, format='(i2.2)')+'.'+string(day, format='(i2.2)')+'/1d]['+string(focus,format='(i2.2)')+']" -P', result
if n_elements(result) ne 1 then begin & print, "no proper flatfield update for this day" & return & endif



ff=fitsio_read_image(result[0]+'/offpoint_flatfield.fits')


spawn, 'show_info -q ds="hmi.offpoint_flatfield['+string(cam,format='(i1.1)')+']['+string(year, format='(i4.4)')+'.'+string(month, format='(i2.2)')+'.'+string(day, format='(i2.2)')+'/1d]['+string(focus,format='(i2.2)')+']" key=T_OBS', result
time_stamp=result[0]


; create bad pixel list
dst=shift(dist(4096),2048,2048)

bad_pix=where(ff le 0.5 and dst le 2200.0)



if n_elements(*hard_bad_list[cam-1]) gt 0 then begin
add_bad=*hard_bad_list[cam-1]
bad_pix=[bad_pix, add_bad]
endif

nbad=n_elements(bad_pix)


filename='badpix_'+camstr[cam-1]+'_'+string(year, format='(i4.4)')+'.'+string(month, format='(i2.2)')+'.'+string(day, format='(i2.2)')+'.bin'
openw,1,filename
writeu,1,bad_pix
close,1

;;

filen='ingest_command_badpix_'+camstr[cam-1]+'_'+string(year, format='(i4.4)')+'.'+string(month, format='(i2.2)')+'.'+string(day, format='(i2.2)')+'.csh'

string='write_badpix instrument="HMI" file_badpix="'+filename+'" series_badpix="hmi.bad_pixel_list" camera='+string(cam,format='(i1.1)')+' t_obs="'+time_stamp+'" nbad='+strcompress(string(nbad),/rem)+';echo $status' 

openw,1,filen
printf,1,string
close,1

if keyword_set(write) then begin
print, "write bad pixel list into data base"
spawn, string, result
stat=result[n_elements(result)-1]
print, 'status: ', stat
if stat ne '0' then print, 'WARNING: Could not ingest bad pixel list into database' else print, "COMPLETED"
endif



end

