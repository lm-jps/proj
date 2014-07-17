pro make_dark, year, month, day, cam, ftsid, time_stamp, darks, write=write



spawn, 'show_info -q ds="hmi.lev1_nrt['+string(year, format='(i4.4)')+'.'+string(month, format='(i2.2)')+'.'+string(day, format='(i2.2)')+'/1d][?HFTSACID = '+strcompress(string(ftsid),/rem)+'?]" key=FSN', result
n=n_elements(result)

if n le 1 then begin
print, "use lev1"
spawn, 'show_info -q ds="hmi.lev1['+string(year, format='(i4.4)')+'.'+string(month, format='(i2.2)')+'.'+string(day, format='(i2.2)')+'/1d][?HFTSACID = '+strcompress(string(ftsid),/rem)+'?]" key=FSN', result

series="hmi.lev1"
endif else begin
series="hmi.lev1_nrt"
endelse


fsn_min=long(result[0])
fsn_max=long(result[n-1])

fsn=lonarr(n)
for i=0, n-1 do fsn[i]=long(result[i])

print, n, fsn_min, fsn_max

selectnames_fsn, nm, "hmi.lev0a", fsn_min, fsn_max, nfiles

if (fsn_max-fsn_min ne (nfiles-1)) then begin print, "incorrect number of files" & stop  & endif 

spawn, 'show_info -q ds="'+series+'['+string(year, format='(i4.4)')+'.'+string(month, format='(i2.2)')+'.'+string(day, format='(i2.2)')+'/1d][?HFTSACID = '+strcompress(string(ftsid),/rem)+'?][?HCAMID='+string(cam-1,format='(i1.1)')+'?]" key=FSN', result

nn_ds=n_elements(result)
fsn_ds=lonarr(nn_ds)
for i=0, nn_ds-1 do fsn_ds[i]=long(result([i]))
darkinds=fsn_ds-fsn_min



if not keyword_set(time_stamp) then begin
spawn, 'show_info -q ds="'+series+'['+string(year, format='(i4.4)')+'.'+string(month, format='(i2.2)')+'.'+string(day, format='(i2.2)')+'/1d][?HFTSACID = '+strcompress(string(ftsid),/rem)+'?][?HCAMID='+string(cam-1,format='(i1.1)')+'?]" key=T_OBS', result_tm

time_stamp=result_tm[n_elements(result_tm)-1]
endif





ndk_s=n_elements(darkinds)/5


darks=fltarr(4096,4096)


for k=0, ndk_s-1 do begin
print, k, ' of ', ndk_s
readimages, nm[darkinds[k*5:k*5+4]], imx, head, /noshow
for j=0, 4095 do for i=0, 4095 do darks[i,j]=darks[i,j]+median(imx[i,j,*])

endfor


darks=darks/float(ndk_s)

camstr=['side', 'front']

openw,1,'dark_'+camstr[cam-1]+'_'+strcompress(string(fsn[0]), /rem)+'.bin'
writeu,1,darks
close,1




string_side='write_dark instrument="HMI" file_dark="dark_'+camstr[cam-1]+'_'+strcompress(string(fsn[0]), /rem)+'.bin"  series_dark="hmi.inspect_dark" camera='+string(cam,format='(i1.1)')+' t_obs="'+time_stamp+'" fsn_list_dark='+strcompress(string(fsn[darkinds[0]]), /rem)+','+strcompress(string(fsn[darkinds[ndk_s-1]]),/rem)+' nx=4096 ny=4096'


openw,1,'ingest_command_dark_'+camstr[cam-1]+'_'+strcompress(string(fsn[0]), /rem)+'.csh'
printf,1,string_side
close,1


print, string_side

if keyword_set(write) then begin
print, "writing dark into database"
spawn, string_side

endif

print, "DARK COMPLETED!"

end

;write_dark instrument="HMI" file_dark="/scr21/richard/hmi/dark_front_12835504.bin"  series_dark="hmi.dark" camera=2 t_obs="2010.11.06_00:00:00_UTC" fsn_list_dark=12835504,12837244 nx=4096 ny=4096
;write_dark instrument="HMI" file_dark="/scr21/richard/hmi/dark_side_12835504.bin"  series_dark="hmi.dark" camera=1 t_obs="2010.11.06_00:00:00_UTC" fsn_list_dark=12835505,12837245 nx=4096 ny=4096
