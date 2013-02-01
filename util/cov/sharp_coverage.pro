;=====================================================================================
; This program requires use of SSWIDL, and uses routines that
; include, but are not limited to, trim.pro, secstr.pro, rd_tfile.pro,
; fndwrd.pro, getwrd.pro, ydn2md.pro, and stress.pro.
;
; The following describes the necessary inputs. 
;
; ------------------------file1------------------------
; show_coverage -qi ds="hmi.mharp_720s[][][? 400000>NPIX>2000 ?]" key=T_REC low=2012.11.01_00_TAI high=2012.12.01_00_TAI | grep 'UNK' | ./reform_coverage_720.csh > mharp_nov2012.txt
;
; ------------------------file2------------------------
; show_coverage -qi ds="hmi.sharp_720s[][][? 400000>NPIX>2000 ?]" key=T_REC low=2012.11.01_00_TAI high=2012.12.01_00_TAI | grep 'UNK' | ./reform_coverage_720.csh > sharp_nov2012.txt
; 
; ------------------------Lookup Table for Start and End Times------------------------
; Times for the beginning and end of each month are calculated
; using a Unix timestamp converter (an online version is available at:
; onlineconversion.com/unix_time). The following are the beginning 
; and end times normalized to 1/1/10,
; i.e. result= (Seconds between n/1/2010 and 1/1/1970)-(Seconds between 1/1/1970 and 1/1/2010)
;      result= (Seconds between n/1/2010 and 1/1/1970)-1262304000L
; 
; maystart      = 10368000L
; mayend        = 13046400L
; junstart      = 13046400L
; junend        = 15638400L
; julystart     = 15638400L
; julyend       = 18316800L
; auguststart   = 18316800L
; augustend     = 20995200L
; septemberstart= 20995200L
; septemberend  = 23587200L
; octoberstart  = 23587200L
; octoberend    = 26265600L
; novemberstart = 26265600L
; novemberend   = 28857600L
; decemberend   = 31536000L
; januarystart  = 31536000L
; januaryend    = 34214400L
;
; n.b.: The number 536457600L is the number of seconds between January
; 1 2010 and January 1 1993
;
; ------------------------Sample Call --------------------------
;coveragemap, monthnumber=11, monthname='November 2012', junstart=89424000L, junend=92016000L, file1='mharp_nov2012.txt', file5='sharp_nov2012.txt',outfile='/home/mbobra/Coverage_Maps/sharp/2012_11_sharp_def.ps', ds='hmi.sharp_720s', cad=720L
;
;-------------------------------------------------------------------------------------
;=====================================================================================

pro coveragemap, monthnumber=monthnumber, monthname=monthname, junstart=junstart, junend=junend, file1=file1, outfile=outfile, ds=ds, cad=cad

text1  = rd_tfile(file1,3)
lines1 = file_lines(file1)

; FOR COMPLETE; these are numbers equal to 0
if lines1 gt 0 then begin

  Start_time_comp = lonarr(n_elements(text1[0,*])) 
  end_time_comp   = lonarr(n_elements(text1[0,*])) 
  for ii=0L, lines1-1 do start_time_comp[ii]=((long(text1[0,ii]))*cad) - 536457600L 
  for jj=0L, lines1-1 do end_time_comp[jj]=cad+start_time_comp[jj]

  value      = lonarr(n_elements(text1[0,*]))  
  for kk=0, lines1-1 do value[kk]=abs( long(text1[1,kk]) - long(text1[2,kk]) )

  value_zero_pointers=where(value eq 0, count)
  if count eq 0 then begin
     start_time=[0]
     end_time=[0]
  endif else begin 
     start_time = start_time_comp[value_zero_pointers]
     end_time   = end_time_comp[value_zero_pointers]
  endelse


  ;; Populate juntime with start and end times, define junstart and junend statically via lookup table
  ;juntime_start=lonarr(lines1)
  ;for ii=0L, lines1-1 do if start_time[ii] ge junstart and start_time[ii] le junend then juntime_start[ii]=start_time[ii]
  ;if lines1 gt 0 then juntime_end=lonarr(lines1) 
  ;for ii=0L, lines1-1 do if end_time[ii] ge junstart and end_time[ii] le junend then juntime_end[ii]=end_time[ii]

  ; Populate juntime with start and end times, define junstart and junend statically via lookup table
  juntime_start=lonarr(lines1)
  for ii=0L, count-1 do if start_time[ii] ge junstart and start_time[ii] le junend then juntime_start[ii]=start_time[ii]
  if lines1 gt 0 then juntime_end=lonarr(lines1) 
  for ii=0L, count-1 do if end_time[ii] ge junstart and end_time[ii] le junend then juntime_end[ii]=end_time[ii]

  ; Try to get juntime into an array without any 0's  
  ;result0  = where(start_time ge junstart and start_time le junend, count)
  ;Result00 = where(end_time   ge junstart and end_time   le junend, count)
  ; These are the correctly sized arrays normalized to seconds starting
  ; on the first day of the month
  ; juntime_start_sized = (juntime_start[result0])- junstart
  ; juntime_end_sized   = (juntime_end[result00]) - junstart
  juntime_start_sized = juntime_start - junstart
  juntime_end_sized = juntime_end - junstart
  ss=n_elements(juntime_start_sized)
endif


; FOR MISSING; these are numbers greater than 0
if lines1 gt 0 then begin

  Start_time_miss = lonarr(n_elements(text1[0,*])) 
  end_time_miss   = lonarr(n_elements(text1[0,*])) 
  for ii=0L, lines1-1 do start_time_miss[ii]= ((long(text1[0,ii]))*cad) - 536457600L 
  for jj=0L, lines1-1 do end_time_miss[jj]  = cad+start_time_miss[jj]

  value_nonzero_pointers=where(value ne 0, count)
  if count eq 0 then begin
     stime_unk=[0]
     etime_unk=[0]
  endif else begin 
     stime_unk = start_time_miss[value_nonzero_pointers]
     etime_unk = end_time_miss[value_nonzero_pointers]
  endelse

  ; Populate juntime with start and end times, define junstart and junend statically via 
  ; lookup table
  unk_start=lonarr(lines1)
  for ii=0L, count-1 do if stime_unk[ii] ge junstart and stime_unk[ii] lt junend then unk_start[ii]=stime_unk[ii]
  unk_end=lonarr(lines1)
  for ii=0L, count-1 do if etime_unk[ii] ge junstart and etime_unk[ii] lt junend then unk_end[ii]=etime_unk[ii]
  
  ; Try to get juntime into an array without any 0's  
  ;result0_unk=(where(stime_unk ge junstart and stime_unk lt junend, count))
  ;Result00_unk=(where(etime_unk ge junstart and etime_unk lt junend, count))
  ; These are the correctly sized arrays normalized to seconds starting
  ; on the first day of the month
  ;ssized = (unk_start[result0_unk] ) - junstart
  ;esized = (unk_end  [result00_unk]) - junstart
  ssized = stime_unk - junstart
  esized = etime_unk - junstart
  ww=n_elements(ssized)
endif 


; Combine all the arrays into one array that looks like:
; x=[start_time[1], start_time[1], end_time[1], end_time[1], start_time[2], 
;      start_time[2], end_time[2], end_time[2] etc.]
; y=[0, 1, 1, 0, 0, 1, 1, 0]

; COMPLETE
if lines1 gt 0 then begin 
  n=(n_elements(juntime_start_sized)+n_elements(juntime_end_sized))*2
  x=lonarr(n)
  y=lonarr(n)
     for jj=0L, ss-1 do x[4*jj+3]=juntime_end_sized[jj]
     for jj=0L, ss-1 do x[4*jj+2]=juntime_end_sized[jj]
     for jj=0L, ss-1 do x[4*jj+1]=juntime_start_sized[jj]
     for jj=0L, ss-1 do x[4*jj+0]=juntime_start_sized[jj]
     for jj=0L, ss-1 do y[4*jj+3]=0
     for jj=0L, ss-1 do y[4*jj+2]=1
     for jj=0L, ss-1 do y[4*jj+1]=1
     for jj=0L, ss-1 do y[4*jj+0]=0
endif

; MISSING
; the array x1 is not in increasing order
if lines1 gt 0 then begin
  n1=(n_elements(ssized)+n_elements(esized))*2
  x1=lonarr(n1)
  y1=lonarr(n1)
     for kk=0L, ww-1L do x1[4*kk+3]=esized[kk]
     for kk=0L, ww-1L do x1[4*kk+2]=esized[kk]
     for kk=0L, ww-1L do x1[4*kk+1]=ssized[kk]
     for kk=0L, ww-1L do x1[4*kk+0]=ssized[kk]
     for kk=0L, ww-1L do y1[4*kk+3]=0
     for kk=0L, ww-1L do y1[4*kk+2]=1
     for kk=0L, ww-1L do y1[4*kk+1]=1
     for kk=0L, ww-1L do y1[4*kk+0]=0
endif 

; There are 86400 seconds in an average day.
day  = lonarr(32)
for i=0, 31 do day[i]=86400*i 

; Create the appropriate x elements, y elements, and plot data
; color tables are according to david fanning's pickcolor.pro
set_plot,'ps'
!p.font=0
device,xsize=xsize, ysize=ysize,xoffset=1,bits=8,/color,/landscape,/helvetica,/bold,filename=outfile
tvlct, 255, 105, 180, 1 ; pink 
tvlct, 30,  144, 255, 2 ; dodger blue
tvlct, 255, 255, 255, 3 ; white

; Plot the 1 through 29th elements
test=lonarr(1)
test[0]=-1

; plot MISSING in pink
if lines1 gt 0 then begin
for i=0,30 do begin 
  date1  = where(x1  ge day[i] and x1  le day[i+1], count1)
    if count1 gt 0 then begin 
       date_xunk=double(x1[date1]-i*86400L)/double(3600L) 
       date_yunk=y1[date1] 
          if (date_yunk[count1-1] eq 1) then begin
             date_xunk=[date_xunk, 24.0, 24.0]
             date_yunk=[date_yunk, 1, 0]
          endif
          if (date_yunk[0] eq 1) then begin
             date_xunk=[0.0, 0.0, date_xunk]
             date_yunk=[0, 1, date_yunk]
          endif
    endif else begin
       date_xunk=test 
       date_yunk=test 
    endelse
  plot, date_xunk, date_yunk, YTICKFORMAT="(A1)", xtickformat="(A1)", xrange=[0, 24], Xstyle = 1, yticklen=1e-5, ytitle=trim(i+1), xthick=4.0,ythick=4.0,xticklen=0.15,  charsize=.9,charthick=3.0,Position=[0.15, 0.95-0.03*i, 0.9, 0.98-0.03*i], /Nodata, /Noerase
    if date_xunk[0] ne -1 then begin
       plots, date_xunk, date_yunk, thick=.1, color=3
       polyfill, date_xunk, date_yunk, color=1
    endif 
endfor 
endif 

; plot COMPLETE in blue
if lines1 gt 0 then begin
for i=0,30 do begin 
  date  = where(x  ge day[i] and x  le day[i+1], count)
    if count gt 0 then begin 
       date_x=double(x[date]-i*86400L)/double(3600L) 
       date_y=y[date] 
          if (date_y[count-1] eq 1) then begin
             date_x=[date_x, 24.0, 24.0]
             date_y=[date_y, 1, 0]
          endif
          if (date_y[0] eq 1) then begin
             date_x=[0.0, 0.0, date_x]
             date_y=[0, 1, date_y]
    endif 
    endif else begin
       date_x=test 
       date_y=test 
    endelse
  plot, date_x, date_y, YTICKFORMAT="(A1)", xtickformat="(A1)", xrange=[0, 24], Xstyle = 1, yticklen=1e-5, ytitle=trim(i+1), xthick=4.0,ythick=4.0,xticklen=0.15,  charsize=.9,charthick=3.0,Position=[0.15, 0.95-0.03*i, 0.9, 0.98-0.03*i], /Nodata, /Noerase
    if date_x[0] ne -1 then begin
       plots, date_x, date_y, thick=.1, color=3
       polyfill, date_x, date_y, color=2
    endif 
endfor 
endif

; Make a vector of 16 points, A[i] = 2pi/16:  
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL  


keywordlist=['PARTIALLY MISSING','COMPLETE','TOTALLY MISSING']
legend,keywordlist,/fill,psym=8+intarr(3), color=[1, 2, 3], chars=.45, pos=[21820,16700],/device

; For the device, !D.X_SIZE = 24130 and !D.Y_SIZE=17780

xyouts, .3, .18, 'Plot Made '+!stime, chars=.7,/device

if lines1 eq 0 then xyouts, 3400, 8000, 'NO SHARPS', chars=7.5, /device 

xyouts, 10200, 17520, 'Coverage Map for ' + ds, chars=.7, /device

xyouts, 12000, 400, 'Hours', chars=1, /device 

xyouts, 1600, 8500, text_axes=0, orientation=90, 'Days of '+Monthname, /device

device,/close

set_plot,'x'

end
