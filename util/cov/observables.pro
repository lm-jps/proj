;=====================================================================================
; This program requires use of SSWIDL, and uses routines that
; include, but are not limited to, trim.pro, secstr.pro, rd_tfile.pro,
; fndwrd.pro, getwrd.pro, ydn2md.pro, and stress.pro.
;
; The following describes the necessary inputs. 
;
; ------------------------file1------------------------
; file1 is the result of show_coverage on any given level 1.5 drms series. 
; In particular, file1 contains times when level 1.5 data are missing.
; I obtain the missing level 1.5 data times using the following commands:
; show_coverage -qi ds=hmi_test.V_45s low=2010.09.01_00_TAI high=2010.10.01_00_TAI | grep 'MISS' | reform_coverage.csh > ~/public_html/Coverage_Maps/miss_month.txt
; for example:
; show_coverage -qi ds=hmi.V_45s low=2011.01.01_00_TAI high=2011.02.01_00_TAI | grep 'MISS' | reform_coverage.csh > ~/public_html/Coverage_Maps/miss_january.txt
; ------------------------file2------------------------
; file2 is the result of x_stats_day_HMI.csh run on the hmi.lev0a data.
; In particular, file2 contains the times that the ISS loop is open:
; I obtain the loop open times using the following commands:
; cp ~thailand/stats_lev0/HMI/HMI_LOOPOPEN_LIST.2010.06* ~/public_html/Coverage_Maps/staging
; grep '2010.06.' HMI_* > ~/public_html/Coverage_Maps/staging/EOI_LOOP_Month.txt
;
; ------------------------file3------------------------
; file3 is the result of a mask being applied to show_coverage in
; order to plot observables that have lower quality bits set
; for example:
; show_coverage -qi ds="hmi.V_45s[][2]" low=2010.04.01_00_TAI high=2010.05.01_00_TAI mask=0xffffffff | grep 'MISS' | ./reform_coverage.csh > ~/Coverage_Maps/miss_april_mask.txt
;
; ------------------------file5------------------------
; file5 is the result of show_coverage on any given level 1.5 drms series. 
; In particular, file5 contains times when level 1.5 data are unknown.
; I obtain the missing level 1.5 data times using the following commands:
; show_coverage -qi ds=hmi.V_45s low=2010.09.01_00_TAI high=2010.10.01_00_TAI | grep UNK | reform_coverage.csh > ~/public_html/Coverage_Maps/unk_month.txt
; for example:
; show_coverage -qi ds=hmi.V_45s low=2011.01.01_00_TAI high=2011.02.01_00_TAI | grep 'UNK' | reform_coverage.csh > ~/public_html/Coverage_Maps/unk_january.txt
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
; coveragemap, monthnumber=1,  junstart=31536000L,junend=34214400L,
; monthname='January 2011', file1='miss_january.txt',file2='EOI_LOOP_January.txt', file3='miss_january_mask.txt',file5='unk_january.txt', outfile='~/public_html/Coverage_Maps/2011_01_map.ps', ds='hmi.V_45s', cad=45L
;-------------------------------------------------------------------------------------
;=====================================================================================

pro coveragemap, monthnumber=monthnumber, monthname=monthname, junstart=junstart, junend=junend, file1=file1, file2=file2, file3=file3, file5=file5, outfile=outfile, ds=ds, cad=cad

text1  = rd_tfile(file1,1)
text2  = rd_tfile(file2,1)
text3  = rd_tfile(file3,1)
text5  = rd_tfile(file5,1)
lines1 = file_lines(file1)
lines3 = file_lines(file3)
lines2 = file_lines(file2)
lines5 = file_lines(file5)

; FOR FILE 1
; populate start times and normalize epoch to 1/1/10

; in the end, the strmid function is not the smartest way to do things. if there is a problem, it is likely with the location of the characters that strmid is parsing.

if lines1 gt 0 then begin
  Start_time=lonarr(n_elements(text1)) 
  for ii=0, lines1-1 do start_time[ii]=((long((strmid(text1[ii], 4, 9))))*cad) - 536457600L 
  end_time=lonarr(n_elements(text1))  
    if (cad eq 45) then begin
      for jj=0, lines1-1 do end_time[jj]=((long((strmid(text1[jj], 14))))*cad)+start_time[jj] 
    endif else begin
      for jj=0, lines1-1 do end_time[jj]=((long((strmid(text1[jj], 12))))*cad)+start_time[jj] 
    endelse 
; Populate juntime with start and end times, define junstart and junend statically via lookup table
  juntime_start=lonarr(lines1)
  for ii=0L, lines1-1 do if start_time[ii] ge junstart and start_time[ii] le junend then juntime_start[ii]=start_time[ii]
  if lines1 gt 0 then juntime_end=lonarr(lines1) 
  for ii=0L, lines1-1 do if end_time[ii] ge junstart and end_time[ii] le junend then juntime_end[ii]=end_time[ii]
  ; Try to get juntime into an array without any 0's  
  result0=where(start_time ge junstart and start_time le junend, count)
  Result00=where(end_time ge junstart and end_time le junend, count)
  ; These are the correctly sized arrays normalized to seconds starting
  ; on the first day of the month
  juntime_start_sized=(juntime_start[result0])-junstart
  juntime_end_sized = (juntime_end[result00])- junstart
  ;juntime_start_sized=start_time - junstart
  ;juntime_end_sized  =end_time - junstart
  ss=n_elements(juntime_start_sized)
endif

; FOR FILE 5
; populate start times and normalize epoch to 1/1/10
if lines5 gt 0 then begin 
  Stime_unk=lonarr(n_elements(text5)) 
  for ii=0L, lines5-1 do stime_unk[ii]=((long((strmid(text5[ii], 4, 9))))*cad)-536457600L
  etime_unk=lonarr(n_elements(text5))  
    if (cad eq 45) then begin
      for jj=0L, lines5-1 do etime_unk[jj]=((long((strmid(text5[jj], 13))))*cad)+stime_unk[jj] 
    endif else begin
      for jj=0L, lines5-1 do etime_unk[jj]=((long((strmid(text5[jj], 11))))*cad)+stime_unk[jj] 
    endelse 
  ; Populate juntime with start and end times, define junstart and junend statically via lookup table
  unk_start=lonarr(lines5)
  for ii=0L, lines5-1 do if stime_unk[ii] ge junstart and stime_unk[ii] lt junend then unk_start[ii]=stime_unk[ii]
  unk_end=lonarr(lines5)
  for ii=0L, lines5-1 do if etime_unk[ii] ge junstart and etime_unk[ii] lt junend then unk_end[ii]=etime_unk[ii]
  ; Try to get juntime into an array without any 0's  
  result0_unk=(where(stime_unk ge junstart and stime_unk lt junend, count))
  Result00_unk=(where(etime_unk ge junstart and etime_unk lt junend, count))
  ; These are the correctly sized arrays normalized to seconds starting
  ; on the first day of the month
  ;ssized = (unk_start[result0_unk] ) - junstart
  ;esized = (unk_end  [result00_unk]) - junstart
  ssized = stime_unk - junstart
  esized = etime_unk - junstart
  ww=n_elements(ssized)
endif 

; FOR FILE 2
; populate start times for loop open
if lines2 gt 0 then begin
  result_loop=where(strmid(text2, 54,1) eq 'N', count)
  ;result_loop=where(strmid(text2, 4,1) eq 'L', count)
  Stime_loop=lonarr(count)
  etime_loop=lonarr(count)
  loop_count=count
  ;For jj=0L, count-1 do Stime_loop[jj]=(secstr(strmid(text2(result_loop[jj]),48, 11))+(double(strmid(text2(result_loop[jj]), 26, 2)*86400L)))-86400L
  For jj=0L, count-1 do Stime_loop[jj]=(secstr(strmid(text2(result_loop[jj]),19, 11))+(double(strmid(text2(result_loop[jj]), 17, 2)*86400L)))-86400L
  for jj=0L, count-1 do etime_loop[jj]=stime_loop[jj]+double(cad)
  ff=n_elements(stime_loop)
endif

; FOR FILE 3
; populate start times and normalize epoch to 1/1/10
if lines3 gt 0 then begin
  Start_time_mask=lonarr(n_elements(text3)) 
  for ii=0, lines3-1 do start_time_mask[ii]=((long((strmid(text3[ii], 4, 9))))*cad) - 536457600L 
  end_time_mask=lonarr(n_elements(text3))  
    if (cad eq 45) then begin
      for jj=0, lines3-1 do end_time_mask[jj]=((long((strmid(text3[jj], 14))))*cad)+start_time_mask[jj] 
    endif else begin
      for jj=0, lines3-1 do end_time_mask[jj]=((long((strmid(text3[jj], 12))))*cad)+start_time_mask[jj] 
    endelse 
; Populate juntime with start and end times, define junstart and junend statically via lookup table
  juntime_start_mask=lonarr(lines3)
  for ii=0L, lines3-1 do if start_time_mask[ii] ge junstart and start_time_mask[ii] le junend then juntime_start_mask[ii]=start_time_mask[ii]
  if lines3 gt 0 then juntime_end_mask=lonarr(lines3) 
  for ii=0L, lines3-1 do if end_time_mask[ii] ge junstart and end_time_mask[ii] le junend then juntime_end_mask[ii]=end_time_mask[ii]
  ; Try to get juntime into an array without any 0's  
  result0_mask=where(start_time_mask ge junstart and start_time_mask le junend+1, count)
  Result00_mask=where(end_time_mask ge junstart and end_time_mask le junend+1, count)
  ; These are the correctly sized arrays normalized to seconds starting
  ; on the first day of the month
  juntime_start_sized_mask=(juntime_start_mask[result0_mask])-junstart
  juntime_end_sized_mask = (juntime_end_mask[result00_mask])- junstart
  rr=n_elements(juntime_start_sized_mask)
  mm=n_elements(juntime_end_sized_mask)
endif

; Combine all the arrays into one array that looks like:
; x=[start_time[1], start_time[1], end_time[1], end_time[1], start_time[2], 
;      start_time[2], end_time[2], end_time[2] etc.]
; y=[0, 1, 1, 0, 0, 1, 1, 0]

; MISS
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

; UNK
; the array x1 is not in increasing order
if lines5 gt 0 then begin
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

; LOOP
if lines2 gt 0 then begin
  n2=(loop_count)*4
  x2=lonarr(n2)
  y2=lonarr(n2)
  for jj=0L, ff-1 do x2[4*jj+3]=etime_loop[jj]
  for jj=0L, ff-1 do x2[4*jj+2]=etime_loop[jj]
  for jj=0L, ff-1 do x2[4*jj+1]=stime_loop[jj]
  for jj=0L, ff-1 do x2[4*jj+0]=stime_loop[jj]
  for jj=0L, ff-1 do y2[4*jj+3]=0
  for jj=0L, ff-1 do y2[4*jj+2]=1
  for jj=0L, ff-1 do y2[4*jj+1]=1
  for jj=0L, ff-1 do y2[4*jj+0]=0
endif 

; MASK
if lines3 gt 0 then begin 
  n3=(n_elements(juntime_start_sized_mask)+n_elements(juntime_end_sized_mask))*2
  x3=lonarr(n3)
  y3=lonarr(n3)
  for jj=0L, mm-1 do x3[4*jj+3]=juntime_end_sized_mask[jj]
  for jj=0L, mm-1 do x3[4*jj+2]=juntime_end_sized_mask[jj]
  for jj=0L, mm-1 do x3[4*jj+1]=juntime_start_sized_mask[jj]
  for jj=0L, mm-1 do x3[4*jj+0]=juntime_start_sized_mask[jj]
  for jj=0L, mm-1 do y3[4*jj+3]=0
  for jj=0L, mm-1 do y3[4*jj+2]=1
  for jj=0L, mm-1 do y3[4*jj+1]=1
  for jj=0L, mm-1 do y3[4*jj+0]=0
endif

; There are 86400 seconds in an average day.
day  = lonarr(32)
for i=0, 31 do day[i]=86400*i 

; Create the appropriate x elements, y elements, and plot data
; color tables are according to david fanning's pickcolor.pro
set_plot,'ps'
!p.font=0
device,xsize=xsize, ysize=ysize,xoffset=1,bits=8,/color,/landscape,/helvetica,/bold,filename=outfile
tvlct, 255, 255, 255, 7 ; white
tvlct, 255, 105, 180, 1 ; pink 
tvlct, 30,  144, 255, 2 ; dodger blue
tvlct, 218, 165, 32,  3 ; gold
tvlct,  50, 205, 50,  4 ; lime green
tvlct, 210, 105, 30,  5 ; chocolate
tvlct,   0,   0,  0,  6 ; black
tvlct, 110, 110, 110, 8 ; dark grey

; Plot the 1 through 29th elements
test=lonarr(1)
test[0]=-1

; plot the UNK data from show_coverage 
; UNK data is plotted in pink
if lines5 gt 0 then begin
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
       plots, date_xunk, date_yunk, thick=.1, color=7
       polyfill, date_xunk, date_yunk, color=1
    endif 
endfor 
endif 

; plot the MASK data from show_coverage
; MASK data is plotted in chocolate
if lines3 gt 0 then begin
for i=0,30 do begin 
  date3 = where(x3 ge day[i] and x3 le day[i+1], count3)
    if count3 gt 0 then begin 
       date_x=double(x3[date3]-i*86400L)/double(3600L) 
       date_y=y3[date3] 
          if (date_y[count3-1] eq 1) then begin
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
  plot, date_x, date_y, YTICKFORMAT="(A1)", xtickformat="(A1)", xrange=[0, 24], Xstyle = 1, yticklen=1e-5, ytitle=trim(i+1), xthick=4.0,ythick=4.0,xticklen=0.15, charsize=.9,charthick=3.0,Position=[0.15, 0.95-0.03*i, 0.9, 0.98-0.03*i], /Nodata, /Noerase
    if date_x[0] ne -1 then begin
       plots, date_x, date_y, thick=.1, color=4
       polyfill, date_x, date_y, color=4
    endif 
endfor 
endif

; plot the MISS data from show_coverage 
; MISS data is plotted in blue
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
       plots, date_x, date_y, thick=.1, color=6
       polyfill, date_x, date_y, color=2
    endif 
endfor 
endif

; plot the iss loop open data
; iss loop open data is plotted in gold
if lines2 gt 0 then begin 
 for i=0,30 do begin 
  date2 = where(x2 ge day[i] and x2 le day[i+1], count2)
    if count2 gt 0 then begin 
      date_x=double(x2[date2]-i*86400L)/double(3600L) 
      date_y=y2[date2] 
      if (date_y[count2-1] eq 1) then begin
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
    plots, date_x, date_y, thick=0.01, color=3
    polyfill, date_x, date_y, color=3
  endif
endfor 
endif 


; Make a vector of 16 points, A[i] = 2pi/16:  
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL  


keywordlist=['missing', 'no guiding','unknown','quality','available']
legend,keywordlist,/fill,psym=8+intarr(5), color=[2, 3, 1, 4,7], chars=.75, pos=[21820,16700],/device

; For the device, !D.X_SIZE = 24130 and !D.Y_SIZE=17780

xyouts, .3, .18, 'Plot Made '+!stime, chars=.7,/device

xyouts, 10200, 17520, 'Coverage Map for ' + ds, chars=.7, /device

xyouts, 12000, 400, 'Hours', chars=1, /device 

xyouts, 1600, 8500, text_axes=0, orientation=90, 'Days of '+Monthname, /device

device,/close

set_plot,'x'

end
