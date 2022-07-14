;=====================================================================================
; This program requires use of SSWIDL, and uses routines that
; include, but are not limited to, trim.pro, secstr.pro, rd_tfile.pro,
; fndwrd.pro, getwrd.pro, ydn2md.pro, and stress.pro.
;
; The following describes the necessary file inputs:
;
; ------------------------file1------------------------
; file1 is the result of show_coverage on any given observables series. 
; In particular, file1 contains times when data are missing.
; I obtain the missing data times using the following example command:
; show_coverage -qi ds="hmi.M_720s[][1]" low=2012.07.01_00_TAI high=2012.08.01_00_TAI | grep 'MISS' | ./reform_coverage_720.csh > /home/mbobra/Coverage_Maps/fd10/miss_july2012_s720m.txt
;
; (n.b. the script reform_coverage creates 24 hour chunks of show_coverage output (for ease of plotting)) 
;
; ------------------------file4------------------------
; file4 is a join on the result of two show_info calls: 
; call 1: show_info -q key=t_rec_index,RSUN_OBS,CDELT1 ds=hmi.S_720s"[2012.07.01_00:00:00_TAI-2012.08.01_00:00:00_TAI][? QUALITY >= 0 ?]" > pixels_jul2012.txt
; call 2: show_info -q key=t_rec_index,INVNPRCS ds=hmi.ME_720s_fd10"[2012.07.01_00:00:00_TAI-2012.08.01_00:00:00_TAI]" > invnprcs_jul2012.txt
; then the following unix command is run:
;
; join -t'	' invnprcs_jul2012.txt pixels_jul2012.txt > file4.txt
;
; (n.b. the join must be tab delimited; to type a tab, press control-v-tab).
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
; ------------------------------------------------------------------------------------------
;
; INPUTS: MONTHNUMBER = number of the month 
;         JUNSTART    = start time as indicated in the lookup table
;         JUNEND      = end time as indicated in the lookup table
;         Monthname   = string value with Month and Year
;         file1       = described above
;         file4       = described above
;         outfile     = string value of the output filename, including path
;         ds          = string value of the drms series
;         cad         = cadence in type long
;
; ------------------------Sample Call --------------------------
; coveragemap, monthnumber=7, junstart=78796800L, junend=81475200L, monthname='July 2012', file1='miss_july2012_s720m.txt', file4='join_jul2012.txt', outfile='/home/mbobra/Coverage_Maps/fd10/2012_07_map___fd10.ps', ds='hmi.ME_720s_fd10', cad=720L
;-------------------------------------------------------------------------------------
;=====================================================================================

pro coveragemap, monthnumber=monthnumber, monthname=monthname, junstart=junstart, junend=junend, file1=file1, file4=file4, outfile=outfile, ds=ds, cad=cad

text1  = rd_tfile(file1,1)
lines1 = file_lines(file1)

text4  = rd_tfile(file4,4)
lines4 = file_lines(file4)

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
  ss=n_elements(juntime_start_sized)
endif


; FOR FILE 4
; populate start times and normalize epoch to 1/1/10

if lines4 gt 0 then begin 
  Stime_inv=lonarr(n_elements(text4)) 
  for ii=0L, lines4-1 do stime_inv[ii] = ((long((text4[0,ii])))*cad)-536457600L
  
  stime_values=fltarr(n_elements(text4)) 
  for ii=0L, lines4-1 do begin
  	;print,'1',(long(text4[1,ii]))
  	;print,'2',float(text4[2,ii]) 
  	;print,'3',float(text4[3,ii]),format='(f10.6)'
  	;print,'original',(long((text4[0,ii])))
        ;print,'ii',ii
    stime_values[ii] =(long(text4[1,ii])) / ( !pi *( float(text4[2,ii])  * (1./float(text4[3,ii]) ) )^2) 
    if stime_values[ii] ge 1.0000000 then stime_values[ii]=0.99999999999
  endfor

;  old method = just divide by max(TOTVALS) in hmi.M_720s
;  for ii=0L, lines4-1 do stime_values[ii]=(float(strmid(text4[ii], 7, 8)))/(12387771.)
;for ii=0L, lines4-1 do stime_values[ii]=   (long(text4[1,ii])) /(12387771.)

  etime_inv=lonarr(n_elements(text4))  
    if (cad eq 45) then begin
      for jj=0L, lines4-1 do etime_inv[jj]=((long((text4[0,jj])))*cad)+stime_inv[jj] 
    endif else begin
      for jj=0L, lines4-1 do etime_inv[jj]=cad+stime_inv[jj] 
    endelse 
  ; Populate juntime with start and end times, define junstart and junend statically via lookup table
  inv_start=lonarr(lines4)
  for ii=0L, lines4-1 do if stime_inv[ii] ge junstart and stime_inv[ii] lt junend then inv_start[ii]=stime_inv[ii]
  inv_end=lonarr(lines4)
  for ii=0L, lines4-1 do if etime_inv[ii] ge junstart and etime_inv[ii] lt junend then inv_end[ii]=etime_inv[ii]
  ; Try to get juntime into an array without any 0's  
  result0_inv=(where(stime_inv ge junstart and stime_inv lt junend, count))
  Result00_inv=(where(etime_inv ge junstart and etime_inv lt junend, count))
  ; These are the correctly sized arrays normalized to seconds starting
  ; on the first day of the month
  invssized = stime_inv - junstart
  invesized = etime_inv - junstart
  uu=n_elements(invssized)

endif 


; Combine all the arrays into one array that looks like:
; x=[start_time[1], start_time[1], end_time[1], end_time[1], start_time[2], 
;      start_time[2], end_time[2], end_time[2] etc.]
; y=[0, 1, 1, 0, 0, 1, 1, 0]


; MISS
; the y-elements are deliberately set to zero so that nothing is
; plotted, however, the MISS values array is kept because if the INV
; array is zero then no plot will be created
if lines1 gt 0 then begin 
  n=(n_elements(juntime_start_sized)+n_elements(juntime_end_sized))*2
  x=lonarr(n)
  y=lonarr(n)
  for jj=0L, ss-1 do x[4*jj+3]=juntime_end_sized[jj]
  for jj=0L, ss-1 do x[4*jj+2]=juntime_end_sized[jj]
  for jj=0L, ss-1 do x[4*jj+1]=juntime_start_sized[jj]
  for jj=0L, ss-1 do x[4*jj+0]=juntime_start_sized[jj]
  for jj=0L, ss-1 do y[4*jj+3]=0
  for jj=0L, ss-1 do y[4*jj+2]=0
  for jj=0L, ss-1 do y[4*jj+1]=0
  for jj=0L, ss-1 do y[4*jj+0]=0
endif


; INV
if lines4 gt 0 then begin 
  n=(n_elements(invssized)+n_elements(invesized))*2
  x4=lonarr(n)
  y4=fltarr(n)
  z4=fltarr(n)
  for jj=0L, uu-1 do x4[4*jj+3]=invesized[jj]
  for jj=0L, uu-1 do x4[4*jj+2]=invesized[jj]
  for jj=0L, uu-1 do x4[4*jj+1]=invssized[jj]
  for jj=0L, uu-1 do x4[4*jj+0]=invssized[jj]
  for jj=0L, uu-1 do y4[4*jj+3]=0
  for jj=0L, uu-1 do y4[4*jj+2]=1
  for jj=0L, uu-1 do y4[4*jj+1]=1
  for jj=0L, uu-1 do y4[4*jj+0]=0
  for jj=0L, uu-1 do z4[4*jj+3]=stime_values[jj]
  for jj=0L, uu-1 do z4[4*jj+2]=stime_values[jj]
  for jj=0L, uu-1 do z4[4*jj+1]=stime_values[jj]
  for jj=0L, uu-1 do z4[4*jj+0]=stime_values[jj]
endif

; There are 86400 seconds in an average day.
day  = lonarr(32)
for i=0, 31 do day[i]=86400*i 

; Create the appropriate x elements, y elements, and plot data
; color tables are according to david fanning's pickcolor.pro

set_plot,'ps'
!p.font=0
device,xsize=xsize, ysize=ysize,xoffset=1,bits=8,/color,/landscape,/helvetica,/bold,filename=outfile

; loading the pre-defined IDL color table 13 (rainbow)
; printing the following lines this way:
; loadct, 13
; tvlct, r,g,b, /get
; for i=0,255 do print,'tvlct,',r[i*2],',',g[i*2],',',b[i*2]

tvlct, 255, 255, 255, 7 ; white

tvlct,   0,   0,   0,      10
tvlct,   9,   0,   7,      11
tvlct,  18,   0,  14,      12
tvlct,  27,   0,  23,      13
tvlct,  36,   0,  32,      14
tvlct,  45,   0,  43,      15
tvlct,  54,   0,  53,      16
tvlct,  61,   0,  63,      17
tvlct,  68,   0,  72,      18
tvlct,  72,   0,  81,      19
tvlct,  77,   0,  91,      20
tvlct,  80,   0, 100,      21
tvlct,  83,   0, 109,      22
tvlct,  84,   0, 118,      23
tvlct,  87,   0, 127,      24
tvlct,  86,   0, 136,      25
tvlct,  87,   0, 145,      26
tvlct,  85,   0, 154,      27
tvlct,  84,   0, 163,      28
tvlct,  83,   0, 173,      29
tvlct,  78,   0, 182,      30
tvlct,  76,   0, 191,      31
tvlct,  70,   0, 200,      32
tvlct,  66,   0, 209,      33
tvlct,  58,   0, 218,      34
tvlct,  53,   0, 227,      35
tvlct,  43,   0, 236,      36
tvlct,  36,   0, 245,      37
tvlct,  25,   0, 255,      38
tvlct,  16,   0, 255,      39
tvlct,   4,   0, 255,      40
tvlct,   0,   4, 255,      41
tvlct,   0,  16, 255,      42
tvlct,   0,  25, 255,      43
tvlct,   0,  38, 255,      44
tvlct,   0,  46, 255,      45
tvlct,   0,  55, 255,      46
tvlct,   0,  67, 255,      47
tvlct,   0,  76, 255,      48
tvlct,   0,  89, 255,      49
tvlct,   0,  97, 255,      50
tvlct,   0, 110, 255,      51
tvlct,   0, 119, 255,      52
tvlct,   0, 131, 255,      53
tvlct,   0, 140, 255,      54
tvlct,   0, 152, 255,      55
tvlct,   0, 161, 255,      56
tvlct,   0, 174, 255,      57
tvlct,   0, 182, 255,      58
tvlct,   0, 195, 255,      59
tvlct,   0, 203, 255,      60
tvlct,   0, 216, 255,      61
tvlct,   0, 225, 255,      62
tvlct,   0, 233, 255,      63
tvlct,   0, 246, 255,      64
tvlct,   0, 255, 255,      65
tvlct,   0, 255, 242,      66
tvlct,   0, 255, 233,      67
tvlct,   0, 255, 220,      68
tvlct,   0, 255, 212,      69
tvlct,   0, 255, 199,      70
tvlct,   0, 255, 191,      71
tvlct,   0, 255, 178,      72
tvlct,   0, 255, 170,      73
tvlct,   0, 255, 157,      74
tvlct,   0, 255, 148,      75
tvlct,   0, 255, 135,      76
tvlct,   0, 255, 127,      77
tvlct,   0, 255, 114,      78
tvlct,   0, 255, 106,      79
tvlct,   0, 255,  97,      80
tvlct,   0, 255,  84,      81
tvlct,   0, 255,  76,      82
tvlct,   0, 255,  63,      83
tvlct,   0, 255,  55,      84
tvlct,   0, 255,  42,      85
tvlct,   0, 255,  34,      86
tvlct,   0, 255,  21,      87
tvlct,   0, 255,  12,      88
tvlct,   0, 255,   0,      89
tvlct,   8, 255,   0,      90
tvlct,  21, 255,   0,      91
tvlct,  29, 255,   0,      92
tvlct,  42, 255,   0,      93
tvlct,  51, 255,   0,      94
tvlct,  63, 255,   0,      95
tvlct,  72, 255,   0,      96
tvlct,  80, 255,   0,      97
tvlct,  93, 255,   0,      98
tvlct, 101, 255,   0,      99
tvlct, 114, 255,   0,     100
tvlct, 123, 255,   0,     101
tvlct, 135, 255,   0,     102
tvlct, 144, 255,   0,     103
tvlct, 157, 255,   0,     104
tvlct, 165, 255,   0,     105
tvlct, 178, 255,   0,     106
tvlct, 187, 255,   0,     107
tvlct, 199, 255,   0,     108
tvlct, 208, 255,   0,     109
tvlct, 221, 255,   0,     110
tvlct, 229, 255,   0,     111
tvlct, 242, 255,   0,     112
tvlct, 250, 255,   0,     113
tvlct, 255, 250,   0,     114
tvlct, 255, 238,   0,     115
tvlct, 255, 229,   0,     116
tvlct, 255, 216,   0,     117
tvlct, 255, 208,   0,     118
tvlct, 255, 195,   0,     119
tvlct, 255, 187,   0,     120
tvlct, 255, 174,   0,     121
tvlct, 255, 165,   0,     122
tvlct, 255, 153,   0,     123
tvlct, 255, 144,   0,     124
tvlct, 255, 131,   0,     125
tvlct, 255, 123,   0,     126
tvlct, 255, 110,   0,     127
tvlct, 255, 102,   0,     128
tvlct, 255,  89,   0,     129
tvlct, 255,  80,   0,     130
tvlct, 255,  72,   0,     131
tvlct, 255,  59,   0,     132
tvlct, 255,  51,   0,     133
tvlct, 255,  38,   0,     134
tvlct, 255,  29,   0,     135
tvlct, 255,  17,   0,     136
tvlct, 255,   8,   0,     137

; Plot the 1 through 29th elements
test=lonarr(1)
test[0]=-1

testinv=fltarr(1)
testinv[0]=-1

; plot the MISS data from show_coverage 
; this section is just to get the grid (y values of MISS are always zero)
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
    plot, date_x, date_y, YTICKFORMAT="(A1)", xtickformat="(A1)", xrange=[0, 24], Xstyle = 1, yticklen=1e-5, ytitle=trim(i+1),xthick=4.0,ythick=4.0,xticklen=0.15, charsize=.9,charthick=3.0,Position=[0.15, 0.95-0.03*i, 0.9, 0.98-0.03*i], /Nodata, /Noerase
    if date_x[0] ne -1 then begin
       plots, date_x, date_y, thick=.1, color=7
       polyfill, date_x, date_y, color=7
    endif 
endfor 
endif

; plot the INV data from show_coverage 
; INV data is plotted in many colors
if lines4 gt 0 then begin
 for i=0,30 do begin 
     date4  = where(x4  ge day[i] and x4  le day[i+1], count4)
     print,count4
    if count4 gt 0 then begin 
       date_xinv=double(x4[date4]-i*86400L)/double(3600L) 
       date_yinv=y4[date4]
       date_zinv=z4[date4]
          if (date_yinv[count4-1] eq 1) then begin
             date_xinv=[date_xinv, 24.0, 24.0]
             date_yinv=[date_yinv, 1., 0.]
             date_zinv=[date_zinv,date_zinv[count4-1],date_zinv[count4-1]]
          endif
          if (date_yinv[0] eq 1.) then begin
             date_xinv=[0.0, 0.0, date_xinv]
             date_yinv=[0., 1., date_yinv]
             date_zinv=[date_zinv[0],date_zinv[0],date_zinv]
          endif
    endif else begin
       date_xinv=testinv
       date_yinv=testinv 
       date_zinv=testinv
   endelse
        plot, date_xinv, date_yinv, YTICKFORMAT="(A1)", xtickformat="(A1)", xrange=[0, 24], Xstyle = 1, yticklen=1e-5, ytitle=trim(i+1),xthick=4.0,ythick=4.0,xticklen=0.15, charsize=.9,charthick=3.0, Position=[0.15, 0.95-0.03*i, 0.9, 0.98-0.03*i], /Nodata,/Noerase
    if date_xinv[0] ne -1. then begin
        plots, date_xinv, date_yinv, thick=.01, color=7
        for ii=10,135 do begin
        w = where(date_zinv gt ((ii-10)/126.) and date_zinv le ((ii-9)/126.),count)
        if count gt 0 then polyfill,date_xinv[w],date_yinv[w],color=ii
        endfor
    endif 
 endfor 
endif 

; For the device, !D.X_SIZE = 24130 and !D.Y_SIZE=17780

; Plot the colorbar
loadct,13,ncolors=255
colorbar,ncolors=255,minrange=0,maxrange=100,chars=.7,divisions=4,pos=[0.92, 0.75, 0.95, 0.95],title='inverted pixels / total pixels on disk',/right,/vertical

xyouts, .3, .18, 'Plot Made '+!stime, chars=.7,/device

xyouts, 10200, 17520, 'Coverage Map for ' + ds, chars=.7, /device

xyouts, 21860, 12180, 'white = missing',chars=.65,/device

xyouts, 12000, 400, 'Hours', chars=1, /device 

xyouts, 1600, 8500, text_axes=0, orientation=90, 'Days of '+Monthname, /device

device,/close

set_plot,'x'

end
