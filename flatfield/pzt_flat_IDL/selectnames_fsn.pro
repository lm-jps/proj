pro selectnames_fsn,names,series, first,last,nfiles,keepdup=keepdup

; keepdup ne zero keeps duplicates

if (n_elements(keepdup) eq 0) then keepdup=0


dir='DS'

if (dir eq '') then begin
  spawn,'hostname',host
  res=stregex(host(0),'.edu',/fold_case,/boolean)
  if res then begin
    dir='DRMS'
  endif else begin
    spawn,'ls -1td /hmi0/egse/logs/200* | head -1',help
    dir=help(0)
  end
end

if ((dir eq 'DRMS') or (dir eq 'WD') or (dir eq 'DS') or (dir eq 'lev0d') or (dir eq 'lev0e') or (dir eq 'lev0d_test') or (dir eq 'lev0f') or (dir eq 'lev0') or (dir eq 'lev1a_nrt')) then begin
  if (dir eq 'DS') then begin
    ds=series
  end

;  first=0l
;  read,first,prompt='First FSN: '
;  last=0l
;  read,last,prompt='Last FSN: '
  if (dir eq 'DRMS') then begin
    spawn,'show_info ds=hmi_ground.lev0\['+strcompress(string(first)+'-'+string(last),/remove_all)+'\] key=FSN,EXTRA_INT_B,T_OBS seg=file -p -r -q',list
    spawn,'echo $$',pid
    openw,lun,'/tmp/show_'+pid,/append,/get_lun
    printf,lun,list
    printf,lun
    close,lun
    nlist=n_elements(list)
    rec=lonarr(nlist)
    fsn=lonarr(nlist)
    miss=lonarr(nlist)
    tobs=strarr(nlist)
    names=strarr(nlist)
    for i=0,nlist-1 do begin
      help=strsplit(list(i),/extract)
      rec(i)=help(0)
      fsn(i)=help(1)
      miss(i)=help(2)
      tobs(i)=help(3)
      names(i)=help(4)
    end
  endif else if (dir eq 'lev0d') then begin
    spawn,'show_info ds=hmi.lev0d\['+strcompress(string(first)+'-'+string(last),/remove_all)+'\] key=FSN,MISSVALS,T_OBS -p -r -q',list
    nlist=n_elements(list)
    rec=lonarr(nlist)
    fsn=lonarr(nlist)
    miss=lonarr(nlist)
    tobs=strarr(nlist)
    names=strarr(nlist)
    for i=0,nlist-1 do begin
      help=strsplit(list(i),/extract)
      rec(i)=help(0)
      fsn(i)=help(1)
      miss(i)=help(2)
      tobs(i)=help(3)
      names(i)=help(4)+'/image.fits'
    end
  endif else if (dir eq 'lev0e') then begin
    spawn,'show_info ds=hmi.lev0e\['+strcompress(string(first)+'-'+string(last),/remove_all)+'\] key=FSN,MISSVALS,T_OBS -p -r -q',list
    nlist=n_elements(list)
    rec=lonarr(nlist)
    fsn=lonarr(nlist)
    miss=lonarr(nlist)
    tobs=strarr(nlist)
    names=strarr(nlist)
    for i=0,nlist-1 do begin
      help=strsplit(list(i),/extract)
      rec(i)=help(0)
      fsn(i)=help(1)
      miss(i)=help(2)
      tobs(i)=help(3)
      names(i)=help(4)+'/image.fits'
    end
  endif else if (dir eq 'lev0f') then begin
    spawn,'show_info ds=hmi.lev0f\['+strcompress(string(first)+'-'+string(last),/remove_all)+'\] key=FSN,MISSVALS,T_OBS -p -r -q',list
    nlist=n_elements(list)
    rec=lonarr(nlist)
    fsn=lonarr(nlist)
    miss=lonarr(nlist)
    tobs=strarr(nlist)
    names=strarr(nlist)
    for i=0,nlist-1 do begin
      help=strsplit(list(i),/extract)
      rec(i)=help(0)
      fsn(i)=help(1)
      miss(i)=help(2)
      tobs(i)=help(3)
      names(i)=help(4)+'/image.fits'
    end
  endif else if (dir eq 'lev0') then begin
    spawn,'show_info ds=hmi.lev0\['+strcompress(string(first)+'-'+string(last),/remove_all)+'\] key=FSN,MISSVALS,T_OBS -p -r -q',list
    nlist=n_elements(list)
    rec=lonarr(nlist)
    fsn=lonarr(nlist)
    miss=lonarr(nlist)
    tobs=strarr(nlist)
    names=strarr(nlist)
    for i=0,nlist-1 do begin
      help=strsplit(list(i),/extract)
      rec(i)=help(0)
      fsn(i)=help(1)
      miss(i)=help(2)
      tobs(i)=help(3)
      names(i)=help(4)+'/image.fits'
    end
  endif else if (dir eq 'lev1a_nrt') then begin
    spawn,'show_info ds=hmi.lev1a_nrt\['+strcompress(string(first)+'-'+string(last),/remove_all)+'\] key=FSN,MISSVALS,T_OBS -p -r -q',list
    nlist=n_elements(list)
    rec=lonarr(nlist)
    fsn=lonarr(nlist)
    miss=lonarr(nlist)
    tobs=strarr(nlist)
    names=strarr(nlist)
    for i=0,nlist-1 do begin
      help=strsplit(list(i),/extract)
      rec(i)=help(0)
      fsn(i)=help(1)
      miss(i)=help(2)
      tobs(i)=help(3)
      names(i)=help(4)+'/image.fits'
    end
  endif else if (dir eq 'lev0d_test') then begin
    spawn,'show_info ds=su_production.lev0d_test\['+strcompress(string(first)+'-'+string(last),/remove_all)+'\] key=FSN,MISSVALS,T_OBS -p -r -q',list
    nlist=n_elements(list)
    rec=lonarr(nlist)
    fsn=lonarr(nlist)
    miss=lonarr(nlist)
    tobs=strarr(nlist)
    names=strarr(nlist)
    for i=0,nlist-1 do begin
      help=strsplit(list(i),/extract)
      rec(i)=help(0)
      fsn(i)=help(1)
      miss(i)=help(2)
      tobs(i)=help(3)
      names(i)=help(4)+'/image.fits'
    end
  endif else if (dir eq 'DS') then begin
    spawn,'show_info ds='+ds+'\['+strcompress(string(first)+'-'+string(last),/remove_all)+'\] key=FSN,MISSVALS,T_OBS -p -r -q',list
    nlist=n_elements(list)
    rec=lonarr(nlist)
    fsn=lonarr(nlist)
    miss=lonarr(nlist)
    tobs=strarr(nlist)
    names=strarr(nlist)
    for i=0,nlist-1 do begin
      help=strsplit(list(i),/extract)
      rec(i)=help(0)
      fsn(i)=help(1)
      miss(i)=help(2)
      tobs(i)=help(3)
      names(i)=help(4)+'/image.fits'
    end
  end else begin ; WD
    wdname='qqq'
    read,wdname,prompt='File to read from: '
    if (wdname eq '') then wdname='/hmi0/egse/logs/wdlog'
;   spawn,'wc -l ' +wdname,help
    spawn,'grep -v # '+wdname+' | wc -l ',help
    nlist0=long(help[0])
    rec0=lonarr(nlist0)
    fsn0=lonarr(nlist0)
    names0=strarr(nlist0)
    fname='qqq'
    openr,1,wdname
    for i=0,nlist0-1 do begin
      repeat begin
        readf,1,fname
      endrep until (stregex(fname,'#') eq -1)
      rec0(i)=i
      fsplit=strsplit(fname,'/.',/extract)
      fsn0(i)=long((fsplit)[n_elements(fsplit)-2])
      names0(i)=fname
    end
    close,1
    w=where((fsn0 ge first) and (fsn0 le last),nlist)
    sw=w(sort(fsn0(w)))
    rec=rec0(sw)
    fsn=fsn0(sw)
    miss=lonarr(nlist)
    tobs=strarr(nlist)
    names=names0(sw)
    fname='qqq'
    for i=0,nlist-1 do begin
      fname=names(i)
      q=readfits(fname,hdr,nslice=0,/silent)
      miss(i)=sxpar(hdr,'MISSVALS')
      tobs(i)=sxpar(hdr,'DATE')
    end
  end
  hfsn=histogram(fsn,min=0)
  w=where(hfsn ne 0,nfiles)
  if ((nfiles ne nlist) and (keepdup eq 0)) then begin
    print,'Duplicates found'
    ix=lonarr(nfiles)
    for i=0,nfiles-1 do begin
      ii=w(i)
      w1=where(fsn eq ii,n1)
      if (n1 eq 1) then begin ; Only one record number with this FSN
        ix(i)=w1
      endif else begin
;       w0=w1(where(miss(w1) eq 0,n0))
        ww=where(miss(w1) eq 0,n0)
        if (n0 ne 0) then w0=w1(ww)
        if (n0 eq 0) then begin ; No good ones.
          print,'All bad. Taking latest for FSN',ii
          ix(i)=where(rec eq max(rec(w1))) ; Take latest file
        end
        if (n0 eq 1) then begin ; Only one good one.
          print,'Taking only complete image for FSN',ii
          ix(i)=w1(where(miss(w1) eq 0))
        end
        if (n0 gt 1) then begin ; Multiple good ones.
          print,'Several good images. Taking latest for FSN',ii
          ix(i)=w0(where(rec eq max(rec(w0)))) ; Take latest file
        end
      end
    end
    names=names(ix)
    miss=miss(ix)
  end
  nfiles=nlist
  print,'Found',nfiles,' files'
end else begin ; regular directory
  dir=dir+'/'
  spawn,'ls -1 '+dir,list0
  repeat begin
    name='qqq'
    read,name,prompt='First file: '
    if (name eq '') then name='i_.*\.fit'
    whit=where(stregex(list0,name,/boolean) ne 0,nhit)
    if (nhit eq 0) then begin
      print,'No matches'
    endif else begin
      list=list0(whit)
      if (nhit eq 1) then print,'Found file: '+list(0)
      if (nhit ge 2) then begin
        print,'Found multiple files:'
        print,list
      end
    end
  endrep until (nhit eq 1)
  ifirst=whit(0)
  
  repeat begin
    name='qqq'
    read,name,prompt='Last file: '
    if (name eq '') then name='i_.*\.fit'
    whit=where(stregex(list0,name,/boolean) ne 0,nhit)
    if (nhit eq 0) then begin
      print,'No matches'
    endif else begin
      list=list0(whit)
      if (nhit eq 1) then print,'Found file: '+list(0)
      if (nhit ge 2) then begin
        print,'Found multiple files:'
        print,list
      end
    end
  endrep until (nhit eq 1)
  ilast=whit(0)
  
  names=dir+list0(ifirst:ilast)
  nfiles=ilast-ifirst+1
  print,'Found',nfiles,' files'
end

end

