function getshs,headers

time=getpar(headers, 'SHS')+getpar(headers, 'SHSS')/2D^32

if (total(time) eq 0) then begin
  nh=n_elements(headers)
  time=dblarr(nh)
  fpt=getpar(headers,'IMGFPT')
  for i=0,nh-1 do begin
    spawn,'time_index -r in='+fpt(i),help
    time(i)=double(help[0])
  end
end

return,time

end
