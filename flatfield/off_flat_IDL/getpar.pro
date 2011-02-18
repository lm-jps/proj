function getpar,headers,name

; returns the values of the fits keyword given by name from the headers

nim=n_elements(headers)

; Get correct type from first header
first=sxpar(*headers(0),name)

values=make_array(nim,value=first)

for i=0,nim-1 do values(i)=sxpar(*headers(i),name)

return,values

end

