pro readimages,names,images,headers,nbin=nbin,silent=silent,noshow=noshow,raw=raw,getnames=getnames,iover=iover,pad=pad,fixcenter=fixcenter

; names contains the list of images read or to be read
; images contain the images
; headers contain the headers
; nbin gives the desired image size. Must divide into the original image
; size (4200 for raw images, 4096 for ones with the cross removed. No binning
; is done if not set.
; /silent makes the program shut up during reading. Otherwise the image
; number (counting from zero), the time part of the image name, the FSN,
; the exposure time, the focus, wavelength and polarization table indices
; and the mean dark value is printed.
; /noshow makes the program not show thumbnails as it goes
; /raw makes the program return raw 4200^2 images for RAL
; /raw makes the program return raw 4096^2 images for CIF
; /getnames makes the program query for the image name instead of using
; those in names
; iover gives the intensity in the overscan pixels
; pad will pad images if CIF and overscan (imcfg>2)
; fixcenter will interpolate center pixel if CIF and imcfg=2


if (n_elements(silent) eq 0) then silent=0
if (n_elements(noshow) eq 0) then noshow=0
if (n_elements(raw) eq 0) then raw=0
if (n_elements(getnames) eq 0) then getnames=0
if (n_elements(pad) eq 0) then pad=1
if (n_elements(fixcenter) eq 0) then fixcenter=1

if (getnames ne 0) then begin
  selectnames,names,nim
end

if (silent eq 0) then spawn,'date'

cfglist=intarr(300)-1
cfglist([2,5,8,11])=0 ; No overscan
cfglist([3,6,9,12])=1 ; One overscan
cfglist([4,7,10,13])=2 ; Two overscan
cfglist(80)=1
cfglist(81)=2
cfglist(84)=0
cfglist(91:107)=1
cfglist(160:184)=1
cfglist(190:212)=1
cfglist(116:118)=2
cfglisth=cfglist
cfglistv=cfglist
cfglisth(123)=1
cfglistv(123)=0
cfglisth(86)=1
cfglistv(86)=0
cfglisth(89)=0
cfglistv(89)=1
cfglisth(121)=0
cfglistv(121)=0
cfglist=0

; Use first image to find type
im=readfits(names(0),header,/silent)
naxis=sxpar(header,'NAXIS')
if (naxis eq 0) then begin
  im=fitsio_read_image(names(0),header)
  nx=sxpar(header,'ZNAXIS1')
  ny=sxpar(header,'ZNAXIS2')
endif else begin
  nx=sxpar(header,'NAXIS1')
  ny=sxpar(header,'NAXIS2')
end
config=sxpar(header,'CONFIG')
sz=size(config,/type)
type=-1
if (sz eq 7) then begin ; Keyword exists and is of type string
  if (config eq 'RAL') then type=0
  if (config eq 'CIF') then type=1
endif else begin ; Otherwise guess based on size
  if (nx eq 4200) then type=0
  if (nx eq 4096) then type=1
end
if (type lt 0) then begin
  print,'Unknown image type'
  stop
end
print,type

nim=n_elements(names)
headers=ptrarr(nim)
iover=fltarr(nim)

if (type eq 0) then begin ; RAL
  nn=4096
  if (raw ne 0) then nn=4200
  if (n_elements(nbin) eq 0) then nbin=nn
  nsmall=64
  if (raw ne 0) then nsmall=60
  
  if (noshow eq 0) then nshow=(!d.x_size/nsmall)*(!d.y_size/nsmall)
  
  images=fltarr(nbin,nbin,nim)
  ix=[2098-2047+indgen(2048),2101+indgen(2048)]
  iy=[2067-2047+indgen(2048),2132+indgen(2048)]
  iy0=2068+indgen(64)
  
  if (silent eq 0) then print,'   #   Name    FSN   EXP CAL  WL POL   Dark     Over'
  for i=0,nim-1 do begin
    name=names(i)
    im=readfits(name,header,/silent)
    q=im(ix,*)
    io=rebin(q(*,iy0)+0.0,1)
    if (raw eq 0) then begin
      if (nn eq nbin) then images(*,*,i)=q(*,iy)+0.0 else images(*,*,i)=rebin(q(*,iy)+0.0,nbin,nbin)
    endif else begin
      if (nn eq nbin) then images(*,*,i)=im+0.0 else images(*,*,i)=rebin(im+0.0,nbin,nbin)
    end
    ln=strlen(names(0))
    sname=strmid(name,ln-10,6)
    iover(i)=io
    if (silent eq 0) then print,i,sname,sxpar(header,'HSQFGSN'),sxpar(header,'HSHIEXP'),sxpar(header,'HCFTID'),sxpar(header,'HWLTID'),sxpar(header,'HPLTID'),rebin(q(0:63,iy(0:63))+0.0,1),io,format='(i4,a8,i6,i6,3i4,2f9.2)'
    headers(i)=ptr_new(header)
    if (noshow eq 0) then tvscl,rebin(images(*,*,i),nsmall,nsmall),i mod nshow
  end
endif else begin ; CIF
  nn=4096
  if (n_elements(nbin) eq 0) then nbin=nn
  nsmall=64
  
  if (noshow eq 0) then nshow=(!d.x_size/nsmall)*(!d.y_size/nsmall)
  
  images=fltarr(nbin,nbin,nim)
  ix=[2098-2047+indgen(2048),2101+indgen(2048)]
  iy=[2067-2047+indgen(2048),2132+indgen(2048)]
  iy0=2068+indgen(64)
  
  if (silent eq 0) then print,'   #   FSN     SEC   EXP CAL  WL POL   Dark     Over'
  for i=0,nim-1 do begin
    name=names(i)
;   im=readfits(name,header,/silent)
    im=fitsio_read_image(name,header)
    if (sxpar(header,'HCFTID') eq 0) then begin ; DRMS
      namesm=strmid(name,0,strlen(name)-5)+'_sm.fits'
;     q=readfits(namesm,header,/silent)
      q=fitsio_read_image(namesm,header)
      sxaddpar,header,'NAXIS1',4096
      sxaddpar,header,'NAXIS2',4096
    end
    if (fixcenter ne 0) then im(2047,2048)=im(2047,2047)
    imcfg=sxpar(header,'H0149')
    if (imcfg eq 0) then imcfg=sxpar(header,'HIMGCFID')
;   if ((not raw) and ((imcfg lt 2) or (imcfg gt 13))) then begin
    if ((raw ne 0) and (imcfg lt 0)) then imcfg=0 ; Bad imcfg
    if (imcfg lt 0) then imcfg=0 ; Bad imcfg
    if ((cfglistv(imcfg) lt 0) or (raw ne 0)) then begin
      if (raw) then begin
        im=im and 32767
        images(*,*,i)=im+0.0
        io=0.0
      endif else begin
        im=fltarr(nn,nn)-1e6
        io=-1e6
        print,'WARNING: BAD IMCFG! (',imcfg,')'
      end
    endif else begin
      im=im and 32767
      q=im
      if (cfglisth(imcfg) eq 0) then begin ; No overscan
        iox=im
      end
      if (cfglisth(imcfg) eq 1) then begin ; 1 overscan
        iox=[im(0:2046,*),im(2049:4095,*)]
        q(1:2047,*)=im(0:2046,*)
        q(2048:4094,*)=im(2049:4095,*)
        if (pad eq 0) then begin
          q(0,*)=-32767
          q(4095,*)=-32767
        end
      end
      if (cfglisth(imcfg) eq 2) then begin ; 2 overscan
        iox=[im(0:2045,*),im(2050:4095,*)]
        q(2:2047,*)=im(0:2045,*)
        q(2048:4093,*)=im(2050:4095,*)
        if (pad eq 0) then begin
          q(0:1,*)=-32767
          q(4094:4095,*)=-32767
        end
      end
      im=q
      if (cfglistv(imcfg) eq 0) then begin
        io=0.0
      end
      if (cfglistv(imcfg) eq 1) then begin ; 1 overscan
        io=mean(iox(*,2047:2048))
        im(*,1:2047)=q(*,0:2046)
        im(*,2048:4094)=q(*,2049:4095)
        if (pad eq 0) then begin
          im(*,0)=-32767
          im(*,4095)=-32767
        end
      end
      if (cfglistv(imcfg) eq 2) then begin ; 2 overscan
        io=mean(iox(*,2046:2049))
        im(*,2:2047)=q(*,0:2045)
        im(*,2048:4093)=q(*,2050:4095)
        if (pad eq 0) then begin
          im(*,0:1)=-32767
          im(*,4094:4095)=-32767
        end
      end
      if (nn eq nbin) then begin
        images(*,*,i)=im+0.0
      endif else begin
        images(*,*,i)=rebin(im+0.0,nbin,nbin)
      end
    end
    iover(i)=io
    ln=strlen(names(0))
    sname=strmid(name,ln-10,6)
    if (silent eq 0) then print,i,sxpar(header,'HSQFGSN'),sxpar(header,'SHS') mod 86400,sxpar(header,'HSHIEXP'),sxpar(header,'HCFTID'),sxpar(header,'HWLTID'),sxpar(header,'HPLTID'),rebin(im(0:63,iy(0:63))+0.0,1),io,format='(i4,i8,i6,i6,3i4,2f9.2)'
    headers(i)=ptr_new(header)
    if (noshow eq 0) then tvscl,rebin(images(*,*,i),nsmall,nsmall),i mod nshow
  end
end

if (silent eq 0) then spawn,'date'

end

