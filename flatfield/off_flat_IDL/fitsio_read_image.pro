;+
; NAME:
;       FITSIO_READ_IMAGE
;
; PURPOSE:
;       Read FITS image and header using external FITSIO library routines
;
; CALLING SEQUENCE:
;       img = FITSIO_READ_IMAGE( filename [, hd] [, /double | /single] $
;                                [, /chksum] )
;
; INPUTS:
;       filename: string containing name of FITS file to read.  By default
;                 the FIRST image in the FITS file will be read.  To read
;                 an image after the first one, one must give the filename
;                 in the FITSIO notation, e.g., 'foo.fits[2]' for the second
;                 extension, or the third HDU, of foo.fits.
;       /single:  force output image to be of type single-precision float
;       /double:  force output image to be of type double-precision float
;       /chksum:  verify FITS checksum
;
; OUTPUTS:
;       img:      image data array
;       hd:       array of strings containing FITS header "cards" which
;                 can be parsed with SXPAR and similar procedures (optional)
;
; MODIFICATION HISTORY:
;       2008.09.08	Keh-Cheng Chu
;	  Initial release
;	2010.05.15	Keh-Cheng Chu
;	  Make second parameter optional
;	  Added checksum verification in external routine
;	2010.06.01	Keh-Cheng Chu
;	  Added /single and /double keywords 
;	2010.06.11	Keh-Cheng Chu
;         Make checksum verification selectable by keyword
;	2010.12.17	Keh-Cheng Chu
;         Instead of calling fits_get_img_equivtype() to determine
;         output data type, always choose float (double) type for 
;         BITPIX = 8 or 16 (32 or 64) when there is non-default
;         BSCALE and BZERO.
;
;-

function FITSIO_READ_IMAGE, filename, hd, single=single, double=double, $
                            chksum=chksum

compile_opt IDL2, logical_predicate, strictarrsubs
on_error, 2

if (n_params() gt 2 or n_params() lt 1) then message, $
    'Usage: img = FITSIO_READ_IMAGE( filename [, hd] [, /single | /double] )'

forcesingle =  keyword_set(single)
forcedouble =  keyword_set(double)
if (forcesingle and forcedouble) then message, $
    'Error: both /single and /double specified!'

;LIB = '/home/kehcheng/idl/fitsio/fitsio.so'
LIB = '/home/jsoc/idl/fitsio.so'

extnum = 999l
dtype = 999l
naxis = 999l
naxes = lon64arr(8)
do_chksum = keyword_set(chksum)
errmsg = call_external(LIB,'get_info',filename,extnum,dtype,naxis,naxes, $
                       do_chksum, /s_value, value=[1,0,0,0,0,0])
if (strlen(errmsg) gt 0) then begin
    print,errmsg
    return, -1
endif

if (n_params() eq 2) then begin
    s = strpos(filename, '[')
    if (s eq -1) then begin
	hd = headfits(filename, exten=extnum) 
    endif else begin
	fn = strmid(filename, 0, s)
	hd = headfits(fn, exten=extnum)
    endelse
endif

naxes = naxes[0:naxis-1]
if (forcesingle) then dtype = 4
if (forcedouble) then dtype = 5
img = make_array(dimension=naxes, type=dtype)
errmsg = call_external(LIB,'get_data',filename,dtype,img, /s_value, value=[1,0,0])
if (strlen(errmsg) gt 0) then begin
    print,errmsg
    return, -1
endif

return, img
end
