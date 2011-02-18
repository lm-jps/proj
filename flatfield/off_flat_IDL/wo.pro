pro WO, reihe, nullen, sv=sv
;+
; NAME:
;       WO
; PURPOSE:
;       Finding invalid points in a time series
; EXPLANATION:
;       giving the first and last point of a sqequence of invalid points in a time series
;
;
; CALLING SEQUENCE:
;       WO, reihe, nullen [,sv=]
; INPUTS:
;       reihe: time series ((n) vector)
; OPTIONAL KEYWORS: sv: value of the invalid points (default criterium: finite(..)=0)
; OUTPUTS:
;       nullen: see explanaition ((2, m) vector)
; PROCEDURES USED:
;None
;
; MODIFICATION HISTORY:
;       Package 'SPM level 1.5' written between 2000 and 2003, Richard Wachter, PMOD/WRC
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 if n_elements(sv) eq 0 then w=where(finite(reihe) eq 0) else w=where(reihe eq sv)



 lw=(size(w))[1]

 seq=lonarr(2)
 nullen=[0L,0L]


 seq[0]=w[0]

 for k=1L, lw-1 do begin

  if ((w[k] - w[k-1]) gt 1) then begin

   seq[1]=w[k-1]
   nullen=[[nullen],[seq]]
   seq[0]=w[k]

  endif

 endfor

 seq[1]=w[lw-1]
 nullen=[[nullen],[seq]]
 nullen=nullen[*,1:*]

end



