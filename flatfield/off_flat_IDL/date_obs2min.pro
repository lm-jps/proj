function date_obs2min, dobs, juld=juld, T_OBS=T_OBS

year=long(strmid(dobs, 0, 4))
month=long(strmid(dobs, 5,2))
day=long(strmid(dobs, 8,2))

hour=long(strmid(dobs, 11, 2))
minute=long(strmid(dobs, 14,2))
if not keyword_set(T_OBS) then sec=float(strmid(dobs, 17, 6)) else sec=float(strmid(dobs, 17, 2))

if keyword_set(juld) then begin & offs=0 & fac=1.0 & taic=0.0 & endif else begin & offs=julday(12,1,1995) & fac=24.*60. & taic=0.5 & endelse

return, (julday(month, day, year, hour, minute, sec)- offs)*fac - taic

end


