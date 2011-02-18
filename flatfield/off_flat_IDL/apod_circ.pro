function apod_circ, nn , rad, nb, offx=offx, offy=offy

pat=fltarr(nn, nn)


rarr=fltarr(nn, nn)
if not keyword_set(offx) then offx=0.0
if not keyword_set(offy) then offy=0.0

for i=0.0, nn-1 do for j=0.0, nn-1 do rarr[i,j]=sqrt((i-(nn/2+offx))^2+(j-(nn/2+offy))^2)


pat[where(rarr lt rad)]=1.0

ixm=where(rarr ge rad and rarr lt rad+nb)
if ixm[0] ne -1 then pat[ixm]=0.5*cos(!pi/nb*(rarr[where(rarr ge rad and rarr lt rad+nb)]-rad))+0.5

ixo=where(rarr ge rad+nb)
if ixo[0] ne -1 then pat[ixo]=0.0

return, pat

end
