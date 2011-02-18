function limb_darkening, n, radius, offset, a

order=5
rad=fltarr(n,n)

for i=0, n-1 do for j=0, n-1 do rad[i,j]=sqrt((i-n/2-offset[0])^2+(j-n/2-offset[1])^2)/radius
  
mu=sqrt(1.0-rad^2)
mu[where(finite(mu) eq 0)]=1.0


; Neckel & Labs 1994 (MDI wavelength)


ldark=poly(mu, a)

return, ldark

end

  

