
; 
pro calc_nn, ra, dec, nth=nth, rarcs=rarcs, dn=dn, nR=nR

; given a list of object ra and dec, calculate the nth nearest
; neighbor
;
; you set nth= any number (3, 5, 7) for the 3rd, 5th, 7th nearest
; neighbor.  nth = 3 by default
; 
; you can also make nth=[3,5,7] 
;
; returns Dn = fltarr(# of elements of nth, # of elements of ra) which
; is the distance to the nth nearest neighbor in arcsec
  
; rarcs = radius (in arcs) of an aperture in which to count neighbors.
; Could be rarcs = 30 (arcs) or rkpc = [20,30,40].
; 
; returns nInR = fltarr(# of of elements in rarcs, # of elements 
; 

if not keyword_set(nth) then nth=3

if size(/dim, ra) ne size(/dim, dec) then begin
   print,'% CALC_NN: error, ra and dec arrays have different sizes'
   return
endif

szN = size(/dim, nth)
Dn = fltarr(szN, size(/dim, ra))

if keyword_set(rarcs) then begin
   szR = size(/dim, rarcs)
   
   nR = intarr(szR, size(/dim, ra))
endif

for i=0,n_elements(ra)-1 do begin

   ;; compute distance to every galaxy and sort: 
   gcirc, 2, ra[i], dec[i], ra, dec, dis
   s = sort(dis)
   ;; take the Nth one (skip the first, which is the galaxy itself)
   Dn[*,i] = dis[s[nth]]

   if keyword_set(rarcs) then begin
      for j=0,n_elements(rarcs)-1 do begin
         xxx=where( dis lt rarcs[j], count)
         nR[j,i] = count-1 ;; don't count the galaxy itself
      endfor
   endif

endfor



end
