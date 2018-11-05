
; Run calc_nn on all galaxies in cdfs
;Nancy edits to allow it to use randomly assigned ra and dec for calculation.
pro calc_4field_randompos, cdfs, eazy, Dn=Dn, Nr=Nr, nth=nth, rkpc=rkpc, $
                 dz=dz, da=da, verbose=verbose, minZred=minZred, $
                 SAM=SAM

; if SAM is set then use format of Henriques et al. 2014 SAM in the
; file 'stripped_MstarGt1e8...fits'
; 
; if using SAM then cdfs should be the h14 struct and nothing is
; needed for eazy

  if not keyword_set(minZred) then minZred=0.05

  if not keyword_set(SAM) then useZ = eazy.z_peak else useZ = cdfs.z_app

  if not keyword_set( da) then begin
     da = angular_diameter_distance(useZ, h=0.7, omega=0.3, lambda=0.7,/arcsec)
  endif

  if not keyword_set( nth) then nth=[2,3,5,7]
  if not keyword_set( rkpc) then rkpc=[200.,400.,500.,1000.]
  if not keyword_set(dz) then dz=0.02*2.5*(useZ+1)

  szN = size(/dim, nth)
  Dn = fltarr(szN, size(/dim, cdfs))
  help,Dn
  help,cdfs

  szR = size(/dim, rkpc)
  nR = intarr(szR, size(/dim, cdfs))

  if not keyword_set(SAM) then td = selStandard( cdfs) else $
     td = where( cdfs.rand_ra le 1./sqrt(2) and cdfs.rand_dec le 1./sqrt(2) and $
                 alog10(cdfs.stellarMass) + 10 ge 9.0)

  for i=0,n_elements(td)-1 do begin
     if keyword_set(verbose) then if i mod 100 eq 0 then print, 'working on ',i,' of ',n_elements(td)
   
     if useZ[td[i]] ge minZred then begin
        if size(/dim, dz) eq size(/dim, useZ) then $
           zd = where( abs( useZ[td[i]] - useZ[td]) le dz[td[i]])$
        else $
           zd = where( abs( useZ[td[i]] - useZ[td]) le dz)
        zd=td[zd]
        
        kpc2arcs = (da[td[i]])[0]
        if n_elements(zd) ge max(nth) then begin
           calc_nn, cdfs[zd].rand_ra, cdfs[zd].rand_dec, nth=nth, Dn=tDn, $
                    rarcs=rkpc/kpc2arcs, nR=tnR

     ;; find the galaxy:
           xx = (where( cdfs[zd].rand_ra eq cdfs[td[i]].rand_ra and cdfs[zd].rand_dec eq cdfs[td[i]].rand_dec))[0]
           nR[*,td[i]] = tnR[*,xx]
           Dn[*,td[i]] = tDn[*,xx]
        endif
     endif 
  endfor

end

;------------------------------------------------------------
;------------------------------------------------------------
;------------------------------------------------------------
; 
; pro calc_4field, cdfs, eazy, Dn=Dn, Nr=Nr, nth=nth, rkpc=rkpc, $
;                  dz=dz, da=da, verbose=verbose, minZred=minZred
; 
;   if not keyword_set(minZred) then minZred=0.05
; 
;   if not keyword_set( da) then begin
;      da = angular_diameter_distance(eazy.z_peak, h=0.7, omega=0.3, lambda=0.7,/arcsec)
;   endif
; 
;   if not keyword_set( nth) then nth=[2,3,5,7]
;   if not keyword_set( rkpc) then rkpc=[200.,400.,500.,1000.]
;   if not keyword_set(dz) then dz=0.02*2.5*(eazy.z_peak+1)
; 
;   szN = size(/dim, nth)
;   Dn = fltarr(szN, size(/dim, cdfs))
; 
;   szR = size(/dim, rkpc)
;   nR = intarr(szR, size(/dim, cdfs))
; 
;   td = selStandard( cdfs)
; 
;   for i=0,n_elements(td)-1 do begin
;      if keyword_set(verbose) then if i mod 100 eq 0 then print, 'working on ',i,' of ',n_elements(td)
;    
;      if eazy[td[i]].z_peak ge minZred then begin
;         if size(/dim, dz) eq size(/dim, eazy) then $
;            zd = where( abs( eazy[td[i]].z_peak - eazy[td].z_peak) le dz[td[i]])$
;         else $
;            zd = where( abs( eazy[td[i]].z_peak - eazy[td].z_peak) le dz)
;         zd=td[zd]
;         
;         kpc2arcs = (da[td[i]])[0]
;         if n_elements(zd) ge max(nth) then begin
;            calc_nn, cdfs[zd].rand_ra, cdfs[zd].rand_dec, nth=nth, Dn=tDn, $
;                     rarcs=rkpc/kpc2arcs, nR=tnR
; 
;      ;; find the galaxy:
;            xx = (where( cdfs[zd].rand_ra eq cdfs[td[i]].rand_ra and cdfs[zd].rand_dec eq cdfs[td[i]].rand_dec))[0]
;            nR[*,td[i]] = tnR[*,xx]
;            Dn[*,td[i]] = tDn[*,xx]
;         endif
;      endif
;   endfor
; 
; end

;------------------------------------------------------------
;------------------------------------------------------------
;------------------------------------------------------------
