;require  calc_4field_randompos.pro
pro calc_4field_zf_randompos
;Read all catalogs
;cdfs = mrdfits('../cdfs.v1.5.5.fits',1)
;cos = mrdfits('../cosmos.v1.2.9.fits',1)
;uds = mrdfits('../uds.v1.4.3.fits',1)

restore,'zfourge2016_ksselectedcat.sav'

;(Uniform) Randomly assign ra and dec to each galaxy
cdfs_randra =min(cdfs.ra)+ (max(cdfs.ra)-min(cdfs.ra))*randomu(seed,n_elements(cdfs.ra))
cdfs_randdec =min(cdfs.dec)+ (max(cdfs.dec)-min(cdfs.dec))*randomu(seed,n_elements(cdfs.dec))
cos_randra =min(cos.ra)+ (max(cos.ra)-min(cos.ra))*randomu(seed,n_elements(cos.ra))
cos_randdec =min(cos.dec)+ (max(cos.dec)-min(cos.dec))*randomu(seed,n_elements(cos.dec))
uds_randra =min(uds.ra)+ (max(uds.ra)-min(uds.ra))*randomu(seed,n_elements(uds.ra))
uds_randdec =min(uds.dec)+ (max(uds.dec)-min(uds.dec))*randomu(seed,n_elements(uds.dec))
;add to the existing structures
struct_add_field,cdfs,'rand_ra',cdfs_randra
struct_add_field,cdfs,'rand_dec',cdfs_randdec
struct_add_field,cos,'rand_ra',cos_randra
struct_add_field,cos,'rand_dec',cos_randdec
struct_add_field,uds,'rand_ra',uds_randra
struct_add_field,uds,'rand_dec',uds_randdec

;td = selStandard(cdfs)
;Read EAZY catalogs
;cdfseazy = mrdfits('../cdfs.v1.5.3.zout.fits',1)
;coseazy=mrdfits('../cosmos.v1.2.7.zout.fits',1)
;udseazy=mrdfits('../uds.v1.4.1.zout.fits',1)  
if 1 then begin
   td=selstandard(cdfs)
   ttd=where( cdfseazy[td].z_peak ge 0.3)
   tda_cdfs = angular_diameter_distance(cdfseazy[td[ttd]].z_peak, h=0.7, omega=0.3, lambda=0.7,/arcsec)
   da_cdfs = replicate(0.0, n_elements(cdfs))
   da_cdfs[td[ttd]] = tda_cdfs
   
   tic
   calc_4field_randompos, cdfs, cdfseazy, Dn=cdfs_dn, nR=cdfs_nR,nth=[1,2,3,4,5,6,7],dz=0.02*2.5*(1+cdfseazy.z_peak), da=da_cdfs, /verb, minZred=0.3
   toc
endif

;------------------------------------------------------------

tu=selstandard(uds)
ttu=where( udseazy[tu].z_peak ge 0.3)
tda_uds = angular_diameter_distance(udseazy[tu[ttu]].z_peak, h=0.7, omega=0.3, lambda=0.7,/arcsec)
da_uds = replicate(0.0, n_elements(uds))
da_uds[tu[ttu]] = tda_uds

tic
calc_4field_randompos, uds, udseazy, Dn=uds_dn, nR=uds_nR,nth=[1,2,3,4,5,6,7], dz=0.02*2.5*(1+udseazy.z_peak), da=da_uds, /verb, minZred=0.3
toc

;------------------------------------------------------------


to=selstandard(cos)
tto=where( coseazy[to].z_peak ge 0.3)
tda_cos = angular_diameter_distance(coseazy[to[tto]].z_peak, h=0.7, omega=0.3, lambda=0.7,/arcsec)
da_cos = replicate(0.0, n_elements(cos))
da_cos[to[tto]] = tda_cos

tic
calc_4field_randompos, cos, coseazy, Dn=cos_dn, nR=cos_nR,nth=[1,2,3,4,5,6,7],dz=0.02*2.5*(1+coseazy.z_peak), da=da_cos, /verb, minZred=0.3
toc

;------------------------------------------------------------

nn = [1,2,3,4,5,6,7]
Rkpc = [200, 400, 500,1000]
save,file='calc_4field_randompos_091216.sav', cos_dn, cos_NR, da_cos, cdfs_dn, cdfs_NR, da_cdfs, uds_dn, uds_nr, da_uds, nn, rkpc,cdfs_randra,cdfs_randdec,cos_randra,cos_randdec,uds_randra,uds_randdec
;save,file='krap.sav', cos_dn, cos_NR, da_cos, cdfs_dn, cdfs_NR, da_cdfs, uds_dn, uds_nr, da_uds, nn, rkpc

; 2015/8/9:  I'm pretty sure NN are measured in arcsecs.  That's what
; calc_nn does and no were do I scale it...

;------------------------------------------------------------

; Do some check stuff -- make maps at z=1.6 for CDFS and UDS and see
;                        if they agree with older maps you have, and
;                        for clusters, etc. 

if 0 then begin

to=selstandard(cos)
tto=where( coseazy[to].z_peak ge 1.95 and coseazy[to].z_peak le 2.15 and cosFast[to].lmass ge 9.0)
tto=to[tto]
;write_reg, cos[tto].ra, cos[tto].dec, 1.0, 'cos_z1p9to2.2_density.reg', text='N2,3,5='+strtrim(string(format='(f4.1)',cos_dn[0,tto]),2)+','+strtrim(string(format='(f4.1)',cos_dn[1,tto]),2)+','+strtrim(string(format='(f4.1)',cos_dn[2,tto]),2)

irkpc=1 ; 
inn=1

write_reg, cos[tto].ra, cos[tto].dec, 1.0, 'cos_z1p95to2.15_N'+strtrim(string(format='(i1)',nn[inn]),2)+'.reg', text='z='+string(format='(f4.2)',cosEazy[tto].z_peak)+' N'+strtrim(string(format='(i1)',nn[inn]),2)+'='+strtrim(string(format='(f4.1)',cos_dn[inn,tto]),2)
write_reg, cos[tto].ra, cos[tto].dec, 5.0, 'cos_z1p95to2.15_R'+strtrim(string(format='(i4)',rkpc[irkpc]),2)+'.reg', text='R'+strtrim(string(format='(i4)',rkpc[irkpc]),2)+'='+strtrim(string(format='(f4.1)',cos_nr[irkpc,tto]),2), color='cyan'

plot, cos_dn[inn,tto], cos_nr[irkpc,tto], psym=4, xr=[0,100],/xst, xtit='D!b'+string(format='(i1)',nn[inn])+'!n [arcsec?]', ytit='N(R<'+string(format='(i4)',rkpc[irkpc])+')'

medNR=findgen(max(cos_nr[irkpc,tto])) 
for i=0,n_elements(medNR)-1 do medNR[i]=median( cos_dn[inn,tto[where(/null,fix(cos_nr[irkpc,tto]) eq i)]])

iarr = findgen(n_elements(medNR))
oplot, medNR, iarr, psym=4, color=cjp_icolor('red'), symsize=3

endif



end  
