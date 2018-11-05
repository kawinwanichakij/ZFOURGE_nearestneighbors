pro convertcat



 
;read the public version of zf catalog and create structure
  rcat,'cosmos.v1.3.6.cat',/unpack
  cos = replicate({catalog,id:0L,ra:0.0,dec:0.0,wmin_optical:0.0,wmin_fs:0.0,use:0.0,star:0.0,f_F160W:0.0,f_F125W:0.0},n_elements(id))
  cos.id = reform(id)
  cos.ra=  reform(ra)
  cos.dec=  reform(dec)
  cos.wmin_optical=  reform(wmin_optical)
  cos.wmin_fs=  reform(wmin_fs)
  cos.use=  reform(use)
  cos.star=  reform(star)
  cos.f_F160W = reform(f_F160W)
  cos.f_F125W = reform(f_F125W)
  rcat,'cosmos.v1.3.6.zout',/unpack
  rcat,'cosmos.v1.3.6.rest.v0.9.cat',/unpack
  Umag = 23.9 - 2.5*alog10(reform(f_U))
  Vmag = 23.9 - 2.5*alog10(reform(f_V))
  Jmag = 23.9 - 2.5*alog10(reform(f_J))
  umv = Umag - Vmag
  vmj = Vmag - Jmag
  coseazy=replicate({eazy,id:0L,z_peak:0.0,umv:0.0,vmj:0.0},n_elements(id))
  coseazy.id = reform(id)
  coseazy.z_peak = reform(z_peak)
  coseazy.umv = umv
  coseazy.vmj = vmj
  ;rcat,'cosmos.v1.3.6.fout',/unpack

  rcat,'cdfs.v1.6.9.cat',/unpack
  cdfs = replicate({catalog,id:0L,ra:0.0,dec:0.0,wmin_optical:0.0,wmin_fs:0.0,use:0.0,star:0.0,f_F160W:0.0,f_F125W:0.0},n_elements(id))
  cdfs.id = reform(id)
  cdfs.ra=  reform(ra)
  cdfs.dec=  reform(dec)
  cdfs.wmin_optical=  reform(wmin_optical)
  cdfs.wmin_fs=  reform(wmin_fs)
  cdfs.use=  reform(use)
  cdfs.star=  reform(star)
  cdfs.f_F160W = reform(f_F160W)
  cdfs.f_F125W = reform(f_F125W)
   
  rcat,'cdfs.v1.6.9.zout',/unpack
  rcat,'cdfs.v1.6.9.rest.v0.9.cat',/unpack
  Umag = 23.9 - 2.5*alog10(reform(f_U))
  Vmag = 23.9 - 2.5*alog10(reform(f_V))
  Jmag = 23.9 - 2.5*alog10(reform(f_J))
  umv = Umag - Vmag
  vmj = Vmag - Jmag
  cdfseazy=replicate({eazy,id:0L,z_peak:0.0,umv:0.0,vmj:0.0},n_elements(id))
  cdfseazy.id = reform(id)
  cdfseazy.z_peak = reform(z_peak)
  cdfseazy.umv = umv
  cdfseazy.vmj = vmj

   rcat,'uds.v1.5.8.cat',/unpack
  uds = replicate({catalog,id:0L,ra:0.0,dec:0.0,wmin_optical:0.0,wmin_fs:0.0,use:0.0,star:0.0,f_F160W:0.0,f_F125W:0.0},n_elements(id))
  uds.id = reform(id)
  uds.ra=  reform(ra)
  uds.dec=  reform(dec)
  uds.wmin_optical=  reform(wmin_optical)
  uds.wmin_fs=  reform(wmin_fs)
  uds.use=  reform(use)
  uds.star=  reform(star)
  uds.f_F160W = reform(f_F160W)
  uds.f_F125W = reform(f_F125W)
  rcat,'uds.v1.5.8.zout',/unpack
  rcat,'uds.v1.5.8.rest.v0.9.cat',/unpack
  Umag = 23.9 - 2.5*alog10(reform(f_U))
  Vmag = 23.9 - 2.5*alog10(reform(f_V))
  Jmag = 23.9 - 2.5*alog10(reform(f_J))
  umv = Umag - Vmag
  vmj = Vmag - Jmag
  udseazy=replicate({eazy,id:0L,z_peak:0.0,umv:0.0,vmj:0.0},n_elements(id))
  udseazy.id = reform(id)
  udseazy.z_peak = reform(z_peak)
  udseazy.umv = umv
  udseazy.vmj = vmj


    rcat,'cdfs.v1.6.9.fout',/unpack
  cdfsfast=replicate({fast,id:0L,lmass:0.0,lage:0.0},n_elements(id))
  cdfsfast.id=reform(id)
  cdfsfast.lmass=reform(lmass)
  cdfsfast.lage =reform(lage)
  rcat,'cosmos.v1.3.6.fout',/unpack
  cosfast=replicate({fast,id:0L,lmass:0.0,lage:0.0},n_elements(id))
  cosfast.id=reform(id)
  cosfast.lmass=reform(lmass)
  cosfast.lage =reform(lage)
  rcat,'uds.v1.5.8.fout',/unpack
  udsfast=replicate({fast,id:0L,lmass:0.0,lage:0.0},n_elements(id))
  udsfast.id=reform(id)
  udsfast.lmass=reform(lmass)
  udsfast.lage = reform(lage)
  readcol,'cdfs.v1.6.9.vdw.v0.4.cat',id,wmin,r,f_r,NUMBER,ra,dec,re,dre,n,dn,q,dq,mag,dmag,pa,dpa,sn,f
  cdfsgalfit = replicate({galfit,id:0L,re:0.0,n:0.0,q:0.0,pa:0.0,mag:0.0},n_elements(id))
  cdfsgalfit.id = id
  cdfsgalfit.re = re
  cdfsgalfit.n = n
  cdfsgalfit.q = q
  cdfsgalfit.pa = pa
  cdfsgalfit.mag=mag
  
  readcol,'cosmos.v1.3.6.vdw.v0.4.cat',id,wmin,r,f_r,NUMBER,ra,dec,re,dre,n,dn,q,dq,mag,dmag,pa,dpa,sn,f
  cosgalfit = replicate({galfit,id:0L,re:0.0,n:0.0,q:0.0,pa:0.0,mag:0.0},n_elements(id))
  cosgalfit.id = id
  cosgalfit.re = re
  cosgalfit.n = n
  cosgalfit.q = q
  cosgalfit.pa =pa
  cosgalfit.mag= mag
  readcol,'uds.v1.5.8.vdw.v0.4.cat',id,wmin,r,f_r,NUMBER,ra,dec,re,dre,n,dn,q,dq,mag,dmag,pa,dpa,sn,f
  udsgalfit = replicate({galfit,id:0L,re:0.0,n:0.0,q:0.0,pa:0.0,mag:0.0},n_elements(id))
  udsgalfit.id = id
  udsgalfit.re = re
  udsgalfit.n = n
  udsgalfit.q = q
  udsgalfit.pa = pa
  udsgalfit.mag = mag
  save,cos,cdfs,uds,coseazy,cdfseazy,udseazy,cdfsfast,cosfast,udsfast,$
       cdfsgalfit,cosgalfit,udsgalfit,$
         filename='zfourge2016_ksselectedcat.sav'
end
 ;; rcat,'cdfs.v1.6.9.cat',/unpack
 ;;    cdfs={id:reform(id),ra:reform(ra),dec:reform(dec),wmin_optical:reform(wmin_optical),wmin_fs:reform(wmin_fs),use:reform(use),star:reform(star)}
 ;;    ;{id:id,ra:ra,dec:dec,wmin_optical:wmin_optical,wmin_fs:wmin_fs,use:use,star:star}
 ;;    rcat,'cdfs.v1.6.9.zout',/unpack
 ;;    rcat,'cdfs.v1.6.9.rest.v0.9.cat',/unpack
 ;;    Umag = 23.9 - 2.5*alog10(reform(f_U))
 ;;    Vmag = 23.9 - 2.5*alog10(reform(f_V))
 ;;    Jmag = 23.9 - 2.5*alog10(reform(f_J))
 ;;    umv = Umag - Vmag
 ;;    vmj = Vmag - Jmag
 ;;    cdfseazy={id:reform(id),z_peak:reform(z_peak),umv:umv,vmj:vmj}
 ;;    ;{id:id,z_peak:z_peak,umv:umv,vmj:vmj}


 ;;    rcat,'uds.v1.5.8.cat',/unpack
 ;;    uds={id:reform(id),ra:reform(ra),dec:reform(dec),wmin_optical:reform(wmin_optical),wmin_fs:reform(wmin_fs),use:reform(use),star:reform(star)}
 ;;    rcat,'uds.v1.5.8.zout',/unpack
 ;;    rcat,'uds.v1.5.8.rest.v0.9.cat',/unpack
 ;;    Umag = 23.9 - 2.5*alog10(reform(f_U))
 ;;    Vmag = 23.9 - 2.5*alog10(reform(f_V))
 ;;    Jmag = 23.9 - 2.5*alog10(reform(f_J))
 ;;    umv = Umag - Vmag
 ;;    vmj = Vmag - Jmag
 ;;    udseazy={id:reform(id),z_peak:reform(z_peak),umv:umv,vmj:vmj}
    
