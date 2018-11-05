pro twodim_binning,xaxis=xaxis,yaxis=yaxis,counts=counts,$
                   binsizex=binsizex,binsizey=binsizey,xcenter=obin1cen,ycenter=obin2cen,$
                   min1=min1,max1=max1,min2=min2,max2=max2
;Make 2D hist and calculate median of U-V color in each bin
;binsizex = 0.2
;binsizey=0.2
  ;if not keyword_set(min1) then min1=min(xaxis)
  ;if not keyword_set(max1) then max1=max(xaxis)
  ;if not keyword_set(min2) then min1=min(yaxis)
  ;if not keyword_set(max2) then max1=max(yaxis)
  ;print,size(xaxis),size(yaxis)
  h2d = HIST2D(xaxis,yaxis, binsize1=binsizex, binsize2=binsizey,min1=min1,max1=max1,min2=min2,max2=max2, $
             obin1=obin1left,obin2=obin2left,BINEDGE1=-1,BINEDGE2=-1)
  h2d = HIST2D(xaxis,yaxis, binsize1=binsizex, binsize2=binsizey,min1=min1,max1=max1,min2=min2,max2=max2, $
             obin1=obin1right,obin2=obin2right,BINEDGE1=1,BINEDGE2=1)
h2d = HIST2D(xaxis,yaxis, binsize1=binsizex, binsize2=binsizey,min1=min1,max1=max1,min2=min2,max2=max2, $
             obin1=obin1cen,obin2=obin2cen,BINEDGE1=0,BINEDGE2=0)
counts=h2d*0.0
for i=0,(size(h2d))[1] -1 do begin
  for j=0,(size(h2d))[2] -1 do begin
     t=where((xaxis ge obin1left[i]) and (xaxis lt obin1right[i]) and (yaxis ge obin2left[j]) and (yaxis lt obin2right[j]),ctotal)
     if (ctotal ge 1) then begin
        counts[i,j]=ctotal
     endif 
     ;print,obin1left[i],obin1right[i],obin2left[j],obin2right[j],ctotal,z_median[i,j]
  endfor  
  endfor

end
; compute the Bayesian surface density  (sigma primed in Eq. 3 from Kawinwanichakij+17) 
pro density_calc,NN=NN,dn=dn,rnd_dn=rnd_dn,$
                 real_bayes_nn=real_bayes_nn, rand_bayes_nn=rand_bayes_nn
  NN=3
  sum=fltarr(n_elements(dn[0,*]))
  sum_rand= fltarr(n_elements(rnd_dn[0,*]))
  for j=0,n_elements(dn[0,*])-1 do begin
     for i=0,NN-1 do begin
        sum[j] = sum[j] + reform(dn[i,j]/60.)^2
        sum_rand[j] = sum_rand[j] + reform(rnd_dn[i,j]/60.)^2
     endfor
  end
  real_bayes_nn=1/sum
  rand_bayes_nn  =1 /sum_rand

  
end

pro meandensity,zmin=zmin,zmax=zmax,lmassmin=lmassmin,lmassmax=lmassmax,$
                real_mean_density=real_mean_density,rand_mean_density=rand_mean_density,number_density=number_density

                                ;zmin=0.5
  ;zmax=2.
  ;lmassmin=9.5
  ;lmassmax=12
  
restore,'../zfourge2016_ksselectedcat.sav'
restore,'../calc_4field_randompos_091216.sav' ;uniformly ranomd ra and dec
_cos_dn_rand=cos_dn
_uds_dn_rand=uds_dn
_cdfs_dn_rand=cdfs_dn
restore,'../calc_4field_091216.sav'
_cos_dn=cos_dn
_cdfs_dn = cdfs_dn
_uds_dn = uds_dn

;;;cdfs
cdfs_masscomp = get_masslim(cdfseazy.z_peak,Klim=25.5,quiescent=1,completeness=90)
diffmass =  cdfsfast.lmass - cdfs_masscomp
cdfs_masscomp_flag = intarr(n_elements(cdfseazy.z_peak))
cdfs_masscomp_flag[where(diffmass ge 0)] = 1
cdfs_masscomp_flag[where(diffmass lt 0)] = 0

td=selstandard(cdfs)
ttd=where((cdfseazy[td].z_peak ge zmin) and (cdfseazy[td].z_peak le zmax) and (cdfsfast[td].lmass ge lmassmin) and (cdfsfast[td].lmass le lmassmax) and (cdfs_masscomp_flag[td] eq 1),counttot)
cdfs_dn = _cdfs_dn[*,td[ttd]]
cdfs_dn_rand = _cdfs_dn_rand[*,td[ttd]]

;;;cos
cos_masscomp = get_masslim(coseazy.z_peak,Klim=25.5,quiescent=1,completeness=90)
diffmass =  cosfast.lmass - cos_masscomp
cos_masscomp_flag = intarr(n_elements(coseazy.z_peak))
cos_masscomp_flag[where(diffmass ge 0)] = 1
cos_masscomp_flag[where(diffmass lt 0)] = 0
td=selstandard(cos)
ttd=where((coseazy[td].z_peak ge zmin) and (coseazy[td].z_peak le zmax) and (cosfast[td].lmass ge lmassmin) and (cosfast[td].lmass le lmassmax) and (cos_masscomp_flag[td] eq 1),counttot)
cos_dn = _cos_dn[td[*,ttd]]
cos_dn_rand = _cos_dn_rand[*,td[ttd]]

;;;uds
uds_masscomp = get_masslim(udseazy.z_peak,Klim=25.5,quiescent=1,completeness=90)
diffmass =  udsfast.lmass - uds_masscomp
uds_masscomp_flag = intarr(n_elements(udseazy.z_peak))
uds_masscomp_flag[where(diffmass ge 0)] = 1
uds_masscomp_flag[where(diffmass lt 0)] = 0
td=selstandard(uds)
ttd=where((udseazy[td].z_peak ge zmin) and (udseazy[td].z_peak le zmax) and (udsfast[td].lmass ge lmassmin) and (udsfast[td].lmass le lmassmax) and (uds_masscomp_flag[td] eq 1),counttot)
uds_dn = _uds_dn[*,td[ttd]]
uds_dn_rand = _uds_dn_rand[*,td[ttd]]


nall = n_elements(cdfs_dn[0,*])+n_elements(cos_dn[0,*])+n_elements(uds_dn[0,*])
dn = make_array(7,nall,/float,value=0.0)
dn_rand = make_array(7,nall,/float,value=0.0)
for i =0,6 do begin 
   dn[i,*] = [reform(cdfs_dn[i,*]),reform(cos_dn[i,*]),reform(uds_dn[i,*])]
   dn_rand[i,*] = [reform(cdfs_dn_rand[i,*]),reform(cos_dn_rand[i,*]),reform(uds_dn_rand[i,*])]
endfor

help,dn
help,dn_rand
density_calc,dn=dn,rnd_dn=dn_rand,real_bayes_nn=bayes_nn,rand_bayes_nn=bayes_nn_rand

number_density= n_elements(dn[0,*])/(11.*11.*3)
real_mean_density=mean(bayes_nn[where(finite(bayes_nn) eq 1)]);n_elements(nnbay)/(11.*11.*3)
print,'mean_density=',real_mean_density,' arcmin^-2'


rand_mean_density = mean(bayes_nn_rand)
print,'rand_mean_density =',rand_mean_density 
C= real_mean_density/rand_mean_density
;print,'C=',C
end



pro calcsigma_1pd_ver2
;,zmin=zmin,zmax=zmax,lmassmin=lmassmin,lmassmax=lmassmax,$
restore,'../zfourge2016_ksselectedcat.sav'
restore,'../calc_4field_randompos_091216.sav' ;uniformly ranomd ra and dec

;;;;; Save density and cobs to file for all ZFOURGE galaxies
cdfs_log1pdelta= replicate(-99.0, n_elements(cdfs))
cdfs_25tile=replicate(-99.0, n_elements(cdfs)) ; for 1+delta
cdfs_50tile=replicate(-99.0, n_elements(cdfs)) ; for 1+delta
cdfs_75tile=replicate(-99.0, n_elements(cdfs)) ; for 1+delta
cdfs_logsigman=replicate(-99.0, n_elements(cdfs)) 
cdfs_25tile_logsigman=replicate(-99.0, n_elements(cdfs)) ; for Sigma N
cdfs_50tile_logsigman=replicate(-99.0, n_elements(cdfs)) ; for Sigma N
cdfs_75tile_logsigman=replicate(-99.0, n_elements(cdfs)) ; for Sigma N


cos_log1pdelta= replicate(-99.0, n_elements(cos))
cos_logsigman=replicate(-99.0, n_elements(cos))
cos_25tile=replicate(-99.0, n_elements(cos)) ; for 1+delta
cos_50tile=replicate(-99.0, n_elements(cos)) ; for 1+delta
cos_75tile=replicate(-99.0, n_elements(cos)) ; for 1+delta
cos_25tile_logsigman=replicate(-99.0, n_elements(cos)) ; for Sigma N
cos_50tile_logsigman=replicate(-99.0, n_elements(cos)) ; for Sigma N
cos_75tile_logsigman=replicate(-99.0, n_elements(cos)) ; for Sigma N

uds_log1pdelta= replicate(-99.0, n_elements(uds))
uds_logsigman=replicate(-99.0, n_elements(uds))
uds_25tile=replicate(-99.0, n_elements(uds)) ; for 1+delta
uds_50tile=replicate(-99.0, n_elements(uds)) ; for 1+delta
uds_75tile=replicate(-99.0, n_elements(uds)) ; for 1+delta
uds_25tile_logsigman=replicate(-99.0, n_elements(uds)) ; for Sigma N
uds_50tile_logsigman=replicate(-99.0, n_elements(uds)) ; for Sigma N
uds_75tile_logsigman=replicate(-99.0, n_elements(uds)) ; for Sigma N

_cos_dn_rand=cos_dn
_uds_dn_rand=uds_dn
_cdfs_dn_rand=cdfs_dn
restore,'../calc_4field_091216.sav'
_cos_dn=cos_dn
_cdfs_dn = cdfs_dn
_uds_dn = uds_dn
;zbin=[0.5,1.0,1.5,2.0]
;lmassbin = [8.8,9.8,10.5,11.5]
;for ii = 0,n_elements(zbin)-2 do begin
;   zmin = zbin[ii]
;   zmax = zbin[ii+1]
;   for jj=0,n_elements(lmassbin)-2 do begin
;      lmassmin = lmassbin[jj]
;      lmassmax = lmassbin[jj+1]
;;;cdfs
zmin = 0.3
zmax=  2.3
lmassmin=8.1;8.8
lmassmax=12

cdfs_masscomp = get_masslim(cdfseazy.z_peak,Klim=25.5,quiescent=1,completeness=90)
diffmass =  cdfsfast.lmass - cdfs_masscomp
cdfs_masscomp_flag = intarr(n_elements(cdfseazy.z_peak))
cdfs_masscomp_flag[where(diffmass ge 0)] = 1
cdfs_masscomp_flag[where(diffmass lt 0)] = 0

td_cdfs=selstandard(cdfs)
ttd_cdfs=where((cdfseazy[td_cdfs].z_peak ge zmin) and (cdfseazy[td_cdfs].z_peak lt zmax) and (cdfsfast[td_cdfs].lmass ge lmassmin) and (cdfsfast[td_cdfs].lmass lt lmassmax) and (cdfs_masscomp_flag[td_cdfs] eq 1),counttot)
cdfs_zpeak=cdfseazy[td_cdfs[ttd_cdfs]].z_peak
cdfs_umv =cdfseazy[td_cdfs[ttd_cdfs]].umv
cdfs_vmj =cdfseazy[td_cdfs[ttd_cdfs]].vmj
cdfs_lmass= cdfsfast[td_cdfs[ttd_cdfs]].lmass
cdfs_sersic = cdfsgalfit[td_cdfs[ttd_cdfs]].n
cdfs_re = cdfsgalfit[td_cdfs[ttd_cdfs]].re
cdfs_dn = _cdfs_dn[*,td_cdfs[ttd_cdfs]]
cdfs_dn_rand = _cdfs_dn_rand[*,td_cdfs[ttd_cdfs]]
;cdfs_masscomp = get_masslim(cdfs_zpeak,Klim=25.5,quiescent=1,completeness=90)
;diffmass =  cdfs_lmass - cdfs_masscomp
;cdfs_masscomp_flag = intarr(n_elements(cdfs_zpeak))
;cdfs_masscomp_flag[where(diffmass ge 0)] = 1
;cdfs_masscomp_flag[where(diffmass lt 0)] = 0
density_calc,NN=3,dn=cdfs_dn,rnd_dn=cdfs_dn_rand,$
             real_bayes_nn=cdfs_nn, rand_bayes_nn=cdfs_nn_rand

;;;;;cosmos
cos_masscomp = get_masslim(coseazy.z_peak,Klim=25.5,quiescent=1,completeness=90)
diffmass =  cosfast.lmass - cos_masscomp
cos_masscomp_flag = intarr(n_elements(coseazy.z_peak))
cos_masscomp_flag[where(diffmass ge 0)] = 1
cos_masscomp_flag[where(diffmass lt 0)] = 0

td_cos=selstandard(cos)
ttd_cos=where((coseazy[td_cos].z_peak ge zmin) and (coseazy[td_cos].z_peak lt zmax) and (cosfast[td_cos].lmass ge lmassmin) and (cosfast[td_cos].lmass lt lmassmax) and (cos_masscomp_flag[td_cos] eq 1),counttot)
cos_zpeak=coseazy[td_cos[ttd_cos]].z_peak
cos_umv =coseazy[td_cos[ttd_cos]].umv
cos_vmj =coseazy[td_cos[ttd_cos]].vmj
cos_lmass= cosfast[td_cos[ttd_cos]].lmass
cos_sersic = cosgalfit[td_cos[ttd_cos]].n
cos_re = cosgalfit[td_cos[ttd_cos]].re
cos_dn = _cos_dn[*,td_cos[ttd_cos]]
cos_dn_rand = _cos_dn_rand[*,td_cos[ttd_cos]]
;cos_masscomp = get_masslim(cos_zpeak,Klim=25.5,quiescent=1,completeness=90)
;diffmass =  cos_lmass - cos_masscomp
;cos_masscomp_flag = intarr(n_elements(cos_zpeak))
;cos_masscomp_flag[where(diffmass ge 0)] = 1
;cos_masscomp_flag[where(diffmass lt 0)] = 0

density_calc,NN=3,dn=cos_dn,rnd_dn=cos_dn_rand,$
             real_bayes_nn=cos_nn, rand_bayes_nn=cos_nn_rand

;;;;;; UDS
uds_masscomp = get_masslim(udseazy.z_peak,Klim=25.5,quiescent=1,completeness=90)
diffmass =  udsfast.lmass - uds_masscomp
uds_masscomp_flag = intarr(n_elements(udseazy.z_peak))
uds_masscomp_flag[where(diffmass ge 0)] = 1
uds_masscomp_flag[where(diffmass lt 0)] = 0
td_uds=selstandard(uds)
ttd_uds=where((udseazy[td_uds].z_peak ge zmin) and (udseazy[td_uds].z_peak lt zmax) and (udsfast[td_uds].lmass ge lmassmin) and (udsfast[td_uds].lmass lt lmassmax) and (uds_masscomp_flag[td_uds] eq 1),counttot)
uds_zpeak=udseazy[td_uds[ttd_uds]].z_peak
uds_umv =udseazy[td_uds[ttd_uds]].umv
uds_vmj =udseazy[td_uds[ttd_uds]].vmj
uds_lmass= udsfast[td_uds[ttd_uds]].lmass
uds_sersic = udsgalfit[td_uds[ttd_uds]].n
uds_re = udsgalfit[td_uds[ttd_uds]].re
uds_dn = _uds_dn[*,td_uds[ttd_uds]]
uds_dn_rand = _uds_dn_rand[*,td_uds[ttd_uds]]
;uds_masscomp = get_masslim(uds_zpeak,Klim=25.5,quiescent=1,completeness=90)
;diffmass =  uds_lmass - uds_masscomp
;uds_masscomp_flag = intarr(n_elements(uds_zpeak))
;uds_masscomp_flag[where(diffmass ge 0)] = 1
;uds_masscomp_flag[where(diffmass lt 0)] = 0

density_calc,NN=3,dn=uds_dn,rnd_dn=uds_dn_rand,$
             real_bayes_nn=uds_nn, rand_bayes_nn=uds_nn_rand 

redshift=[cdfs_zpeak,cos_zpeak,uds_zpeak]
fields = [1+intarr(n_elements(cdfs_zpeak)),2+intarr(n_elements(cos_zpeak)),3+intarr(n_elements(uds_zpeak))]

nall = n_elements(cdfs_dn[0,*])+n_elements(cos_dn[0,*])+n_elements(uds_dn[0,*])
dn_all = make_array(7,nall,/float,value=0.0)
dn_rand_all = make_array(7,nall,/float,value=0.0)
for i =0,6 do begin 
   dn_all[i,*] = [reform(cdfs_dn[i,*]),reform(cos_dn[i,*]),reform(uds_dn[i,*])]
   dn_rand_all[i,*] = [reform(cdfs_dn_rand[i,*]),reform(cos_dn_rand[i,*]),reform(uds_dn_rand[i,*])]
endfor


density_calc,NN=3,dn=dn_all,rnd_dn=dn_rand_all,$
             real_bayes_nn=bayes_nn_all, rand_bayes_nn= bayes_nn_rand_all

;mlim=get_masslim(zmax,Klim=25.5,quiescent=1,completeness=90)
meandensity,zmin=zmin,zmax=zmax,lmassmin=lmassmin,lmassmax=lmassmax,real_mean_density=mean_density,rand_mean_density=mean_density_rand,number_density=number_density

C=  mean_density/mean_density_rand
;print,'C=',C
sigman_all = C*bayes_nn_all
onepdelta_all = 1+(sigman_all - number_density)/number_density
log1pd_all = alog10(onepdelta_all)
;;;;;;;;;;Quantile regression for 1+delta
filename='overdensity_redshift.txt'
openw,lun,filename,/get_lun
for i =0,n_elements(redshift)-1 do begin
   ;printf,lun,redshift[i],sigman_all[i],onepdelta_all[i],format='(f,x,f,x,f)'
   printf,lun,redshift[i],onepdelta_all[i],fields[i],format='(f,x,f,x,i)'
endfor
close,lun
free_lun,lun
;plot,redshift,alog10(sigman_all),psym=1
;Run Cobs in R
spawn,'Rscript runcobs.r'
readcol,'cobs_output.txt',id,z,q25,q50,q75,field 

diffq25  = log1pd_all - q25
diffq75 = log1pd_all - q75

t=where((diffq25 lt 0) ,clow)

tt=where(diffq75 gt 0,chigh)

print,zmin,zmax,lmassmin,lmassmax,clow,chigh
;;;;; Calculation overdensity for each field
onepdeltaN_cdfs = 1+(C*cdfs_nn - number_density)/number_density
onepdeltaN_cos = 1+(C*cos_nn - number_density)/number_density
onepdeltaN_uds = 1+(C*uds_nn - number_density)/number_density
cdfs_q25=q25[where(field eq 1)]
cos_q25=q25[where(field eq 2)]
uds_q25 = q25[where(field eq 3)]

cdfs_q50=q50[where(field eq 1)]
cos_q50=q50[where(field eq 2)]
uds_q50 = q50[where(field eq 3)]

cdfs_q75=q75[where(field eq 1)]
cos_q75=q75[where(field eq 2)]
uds_q75 = q75[where(field eq 3)]

cdfs_log1pdelta[td_cdfs[ttd_cdfs]] = alog10(onepdeltaN_cdfs)
cdfs_logsigman[td_cdfs[ttd_cdfs]]=alog10(C*cdfs_nn)
cdfs_25tile[td_cdfs[ttd_cdfs]] = cdfs_q25
cdfs_50tile[td_cdfs[ttd_cdfs]] = cdfs_q50
cdfs_75tile[td_cdfs[ttd_cdfs]] = cdfs_q75

cos_log1pdelta[td_cos[ttd_cos]] = alog10(onepdeltaN_cos)
cos_logsigman[td_cos[ttd_cos]]=alog10(C*cos_nn)
cos_25tile[td_cos[ttd_cos]] = cos_q25
cos_50tile[td_cos[ttd_cos]] = cos_q50
cos_75tile[td_cos[ttd_cos]] = cos_q75

uds_log1pdelta[td_uds[ttd_uds]] = alog10(onepdeltaN_uds)
uds_logsigman[td_uds[ttd_uds]]=alog10(C*uds_nn)
uds_25tile[td_uds[ttd_uds]] = uds_q25
uds_50tile[td_uds[ttd_uds]] = uds_q50
uds_75tile[td_uds[ttd_uds]] = uds_q75

;;;;; Quantile regression for SIGMA
filename='sigma_redshift.txt'
openw,lun,filename,/get_lun
for i =0,n_elements(redshift)-1 do begin
   printf,lun,redshift[i],sigman_all[i],fields[i],format='(f,x,f,x,i)'
endfor
close,lun
free_lun,lun
;Run Cobs in R
spawn,'Rscript runcobs_sigma.r'
readcol,'cobs_sigma_output.txt',id,z,q25_logsigma,q50_logsigma,q75_logsigma,field

cdfs_q25_logsigma=q25_logsigma[where(field eq 1)]
cos_q25_logsigma=q25_logsigma[where(field eq 2)]
uds_q25_logsigma = q25_logsigma[where(field eq 3)]

cdfs_q50_logsigma=q50_logsigma[where(field eq 1)]
cos_q50_logsigma=q50_logsigma[where(field eq 2)]
uds_q50_logsigma = q50_logsigma[where(field eq 3)]

cdfs_q75_logsigma=q75_logsigma[where(field eq 1)]
cos_q75_logsigma=q75_logsigma[where(field eq 2)]
uds_q75_logsigma = q75_logsigma[where(field eq 3)]

cdfs_25tile_logsigman[td_cdfs[ttd_cdfs]] = cdfs_q25_logsigma
cdfs_50tile_logsigman[td_cdfs[ttd_cdfs]] = cdfs_q50_logsigma
cdfs_75tile_logsigman[td_cdfs[ttd_cdfs]] = cdfs_q75_logsigma

cos_25tile_logsigman[td_cos[ttd_cos]] = cos_q25_logsigma
cos_50tile_logsigman[td_cos[ttd_cos]] = cos_q50_logsigma
cos_75tile_logsigman[td_cos[ttd_cos]] = cos_q75_logsigma

uds_25tile_logsigman[td_uds[ttd_uds]] = uds_q25_logsigma
uds_50tile_logsigman[td_uds[ttd_uds]] = uds_q50_logsigma
uds_75tile_logsigman[td_uds[ttd_uds]] = uds_q75_logsigma




save,cdfs_25tile,cdfs_50tile,cdfs_75tile,$
     cos_25tile,cos_50tile,cos_75tile,$
     uds_25tile,uds_50tile,uds_75tile,$
     cdfs_log1pdelta,cos_log1pdelta,uds_log1pdelta,$
     cdfs_25tile_logsigman,cdfs_50tile_logsigman,cdfs_75tile_logsigman,$
     cos_25tile_logsigman,cos_50tile_logsigman,cos_75tile_logsigman,$
     uds_25tile_logsigman,uds_50tile_logsigman,uds_75tile_logsigman,$
     cdfs_logsigman,cos_logsigman,uds_logsigman,$
     
     ;filename = 'nnz_N3_bayesian_nancy_cobs_v3.sav'
     filename = 'nnz_N3_bayesian_nancy_cobs_v4.1.sav'


end
pro plotdensity_cobs
  zmin = 0.5
  zmax = 2.
  lmassmin = 8.1
  lmassmax = 12
  restore,'nnz_N3_bayesian_nancy_cobs_v4.sav'
  restore,'../zfourge2016_ksselectedcat.sav'
  ;;;cdfs
cdfs_masscomp = get_masslim(cdfseazy.z_peak,Klim=25.5,quiescent=1,completeness=90)
diffmass =  cdfsfast.lmass - cdfs_masscomp
cdfs_masscomp_flag = intarr(n_elements(cdfseazy.z_peak))
cdfs_masscomp_flag[where(diffmass ge 0)] = 1
cdfs_masscomp_flag[where(diffmass lt 0)] = 0

td=selstandard(cdfs)
ttd=where((cdfseazy[td].z_peak ge zmin) and (cdfseazy[td].z_peak le zmax) and (cdfsfast[td].lmass ge lmassmin) and (cdfsfast[td].lmass le lmassmax) and (cdfs_masscomp_flag[td] eq 1),counttot)
cdfs_zpeak=cdfseazy[td[ttd]].z_peak
cdfs_umv =cdfseazy[td[ttd]].umv
cdfs_vmj =cdfseazy[td[ttd]].vmj
cdfs_lmass= cdfsfast[td[ttd]].lmass
cdfs_sersic = cdfsgalfit[td[ttd]].n
cdfs_re = cdfsgalfit[td[ttd]].re
cdfs_log1pd = cdfs_log1pdelta[td[ttd]]
cdfs_logsigman = cdfs_logsigman[td[ttd]]
cdfs_q25 = cdfs_25tile[td[ttd]]
cdfs_q50 = cdfs_50tile[td[ttd]]
cdfs_q75 = cdfs_75tile[td[ttd]]
cdfs_q25_logsigman = cdfs_25tile_logsigman[td[ttd]]
cdfs_q50_logsigman = cdfs_50tile_logsigman[td[ttd]]
cdfs_q75_logsigman = cdfs_75tile_logsigman[td[ttd]]


cos_masscomp = get_masslim(coseazy.z_peak,Klim=25.5,quiescent=1,completeness=90)
diffmass =  cosfast.lmass - cos_masscomp
cos_masscomp_flag = intarr(n_elements(coseazy.z_peak))
cos_masscomp_flag[where(diffmass ge 0)] = 1
cos_masscomp_flag[where(diffmass lt 0)] = 0
td=selstandard(cos)
ttd=where((coseazy[td].z_peak ge zmin) and (coseazy[td].z_peak le zmax) and (cosfast[td].lmass ge lmassmin) and (cosfast[td].lmass le lmassmax) and (cos_masscomp_flag[td] eq 1),counttot)
cos_zpeak=coseazy[td[ttd]].z_peak
cos_umv =coseazy[td[ttd]].umv
cos_vmj =coseazy[td[ttd]].vmj
cos_lmass= cosfast[td[ttd]].lmass
cos_sersic = cosgalfit[td[ttd]].n
cos_re = cosgalfit[td[ttd]].re
cos_log1pd = cos_log1pdelta[td[ttd]]
cos_logsigman = cos_logsigman[td[ttd]]
cos_q25 = cos_25tile[td[ttd]]
cos_q50 = cos_50tile[td[ttd]]
cos_q75 = cos_75tile[td[ttd]]
cos_q25_logsigman = cos_25tile_logsigman[td[ttd]]
cos_q50_logsigman = cos_50tile_logsigman[td[ttd]]
cos_q75_logsigman = cos_75tile_logsigman[td[ttd]]

uds_masscomp = get_masslim(udseazy.z_peak,Klim=25.5,quiescent=1,completeness=90)
diffmass =  udsfast.lmass - uds_masscomp
uds_masscomp_flag = intarr(n_elements(udseazy.z_peak))
uds_masscomp_flag[where(diffmass ge 0)] = 1
uds_masscomp_flag[where(diffmass lt 0)] = 0
td=selstandard(uds)
ttd=where((udseazy[td].z_peak ge zmin) and (udseazy[td].z_peak le zmax) and (udsfast[td].lmass ge lmassmin) and (udsfast[td].lmass le lmassmax) and (uds_masscomp_flag[td] eq 1),counttot)
uds_zpeak=udseazy[td[ttd]].z_peak
uds_umv =udseazy[td[ttd]].umv
uds_vmj =udseazy[td[ttd]].vmj
uds_lmass= udsfast[td[ttd]].lmass
uds_sersic = udsgalfit[td[ttd]].n
uds_re = udsgalfit[td[ttd]].re
uds_log1pd = uds_log1pdelta[td[ttd]]
uds_logsigman = uds_logsigman[td[ttd]]
uds_q25 = uds_25tile[td[ttd]]
uds_q50 = uds_50tile[td[ttd]]
uds_q75 = uds_75tile[td[ttd]]
uds_q25_logsigman = uds_25tile_logsigman[td[ttd]]
uds_q50_logsigman = uds_50tile_logsigman[td[ttd]]
uds_q75_logsigman = uds_75tile_logsigman[td[ttd]]

q25= [cdfs_q25,cos_q25,uds_q25]
q50= [cdfs_q50,cos_q50,uds_q50]
q75= [cdfs_q75,cos_q75,uds_q75]
redshift = [cdfs_zpeak,cos_zpeak,uds_zpeak]
log1pd  =  [cdfs_log1pd,cos_log1pd,uds_log1pd]
logsigman = [cdfs_logsigman,cos_logsigman,uds_logsigman]
q25_logsigman= [cdfs_q25_logsigman,cos_q25_logsigman,uds_q25_logsigman]
q50_logsigman= [cdfs_q50_logsigman,cos_q50_logsigman,uds_q50_logsigman]
q75_logsigman= [cdfs_q75_logsigman,cos_q75_logsigman,uds_q75_logsigman]

N=3
;; overplot the lines of HDE, LDE, MDE from the previous version (not
;; using COBs) on 240517
restore,'~/ZF_environ/nearest_neighbor/bayesian_n3_z_240517_forplot.sav'

;psplot, file='plots/bayesian_n'+StrTrim(String(N,format='(i)'),2)+'_z_040817_with240517.ps', xsize=7, ysize=10.5,/color
psplot, file='plots/bayesian_n'+StrTrim(String(N,format='(i)'),2)+'_z_110817_with240517_largelabel.ps', xsize=7, ysize=10.5,/color
!P.multi=[0,1,2]
;;;;plot Sigma N

binsizex=0.05
binsizey=0.1
minx=0.45
maxx=2.05
miny=0.
maxy=3.0
xr=[0.4,2.1]
yr=[0.,2.5]
plot,redshift,logsigman,psym=3,xr=xr,yr=yr,xtit='z',$
     ytit='Log('+cgGreek('Sigma')+cgSymbol("prime")+'!B'+StrTrim(String(N,format='(i)'),2)+'!N/arcmin!E-2!N)',$
     color=cgColor('Black'),/nodata,/xst,/yst,charsize=1.7,charthick=6

twodim_binning,xaxis=redshift,yaxis=logsigman,counts=hist,binsizex=binsizex,$
               binsizey=binsizey,xcenter=xcenbin,ycenter=ycenbin,min1=minx,max1=maxx,min2=miny,max2=maxy
scaledDensity = BytScl(hist,min=min(hist),max=max(hist),top=255)
loadct,63,/silent
plotsym,8,0.5,/fill
for k=0,n_elements(xcenbin)-1 do begin
   for m=0,n_elements(ycenbin)-1 do begin
      if (scaledDensity[k,m] gt 0) then begin
         boxx=(binsizex/2.)
         boxy=(binsizey/2.)
         if  ( (ycenbin[m]+boxy lt yr[1]) and (xcenbin[k]+boxx lt xr[1]) and (ycenbin[m]-boxy gt yr[0]+0.05)  and (xcenbin[k]-boxx gt xr[0]) ) then begin
          Polyfill,[xcenbin[k]-boxx,xcenbin[k]+boxx,xcenbin[k]+boxx,xcenbin[k]-boxx],$
                  [ycenbin[m]-boxy,ycenbin[m]-boxy,ycenbin[m]+boxy,ycenbin[m]+boxy],Color=scaledDensity[k,m]
         endif
      ;oplot,[xcenbin[k]],[ycenbin[m]],psym=8,color=scaledDensity[k,m]
   endif 
   endfor
endfor 
;readcol,'sigman_cobs_output_0.5z2.txt',id,z,q25,q50,q75
oplot,redshift[sort(redshift)],q25_logsigman[sort(redshift)],linestyle=2,color=cgColor('Salmon'),thick=8
oplot,redshift[sort(redshift)],q50_logsigman[sort(redshift)],linestyle=2,color=cgColor('Black'),thick=8
oplot,redshift[sort(redshift)],q75_logsigman[sort(redshift)],linestyle=2,color=cgColor('red'),thick=8

;;Not COBs
;oplot,zcenter,alog10(sigma_MDE),linestyle=3,color=cgColor('black')
;oplot,zcenter,alog10(sigma_LDE),linestyle=3,color=cgColor('Salmon')
;oplot,zcenter,alog10(sigma_HDE),linestyle=3,color=cgColor('red')

;;;; Plot log(1+delta)
binsizex=0.05
binsizey=0.1
minx=0.45
maxx=2.05
miny=-1.5
maxy=2.0
xr=[0.4,2.1]
yr=[-1.5,1.5]
plot, redshift,log1pd,$
      xr=xr,yr=yr,xtit='z',$
     ytit='log(1+'+cgGreek("delta")+cgSymbol("prime")+')!B'+StrTrim(String(N,format='(i)'),2)+'!N',$
     color=cgColor('Black'),/nodata,/xst,/yst,charsize=1.7,charthick=6
twodim_binning,xaxis=redshift,yaxis=log1pd,counts=hist,binsizex=binsizex,$
               binsizey=binsizey,xcenter=xcenbin,ycenter=ycenbin,min1=minx,max1=maxx,min2=miny,max2=maxy
scaledDensity = BytScl(hist,min=min(hist),max=max(hist),top=255)
loadct,63,/silent
plotsym,8,0.5,/fill
for k=0,n_elements(xcenbin)-1 do begin
   for m=0,n_elements(ycenbin)-1 do begin
      if (scaledDensity[k,m] gt 0) then begin
         boxx=(binsizex/2.)
         boxy=(binsizey/2.)
         if  ( (ycenbin[m]+boxy lt yr[1]) and (xcenbin[k]+boxx lt xr[1]) and (ycenbin[m]-boxy gt yr[0]+0.05)  and (xcenbin[k]-boxx gt xr[0]) ) then begin
          Polyfill,[xcenbin[k]-boxx,xcenbin[k]+boxx,xcenbin[k]+boxx,xcenbin[k]-boxx],$
                  [ycenbin[m]-boxy,ycenbin[m]-boxy,ycenbin[m]+boxy,ycenbin[m]+boxy],Color=scaledDensity[k,m]
         endif
      ;oplot,[xcenbin[k]],[ycenbin[m]],psym=8,color=scaledDensity[k,m]
   endif 
   endfor
endfor 


;oplot,cdfs_zpeak, cdfs_log1pd,psym=2,color=cgColor('red')
;oplot,cos_zpeak, cos_log1pd,psym=2,color=cgColor('blue')
;oplot,uds_zpeak, uds_log1pd,psym=2,color=cgColor('green')
;quar25= [cdfs_q25,cos_q25,uds_q25]
;quar50= [cdfs_q50,cos_q50,uds_q50]
;quar75= [cdfs_q75,cos_q75,uds_q75]

;readcol,'1pd_cobs_output_0.5z2.txt',id,z,q25,q50,q75
oplot,redshift[sort(redshift)],q25[sort(redshift)],linestyle=2,color=cgColor('Salmon'),thick=8
oplot,redshift[sort(redshift)],q50[sort(redshift)],linestyle=2,color=cgColor('Black'),thick=8
oplot,redshift[sort(redshift)],q75[sort(redshift)],linestyle=2,color=cgColor('red'),thick=8

;oplot,z[sort(z)],q25[sort(z)],linestyle=2,color=cgColor('Salmon'),thick=8
;oplot,z[sort(z)],q50[sort(z)],linestyle=2,color=cgColor('Black'),thick=8
;oplot,z[sort(z)],q75[sort(z)],linestyle=2,color=cgColor('red'),thick=8

;xyouts,0.8,0.6,'Highest density quartile',charsize=1.9,charthick=20,color=cgColor('White')
;xyouts,0.8,0.6,'Highest density quartile',charsize=1.9,charthick=6,color=cgColor('Red')
;xyouts,0.8,-0.6,'Lowest density quartile',charsize=1.9,charthick=20,color=cgColor('White')
;xyouts,0.8,-0.6,'Lowest density quartile',charsize=1.9,charthick=6,color=cgColor('Salmon')

xyouts,0.6,0.6,'Highest density quartile',charsize=2.5,charthick=20,color=cgColor('White')
xyouts,0.6,0.6,'Highest density quartile',charsize=2.5,charthick=6,color=cgColor('Red')
xyouts,0.6,-0.65,'Lowest density quartile',charsize=2.5,charthick=20,color=cgColor('White')
xyouts,0.6,-0.65,'Lowest density quartile',charsize=2.5,charthick=6,color=cgColor('Salmon')
;Not COBS
;oplot,zcenter,alog10(onepd_MDE),linestyle=3,color=cgColor('black')
;oplot,zcenter,alog10(onepd_LDE),linestyle=3,color=cgColor('Salmon')
;oplot,zcenter,alog10(onepd_HDE),linestyle=3,color=cgColor('red')

psplot,/stop

end

;Assign galaxy into low and high density quartile
;and see if they are equal number of galaxies in D1 and D4
pro lowhighq
filename='sigma_overdensity_redshift_0.3z2.txt'
readcol,filename,redshift,sigman,onepdelta
log1pd= alog10(onepdelta)
readcol,'cobs_output3.txt',z,q25,q50,q75;,format='a,x,f,x,f,x,f'
Dq=intarr(n_elements(redshift))
diffq25  = log1pd - q25
diffq75 = log1pd - q75

t=where((diffq25 lt 0) ,clow)
Dq[t] = 25
tt=where(diffq75 gt 0,chigh)
Dq[tt] = 75
print,clow,chigh
;for i =0,n_elements(redshift)-1 do begin
;   if (log1pd[i] lt q25[i]) then  Dq[i] = 25
;   if (log1pd[i] gt q75[i]) then  Dq[i] = 75
;endfor
;t=where(Dq eq 25,clow)
;tt=where(Dq eq 75,chigh)
;print,clow,chigh
end



