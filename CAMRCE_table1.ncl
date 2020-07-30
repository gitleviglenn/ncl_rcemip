;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; CAMRCE_table1.ncl
;
; calculate the global mean values for table 1 in RCEMIP papers 
; for CAM5 and CAM6 at 295, 300, and 305K
;
; levi silvers                                                            july 2020
;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

path_cheyenne="/glade/scratch/levi/RCEMIP/"
path_C5="CAM5_GCM/";RCE_large305/2D/"
path_C6="CAM6_GCM/";RCE_large305/2D/"

; cl
file_cl_5_295=path_cheyenne+path_C5+"RCE_large295/2D/CAM5_GCM_RCE_large295_2D_cl.nc"
file_cl_5_300=path_cheyenne+path_C5+"RCE_large300/2D/CAM5_GCM_RCE_large300_2D_cl.nc"
file_cl_5_305=path_cheyenne+path_C5+"RCE_large305/2D/CAM5_GCM_RCE_large305_2D_cl.nc"

file_cl_6_295=path_cheyenne+path_C6+"RCE_large295/2D/CAM6_GCM_RCE_large295_2D_cl.nc"
file_cl_6_300=path_cheyenne+path_C6+"RCE_large300/2D/CAM6_GCM_RCE_large300_2D_cl.nc"
file_cl_6_305=path_cheyenne+path_C6+"RCE_large305/2D/CAM6_GCM_RCE_large305_2D_cl.nc"

; pr
pscale=86400; converts precipitation into mm/day
file_pr_avg_5_295=path_cheyenne+path_C5+"RCE_large295/0D/CAM5_GCM_RCE_large295_0D_pr_avg.nc"
file_pr_avg_5_300=path_cheyenne+path_C5+"RCE_large300/0D/CAM5_GCM_RCE_large300_0D_pr_avg.nc"
file_pr_avg_5_305=path_cheyenne+path_C5+"RCE_large305/0D/CAM5_GCM_RCE_large305_0D_pr_avg.nc"

file_pr_avg_6_295=path_cheyenne+path_C6+"RCE_large295/0D/CAM6_GCM_RCE_large295_0D_pr_avg.nc"
file_pr_avg_6_300=path_cheyenne+path_C6+"RCE_large300/0D/CAM6_GCM_RCE_large300_0D_pr_avg.nc"
file_pr_avg_6_305=path_cheyenne+path_C6+"RCE_large305/0D/CAM6_GCM_RCE_large305_0D_pr_avg.nc"

; prw
file_prw_avg_5_295=path_cheyenne+path_C5+"RCE_large295/0D/CAM5_GCM_RCE_large295_0D_prw_avg.nc"
file_prw_avg_5_300=path_cheyenne+path_C5+"RCE_large300/0D/CAM5_GCM_RCE_large300_0D_prw_avg.nc"
file_prw_avg_5_305=path_cheyenne+path_C5+"RCE_large305/0D/CAM5_GCM_RCE_large305_0D_prw_avg.nc"

file_prw_avg_6_295=path_cheyenne+path_C6+"RCE_large295/0D/CAM6_GCM_RCE_large295_0D_prw_avg.nc"
file_prw_avg_6_300=path_cheyenne+path_C6+"RCE_large300/0D/CAM6_GCM_RCE_large300_0D_prw_avg.nc"
file_prw_avg_6_305=path_cheyenne+path_C6+"RCE_large305/0D/CAM6_GCM_RCE_large305_0D_prw_avg.nc"

; rlut
file_rlut_avg_5_295=path_cheyenne+path_C5+"RCE_large295/0D/CAM5_GCM_RCE_large295_0D_rlut_avg.nc"
file_rlut_avg_5_300=path_cheyenne+path_C5+"RCE_large300/0D/CAM5_GCM_RCE_large300_0D_rlut_avg.nc"
file_rlut_avg_5_305=path_cheyenne+path_C5+"RCE_large305/0D/CAM5_GCM_RCE_large305_0D_rlut_avg.nc"

file_rlut_avg_6_295=path_cheyenne+path_C6+"RCE_large295/0D/CAM6_GCM_RCE_large295_0D_rlut_avg.nc"
file_rlut_avg_6_300=path_cheyenne+path_C6+"RCE_large300/0D/CAM6_GCM_RCE_large300_0D_rlut_avg.nc"
file_rlut_avg_6_305=path_cheyenne+path_C6+"RCE_large305/0D/CAM6_GCM_RCE_large305_0D_rlut_avg.nc"

;----

print("begin adding files...")
timelength=8760 ; 1 year 
endtime=26200 ; ICON GCM has fewer timesteps that CAM or ECHAM
startime=endtime-timelength+1

; below are files from CAM5 and CAM6
; cl
fc5_295_cl=addfile(file_cl_5_295,"r")
c5a=fc5_295_cl->cl(startime:endtime,:,:)
fc5_300=addfile(file_cl_5_300,"r")
c5b=fc5_300->cl(startime:endtime,:,:)
fc5_305=addfile(file_cl_5_305,"r")
c5c=fc5_305->cl(startime:endtime,:,:)

fc6_295=addfile(file_cl_6_295,"r")
c6a=fc6_295->cl(startime:endtime,:,:)
fc6_300=addfile(file_cl_6_300,"r")
c6b=fc6_300->cl(startime:endtime,:,:)
fc6_305=addfile(file_cl_6_305,"r")
c6c=fc6_305->cl(startime:endtime,:,:)

; pr 
fc5_295=addfile(file_pr_avg_5_295,"r")
c5a_pr=fc5_295->pr_avg(startime:endtime)
fc5_300=addfile(file_pr_avg_5_300,"r")
c5b_pr=fc5_300->pr_avg(startime:endtime)
fc5_305=addfile(file_pr_avg_5_305,"r")
c5c_pr=fc5_305->pr_avg(startime:endtime)

fc6_295=addfile(file_pr_avg_6_295,"r")
c6a_pr=fc6_295->pr_avg(startime:endtime)
fc6_300=addfile(file_pr_avg_6_300,"r")
c6b_pr=fc6_300->pr_avg(startime:endtime)
fc6_305=addfile(file_pr_avg_6_305,"r")
c6c_pr=fc6_305->pr_avg(startime:endtime)

; prw 
fc5_295=addfile(file_prw_avg_5_295,"r")
c5a_prw=fc5_295->prw_avg(startime:endtime)
fc5_300=addfile(file_prw_avg_5_300,"r")
c5b_prw=fc5_300->prw_avg(startime:endtime)
fc5_305=addfile(file_prw_avg_5_305,"r")
c5c_prw=fc5_305->prw_avg(startime:endtime)

fc6_295=addfile(file_prw_avg_6_295,"r")
c6a_prw=fc6_295->prw_avg(startime:endtime)
fc6_300=addfile(file_prw_avg_6_300,"r")
c6b_prw=fc6_300->prw_avg(startime:endtime)
fc6_305=addfile(file_prw_avg_6_305,"r")
c6c_prw=fc6_305->prw_avg(startime:endtime)

; rlut 
fc5_295=addfile(file_rlut_avg_5_295,"r")
c5a_rlut=fc5_295->rlut_avg(startime:endtime)
fc5_300=addfile(file_rlut_avg_5_300,"r")
c5b_rlut=fc5_300->rlut_avg(startime:endtime)
fc5_305=addfile(file_rlut_avg_5_305,"r")
c5c_rlut=fc5_305->rlut_avg(startime:endtime)

fc6_295=addfile(file_rlut_avg_6_295,"r")
c6a_rlut=fc6_295->rlut_avg(startime:endtime)
fc6_300=addfile(file_rlut_avg_6_300,"r")
c6b_rlut=fc6_300->rlut_avg(startime:endtime)
fc6_305=addfile(file_rlut_avg_6_305,"r")
c6c_rlut=fc6_305->rlut_avg(startime:endtime)

  printVarSummary(c5a)
  printVarSummary(c6a)

lat = fc5_295_cl->lat
rad = 4.0*atan(1.0)/180.0
clat = cos(lat*rad)

field_c5_a_weighted=new(timelength,float)
field_c5_b_weighted=new(timelength,float)
field_c5_c_weighted=new(timelength,float)

field_c6_a_weighted=new(timelength,float)
field_c6_b_weighted=new(timelength,float)
field_c6_c_weighted=new(timelength,float)

do ti=0,timelength-1
  field_c5_a_weighted(ti)     = wgt_areaave(c5a(ti,:,:),clat,1.0,0)
  field_c5_b_weighted(ti)     = wgt_areaave(c5b(ti,:,:),clat,1.0,0)
  field_c5_c_weighted(ti)     = wgt_areaave(c5c(ti,:,:),clat,1.0,0)
  ;field_c5_pr_a_weighted(ti)     = wgt_areaave(c5a_pr(ti,:,:),clat,1.0,0)
  ;field_c5_pr_b_weighted(ti)     = wgt_areaave(c5b_pr(ti,:,:),clat,1.0,0)
  ;field_c5_pr_c_weighted(ti)     = wgt_areaave(c5c_pr(ti,:,:),clat,1.0,0)
  ;field_c5_prw_a_weighted(ti)     = wgt_areaave(c5a_prw(ti,:,:),clat,1.0,0)
  ;field_c5_prw_b_weighted(ti)     = wgt_areaave(c5b_prw(ti,:,:),clat,1.0,0)
  ;field_c5_prw_c_weighted(ti)     = wgt_areaave(c5c_prw(ti,:,:),clat,1.0,0)
  ;field_c5_rlut_a_weighted(ti)     = wgt_areaave(c5a_rlut(ti,:,:),clat,1.0,0)
  ;field_c5_rlut_b_weighted(ti)     = wgt_areaave(c5b_rlut(ti,:,:),clat,1.0,0)
  ;field_c5_rlut_c_weighted(ti)     = wgt_areaave(c5c_rlut(ti,:,:),clat,1.0,0)
;
;
  field_c6_a_weighted(ti)     = wgt_areaave(c6a(ti,:,:),clat,1.0,0)
  field_c6_b_weighted(ti)     = wgt_areaave(c6b(ti,:,:),clat,1.0,0)
  field_c6_c_weighted(ti)     = wgt_areaave(c6c(ti,:,:),clat,1.0,0)
  ;field_c6_pr_a_weighted(ti)     = wgt_areaave(c6a_pr(ti,:,:),clat,1.0,0)
  ;field_c6_pr_b_weighted(ti)     = wgt_areaave(c6b_pr(ti,:,:),clat,1.0,0)
  ;field_c6_pr_c_weighted(ti)     = wgt_areaave(c6c_pr(ti,:,:),clat,1.0,0)
  ;field_c6_prw_a_weighted(ti)     = wgt_areaave(c6a_prw(ti,:,:),clat,1.0,0)
  ;field_c6_prw_b_weighted(ti)     = wgt_areaave(c6b_prw(ti,:,:),clat,1.0,0)
  ;field_c6_prw_c_weighted(ti)     = wgt_areaave(c6c_prw(ti,:,:),clat,1.0,0)
  ;field_c6_rlut_a_weighted(ti)     = wgt_areaave(c6a_rlut(ti,:,:),clat,1.0,0)
  ;field_c6_rlut_b_weighted(ti)     = wgt_areaave(c6b_rlut(ti,:,:),clat,1.0,0)
  ;field_c6_rlut_c_weighted(ti)     = wgt_areaave(c6c_rlut(ti,:,:),clat,1.0,0)
end do

;  print("dimensions of weighted mean field for 295K: "+dimsizes(field_c5_a_weighted))

  print("mean values of cl for CAM5 295K: "+dim_avg_n(field_c5_a_weighted,0)+" and 300K: "+dim_avg_n(field_c5_b_weighted,0)+" and 305K: "+dim_avg_n(field_c5_c_weighted,0))
  print("mean values of pr for CAM5 295K: "+pscale*dim_avg_n(c5a_pr,0)+" and 300K: "+pscale*dim_avg_n(c5b_pr,0)+" and 305K: "+pscale*dim_avg_n(c5c_pr,0))
  print("mean values of prw for CAM5 295K: "+dim_avg_n(c5a_prw,0)+" and 300K: "+dim_avg_n(c5b_prw,0)+" and 305K: "+dim_avg_n(c5c_prw,0))
  print("mean values of rlut for CAM5 295K: "+dim_avg_n(c5a_rlut,0)+" and 300K: "+dim_avg_n(c5b_rlut,0)+" and 305K: "+dim_avg_n(c5c_rlut,0))
;
  ;print("mean values of cl for CAM6 295K: "+dim_avg_n(field_c6_a_weighted,0))
  print("mean values of cl for CAM6 295K: "+dim_avg_n(field_c6_a_weighted,0)+" and 300K: "+dim_avg_n(field_c6_b_weighted,0)+" and 305K: "+dim_avg_n(field_c6_c_weighted,0))
  print("mean values of pr for CAM6 295K: "+pscale*dim_avg_n(c6a_pr,0)+" and 300K: "+pscale*dim_avg_n(c6b_pr,0)+" and 305K: "+pscale*dim_avg_n(c6c_pr,0))
  print("mean values of prw for CAM6 295K: "+dim_avg_n(c6a_prw,0)+" and 300K: "+dim_avg_n(c6b_prw,0)+" and 305K: "+dim_avg_n(c6c_prw,0))
  print("mean values of rlut for CAM6 295K: "+dim_avg_n(c6a_rlut,0)+" and 300K: "+dim_avg_n(c6b_rlut,0)+" and 305K: "+dim_avg_n(c6c_rlut,0))
;
end
