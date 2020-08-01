;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; clouds_CAMRCE.ncl
;
; calculate the global mean of the cloud fraction for CAM5 and CAM6 at 295, 300, and 305K
;
; levi silvers                                                            july 2020
;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

path_cheyenne="/glade/scratch/levi/RCEMIP/"
path_C5="CAM5_GCM/";RCE_large305/2D/"
path_C6="CAM6_GCM/";RCE_large305/2D/"

file_cl_5_295=path_cheyenne+path_C5+"RCE_large295/2D/CAM5_GCM_RCE_large295_2D_cl.nc"
file_cl_5_300=path_cheyenne+path_C5+"RCE_large300/2D/CAM5_GCM_RCE_large300_2D_cl.nc"
file_cl_5_305=path_cheyenne+path_C5+"RCE_large305/2D/CAM5_GCM_RCE_large305_2D_cl.nc"

file_cl_6_295=path_cheyenne+path_C6+"RCE_large295/2D/CAM6_GCM_RCE_large295_2D_cl.nc"
file_cl_6_300=path_cheyenne+path_C6+"RCE_large300/2D/CAM6_GCM_RCE_large300_2D_cl.nc"
file_cl_6_305=path_cheyenne+path_C6+"RCE_large305/2D/CAM6_GCM_RCE_large305_2D_cl.nc"

print("begin adding files...")
timelength=8760 ; 1 year 
endtime=26000 ; ICON GCM has fewer timesteps that CAM or ECHAM
startime=endtime-timelength+1

; below are files from CAM5 and CAM6
fc5_295=addfile(file_cl_5_295,"r")
c5a=fc5_295->cl(startime:endtime,:,:)
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

  printVarSummary(c5a)
  printVarSummary(c6a)

lat = fc5_295->lat
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
;
  field_c6_a_weighted(ti)     = wgt_areaave(c6a(ti,:,:),clat,1.0,0)
  field_c6_b_weighted(ti)     = wgt_areaave(c6b(ti,:,:),clat,1.0,0)
  field_c6_c_weighted(ti)     = wgt_areaave(c6c(ti,:,:),clat,1.0,0)
end do

;  print("dimensions of weighted mean field for 295K: "+dimsizes(field_c5_a_weighted))

  print("mean values of cl for CAM5 295K: "+dim_avg_n(field_c5_a_weighted,0)+" and 300K: "+dim_avg_n(field_c5_b_weighted,0)+" and 305K: "+dim_avg_n(field_c5_c_weighted,0))
;
  ;print("mean values of cl for CAM6 295K: "+dim_avg_n(field_c6_a_weighted,0))
  print("mean values of cl for CAM6 295K: "+dim_avg_n(field_c6_a_weighted,0)+" and 300K: "+dim_avg_n(field_c6_b_weighted,0)+" and 305K: "+dim_avg_n(field_c6_c_weighted,0))
;
end
