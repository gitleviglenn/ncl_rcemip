;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; omega500_CAMRCE.ncl
;
; analyze characteristics of the circulation using omega at 500 hPa
;
; calculates the intensity of circulation as measured by the difference 
; between the upwelling and downwelling regions.
;
; data is from CAM5 GCM RCE experiments and is at hourly frequency over 3 years
; dimensions of omega:	[time | 26281] x [lat | 192] x [lon | 288]
;
; timing: 
; 1 year:         8760 time steps
; 1 30 day month: 720 time steps
; 75 days:        1800 time steps
;
; specify the block of time over which values are averaged before calculating the pdf.
; Cronin and Wing used 5 day chunks.  Kevins tendency seems to be to not chunk at all
; but compute the pdf on as many samples as possible.  
;
; check how much of a difference the native grid makes to the PDFs.
;
; levi silvers                                                           may 2020
;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

; define indices for time steps read and used...
;timelength=24479 ; all but the first 75 days
;timelength=2880 ; 4 months
;timelength=240 ; 100 days
timelength=8760 ; 1 year 
;timelength=17520 ; 2 years
endtime=26200;
;endtime=25000 ; ICON GCM has fewer timesteps that CAM or ECHAM
;echam has 26280
;endtime=26279 ; this is the end time for CAM5 and CAM6
startime=endtime-timelength+1

; number of bins used for the PDFs
;nbins=250
nbins=50

; define paths to data
path_laptop="/Users/silvers/data/CAM/"
path_C5_lt="CAM5_GCM_RCE/"
path_C6_lt="CAM6_GCM_RCE/"

path_cheyenne="/glade/scratch/levi/RCEMIP/"
path_2_native="/glade/scratch/kareed/RCEMIP/"
path_C5="CAM5_GCM/";RCE_large305/2D/"
path_C6="CAM6_GCM/";RCE_large305/2D/"

;file_omega_5_295="/Users/silvers/data/CAM/CAM5_GCM_RCE/CAM5_GCM_RCE_large295_2D_wap500.nc"
;file_omega_6_295="/Users/silvers/data/CAM/CAM6_GCM_RCE/CAM6_GCM_RCE_large295_2D_wap500.nc"
;file_omega_5_300="/Users/silvers/data/CAM/CAM5_GCM_RCE/CAM5_GCM_RCE_large300_2D_wap500.nc"
;file_omega_6_300="/Users/silvers/data/CAM/CAM6_GCM_RCE/CAM6_GCM_RCE_large300_2D_wap500.nc"
;file_omega_5_305="/Users/silvers/data/CAM/CAM5_GCM_RCE/CAM5_GCM_RCE_large305_2D_wap500.nc"
;file_omega_6_305="/Users/silvers/data/CAM/CAM6_GCM_RCE/CAM6_GCM_RCE_large305_2D_wap500.nc"

file_omega_5_295=path_cheyenne+path_C5+"RCE_large295/2D/CAM5_GCM_RCE_large295_2D_wap500.nc"
file_omega_6_295=path_cheyenne+path_C6+"RCE_large295/2D/CAM6_GCM_RCE_large295_2D_wap500.nc"
file_omega_5_300=path_cheyenne+path_C5+"RCE_large300/2D/CAM5_GCM_RCE_large300_2D_wap500.nc"
file_omega_6_300=path_cheyenne+path_C6+"RCE_large300/2D/CAM6_GCM_RCE_large300_2D_wap500.nc"
file_omega_5_305=path_cheyenne+path_C5+"RCE_large305/2D/CAM5_GCM_RCE_large305_2D_wap500.nc"
file_omega_6_305=path_cheyenne+path_C6+"RCE_large305/2D/CAM6_GCM_RCE_large305_2D_wap500.nc"

  print("begin adding files...")

; below are files from CAM5 and CAM6
fc5_295=addfile(file_omega_5_295,"r")
c5a=fc5_295->wap500(startime:endtime,:,:)
fc6_295=addfile(file_omega_6_295,"r")
c6a=fc6_295->wap500(startime:endtime,:,:)

fc5_300=addfile(file_omega_5_300,"r")
c5b=fc5_300->wap500(startime:endtime,:,:)
fc6_300=addfile(file_omega_6_300,"r")
c6b=fc6_300->wap500(startime:endtime,:,:)
;
fc5_305=addfile(file_omega_5_305,"r")
c5c=fc5_305->wap500(startime:endtime,:,:)
fc6_305=addfile(file_omega_6_305,"r")
c6c=fc6_305->wap500(startime:endtime,:,:)
  
  print("convert to hPa/d")

c5a=864*c5a  ; convert to hPa/d
c6a=864*c6a  ; convert to hPa/d
c5b=864*c5b  ; convert to hPa/d
c6b=864*c6b  ; convert to hPa/d
c5c=864*c5c  ; convert to hPa/d
c6c=864*c6c  ; convert to hPa/d

;printVarSummary(c6a)

ntim=fc5_295->time

; look at 5 day chunks (for files that have hourly output)
;hourChunk=24
hourChunk=120
n5d=ntim/hourChunk
numperiods=timelength/hourChunk
print("numperiods ="+numperiods)


;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; import the data that is on the native grid...

ncols = 48602
sstarray=(/"295","300","305"/)
blab = new((/3,nbins/),float)
blab_a = new((/1,nbins/),float)
test_reshape=new((/hourChunk,numperiods,ncols/),float)
in = 0
; CAM5 on the native grid
do while(in .le. 2) 
  sststring=sstarray(in)
  files_native  = systemfunc ("ls "+path_2_native+"/CAM5_GCM_RCE_large"+sststring+"/run/*cam.h1.*nc")
  natives       = addfiles (files_native,"r")
  ListSetType(natives,"cat")
    print("natives names are: "+files_native)
  
  date = natives[:]->date
  
  c5a_test=natives[:]->OMEGA500(startime:endtime,:)
  c5a_test=864*c5a_test  ; convert to hPa/d
  
  test_reshape=reshape(c5a_test,(/hourChunk,numperiods,ncols/))
  
  ; compute time average...
  test_dmn=dim_avg_n(test_reshape,0)
  pdf_test=pdfx(test_dmn,nbins,False)

; Calculate PDFs for CAM5 
  if (in .lt. 1) then
    pdf_a=pdf_test
  else if (in .lt. 2) then
    pdf_b=pdf_test 
  else
    pdf_c=pdf_test
  end if
  end if

; CAM6 on the native grid
  files_native  = systemfunc ("ls "+path_2_native+"/CAM6_GCM_RCE_large"+sststring+"/run/*cam.h1.*nc")
  natives       = addfiles (files_native,"r")
  ListSetType(natives,"cat")
    print("natives names are: "+files_native)
  
  c_test=natives[:]->OMEGA500(startime:endtime,:)
  c_test=864*c_test  ; convert to hPa/d
  
  test_reshape=reshape(c_test,(/hourChunk,numperiods,ncols/))
  
  ; compute time average...
  test_dmn=dim_avg_n(test_reshape,0)
  pdf_test=pdfx(test_dmn,nbins,False)

; Calculate PDFs for CAM6 
  if (in .lt. 1) then
    pdf_d=pdf_test
  else if (in .lt. 2) then
    pdf_e=pdf_test 
  else
    pdf_f=pdf_test
  end if
  end if
  in=in+1  
end do

;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

print("begin to reshape arrays")

;c5a_reshape=new((/120,219,192,288/),float)
;c5a_reshape=reshape(play5a,(/120,219,192,288/))
c5a_reshape=new((/hourChunk,numperiods,192,288/),float)
c5b_reshape=new((/hourChunk,numperiods,192,288/),float)
c5c_reshape=new((/hourChunk,numperiods,192,288/),float)
c6c_reshape=new((/hourChunk,numperiods,192,288/),float)
c6b_reshape=new((/hourChunk,numperiods,192,288/),float)
c6a_reshape=new((/hourChunk,numperiods,192,288/),float)

c5a_reshape=reshape(c5a,(/hourChunk,numperiods,192,288/))
c5b_reshape=reshape(c5b,(/hourChunk,numperiods,192,288/))
c5c_reshape=reshape(c5c,(/hourChunk,numperiods,192,288/))
c6a_reshape=reshape(c6a,(/hourChunk,numperiods,192,288/))
c6b_reshape=reshape(c6b,(/hourChunk,numperiods,192,288/))
c6c_reshape=reshape(c6c,(/hourChunk,numperiods,192,288/))

;printVarSummary(c5a_reshape)

print("average arrays")

; compute time average.  The resulting arrays are then binned
; into pdfs for plotting....
c5a_5dmn=dim_avg_n(c5a_reshape,0)
c5b_5dmn=dim_avg_n(c5b_reshape,0)
c5c_5dmn=dim_avg_n(c5c_reshape,0)
c6a_5dmn=dim_avg_n(c6a_reshape,0)
c6b_5dmn=dim_avg_n(c6b_reshape,0)
c6c_5dmn=dim_avg_n(c6c_reshape,0)

printVarSummary(c5a_5dmn)

;-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-===-=-=-=-=--=-=-=-=-=
print("compute omega up and omega dn")

; where(cond_expr,true,false)
omega_u_5a=where((c5a_5dmn .le. 0),c5a_5dmn,c5a_5dmn@_FillValue)
omega_d_5a=where((c5a_5dmn .ge. 0),c5a_5dmn,c5a_5dmn@_FillValue)
omega_u_5b=where((c5b_5dmn .le. 0),c5b_5dmn,c5b_5dmn@_FillValue)
omega_d_5b=where((c5b_5dmn .ge. 0),c5b_5dmn,c5b_5dmn@_FillValue)
omega_u_5c=where((c5c_5dmn .le. 0),c5c_5dmn,c5c_5dmn@_FillValue)
omega_d_5c=where((c5c_5dmn .ge. 0),c5c_5dmn,c5c_5dmn@_FillValue)

;omega_u_a=dim_avg_n(omega_up,0)
;;print("dimsizes of a are: "+dimsizes(omega_u_a))
;omega_u_b=dim_avg_n(omega_u_a,0)
;;print("dimsizes of b are: "+dimsizes(omega_u_b))
;omega_u_c=dim_avg_n(omega_u_b,0)
;;print("dimsizes of c are: "+dimsizes(omega_u_c))
;print("mean omega up for CAM5 295 is: "+omega_u_c)
;omega_u_c5_295=omega_u_c

;omega_d_a=dim_avg_n(omega_dn,0)
;omega_d_b=dim_avg_n(omega_d_a,0)
;omega_d_c=dim_avg_n(omega_d_b,0)
;print("mean omega down for CAM5 295 is: "+omega_d_c)

omega_u_6a=where((c6a_5dmn .le. 0),c6a_5dmn,c6a_5dmn@_FillValue)
omega_d_6a=where((c6a_5dmn .ge. 0),c6a_5dmn,c6a_5dmn@_FillValue)
omega_u_6b=where((c6b_5dmn .le. 0),c6b_5dmn,c6b_5dmn@_FillValue)
omega_d_6b=where((c6b_5dmn .ge. 0),c6b_5dmn,c6b_5dmn@_FillValue)
omega_u_6c=where((c6c_5dmn .le. 0),c6c_5dmn,c6c_5dmn@_FillValue)
omega_d_6c=where((c6c_5dmn .ge. 0),c6c_5dmn,c6c_5dmn@_FillValue)

;omega_u_6_a=dim_avg_n(omega_6_up,0)
;omega_u_6_b=dim_avg_n(omega_u_6_a,0)
;omega_u_6_c=dim_avg_n(omega_u_6_b,0)
;print("mean omega up for CAM6 295 is: "+omega_u_6_c)
;
;omega_d_6_a=dim_avg_n(omega_6_dn,0)
;omega_d_6_b=dim_avg_n(omega_d_6_a,0)
;omega_d_6_c=dim_avg_n(omega_d_6_b,0)
;print("mean omega down for CAM6 295 is: "+omega_d_6_c)

omega_up   =new((/6,numperiods,192,288/),float)
omega_down =new((/6,numperiods,192,288/),float)
omega_u    =new((/6/),float)
omega_d    =new((/6/),float)


omega_up(0,:,:,:)=omega_u_5a
omega_up(1,:,:,:)=omega_u_5b
omega_up(2,:,:,:)=omega_u_5c
omega_up(3,:,:,:)=omega_u_6a
omega_up(4,:,:,:)=omega_u_6b
omega_up(5,:,:,:)=omega_u_6c

omega_down(0,:,:,:)=omega_d_5a
omega_down(1,:,:,:)=omega_d_5b
omega_down(2,:,:,:)=omega_d_5c
omega_down(3,:,:,:)=omega_d_6a
omega_down(4,:,:,:)=omega_d_6b
omega_down(5,:,:,:)=omega_d_6c

blue=omega_down(1,:,:,:)

;print("dimsizes of omega_down are: "+dimsizes(omega_down))
;print("dimsizes of blue are: "+dimsizes(blue))

do i=0,5
  omega_dn= omega_down(i,:,:,:)
  av_a=dim_avg_n(omega_dn,0)
  av_b=dim_avg_n(av_a,0)
  omega_d(i)=dim_avg_n(av_b,0)
  omega_uu= omega_up(i,:,:,:)
  av_a=dim_avg_n(omega_uu,0)
  av_b=dim_avg_n(av_a,0)
  omega_u(i)=dim_avg_n(av_b,0)
end do

down_min_up=omega_d-omega_u
down_plu_up=omega_d+omega_u

print("!@#$%@%^@#%$")
print("omega up: "+omega_u)
print("omega down: "+omega_d)
print("!@#$%@%^@#%$")
print("omega down for CAM5 300: "+omega_d(1))
print("Circulation Strength: omega down minus omega up")
print(down_min_up)
print("Circulation Strength: omega down plus omega up")
print(down_plu_up)
print("!@#$%@%^@#%$")

;;---------------------
;; make plots
;;---------------------
print("begin to make plots")

plot = new(1,"graphic")

pdf_c5a=pdfx(c5a_5dmn,nbins,False)
pdf_c5b=pdfx(c5b_5dmn,nbins,False)
pdf_c5c=pdfx(c5c_5dmn,nbins,False)
pdf_c5a_morebins=pdfx(c5a_5dmn,50,False)
;
;pdf_test=pdfx(test_dmn,nbins,False)
;
pdf_c6a=pdfx(c6a_5dmn,nbins,False)
pdf_c6b=pdfx(c6b_5dmn,nbins,False)
pdf_c6c=pdfx(c6c_5dmn,nbins,False)
pdf_c6a_morebins=pdfx(c6a_5dmn,50,False)

printVarSummary(pdf_c6a)
printVarSummary(pdf_test)
;
; colors = (/"black","darkgoldenrod","darkolivegreen3","steelblue","mediumorchid3","darkslategray4","cyan3","red"/)
  ;colors = (/"black","orangered","orangered3","orangered4"/)
  colors = (/"chocolate2","aquamarine3"/)

;
wks = gsn_open_wks("eps","testplot")

res                    = True
res@gsnDraw            = False
res@gsnFrame           = False
res@xyLineThicknessF   = 5.
res@tiYAxisString      = "PDF (%)"
res@tiXAxisString      = "circulation regime (hPa/d)"

res@trXMinF = -100.0
res@trXMaxF = 75.0
res@trYMinF = 0.0
res@trYMaxF = 18.0

res@tmXTOn              = False
res@tmXTBorderOn        = False
res@tmYROn              = False
res@tmYRBorderOn        = False
res@vpWidthF            = 0.9
res@vpHeightF           = 0.5
;res@vpXF                = .9
;res@vpYF                =.45
;res@gsnCenterString    = "default bins"
res@xyLineColors = colors(0)
res@xyDashPattern =  0
res@xyLineThicknessF = 2.
plot0_a    =  gsn_csm_xy (wks,pdf_c5a@bin_center, pdf_c5a, res)
res@xyLineColors = colors(0)
;res@xyDashPattern =  3
res@xyLineThicknessF = 4.
plot0_b    =  gsn_csm_xy (wks,pdf_c5b@bin_center, pdf_c5b, res)
res@xyLineColors = colors(0)
;res@xyDashPattern =  16
res@xyLineThicknessF = 6.
plot0_c    =  gsn_csm_xy (wks,pdf_c5c@bin_center, pdf_c5c, res)
res@xyLineColors = colors(1)
;res@xyDashPattern = 0 
res@xyLineThicknessF = 2.
plot0_d    =  gsn_csm_xy (wks,pdf_c6a@bin_center, pdf_c6a, res)
res@xyLineColors = colors(1)
;res@xyDashPattern = 3 
res@xyLineThicknessF = 4.
plot0_e    =  gsn_csm_xy (wks,pdf_c6b@bin_center, pdf_c6b, res)
res@xyLineColors = colors(1)
;res@xyDashPattern =  16
res@xyLineThicknessF = 6.
plot0_f    =  gsn_csm_xy (wks,pdf_c6c@bin_center, pdf_c6c, res)
;test native grid...
; pdf_a, pdf_b, and pdf_c are CAM5
res@xyLineColors = colors(0)
res@xyLineThicknessF = 2.
plot0_ng_a    =  gsn_csm_xy (wks,pdf_a@bin_center, pdf_a, res)
res@xyLineThicknessF = 4.
plot0_ng_b    =  gsn_csm_xy (wks,pdf_b@bin_center, pdf_b, res)
res@xyLineThicknessF = 6.
plot0_ng_c    =  gsn_csm_xy (wks,pdf_c@bin_center, pdf_c, res)
; pdf_d, pdf_e, and pdf_f are CAM6
res@xyLineColors = colors(1)
res@xyLineThicknessF = 2.
plot0_ng_d    =  gsn_csm_xy (wks,pdf_d@bin_center, pdf_d, res)
res@xyLineThicknessF = 4.
plot0_ng_e    =  gsn_csm_xy (wks,pdf_e@bin_center, pdf_e, res)
res@xyLineThicknessF = 6.
plot0_ng_f    =  gsn_csm_xy (wks,pdf_f@bin_center, pdf_f, res)

plotnative=1
; plot all lines
if (plotnative .lt. 1) then
  overlay(plot0_ng_a,plot0_a)
  overlay(plot0_ng_a,plot0_b)
  overlay(plot0_ng_a,plot0_c)
  overlay(plot0_ng_a,plot0_d)
  overlay(plot0_ng_a,plot0_e)
  overlay(plot0_ng_a,plot0_f)
  overlay(plot0_ng_a,plot0_ng_b)
  overlay(plot0_ng_a,plot0_ng_c)
  overlay(plot0_ng_a,plot0_ng_d)
  overlay(plot0_ng_a,plot0_ng_e)
  overlay(plot0_ng_a,plot0_ng_f)
else ; plot only data from the native grid
  overlay(plot0_ng_a,plot0_ng_b)
  overlay(plot0_ng_a,plot0_ng_c)
  overlay(plot0_ng_a,plot0_ng_d)
  overlay(plot0_ng_a,plot0_ng_e)
  overlay(plot0_ng_a,plot0_ng_f)
end if
;
plot(0)=plot0_ng_a

;
;resP    = True
;resP@gsnPanelMainString  = "Circulation Regime over last year"
;gsn_panel(wks,plot,(/1,1/),resP)

lgres                    = True
lgres@lgLineColors       = (/"chocolate2","aquamarine3"/)    ; colors for legend lines
lgres@vpWidthF           = 0.12                   ; width of legend (NDC)
lgres@vpHeightF          = 0.08                   ; height of legend (NDC)
lgres@lgPerimOn          = False                 ; turn off perimeter
lgres@lgDashIndexes      = (/0,0/)              ; dash indexes
lgres@lgLineLabelStrings = (/"",""/)             ; no labels in lines
lgres@lgLineThicknessF   = 8.0
lgres@lgLabelFontHeightF = 0.1
lgres@lgLabelFont        = 22
lgres@lgBoxMinorExtentF  = .2
gsn_legend_ndc(wks,2,(/" CAM5"," CAM6"/),.25,.7,lgres)

draw(plot(0))
 
frame(wks)

end
