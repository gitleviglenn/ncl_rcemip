;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; omega500_circPDF_1SST.ncl
;
; levi silvers                                                           sep 2020
;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 
begin

; define indices for time steps read and used...
; unites for all of these are hours...
;timelength=24479 ; all but the first 75 days
;timelength=2880 ; 4 months
;timelength=2400 ; 100 days
;timelength=8760 ; 1 year 
;timelength=17520 ; 2 years
endtime=26200;
;endtime=25000 ; ICON GCM has fewer timesteps that CAM or ECHAM
;echam has 26280
;endtime=26279 ; this is the end time for CAM5 and CAM6

;startime=endtime-timelength+1

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

; look at 5 day chunks (for files that have hourly output)
hourChunk=120
;hourChunk=24

;numperiods=timelength/hourChunk
;print("numperiods ="+numperiods)

;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; import the data that is on the native grid...

ncols = 48602
sstarray=(/"295","300","305"/)
; 100days, 1 year, 2 years
;firsttime_array=(/25961,17441,8681/)

; 2400 -> 100 days
; 8760 -> 365 days
; 17520 -> 730 days

timelength_array=(/2400,8760,17520/)  ; probably needs to be exactly divisible by hourChunk

stime = new((/3/),"integer")
numperiods = new((/3/),"integer")
;omega_test = new((/2400,ncols/),"float")

in = 1
;do while(in .le. 2) ; loop over SST
  sststring=sstarray(in)
  files_native  = systemfunc ("ls "+path_2_native+"/CAM5_GCM_RCE_large"+sststring+"/run/*cam.h1.*nc")
  natives       = addfiles (files_native,"r")
  ListSetType(natives,"cat")
   ; print("natives names are: "+files_native)
  
  date = natives[:]->date
;end do

index=0 
do while(index .le. 2) ; loop over start time
  
  omega_test = new((/timelength_array(index),ncols/),"float")
  
  numperiods(index)=timelength_array(index)/hourChunk

  stime(index)=endtime-timelength_array(index)+1

  ;omega_test=natives[:]->OMEGA500(stime(index):endtime,:)
  omega_in=natives[:]->OMEGA500(:,:)
    print("dimensions of incoming omega are: "+dimsizes(omega_in))
  omega_test=omega_in(stime(index):endtime,:)
  omega_test=864*omega_test  ; convert to hPa/d

  ;numperiods(index)=timelength(index)/hourChunk
    print("numperiods ="+numperiods)
    print("stime ="+stime(index))
    print("dimensions of omega_test are: "+dimsizes(omega_test))
    print("===========================================")
    print("numperiod array is ="+numperiods(0)+" and "+numperiods(1)+" and "+numperiods(2))
    print("===========================================")

  perioddim=numperiods(index)
  test_reshape=new((/hourChunk,perioddim,ncols/),float)
    print("dimensions of test_reshape are: "+dimsizes(test_reshape))
  
  test_reshape=reshape(omega_test,(/hourChunk,perioddim,ncols/))
  
  ; compute time average...
  test_dmn=dim_avg_n(test_reshape,0)
  pdf_test=pdfx(test_dmn,nbins,False)

; Calculate PDFs 
  if (index .lt. 1) then    ; 100 days
    print("index is: "+index)
    pdf_a=pdf_test
  else if (index .lt. 2) then  ; 365 days
    print("index is: "+index)
    pdf_b=pdf_test 
  else                         ; 730 days
    print("index is: "+index)
    pdf_c=pdf_test
  end if
  end if
  index=index+1
  delete(test_dmn)
  delete(omega_test)
  delete(test_reshape)
end do

;;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

;;---------------------
;; make plots
;;---------------------
print("begin to make plots")

plot = new(1,"graphic")

printVarSummary(pdf_test)
;
 colors = (/"black","darkgoldenrod","darkolivegreen3","steelblue","mediumorchid3","darkslategray4","cyan3","red"/)
  ;colors = (/"black","orangered","orangered3","orangered4"/)
;  colors = (/"chocolate2","aquamarine3"/)

;
wks = gsn_open_wks("eps","testplot_circ")

res                    = True
res@gsnDraw            = False
res@gsnFrame           = False
res@xyLineThicknessF   = 5.
res@tiMainString       = "5 days Chunks, SST = 300K"
res@tiYAxisString      = "Probability (%)"
res@tiXAxisString      = "Circulation Regime (hPa/day)"

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
res@xyLineColors = colors(0)
res@xyDashPattern =  0
res@xyLineThicknessF = 2.
res@xyLineColors = colors(0)
;;res@xyDashPattern =  3
;res@xyLineThicknessF = 4.
;res@xyLineColors = colors(0)
;;res@xyDashPattern =  16
;res@xyLineThicknessF = 6.
;res@xyLineColors = colors(1)
;;res@xyDashPattern = 0 
;res@xyLineThicknessF = 2.
;res@xyLineColors = colors(1)
;res@xyLineThicknessF = 4.
;res@xyLineColors = colors(1)
;res@xyLineThicknessF = 6.
;native grid...
res@xyLineColors = colors(0)
res@xyLineThicknessF = 4.
plot0_ng_a    =  gsn_csm_xy (wks,pdf_a@bin_center, pdf_a, res)
res@xyLineColors = colors(1)
plot0_ng_b    =  gsn_csm_xy (wks,pdf_b@bin_center, pdf_b, res)
res@xyLineColors = colors(2)
plot0_ng_c    =  gsn_csm_xy (wks,pdf_c@bin_center, pdf_c, res)

plotnative=1
if (plotnative .lt. 1) then
else ; plot only data from the native grid
  overlay(plot0_ng_a,plot0_ng_b)
  overlay(plot0_ng_a,plot0_ng_c)
end if
;
plot(0)=plot0_ng_a

;
;resP    = True
;resP@gsnPanelMainString  = "5 days Chunks, SST = 295K"
;gsn_panel(wks,plot,(/1,1/),resP)

lgres                    = True
lgres@lgLineColors       = (/"black","darkgoldenrod","darkolivegreen3"/)    ; colors for legend lines
lgres@vpWidthF           = 0.12                   ; width of legend (NDC)
lgres@vpHeightF          = 0.08                   ; height of legend (NDC)
lgres@lgPerimOn          = False                 ; turn off perimeter
lgres@lgDashIndexes      = (/0,0,0/)              ; dash indexes
lgres@lgLineLabelStrings = (/"","",""/)             ; no labels in lines
lgres@lgLineThicknessF   = 8.0
lgres@lgLabelFontHeightF = 0.2
lgres@lgLabelFont        = 22
lgres@lgBoxMinorExtentF  = .2
gsn_legend_ndc(wks,3,(/" 100 days "," 365 days"," 730 days"/),.25,.7,lgres)

draw(plot(0))

frame(wks)

end
