;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; cam_2D_omega.ncl
;
; plot 2D maps of omega for CAM5 and CAM6
; the period over which the omega fields are averaged should be specified
;
; levi silvers                                                                 sep 2020
;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
;******************************************************
;  Open a netCDF file containing the grid and data from the HOMME
;  (High-Order Multiscale Modeling Environment) model, formerly
;  called SEAM (Spectral Element Atmosphere Model).
;******************************************************
  wks_type = "eps"
  wks_type@wkPaperWidthF  =  8.5  ; in inches
  wks_type@wkPaperHeightF = 11.0  ; in inches
  wks = gsn_open_wks(wks_type,"omega_comparison")
  ;gsn_define_colormap(wks,"gsltod")   ; choose colormap
  gsn_define_colormap(wks,"NCV_jaisnd")   ; choose colormap

  RCEMIP_dir = "/glade/scratch/levi/RCEMIP/"
;  plot_time = 1*24*365*3   ; this will need to be changed evenutally...
  plot_time = 26200; 
; set SST used with i: 
  temp  = (/295,300,305/)
  tempindex = 0
  model = (/"CAM5","CAM6"/)
  ;chunklength=(/"1 d chunk","5 d chunks","10 d chunks"/) 
  chunklength=(/"1 d chunk","5 d chunks","30 d chunks"/) 

  ntemp = dimsizes(temp)
  nmodel= dimsizes(model)
  plot1 = new(ntemp*nmodel,graphic)
  plot2 = new(ntemp*nmodel,graphic)

;  filein_c5 = addfile(RCEMIP_dir+model(0)+"_GCM/RCE_large"+temp(tempindex)+"/2D/"+model(0)+"_GCM_RCE_large"+temp(tempindex)+"_2D_wap500.nc","r")
;  filein_c6 = addfile(RCEMIP_dir+model(1)+"_GCM/RCE_large"+temp(tempindex)+"/2D/"+model(1)+"_GCM_RCE_large"+temp(tempindex)+"_2D_wap500.nc","r")
;
;  testdims=filein_c5->wap500(plot_time,:,:)
;  testdims2=filein_c6->wap500(plot_time,:,:)
;  ;printVarSummary(testdims)
;  ;printVarSummary(testdims2)
;  nlat = filein_c5->lat
;  nlon = filein_c5->lon
;  ;omega_array = new((/3,192,288/),"float") 
 
;  do i = 0,ntemp-1 
    do j = 0,nmodel-1
  omega_array = new((/3,192,288/),"float") 
    print("Creating plots for temperature = "+temp(tempindex))
      ;file_omega = addfile(RCEMIP_dir+model(j)+"_GCM/RCE_large"+temp(i)+"/2D/"+model(j)+"_GCM_RCE_large"+temp(i)+"_2D_wap500.nc","r")

      file_omega = addfile(RCEMIP_dir+model(j)+"_GCM/RCE_large"+temp(tempindex)+"/2D/"+model(j)+"_GCM_RCE_large"+temp(tempindex)+"_2D_wap500.nc","r")
  
      omega      = file_omega->wap500(plot_time,:,:)
      omega_full = file_omega->wap500(:,:,:)
      omega_5d   = omega_full(26081:26200,:,:)
      omega_10d  = omega_full(25960:26200,:,:)
      omega_30d  = omega_full(25480:26200,:,:)
      time       = file_omega->time(plot_time)
    
      omega_5dmn  = dim_avg_n(omega_5d,0)
      omega_10dmn = dim_avg_n(omega_10d,0)
      omega_30dmn = dim_avg_n(omega_30d,0)

      copy_VarMeta(omega,omega_5dmn)
      copy_VarMeta(omega,omega_10dmn)
      copy_VarMeta(omega,omega_30dmn)

    print("dimensions of omega_full are: "+dimsizes(omega_full))
    print("dimensions of omega are: "+dimsizes(omega))
    print("dimensions of omega_5d are: "+dimsizes(omega_5d))
    print("dimensions of omega_5dmn are: "+dimsizes(omega_5dmn))
    print("dimensions of omega_10dmn are: "+dimsizes(omega_10dmn))
    print("dimensions of omega_array are: "+dimsizes(omega_array))
    print("Plotting for time = "+time)

;    printVarSummary(omega)
;    printVarSummary(omega_array)

      omega_array(0,:,:)=omega(:,:)
      omega_array(1,:,:)=omega_5dmn(:,:)  
      omega_array(2,:,:)=omega_30dmn(:,:)  

      omega_array!0 = "ChunkPeriod"
      omega_array!1 = "lat"
      omega_array!2 = "lon"
      omega_array&lat = omega&lat    
      omega_array&lon = omega&lon    
;      copy_VarMeta(omega,omega_array)

      ;prec  = prec*3.6e6/1000
      ;prec@units = "mm/hr"
  
      omega_array = 864*omega_array ; should scale Pa/s to hPa/day
  
    print("min omega: "+min(omega))
    print("max omega: "+max(omega))
    
      res1                      = True         ; plot modifications desired
      res1@gsnDraw              = False        ; don't draw plot
      res1@gsnFrame             = False        ; don't advance frame
  
      res1@cnLevelSelectionMode = "manualLevels"    ; set contour intervals
      res1@cnLevelSpacingF      =  10.                ;
      res1@cnMinLevelValF       =  -150.0           ;
      res1@cnMaxLevelValF       =  50.0 ;
      res1@cnFillOn             = True         ; Turn on contour fill
      res1@cnFillMode           = "RasterFill"
      res1@cnLinesOn            = False        ; Turn off contour lines
      res1@cnLineLabelsOn       = False        ; Turn off contour line labels
      res1@lbLabelBarOn         = False
  
      res1@mpGridAndLimbOn   = True           ; turn on lat/lon lines
      res1@mpGridLatSpacingF = 30.
      res1@mpGridLonSpacingF = 30.
      res1@mpPerimOn         = False               ; turn off perimeter
      res1@mpFillOn         = False
      res1@mpOutlineBoundarySets     = "NoBoundaries"
      res1@mpProjection = "Mollweide"
  
      res1@gsnRightString  = ""
      res1@gsnLeftString   = ""
  
      res1@tiYAxisFontHeightF    = 0.025
      res1@tiYAxisFont           = 22   
      res1@tiYAxisSide           = "Left"
 
  ;chunklength=(/"1 d chunk","5 d chunks", 10 d chunks"/) 
  do ii = 0,ntemp-1 
      if (j .eq. 0) then
        ;res1@tiYAxisString       = temp(i)+" K"
        res1@tiYAxisString       = chunklength(ii)
      else
        res1@tiYAxisString       = ""  
      end if
  
      if (ii .eq. 0) then
        res1@gsnCenterString        = model(j) 
        res1@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
        res1@gsnCenterStringFontHeightF = 0.025              ; ditto
      else
        res1@gsnCenterString        = ""
      end if
  
      ;plot1(2*i+j)   = gsn_csm_contour_map(wks,omega(:,:),res1)
      plot1(2*ii+j)   = gsn_csm_contour_map(wks,omega_array(ii,:,:),res1)
;      plot1(j)     = gsn_csm_contour_map(wks,omega(:,:),res1)
;      plot1(2+j)   = gsn_csm_contour_map(wks,omega_array(i,:,:),res1)
;      plot1(4+j)   = gsn_csm_contour_map(wks,omega_array(i,:,:),res1)
  
    end do 
    delete(omega_array)
    delete(res1)
  end do

 panelres1                            = True
 panelres1@gsnFrame                   = False    ; save panel until both ready 
 panelres1@gsnPanelBottom             = 0.2			 ; shrink panel plot by setting bottom edge of plot
 panelres1@gsnPanelTop                = 0.9			 ; shrink panel plot by setting top edge of plot
 panelres1@gsnPanelLabelBar           = True
 panelres1@vpXF                       = 0.518-.22         ; Position labelbar at
 panelres1@vpHeightF                  = 0.095     ; Height of labelbar
 panelres1@vpWidthF                   = 0.44      ; Width of labelbar
 panelres1@lbLabelFontHeightF         = 0.01
 panelres1@lbLabelStride              = 2
 panelres1@txPosXF                    = 0.72
 panelres1@txPosYF                    = 0.925
 panelres1@txFont                     = 22 
 ;panelres1@lbTitleString              = "Outgoing Longwave Radiation (W/m~S~2~N~)"
 panelres1@lbTitleString              = "Instantaneous Pressure Velocity (hPa/day) "
 panelres1@lbTitlePosition            = "Top"
 panelres1@lbTitleFontHeightF         = 16
 panelres1@lbTitleFont                = 22
 panelres1@pmLabelBarOrthogonalPosF   = -.01
 gsn_panel(wks,plot1, (/ntemp,nmodel/),panelres1)


  frame(wks)

end

