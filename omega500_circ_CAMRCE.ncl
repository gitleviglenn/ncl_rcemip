;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; omega500_circ_CAMRCE.ncl
;
; script to compute the intensity of the overturning circulation
;
; levi silvers                                                           may 2020
;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

; define indices for time steps read and used...
;timelength=240 ; 100 days
timelength=8760 ; 1 year 
;timelength=17520 ; 2 years
endtime=26200;
startime=endtime-timelength+1

; number of bins used for the PDFs
nbins=50

; period over which averaging occurs before pdf is computed
hourChunk=120
numperiods=timelength/hourChunk
print("numperiods ="+numperiods)

; define paths to data
path_laptop="/Users/silvers/data/CAM/"
path_C5_lt="CAM5_GCM_RCE/"
path_C6_lt="CAM6_GCM_RCE/"

path_cheyenne="/glade/scratch/levi/RCEMIP/"
path_2_native="/glade/scratch/kareed/RCEMIP/"
path_C5="CAM5_GCM/";RCE_large305/2D/"
path_C6="CAM6_GCM/";RCE_large305/2D/"

file_omega_5_295=path_cheyenne+path_C5+"RCE_large295/2D/CAM5_GCM_RCE_large295_2D_wap500.nc"
file_omega_6_295=path_cheyenne+path_C6+"RCE_large295/2D/CAM6_GCM_RCE_large295_2D_wap500.nc"
file_omega_5_300=path_cheyenne+path_C5+"RCE_large300/2D/CAM5_GCM_RCE_large300_2D_wap500.nc"
file_omega_6_300=path_cheyenne+path_C6+"RCE_large300/2D/CAM6_GCM_RCE_large300_2D_wap500.nc"
file_omega_5_305=path_cheyenne+path_C5+"RCE_large305/2D/CAM5_GCM_RCE_large305_2D_wap500.nc"
file_omega_6_305=path_cheyenne+path_C6+"RCE_large305/2D/CAM6_GCM_RCE_large305_2D_wap500.nc"


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
;  ntim = natives[:]->time
  
  c5a_test=natives[:]->OMEGA500(startime:endtime,:)
  c5a_test=864*c5a_test  ; convert to hPa/d
  
  test_reshape=reshape(c5a_test,(/hourChunk,numperiods,ncols/))
  
  ; compute time average...
  test_dmn=dim_avg_n(test_reshape,0)
;  pdf_test=pdfx(test_dmn,nbins,False)

; Calculate PDFs for CAM5 
  if (in .lt. 1) then
    print("files being used: "+files_native)
;    pdf_a=pdf_test
  else if (in .lt. 2) then
    print("files being used: "+files_native)
;    pdf_b=pdf_test 
  else
    print("files being used: "+files_native)
;    pdf_c=pdf_test
  end if
  end if

; CAM6 on the native grid
  files_native  = systemfunc ("ls "+path_2_native+"/CAM6_GCM_RCE_large"+sststring+"/run/*cam.h1.*nc")
  natives       = addfiles (files_native,"r")
  ListSetType(natives,"cat")
    print("natives names are: "+files_native)
  
  data_c6=natives[:]->OMEGA500(startime:endtime,:)
  data_c6=864*data_c6  ; convert to hPa/d
  
  test_reshape=reshape(data_c6,(/hourChunk,numperiods,ncols/))
  
  ; compute time average...
;  test_dmn=dim_avg_n(test_reshape,0)
;  pdf_test=pdfx(test_dmn,nbins,False)

; Calculate PDFs for CAM6 
  if (in .lt. 1) then
    print("files being used: "+files_native)
;    pdf_d=pdf_test
  else if (in .lt. 2) then
    print("files being used: "+files_native)
;    pdf_e=pdf_test 
  else
    print("files being used: "+files_native)
;    pdf_f=pdf_test
  end if
  end if
  in=in+1  
end do

  print("done with loops for now")

;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

print("begin to calculate circulation strength")

; read omega
; calculate global weighting if not using the native grid
; separate omega into regions of upwelling and downwelling
; compute average of these regions
; compute difference of the averages

;c5a_reshape=new((/120,219,192,288/),float)

c5a_reshape=new((/120,219,ncols/),float)
c5a_reshape=reshape(play5a,(/hourChunk,numperiods,ncols/))
c5a_5dmn=dim_avg_n(c5a_reshape,0)
printVarSummary(c5a_5dmn)
; where(cond_expr,true,false)
omega_u_5a=where((c5a_5dmn .le. 0),c5a_5dmn,c5a_5dmn@_FillValue)
omega_d_5a=where((c5a_5dmn .ge. 0),c5a_5dmn,c5a_5dmn@_FillValue)

omega_up   =new((/6,numperiods,192,288/),float)
omega_down =new((/6,numperiods,192,288/),float)
omega_u    =new((/6/),float)
omega_d    =new((/6/),float)

omega_up(0,:,:,:)=omega_u_5a
omega_down(0,:,:,:)=omega_d_5a

data_c6_reshape=new((/hourChunk,numperiods,ncols/),float)
data_c6_reshape=reshape(data_c6,(/120,219,ncols/))

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




end


