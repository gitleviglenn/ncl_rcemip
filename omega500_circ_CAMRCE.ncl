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
timelength=240 ; 100 days
;timelength=8760 ; 1 year 
;timelength=17520 ; 2 years
endtime=26200;
startime=endtime-timelength+1

; number of bins used for the PDFs
;nbins=50

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

;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; import the data that is on the native grid...

ncols = 48602
sstarray=(/"295","300","305"/)

omega_up   =new((/6,numperiods,ncols/),float)
omega_down =new((/6,numperiods,ncols/),float)
omega_u    =new((/6/),float)
omega_d    =new((/6/),float)

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
  
  data_c5=natives[:]->OMEGA500(startime:endtime,:)
  data_c5=864*data_c5 ; convert to hPa/d
  
  data_c5_generic_reshape=new((/hourChunk,numperiods,ncols/),float)
  data_c5_generic_reshape=reshape(data_c5,(/hourChunk,numperiods,ncols/))

  c5_chunk_mn=dim_avg_n(data_c5_generic_reshape,0) ; average over chunks
  omega_u_5=where((c5_chunk_mn .le. 0),c5_chunk_mn,c5_chunk_mn@_FillValue)
  omega_d_5=where((c5_chunk_mn .ge. 0),c5_chunk_mn,c5_chunk_mn@_FillValue)

  if (in .lt. 1) then
      print("in is: "+in+" and sststring for CAM5 is: "+sststring)
    omega_up(0,:,:)    =omega_u_5
    omega_down(0,:,:)  =omega_d_5
  else if (in .lt. 2) then
      print("in is: "+in+" and sststring for CAM5 is: "+sststring)
    omega_up(1,:,:)    =omega_u_5
    omega_down(1,:,:)  =omega_d_5
  else
      print("in is: "+in+" and sststring for CAM5 is: "+sststring)
    omega_up(2,:,:)    =omega_u_5
    omega_down(2,:,:)  =omega_d_5
  end if
  end if

; CAM6 on the native grid
    print("begin work with CAM6")
  files_native  = systemfunc ("ls "+path_2_native+"/CAM6_GCM_RCE_large"+sststring+"/run/*cam.h1.*nc")
  natives       = addfiles (files_native,"r")
  ListSetType(natives,"cat")
    print("natives names are: "+files_native)
  
  data_c6=natives[:]->OMEGA500(startime:endtime,:)
  data_c6=864*data_c6  ; convert to hPa/d
  
  data_c6_generic_reshape=new((/hourChunk,numperiods,ncols/),float)
  data_c6_generic_reshape=reshape(data_c6,(/hourChunk,numperiods,ncols/))

;data_c6_generic_reshape=new((/hourChunk,numperiods,ncols/),float)
;data_c6_generic_reshape=reshape(data_c6a,(/hourChunk,numperiods,ncols/))
;; average over chunks....  should I average over chunks before calculating up down?
;c6a_chunk_mn=dim_avg_n(data_c6_reshape,0)
;printVarSummary(c6a_chunk_mn)
  
  ; compute time average...
;  test_dmn=dim_avg_n(test_reshape,0)
;  pdf_test=pdfx(test_dmn,nbins,False)

; Calculate PDFs for CAM6 
  if (in .lt. 1) then
      print("in is: "+in+" and sststring for CAM6 is: "+sststring)
    c6_chunk_mn=dim_avg_n(data_c6_generic_reshape,0) ; average over chunks
    omega_u_6=where((c6_chunk_mn .le. 0),c6_chunk_mn,c6_chunk_mn@_FillValue)
    omega_d_6=where((c6_chunk_mn .ge. 0),c6_chunk_mn,c6_chunk_mn@_FillValue)
    omega_up(3,:,:)    =omega_u_6
    omega_down(3,:,:)  =omega_d_6
  else if (in .lt. 2) then
      print("in is: "+in+" and sststring for CAM6 is: "+sststring)
    c6_chunk_mn=dim_avg_n(data_c6_generic_reshape,0) ; average over chunks
    omega_u_6=where((c6_chunk_mn .le. 0),c6_chunk_mn,c6_chunk_mn@_FillValue)
    omega_d_6=where((c6_chunk_mn .ge. 0),c6_chunk_mn,c6_chunk_mn@_FillValue)
    omega_up(4,:,:)    =omega_u_6
    omega_down(4,:,:)  =omega_d_6
  else
      print("in is: "+in+" and sststring for CAM6 is: "+sststring)
    c6_chunk_mn=dim_avg_n(data_c6_generic_reshape,0) ; average over chunks
    omega_u_6=where((c6_chunk_mn .le. 0),c6_chunk_mn,c6_chunk_mn@_FillValue)
    omega_d_6=where((c6_chunk_mn .ge. 0),c6_chunk_mn,c6_chunk_mn@_FillValue)
    omega_up(5,:,:)    =omega_u_6
    omega_down(5,:,:)  =omega_d_6
  end if
  end if
  in=in+1  
end do

  print("done with loops for now")

do i=0,5
  omega_dn= omega_down(i,:,:)
  av_a=dim_avg_n(omega_dn,0)    ; average in time
  av_b=dim_avg_n(av_a,0)        ; average over latitude
  omega_d(i)=dim_avg_n(av_b,0)  ; average over longitude
  omega_uu= omega_up(i,:,:)
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

;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

print("begin to calculate circulation strength")

;-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-===-=-=-=-=--=-=-=-=-=
print("compute omega up and omega dn")



end


