PRO stuff
;
; Stuff for getting OH temperatures for Matthew Cooper (HODI)
;
;

h_setup
read_lut


; Wavelength calibration:

read_tim,'20/12/2022 22:04:30',19.0/60d,mjs0,time,dseq,/nophot

av=reform(total(dseq,1)/double(n_elements(time)),[1,512,512])
spectra,2,mjs0,time[0],av,sp
get_w,mjs0,2,wl

auto_wlcal,mjs0,time[0],av,2,7280d,7430d,p0,wout,/mid,ctype=1


; Specific time of interest:

read_tim,'20/12/2022 09:26:00',5.0/60d,mjs0,time,dseq,/nophot

av=reform(total(dseq,1)/double(n_elements(time)),[1,512,512])
spectra,2,mjs0,time[0],av,sp
get_w,mjs0,2,wl

; Maybe do it like this, but probably not:
; Fitting setup
;dt = 120.
;N2file = '$HDIR/N2spec_hwhm08.idl'
;Tlinelist = ['OH(8-3)P1(2)','OH(8-3)P1(3)','OH(8-3)P1(4)','OH(8-3)P1(5)']
;linefile = '$HDIR/input_fit_lines_Tnpanel_v3.dat'
;pnum=2
;dir='/home/dkw/Work/People/HODI/specfit_output/'
;fit_o_matic_evolv, wl,sp, status,error,params,fitsp,fitwl, .....

; Or use specfit - need to remember how
; See /home/dkw/Work/People/Natalija/specfit_wrap.pro
;specfit,mjs0+(indgen(3)*120),120d,Nspec,pnum,linefile,Tlinelist,N2file=N2file,dir=dir,/correct_banding

; Or retrieve_spec_params:
.r /stp/raid1/sif/idl/hsoft_gather/work_version/run_HiTIES_fit.pro
.r /stp/raid1/sif/idl/hsoft_gather/work_version/get_Einstein_OH.pro
.r /stp/raid1/sif/idl/hsoft_gather/work_version/get_Elevels_OH.pro
.r /stp/raid1/sif/idl/hsoft_gather/work_version/get_OHline_ratio.pro
;.r /home/dkw/Work/People/Natalija/retrieve_spec_params.pro
.r /stp/raid1/sif/idl/hsoft_gather/work_version/retrieve_spec_params.pro
.r /stp/raid1/sif/idl/astron/pro/readcol
.r /stp/raid1/sif/idl/astron/pro/remchar
.r /stp/raid1/sif/idl/astron/pro/gettok
.r /stp/raid1/sif/idl/astron/pro/strnumber
retrieve_spec_params, mjs0,time[0],av,1,pnum,linefile,Tlinelist,N2file=N2file,/plot,TOH,subTOH,TN2,IN2,PWV,params,subTerror,I_Op,errortag,V_P1,V_P2,V_Q

; ***********
; Get OH temperatures all day
; ***********
tt_mjs,2022,12,20,0,0,0,0,mjs
tt_mjs,2022,12,20,23,59,0,0,mjsend
N2file = '$HDIR/N2spec_hwhm08.idl'
Tlinelist = ['OH(8-3)P1(2)','OH(8-3)P1(3)','OH(8-3)P1(4)','OH(8-3)P1(5)','OH(8-3)P2(4)','OH(8-3)P2(2)','OH(8-3)P2(3)','OH(8-3)P2(5)']
linefile = '$HDIR/input_fit_lines_Tnpanel_v4.dat'
pnum=2
step=5*60d ; Spectrum integration time in seconds
outarr=[!values.F_NaN,!values.F_NaN,!values.F_NaN,!values.F_NaN]
.r
while mjs lt mjsend do begin
 datetime=strmid(mjs2str(mjs,/other),0,19)
 print,strmid(mjs2str(mjs),0,19)+': ',form='(a,$)'
 read_tim,datetime,step/3600d,mjs0,time,dseq,/nophot
 if n_elements(dseq) lt 512 then begin
   mjs+=step
   continue
 endif
 av=reform(total(dseq,1)/double(n_elements(time)),[1,512,512])
 retrieve_spec_params, mjs0,time[0],av,1,pnum,linefile,Tlinelist,N2file=N2file,TOH,subTOH,TN2,IN2,PWV,params,subTerror,I_Op,errortag,V_P1,V_P2,V_Q;,/plot
 print,'T = '+string(subTOH,form='(f5.1)')+' +/-'+string(subTerror,form='(f5.1)')+' K, PWV = '+string(max([PWV,0]),form='(f4.1)')+' mm'
 outarr=[[outarr],[mjs,subTOH,subTerror,PWV]]
 mjs+=step
endwhile
end


openw,1,'OH_T_20221220.txt'
for i=0,(n_elements(outarr)/4)-1 do if finite(outarr[0,i]) then printf,1,strmid(mjs2str(outarr[0,i]),0,19)+string(step,form='(i6)')+string(outarr[1,i],form='(f7.1)')+string(outarr[2,i],form='(f6.1)')
close,1


; Make some plots
read_tim,'20/12/2022 09:25:00',5.0/60d,mjs0,time,dseq,/nophot
av=reform(total(dseq,1)/double(n_elements(time)),[1,512,512])
N2file = '$HDIR/N2spec_hwhm08.idl'
Tlinelist = ['OH(8-3)P1(2)','OH(8-3)P1(3)','OH(8-3)P1(4)','OH(8-3)P1(5)','OH(8-3)P2(4)','OH(8-3)P2(2)','OH(8-3)P2(3)','OH(8-3)P2(5)']
linefile = '$HDIR/input_fit_lines_Tnpanel_v4.dat'
pnum=2
retrieve_spec_params, mjs0,time[0],av,1,pnum,linefile,Tlinelist,N2file=N2file,plot='20221220-092500',TOH,subTOH,TN2,IN2,PWV,params,subTerror,I_Op,errortag,V_P1,V_P2,V_Q



end
