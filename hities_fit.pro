pro hities_fitting.pro

h_setup
read_lut

read_tim,'20/12/2022 09:26:00',5.0/60d,mjs0,time,dseq,/nophot

av=reform(total(dseq,1)/double(n_elements(time)),[1,512,512])
spectra,2,mjs0,time[0],av,sp
get_w,mjs0,2,wl

tt_mjs,2022,12,20,0,0,0,0,mjs
tt_mjs,2022,12,20,23,59,0,0,mjsend
N2file = '$HDIR/N2spec_hwhm08.idl'
Tlinelist = ['OH(8-3)P1(2)','OH(8-3)P1(3)','OH(8-3)P1(4)','OH(8-3)P1(5)','OH(8-3)P2(4)','OH(8-3)P2(2)','OH(8-3)P2(3)','OH(8-3)P2(5)']
linefile = '$HDIR/input_fit_lines_Tnpanel_v4.dat'
pnum=2
step=5*60d ; Spectrum integration time in seconds
outarr=[!values.F_NaN,!values.F_NaN,!values.F_NaN,!values.F_NaN]

specfit,mjs0+(indgen(3)*120),120d,Nspec,pnum,linefile,Tlinelist,N2file=N2file,dir=dir,/correct_banding

end
