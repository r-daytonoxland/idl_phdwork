.r air2vac.pro, constant.pro, oh_label.pro, pregen_o2p_1n.pro, prog_n2_1p_v3_exp.pro, prog_o2_1n_v3_exp.pro, synth_n2_1p.pro, synth_o2p_1n.pro, synth_oh.pro, synth_op2p.pro, synth_panel_op.pro, vac2air.pro
.r ../ProtonHeating/proton_heating

; HiTIES wavelength calibration
read_tim, '03/01/2020 20:00:00', 1., mjs0, time, dseq, icount, /nophot
av=reform(total(dseq,1)/double(n_elements(time)),[1,512,512])
spectra,3,mjs0,time[0],av,sp
get_w,mjs0,3,wl
plot,wl,sp
auto_wlcal,mjs0,time[0],av,3,6490d, 6660d,p0,wout,/mid,ctype=1,hydroxyl=2

; Mke sure you ignore the H geocorona!
wout = 1.0000000       201.00000       6556.9836     -0.34350718   0.00020035277  -1.8392371e-06

;Add H rest wavelength
oplot, fltarr(2) + 6563, !y.crange

;Continuum spectrum
;timeframe 55-62 ON
;timeframe 25-32 OFF
read_tim, '03/01/2020 09:10:00', 10/60., mjs0, time, dseq, icount, /nophot
proton_keogram, '03/01/2020 09:10:00', 10/60., keo, 'blah.png'
on = reform(total(keo[55:62,*]/8., 1),[1,402])
off = reform(total(keo[25:32,*]/8., 1),[1,402])
get_w, mjs0, 3, wl
plot, wl, on, yrange=[1000, 1006]
oplot, wl, off

;Do it using spectra
read_tim, '03/01/2020 09:10:00', 10/60., mjs0, time, dseq, icount, /nophot

on = dseq[550:620,*,*]
on_dseq = reform(total(on,1)/80.,[1,512,512])
off = dseq[250:320,*,*]
off_dseq = reform(total(off,1)/80.,[1,512,512])

spectra, 3, mjs0, time, on_dseq, sp_on
spectra, 3, mjs0, time, off_dseq, sp_off
get_w, mjs0, 3, wl