pro calibrate_panel, start_string, len, pnum, p0, wout

if pnum eq 1 then begin
     wls = [7920d, 8040d]
endif
if pnum eq 2 then begin
     wls = [7280d, 7430d] 
endif
if pnum eq 3 then begin
     wls = [6500d, 6620d] 
endif

read_tim, start_string, len/60d, mjs0, time, dseq, /nophot

av = reform(total(dseq, 1) / double(n_elements(time)), [1, 512, 512])
spectra, pnum, mjs0, time[0], av, sp
get_w, mjs0, pnum, wl

window, 10
plot, wl, sp, title='Reference'

auto_wlcal, mjs0, time[0], av, pnum, wls[0], wls[1], p0, wout, /mid, ctype=1, hydroxyl=2

end

; .r air2vac.pro, constant.pro, get_pwv_coeff.pro, hities_fit.pro, hities_stuff.pro, hitran_absorption.pro, hitran_intensity.pro, oh_label.pro, pregen_o2p_1n.pro, prog_n2_1p_v3_exp.pro, prog_o2_1n_v3_exp.pro, pwv_absorption.pro, read_hitran.pro, synth_n2_1p.pro, synth_o2p_1n.pro, synth_oh.pro, synth_op2p.pro, synth_panel_op.pro, vac2air.pro