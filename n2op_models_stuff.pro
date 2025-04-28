; Using synth_n1_1p.pro and synth_op2p.pro

.r synth_n2_1p.pro, synth_op2p.pro, pregen_n2.pro, prog_n2_1p_v3_exp.pro, constant.pro
.r ProtonHeating/proton_heating

; Oplus in panel 2

read_tim, '30/12/2021 08:00:00', 1/60., mjs0, time, dseq, icount, /nophot, tadd=60
get_w, mjs0, 1, wl1
get_w, mjs0, 2, wl2
get_w, mjs0, 3, wl3

pnum=2
get_wlrange, pnum, range

; What on earth is T? Ratio of densities for states? What is the usual one?
T = 0.367 ; From Dan but I should try fitting it
synth_op2p, range, T, wl, int, width, labels, weights=weights
convolve_sp, wl, int, 0.6d, wl2, sp2

save, filename='op_model_30122021.sav', wl2, sp2

; N2 in all panels
T = 300.

pnum=1
get_wlrange, pnum, range
synth_n2_1p, range, T, wl, int, width, labels, weights=weights
convolve_sp, wl, int, 0.6d, wl1, sp1

pnum=2
get_wlrange, pnum, range
synth_n2_1p, range, T, wl, int, width, labels, weights=weights
convolve_sp, wl, int, 0.6d, wl2, sp2

pnum=3
get_wlrange, pnum, range
synth_n2_1p, range, T, wl, int, width, labels, weights=weights
convolve_sp, wl, int, 0.6d, wl3, sp3

save, filename='n21p_model_3012021.sav', wl1, wl2, wl3, sp1, sp2, sp3