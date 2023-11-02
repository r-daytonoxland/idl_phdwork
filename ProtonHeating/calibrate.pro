pro calibrate_panel, start_string, pnum

if pnum = 1 then wls = [7920d, 8040d] endif
if pnum = 2 then wls = [7280d, 7430d] endif
if pnum = 3 then wls = [6500d, 6620d] endif

read_tim, start_string, 19.0/60d, mjs0, time, dseq, /nophot

av = reform(total(dseq, 1) / double(n_elements(time)), [1, 512, 512])
spectra, pnum, mjs0, time[0], av, sp
get_w, mjs0, pnum, wl

auto_wlcal, mjs0, time[0], av, pnum, wls[0], wls[1], p0, wout, /mid, ctype=1, hydroxyl=2

end