pro spaceintegrate, mjs0, dseq, pno, space_integrated

get_w, mjs0, pno, wl
get_p, mjs0, dseq, pno, panel, /percentff
space_integrated = total(panel, 2)/144

end

;---------------------

pro wlintegrate, space_integrated, peakint

peakint = total(space_integrated[*,210:270], 1)

end
