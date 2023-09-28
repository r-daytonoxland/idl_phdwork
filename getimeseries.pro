pro get_time_series, t1, length, peakint

print, 'Reading data...'
read_tim, t1, length, mjs0, time, dseq, icount, /nophot

panel=3
get_p, mjs0, dseq, 3, panel, /percentff

space_integrated = total(panel, 2)/144
peakint = total(space_integrated[*,210:270], 2)
peakint = peakint/120.

end
