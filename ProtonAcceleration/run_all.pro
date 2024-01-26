pro runall

get_timings, 'JUSTprotons.csv', years, months, days, hours, mins, secs, msecs, lens
spectrathing, years, months, days, hours, mins, secs, msecs, lens, spectrum_protons

save, spectrum_protons, filename = 'spectrumprotons.sav'

get_timings, 'pc1ANDprotons.csv', years, months, days, hours, mins, secs, msecs, lens
spectrathing, years, months, days, hours, mins, secs, msecs, lens, spectrum_pc1

save, spectrum_pc1, filename = 'spectrumpc1.sav'

end