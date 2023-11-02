pro temp_space, pnum, wl, temperatures, synth_spectra

;Create empty lists for the temperatures and their associated OH spectra
ts = list()
ss = list()

get_wlrange, pnum, wls

; For a range of physically probably OH temperatures
for T=150d, 250d, 1 do begin
	; Generate the synthetic OH spectrum
	synth_oh, wls, T, ohwl, ohint, width, upperv=upperv
	; Convolve it to the required wavelengths
	convolve_sp, ohwl, ohint, 0.6d, wl, sp
	; Append
	ts.add, T
	ss.add, sp
endfor

temperatures = ts.toarray()
synth_spectra = ss.toarray()

end
