pro temp_space, wl, temperatures, synth_spectra

;Create empty lists for the temperatures and their associated OH spectra
ts = list()
ss = list()

; For a range of physically probably OH temperatures
for T=150d, 250d, 1 do begin
	; Generate the synthetic OH spectrum
	synth_oh, [6500d,6620d], T, ohwl, ohint, width
	; Convolve it to the required wavelengths
	convolve_sp, ohwl, ohint, 0.6d, wl, sp
	; Append
	ts.add, T
	ss.add, sp
endfor

temperatures = ts.toarray()
synth_spectra = ss.toarray()

end
