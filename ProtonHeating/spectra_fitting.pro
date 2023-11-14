pro temp_space, pnum, wl, temperatures, synth_spectra, upperv=upperv

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


pro leastsquaresfit, synth_spectra, data, arange, brange, besta, bestb, bests, halpha=halpha
; Fit OH spectrum to H-alpha panel
; Inputs
;	synth_spectra : Output from temp_space with a big list of possible synthetic spectra
;   data : teh spectrum you're trying to fit
; 	arange : range of test a's, where a is the y background
;	brange : ditto, where b is the signal amplitude
; Outputs
;	besta : Fitting output parameter (background)
;   bestb : Ditto (signal amplitude)
;   bests : Index of the best fit synthetic spectrum
; Keywords
;   halpha : Removes the proton peak from the halpha panel before fitting by multiplying both the data and model by a step function

bestfit = 1

if keywordset(halpha)
	step = fltarr(402) + 1
	step[200:280] = 0
endif

a = size(synth_spectra)
len = a[1]

for a=arange[0], arange[1], arange[2] do begin
	for b=brange[0], brange[1], brange[2] do begin
		for s=0, len-1 do begin

			model = synth_spectra[s, *]

			if keywordset(halpha) then begin
				 oh = (b * model * step)
				 tidy = (data - a)*step
				 lsf = total((tidy - oh)^2)
			endif else begin
				 oh = (b * model)
				 tidy = data - a
			 	 lsf = total((tidy - oh)^2)
			endelse

			if lsf lt bestfit then begin
				bestfit = lsf
				besta = a
				bestb = b
				bests = s
			endif
		endfor
	endfor
endfor

end
