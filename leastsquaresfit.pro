pro leastsquaresfit, synth_spectra, data, arange, brange, besta, bestb, bests
; Fit OH spectrum to H-alpha panel

bestfit = 1

step = fltarr(402) + 1
step[200:280] = 0

a = size(synth_spectra)
len = a[1]

for a=arange[0], arange[1], arange[2] do begin
	for b=brange[0], brange[1], brange[2] do begin
		for s=0, len-1 do begin

			model = synth_spectra[s, *]
			oh = (b * model * step)
			tidy = (data - a)*step
			lsf = total((tidy - oh)^2)

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
