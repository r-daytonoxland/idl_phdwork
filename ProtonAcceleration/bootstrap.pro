pro bootstrap, space_integrated, lags, noperts, max

; Peakint is the intensity of the proton peak at each second -120
peakint = total(space_integrated[*,210:270], 2)
peakint = peakint/120. ; mean

; Set up max amplitides
maxamps = fltarr(noperts, 60)

; Iterate over noperts
for i = 0,(noperts-1) do begin
	y = randomu(dseed, 120) ; Generate random order
	x = peakint[sort(y)]    ; Sort times into that order
	acorr = a_correlate(x, lags) ; Do autocorrelation
	maxamps[i, *] = acorr   ; Save it in maxamps
endfor

max = max(abs(maxamps), dimension=1)

plot, lags/2., max

end

