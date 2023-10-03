pro protonintensity, dseq, mjs0

get_w, mjs0, 3, wl

;Define my panel 3, x and y coordinates from panels.lut (idk how to use it)
x1=294 & y1=60 & nx=144 & ny=402
x2=x1+nx & y2=y1+ny
panel3 = dseq[*, x1:x2, y1:y2]

;Get panel intensity vs wavelength
svw = total(panel3, 2)
svw = svw/403

;Plot spectrum
window, 1
plot, wl, svw[1, *], yran=[358, 364]

;Integrate peak, from 210-270
int_peak = total(svw[*, 210:270], 2)
window, 2
plot, int_peak, yran=[2.2e4, 2.21e4] 

end

;--------------------------------------

pro autocorrtest, int_peak

lags = indgen(20)/2.
autocorr = a_correlate(int_peak, lags)

end



