pro sp_func_h, X, A, F
; Input function for fitting the spectrum in the H panel
; Inputs
;       X is wl (the wavelengths of the panel)
;       A is the guess free parameters A[0] = Temperature, A[1] = Intensity, A[2] = Background
; Outputs
;       A is the updated parameters
;       F is the function y values

;step = fltarr(402) + 1
;step[200:280] = 0

pnum = 3
get_wlrange, pnum, wls

synth_oh, wls, A[0], ohwl, ohint, width, upperv=6
convolve_sp, ohwl, ohint, 0.6d, X, sp

F = (A[1] * sp) + A[2]
;F = F * step

end


pro crop, vector, cropped

cropped = vector[70:401]

end


pro hpanel_fit, wl, sp, A, result
; Fitting function for the spectrum in the H panel
; Inputs
;       wl is wl (the wavelengths of the panel)
;       sp is the data to be fit
; Outputs
;       A is the fit parameters A[0] = Temperature, A[1] = Intensity, A[2] = Background
;       result is the fitted spectrum

A = [200d, 50d, 5d]
fita = [1, 1, 1]
weights = fltarr(402) + 1

;step = fltarr(402) + 1
;step[200:280] = 0

crop, wl, wlc
crop, sp, spc

;spc = spc*step

result = curvefit(wl, sp, weights, A, sigma, function_name='sp_func_h', /noderivative, itmax=100, /double, fita=fita)
;result = curvefit(wlc, spc, weights, A, sigma, function_name='sp_func_h', /noderivative, itmax=100, /double, fita=fita)

end