PRO synth_panel_op, wl, param, spectrum
;
; Make a synthetic spectrum for the "O+ panel" of HiTIES, including:
; - OH (8-3) emission at a single temperature
; - O+ 2P emission with a single density ratio
; - N2 1P emission at a single temperature
; - O2+ 1N emission at a single temperature
; - Water vapour absorption
; - Uniform background
;
; Inputs:
;  wl    - Desired wavelength grid points, in A
;  param - A vector of parameters defining the spectrum, with these elements:
;       0: Background (R/A)
;       1: Instrument function, FWHM of a Gaussian (A)
;       2: PWV (mm)
;       3: OH intensity, sum of all lines within range of wl (R)
;       4: OH temperature (K)
;       5: O+ intensity, sum of all lines within range (R)
;       6: O+ density ratio of 2P1/2 to 2P3/2 states
;       7: N2 intensity, sum of all lines within range (R)
;       8: N2 temperature (K)
;       9: *Reserved for future use
;      10: *Reserved for future use
;      11: O2+ intensity, sum of all lines within range (R)
;      12: O2+ temperature (K)
;    O2+ is optional - if param has less than 12 elements O2+ is ignored.
;
; Outputs:
;  spectrum - Intensity of the spectrum in R/A, at each of the wavelengths in
;             the wl array.
;


if n_elements(param) lt 9 or n_elements(param) gt 13 then message,'Error: incorrect number of parameters'

; The wavelength range
range=[min(wl),max(wl)]

synth_oh, range, param[4], wlOH, intOH
synth_op2p, range, param[6], wlOp, intOp
synth_n2_1p, range, param[8], wlN2, intN2
if n_elements(param) gt 11 then begin
 synth_o2p_1n, range, param[12], wlO2p, intO2p
 intO2p*=param[11]
endif else begin
 wlO2p=min(wl)-100 ; Set something reasonable but outside the wl range,
 intO2p=0          ; intensity is zero anyway
endelse

lines=[wlOH, wlOp, wlN2, wlO2p]
ints=[intOH*param[3], intOp*param[5], intN2*param[7], intO2p]

; Apply water vapour absorption
if (param[2] ne 0) then ints=ints*pwv_absorption(param[2],lines)

; Convolve with intrument function
convolve_sp, lines,ints, param[1]*0.5, wl, spectrum

; Add background
spectrum+=param[0]

;return,spectrum

end
