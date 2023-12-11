FUNCTION pwv_absorption, pwv, wl;, width
;
; Calculates the absorption of an emission line due to water vapour.
; Actually returns the *transmission*, i.e. if there is no absorption
; returns 1, if there is 100% absorption returns 0.
;
; Inputs:
;  pwv   - Amount of precipitable water vapour, in mm.
;          Currently must be scalar.
;  wl    - Wavelength of the emission, in A. Scalar or array.
;  width - The width of the line due to Doppler broadening,
;          as FWHM in A. This is optional, and usually not needed
;          except for light emitting species at high temperature.
;          If width is not present then the PWV absorption is just
;          calculated at the line central wavelength, as in
;          Chadney 2017. If width is included it must have the same
;          number of elements as wl.
;   ** Width currently not implemented! **
;
; Returns:
;  Transmission through the water vapour at each wl.
;

; Store coefficients in a common block to save calculating them every time
common _pwv_absorption, pwv_wl,pwv_coeff

if n_elements(pwv_wl) le 0 || ~array_equal(pwv_wl,wl) then begin
 ; Calculate the PWV coefficients
 pwv_coeff=get_pwv_coeff(wl)
 pwv_wl=wl
endif

return,exp(-pwv_coeff*pwv)


; Old stuff below here, moved to get_pwv_coeff

; ** Put all this stuff in another routine, which determines a Josh-style
;    exponential coefficient from fitting to a few values of PWV.

; Convert to wavenumber
v=1/(wl*1d-8)
nwl=n_elements(v)

; Get temperature (K), pressure (atm), and H2O number density (m^-3) profiles.
; Also dz, altitude step in m. All corresponding to 1mm of PWV.
restore,file='~/lib/NYA_water_vapour.idl'

; Scale H2O number density by the PWV
nH2O*=pwv

; Partial pressure of H2O, in atm
pp=nH2O*constant(/kb)*T/1.01325d5

tim=systime(1)
; Absorption cross section, dimensions (P, wl), converted to m^2
sigma=hitran_absorption(v, T, P, pp)*1d-4
print,'abs cs',systime(1)-tim

; Transmission coefficient - here we do the integration in altitude
; to get optical depth and then convert to transmission coefficient.
; Will have dimensions (nwl).
Tlambda=exp(-total(sigma*((nH2O*dz) # (dblarr(nwl)+1)),1))

return,Tlambda

end
