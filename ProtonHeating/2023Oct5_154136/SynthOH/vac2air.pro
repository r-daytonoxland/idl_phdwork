FUNCTION vac2air, wl
;
; Converts wavelength in vacuum to wavelength in air. This uses the conversion
; from Ciddor 1996, referenced in Chadney 2017. Works with arrays. Wavelength
; units are Angstroms for both input and output.
;

k0=238.0185d
k1=5792105d
k2=57.362d
k3=167917.0d

s2 = 1/(wl*1d-4)^2
n = (((k1/(k0-s2)) + (k3/(k2-s2)))/1d8) + 1

return,wl/n

end
