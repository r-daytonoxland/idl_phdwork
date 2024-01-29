FUNCTION air2vac, wl
;
; Converts wavelength in air to wavelength in vacuum. This is much less
; trivial than the opposite, vacuum to air, since refractive index n depends
; on vacuum wavelenth. This equation comes from Uppsala:
; https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
;
; Works with arrays.
; Wavelength units are Angstroms for both input and output.
;

k0=8.336624212083d-5
k1=2.408926869968d-2
k2=1.301065924522d2
k3=1.599740894897d-4
k4=3.892568793293d1

s2=(1d4/wl)^2
n=1+k0+(k1/(k2-s2))+(k3/(k4-s2))

return,wl*n


end