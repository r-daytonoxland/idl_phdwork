FUNCTION hitran_absorption, v, T, P, pp
;
; Calculates absorption coefficient, i.e. cross secton, at a given
; wavenumber for certain temperature and pressure conditions.
;
; Applies pressure broadening with a Lorentzian profile, but not
; Doppler broadening, so this is suitable for absorption in the lower
; atmosphere e.g. by tropospheric water vapour.
;
; Inputs:
;  v - Wavenumber at which to calculate the absorption cross section in cm^-1.
;  T - Temperature in K.
;  P - Pressure in atm.
;  pp - Partial pressure of the absorper in atm, currently just H2O.
;
; Output:
;  Absorption cross section in cm^2/molecule.
;
; P, T and pp can be vector arrays but they must have the same number of
; elements, e.g. describing a height profile.
; v can also be an array, with a different number of dimensions to P and T,
;  to calculate for multiple lines.
; The output will be an array with dimensions (pressures, wavenumbers).
;
; Based on equations at https://hitran.org/docs/definitions-and-units/
;

if n_elements(P) ne n_elements(T) then message,'Error: T and P must have the same number of elements.'
if n_elements(P) ne n_elements(pp) then message,'Error: P and pp must have the same number of elements.'

read_hitran
Tref=296d0 ; HITRAN reference temperature in K

common _hitran,ht_file,ht_v,ht_S,ht_A,ht_gair,ht_gself,ht_E,ht_nair,ht_dair

if (min(v) lt min(ht_v)) || (max(v) gt max(ht_v)) then message,'Warning: HITRAN data file may not cover the requested wavelengths.',/continue

n=n_elements(P) ; number of pressures and temperatures

; Only worry about lines near the wavelengths we're interested in
vtol=20.0 ; Corresponds to about 10A
s=where(ht_v ge (min(v)-vtol) and ht_v le (max(v)+vtol),m)
; m is now number of absorption lines

; Pressure broadened HWHM, dimensions (m,n)
gam=(((dblarr(m)+1) # (Tref/T))^(ht_nair[s] # (dblarr(n)+1)))*((ht_gair[s] # (P-pp))+(ht_gself[s] # pp))

; Pressure corrected line positions, dimensions (m,n)
vcorr=(ht_v[s] # (dblarr(n)+1))+(ht_dair[s] # P)

; Absorption line intensities, temperature corrected, dimensions (m,n)
Scorr=hitran_intensity(ht_v[s], ht_E[s], ht_S[s], T)

; Now too many dimensions so go in a for-loop over each wavenumber desired...
; Calculate absorption at each T,P due to each line, sum over all lines,
; then integrate the profile.
k=dblarr(n,n_elements(v))
; Lorentz profile multiplied by line intensity, summed over all lines
for i=0L,n_elements(v)-1 do k[*,i]=total((gam/((gam^2)+(v[i]-vcorr)^2))*Scorr,1)/!dpi
 
return,reform(k)

end
