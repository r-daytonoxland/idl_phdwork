PRO synth_n2_1p, range, T, wl, int, width, labels, weights=weights
;
; Function to generate a synthetic N2 1P auroral spectrum.
;
; At the moment this does it the slow way! Have a look at speeding things
; up at some point...
;
; Inputs:
;  range - 2-element vector with wavelength range desired, in A
;  T     - Temperature, either scalar or a vector array, in K
; Outputs:
;  All arrays with the same number of elements.
;  wl    - Wavelengths of the N2 lines, in A
;  int   - Intensities of the N2 lines, normalised so the total is 1.
;  width - The width of each line due to Doppler broadening, FWHM in A.
;  labels - A string array of labels for each line, in the format
;           'OH (8-3) P1(3)'. *NOT YET IMPLEMENTED*
; Keywords:
;  weights - Ignored unless T is an array, see description below.
;
; If T is a single number this routine gives the spectrum for that
; temperature.
; If T is an array, the spectrum at each temperature in T is generated
; and then summed. If the weights keyword is set to an array with the
; same number of elements as T then the weighted sum is calculated,
; useful e.g. for generating an integrated spectrum for a temperature
; profile.
;
; Note water vapour absorption is not considered in this routine.
;

if n_elements(weights) eq 0 then weights=dblarr(n_elements(T))+1
if n_elements(weights) ne n_elements(T) then message,'Error: T and weights are different sizes'

; This will always output 2d arrays, first dimension temperature,
; even for only one temperature. This is slow!
;pregen_n2,T,min(range),max(range),N2_i,N2_w

; Now just generate spectra for a bunch of temperatures and then interpolate,
; until prog_n2_1p... can be made much faster. This speeds things up when
; using this routine in an iterative spectrum fitting process.
Tpoints=[(dindgen(59)+1)*20,(dindgen(20)*40)+1200,(dindgen(11)*100)+2000]
pregen_n2,Tpoints,min(range),max(range),N2_ipoints,N2_w
nlines=n_elements(N2_ipoints)/n_elements(Tpoints)
N2_i=dblarr(n_elements(T),nlines)
; Need to interpolate each individual emission line unfortunately
for i=0,nlines-1 do N2_i[*,i]=interpol(reform(N2_ipoints[*,i]),Tpoints,T,/quad);/spline)

; wl for each T should be the same...
wl=reform(N2_w[0,*])
int=total((weights # (dblarr(n_elements(wl))+1))*N2_i,1,/nan)

; Gives a NaN for some reason so need to get rid of it
s=where(finite(wl) and finite(int))
wl=wl[s]
int=int[s]

if arg_present(width) then begin
 ; Calculate FWHM of Doppler broadening

 ; Width factor - multiply by sqrt(T) and wl to get FWHM due to Doppler broadening
 w_factor=sqrt(8*alog(2)*constant(/kB)/(constant(/amu)*28*(constant(/c)^2)))

 if n_elements(T) eq 1 then width=wl*w_factor*sqrt(T) else $
 ; Need to normalise weights here but not above since intensity is normalised
 ; anyway afterwards, but not the width - it needs to be in proper units
 width=total((wl*w_factor) # (sqrt(T)*weights/total(weights)),2)

endif

int/=total(int)

end
