PRO synth_op2p, range, T, wl, int, width, labels, weights=weights
;
; Function to generate a synthetic O+ 2P auroral spectrum.
; Inputs:
;  range - 2-element vector with wavelength range desired, in A
;  T     - Ratio of densities of the 2P1/2 to 2P3/2 states
;     ** Note T is NOT temperature in this routine, yet **
; Outputs:
;  All arrays with the same number of elements.
;  wl    - Wavelengths of the O+ 2P lines, in A
;  int   - Intensities of the O+ 2P lines, normalised so the total is 1.
;  width - The width of each line due to Doppler broadening, FWHM in A.
;          ** Note we assume a temperature of 500K here **
;  labels - A string array of labels for each line, in the format
;           'O+ 2D3/2 - 2P1/2', with IDL formatting for sub/superscript.
; Keywords:
;  weights - Ignored unless T is an array, see description below.
;
; If T is a single number this routine gives the spectrum for that
; ratio value.
; If T is an array, the spectrum at each ratio in T is generated
; and then summed. If the weights keyword is set to an array with the
; same number of elements as T then the weighted sum is calculated,
; useful e.g. for generating an integrated spectrum for a temperature
; profile.
;
; Note water vapour absorption is not considered in this routine.

if n_elements(weights) eq 0 then weights=dblarr(n_elements(T))+1
if n_elements(weights) ne n_elements(T) then message,'Error: T and weights are different sizes'

; Data, as vectors all with the same number of elements:
; lab - Line label
; A - Einstein coefficients, from Zeippen 1987 (s^-1)
; wl - Wavelengths in air from Sharpee 2004 (A)
st=['!U2!NP!D3/2!N','!U2!NP!D1/2!N','!U2!ND!D5/2!N','!U2!ND!D3/2!N']
lab=[st[2]+'-'+st[1],st[2]+'-'+st[0],st[3]+'-'+st[1],st[3]+'-'+st[0]]
;   D5/2-P1/2,D5/2-P3/2,D3/2-P1/2,D3/2-P3/2
A =[  5.63d-2, 1.074d-1,  9.41d-2,  5.80d-2]
; Was these until 2023/10/25:
;A =[  5.63d-2,  1.07d-1,  9.39d-2,  5.78d-2]
wl=[ 7319.044, 7320.121, 7329.675, 7330.755] ; As in Whiter 2014
;wl=[7318.92,7319.99,7329.67,7330.73] ; From NIST
; 1 if it comes from 2P3/2 state, 0 otherwise:
up=[        0,        1,        0,        1]

; Filter to just the lines we care about
s=where(wl ge min(range) and wl le max(range),c)
if (c le 0) then begin
 ; There are no lines present
 wl=!values.F_NaN
 int=!values.F_NaN
 width=!values.F_NaN
 labels=''
 return
endif
if arg_present(labels) then labels='O+ '+lab[s]
A=A[s]
up=up[s]
wl=double(wl[s])

R=T ; Density ratio of 2P states - hopefully in future derive from temperature
; Assume number density of 2P3/2 state is 1.
; Then number density of 2P1/2 state is R.

if n_elements(R) eq 1 then int=(((1-up)*R)+up)*A else $
int=total( (((1-up)*A) # (R*weights/total(weights))) + ((up*A) # (dblarr(n_elements(R))+(1d/n_elements(R)))),2)

if arg_present(width) then begin
 ; Calculate FWHM of Doppler broadening

 ; Width factor - multiply by sqrt(Temperature) and wl to get FWHM due to Doppler broadening
 w_factor=sqrt(8*alog(2)*constant(/kB)/(constant(/amu)*16*(constant(/c)^2)))

 ; Currently assume a temperature of 500K - could fix later
 width=wl*w_factor*sqrt(500d)

endif

int/=total(int)

end
