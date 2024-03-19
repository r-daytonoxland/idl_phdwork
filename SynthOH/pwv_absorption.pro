FUNCTION pwv_absorption, pwv, wl, recalc=recalc;, width
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
; Keywords:
;  recalc - Force this routine to forget all stored coefficients.
;
; Returns:
;  Transmission through the water vapour at each wl.
;

; Store coefficients in a common block to save calculating them every time
common _pwv_absorption, pwv_wl,pwv_coeff

if n_elements(pwv_wl) le 0 || keyword_set(recalc) || ~array_equal(pwv_wl,wl) then begin
 ; We don't have any yet, so calculate all the PWV coefficients
 pwv_coeff=get_pwv_coeff(wl)
 pwv_wl=wl
 sel=indgen(n_elements(wl))

;;; This bit doesn't work at the moment, need to get it to deal with
; duplicate elements in wl correctly.
endif else if 0 then begin ;if ~array_equal(pwv_wl,wl) then begin
 ; Calculate just PWV coefficients we don't already have
 sel=intarr(n_elements(wl))-1

 ; First find matching indicies for wl in pwv_wl
 srt1=sort(wl)
 unq1=uniq(wl[srt1])
 srt2=sort(pwv_wl)
 got1=array_intersect(wl[srt1[unq1]],pwv_wl[srt2],/index)
 got2=array_intersect(pwv_wl[srt2],wl[srt1],/index)
 ; got1 and got2 correspond to each other
 ; - srt2[got2] are the indicies of pwv_wl which correspond to the srt1[got1] indices of wl
 if finite(got1[0]) and finite(got2[0]) then sel[srt1[got1]]=srt2[got2]

 ; Then find missing coefficients and add them to the stored arrays
 missing=where(sel eq -1,c)
 if c gt 0 then begin
  sel[missing]=indgen(c)+n_elements(pwv_wl)
  pwv_coeff=[pwv_coeff,get_pwv_coeff(wl[missing])]
  pwv_wl=[pwv_wl,wl[missing]]
 endif
endif else sel=indgen(n_elements(wl))

return,exp(-pwv_coeff[sel]*pwv)


end
