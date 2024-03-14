PRO synth_oh, range, T, wl, int, width, labels, weights=weights, file=file,$
          chadney=chadney, upperv=upperv
;
; Function to generate a synthetic OH airglow spectrum.
;
; IMPORTANT NOTE: The equation used to calculate the line intensities
;   will give incorrect ratios between vibrational bands with different
;   upper states, since the populations will be different. If looking
;   at a region with multiple bands, you probably want to call this routine
;   multiple times with the upperv keyword set to each upper vibrational
;   state in turn.
;
; Inputs:
;  range - 2-element vector with wavelength range desired, in A
;  T     - Temperature, either scalar or a vector array, in K
; Outputs:
;  All arrays with the same number of elements.
;  wl    - Wavelengths of the OH lines, in A
;  int   - Intensities of the OH lines, normalised so the total is 1.
;  width - The width of each line due to Doppler broadening, FWHM in A.
;  labels - A string array of labels for each line, in the format
;           'OH (8-3) P1(3)'.
; Keywords:
;  weights - Ignored unless T is an array, see description below.
;  file  - The file with the OH data, the "line list" from the
;          supplementary material of Brooke et al. 2016
;          doi:10.1016/j.jqsrt.2015.07.021
;          Default ~/lib/OH-XX-Line_list_201804001442519260.txt.
;  chadney - Use the coefficients from Mies 1974 and Krassovsky 1962
;          instead of the ones from the Brooke 2016 file. Only relevant
;          for the 8-3 band - coefficients for others not included.
;  upperv - Only include lines from bands originating in the upperv
;          vibrational level. Can be a scalar (e.g. upperv=8 will give
;          the 8-3 band, 8-2 band, etc), or a 2-element array for a
;          range of v, e.g. upperv=[6,8] will give 6-4, 7-2, 8-3, bands etc.
;          * NOT YET IMPLEMENTED *
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

if n_elements(file) lt 1 then file='~/lib/OH-XX-Line_list_201804001442519260.txt'

; Line notation e.g. P1(4)
; delJ = -1, 0, 1 for P, Q, R branches
; s  = +1/2, -1/2 for 1, 2 subscript
; L  = 1, 2, 3..., in brackets

if n_elements(chadney) lt 1 then begin ; Use the newer coefficients

 ; Store OH data in common block to save reading the file everytime.
 ; If the range changes then the file is re-read.
 common _synth_oh, OH_file,OH_range,_lab,_A,_E,_J,_w,_vu

 if (n_elements(OH_file) lt 1) || (OH_file ne file) || ~array_equal(OH_range, range) then begin
  ; Read the data
  _lab=strarr(50000)
  _A=dblarr(50000)
  _E=dblarr(50000)
  _J=dblarr(50000)
  _w=dblarr(50000)
  _vu=lonarr(50000)
  vrange=1d8/air2vac(range)
  vrange=double([floor(min(vrange)),ceil(max(vrange))]) ; Range as wavenumber in cm^-1
  i=0
  openr,lun,file,/get_lun
  aa='blah'
  while strmid(aa,0,7) ne " v' v''" do readf,lun,aa ; Read the header
  while ~eof(lun) do begin
   readf,lun,aa
   dat=strsplit(aa,/extract)
   ; Don't include lines outside our wavelength range
   if (double(dat[9]) lt vrange[0]) || (double(dat[9]) gt vrange[1]) then continue ; Not in our range, don't care about it
   ; Don't include transitions from vibrational levels > 9 or < 6
   ; Not needed now we have upperv keyword
   ;;if (fix(dat[0]) gt 9) || (fix(dat[0]) lt 6) then continue
   ; Don't include transitions with different F-levels
   ;if (dat[4] ne dat[5]) then continue
   ; Don't include transitions with J > 9
   ;if (float(dat[2]) gt 9) then continue
   ; For Q-branch, don't include L>6 ; Note: was 3
   ;if (dat[2] eq dat[3]) && ((float(dat[2])+float(dat[4])-1.5) gt 6) then continue
   _w[i]=vac2air(1d8/double(dat[9]))
   _A[i]=double(dat[12])
   _E[i]=(double(dat[11])+double(dat[9]))*100 ; Need upper state energy, so add wavenumber to lower state energy, then convert to m^-1
   _J[i]=double(dat[2])
   _vu[i]=long(dat[0])
   _lab[i]=oh_label(double(dat[0:5]))+dat[7]
   i+=1
  endwhile
  close,lun
  free_lun, lun
  if i gt 0 then begin
   _lab=_lab[0:i-1]
   _A=_A[0:i-1]
   _E=_E[0:i-1]
   _J=_J[0:i-1]
   _w=_w[0:i-1]
   _vu=_vu[0:i-1]
  endif else begin ; No lines in this wavelength range
   _lab=''
   _A=!values.F_NaN
   _E=!values.F_NaN
   _J=!values.F_NaN
   _w=!values.F_NaN
   _vu=!values.F_NaN
  endelse
  OH_file=file
  OH_range=range
 endif

 ; Filter on upperv
 if (n_elements(upperv) eq 1) then s=where(_vu eq upperv) $
  else if (n_elements(upperv) eq 2) then s=where(_vu ge min(upperv) and _vu le max(upperv)) $
  else s=indgen(n_elements(_vu))
 lab=_lab[s]
 A=_A[s]
 E=_E[s]
 J=_J[s]
 w=_w[s]
 vu=_vu[s]

 if (n_elements(w) le 1) && (w eq !values.F_NaN) then begin ; No lines in this wavelength range
  wl=!values.F_NaN
  int=!values.F_NaN
  width=!values.F_NaN
  labels=''
  return
 endif
 if arg_present(labels) then labels=lab
 wl=w

endif else begin ; Use the old coefficients

 ; Old data, as vectors all with the same number of elements:
 ; lab - Line label
 ; A - Einstein coefficients, from Mies 1974 (s^-1)
 ; E - Energy levels, from Krassovsky 1962 (cm^-1 relative to ground state)
 ; J - Total angular momentum quantum number, delJ + L + s
 ; wl - Wavelengths in air from Chadney 2017 (A)

 ; Only 8-3 band at the moment
 lab=[ 'P1(2)e', 'P1(2)f', 'P1(3)e', 'P1(3)f', 'P1(4)e', 'P1(4)f', 'P1(5)e', 'P1(5)f', 'P1(6)e', 'P1(6)f', 'P1(7)e', 'P1(7)f', 'P1(8)e', 'P1(8)f', 'P1(9)e', 'P1(9)f', 'P2(2)e', 'P2(2)f', 'P2(3)e', 'P2(3)f', 'P2(4)e', 'P2(4)f', 'P2(5)e', 'P2(5)f', 'P2(6)e', 'P2(6)f', 'P2(7)e', 'P2(7)f', 'P2(8)e', 'P2(8)f', 'P2(9)e', 'P2(9)f', 'Q1(1)e', 'Q1(1)f', 'Q1(2)e', 'Q1(2)f', 'Q1(3)e', 'Q1(3)f', 'Q1(4)e', 'Q1(4)f', 'Q1(5)e', 'Q1(5)f', 'Q1(6)e', 'Q1(6)f', 'Q1(7)e', 'Q1(7)f', 'Q1(8)e', 'Q1(8)f', 'Q2(1)e', 'Q2(1)f', 'Q2(2)e', 'Q2(2)f', 'Q2(3)e', 'Q2(3)f', 'Q2(4)e', 'Q2(4)f', 'Q2(5)e', 'Q2(5)f', 'Q2(6)e', 'Q2(6)f', 'Q2(7)e', 'Q2(7)f', 'Q2(8)e', 'Q2(8)f']
 A  =[    0.243,    0.243,    0.298,    0.298,    0.321,    0.321,    0.335,    0.335,    0.346,    0.346,    0.355,    0.355,    0.363,    0.363,    0.370,    0.370,    0.387,    0.387,    0.360,    0.360,    0.354,    0.354,    0.355,    0.355,    0.359,    0.359,    0.364,    0.364,    0.369,    0.369,    0.375,    0.375,    0.329,    0.329,    0.135,    0.135,    0.072,    0.072,    0.044,    0.044,    0.029,    0.029,    0.021,    0.021,    0.015,    0.015,    0.012,    0.012,    0.188,    0.188,    0.041,    0.041,    0.020,    0.020,    0.012,    0.012,    0.009,    0.009,    0.007,    0.007,    0.006,    0.006,    0.005,    0.005]
 E  =[    -47.8,    -47.8,     12.1,     12.1,     96.0,     96.0,    204.6,    204.6,    338.1,    338.1,    494.6,    494.6,    675.6,    675.6,    811.6,    811.6,     83.0,     83.0,    125.2,    125.2,    194.2,    194.2,    290.7,    290.7,    413.2,    413.2,    562.2,    562.2,    737.2,    737.2,    938.2,    938.2,    -47.8,    -47.8,     12.1,     12.1,     96.0,     96.0,    204.6,    204.6,    338.1,    338.1,    494.6,    494.6,    675.6,    675.6,    811.6,    811.6,     83.0,     83.0,    125.2,    125.2,    194.2,    194.2,    290.7,    290.7,    413.2,    413.2,    562.2,    562.2,    737.2,    737.2,    938.2,    938.2]
 J  =[      1.5,      1.5,      2.5,      2.5,      3.5,      3.5,      4.5,      4.5,      5.5,      5.5,      6.5,      6.5,      7.5,      7.5,      8.5,      8.5,      0.5,      0.5,      1.5,      1.5,      2.5,      2.5,      3.5,      3.5,      4.5,      4.5,      5.5,      5.5,      6.5,      6.5,      7.5,      7.5,      1.5,      1.5,      2.5,      2.5,      3.5,      3.5,      4.5,      4.5,      5.5,      5.5,      6.5,      6.5,      7.5,      7.5,      8.5,      8.5,      0.5,      0.5,      1.5,      1.5,      2.5,      2.5,      3.5,      3.5,      4.5,      4.5,      5.5,      5.5,      6.5,      6.5,      7.5,      7.5]
 wl =[ 7316.246, 7316.318, 7340.813, 7340.956, 7369.248, 7369.483, 7401.688, 7402.028,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    , 7303.679, 7303.754, 7329.123, 7329.173, 7358.657, 7358.661, 7392.169, 7392.228,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    , 7276.400, 7276.410, 7284.421, 7284.457, 7295.901, 7295.978,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    , 7275.095, 7275.120,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ]
 ; Note, some wavelengths missing for higher L values, especially in Q-branch

 ; Filter to just the lines we care about
 s=where(wl ge min(range) and wl le max(range),c)
 if (c le 0 || (n_elements(upperv) ne 0 && upperv[0] ne 8)) then begin
 ; ; There are no lines present
  wl=!values.F_NaN
  int=!values.F_NaN
  width=!values.F_NaN
  labels=''
  return
 endif
 if arg_present(labels) then labels='OH (8-3) '+lab[s]
 A=A[s]
 E=E[s]*100 ; Converted to m^-1
 J=J[s]
 wl=double(wl[s])

endelse ; Done sorting coefficients

; For OH, equation for relative intensity of a line is:
; Intensity = A(2J + 1) * exp(-hcE/kT)
; (Won't work across multiple vibrational bands, since populations
;  will be different)

; Constants inside the exponential
exppart=-constant(/h)*constant(/c)/constant(/kb)

; Matrix, first dimension lines, second dimension temperatures,
; and then totaled in the temperature dimension:
int=total((weights # (A*(2*J + 1)))*exp((exppart/T) # E),1)

if arg_present(width) then begin
 ; Calculate FWHM of Doppler broadening

 ; Width factor - multiply by sqrt(T) and wl to get FWHM due to Doppler broadening
 w_factor=sqrt(8*alog(2)*constant(/kB)/(constant(/amu)*17*(constant(/c)^2)))

 ; Need to normalise weights here but not above since intensity is normalised
 ; anyway afterwards, but not the width - it needs to be in proper units
 width=total((sqrt(T)*weights/total(weights)) # (wl*w_factor),1)

endif

int/=total(int)

end
