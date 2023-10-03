pro blueshift, mjs0, dseq, blueshifts
; Calculate centre of mass of peak and blueshift from H-alpha.
pnum = 3

; Get panel and wavelength arrays
get_p, mjs0, dseq, pnum, panel, /percentff
get_w, mjs0, pnum, wl
; Integrate panel to get 1D spectrum
space_integrated = total(panel, 2)/144
; Crop out region of spectrum of H-alpha line
peak = space_integrated[*,210:270]
a = size(peak)
b = a[1]

; pp is an array which holds the wavelength [0,*] and intensity [1,*] of the spectrum at 'time' [i]

pp = fltarr(2, 61)
pp[0,*] = wl[210:270]
blueshifts = fltarr(b)

for i=0,b-1 do begin ; Iterate in 'time'
;	pp[1,*] = peak[i,*]
;	j=0
;	psum=0
 ;	while psum le total(pp[1,*])/2. do begin ; Iterate over wavelength
;		psum += pp[1,j]
;		blueshifts[i] = 6562.81 - pp[0,j] ;H-alpha - centre of peak
;		j+=1
;	endwhile
	
	pp[1,*] = peak[i,*]
	centre = total(pp[1,*] * pp[0,*])/total(pp[1,*])
	blueshifts[i] = 6562.81 - centre
endfor 

end

;
;
;
;
;
