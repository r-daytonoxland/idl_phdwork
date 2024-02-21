pro proton_keogram, t1, length, keogram, fname
; Generates a proton line keogram from HiTIES spectra for a given start time, panel, and duration.
; Created by Rowan 14/12/2021
;
; Inputs:
;          t1 is the start time for the keogram, string of the form 'dd/mm/yyyy hh:mm:ss'
;          length is the duration of the keogram in hours
; Outputs:
;          keogram is an array containing the final keogram data
;-------------------------------------------------------------------------------------------------------
; Run read_tim to get the necessary data and get_panels for the panel boundaries
print, 'read_tim in progress'
read_tim, t1, length, mjs0, time, dseq, /nophot, tadd=5
get_panels, mjs0, np, xx, yy
; Check which panel is the proton line
print, 'creating keogram'
if np eq 4 then n=2   ; For the 4 panel mosaic H-alpha is panel 2.
if np eq 3 then n=3   ; For the 3 panel mosiac H-alpha is panel 3.
; Cut out the correct panel for each slide, this is panel #n.
panel = dseq[*, xx[n-1,0]:xx[n-1,1], yy[n-1,0]:yy[n-1,1]]
dims = size(panel)
; For each time dseq[t, x, y] average over y
summed_dseq = total(panel, 2) ; mean of dseq in that direction
keogram = summed_dseq/dims[2]

black = 1002.
white = 1005.

window, 0, xsize = 256*2, ysize = 512*2
tvin, bytscl(reverse(keogram, 2), min = black, max = white)
write_png, fname, tvrd(/true)

end