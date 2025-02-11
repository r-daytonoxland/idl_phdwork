; tried 23/12/2014 20:00:00, sih-0-1 - too cloudy for good stars
; tried 22/12/2019 02:00:00, sih-0-3 - shit stars
; 01/12/2014 00:00:00 sih-0-4 - oads of aurora aurora
; '22/11/2014 20:00:00' - no stars

read_tim, '25/12/2014 17:00:00', 4., mjs0, time, dseq, icount, /nophot, tadd=10
save, mjs0, time, dseq, filename='slit_calibration-0-6.sav'

restore, filename='slit_calibration-0-6.sav', /verbose

get_p,mjs0,dseq,1,panel1out
slice_sp,panel1out,keoout
keoscale=bytscl(keoout, min=0, max=100)
read_off_sih, mjs0, keoscale, 'sih-0-6.dat'

; Click on the right edge to move forward
; Click on the stars from left to right
; Click on the left edge to end
; Should have written a file 'test.dat' with a time and a bunch of locations

read_sao
associate_slit, mjs0, time, keoscale, 'sih-0-6.dat', pnum=1
; y for yes and n for no
fit_slit, 'sih-0-6.dat'

; 1 is O+, 2 is Hbeta, 3 is O2+, 4 is O
; Do for panels 1/3 and 2e