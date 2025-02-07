read_tim, '22/12/2014 21:30:00', 4., mjs0, time, dseq, icount, /nophot, tadd=10

get_p,mjs0,dseq,1,panel1out
slice_sp,panel1out,keoout
keoscale=bytscl(keoout)

; Click on the right edge to move forward
; Click on the stars from left to right
; Click on the left edge to end
; Should have written a file 'test.dat' with a time and a bunch of locations

read_sao
associate_slit, mjs0, time, ..
; y for yes and n for no
fit_slit

; 1 is O+, 2 is Hbeta, 3 is O2+, 4 is O
; Do for panels 1/3 and 2