; Conditional integration

; Read in FROG
times = ['01/01/2019 20:57:01',
        '02/03/2017 19:44:16', 
        '22/12/2014 06:25:46',
        '22/01/2020 22:01:07',
        '20/12/2014 19:13:40']
event_names=['A', 'B', 'C', 'D', 'E']

i = 0
read_tim, times[i], 20/60., mjs0, time, dseq, icount, /nophot
save, mjs0, time, dseq, filename='FROG' + event_names[i] + 'data.sav'

; Get right mjs start and end from some frogsinfo file

;frog_present
where is mjs start and mjs fin? find indexes
mjz = mjs[start:fin]
time = time[start:fin]
dseq = dseq[start:fin, *, *]

background = ?

if dseq is > e* background then?



spectrae = fltarr(481, 240)

for i=0,480 do begin
    spectra, pnum, mjs0, time[i], dseq[i,*,*], sp
    spectrae[i,*] = sp
endfor

read_tim, '22/12/2014 06:25:00', 2/60., mjs0, time, dseq, icount, /nophot
av=reform(total(dseq,1)/double(n_elements(time)),[1,512,512])