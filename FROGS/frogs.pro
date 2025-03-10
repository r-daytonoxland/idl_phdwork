;read_tim, '22/12/2014 06:25:00', 5/60., mjs0, time, dseq, icount, /nophot
;av=reform(total(dseq,1)/double(n_elements(time)),[1,512,512])

; Use this to produce a spectrum for every time for each panel and export to Python
pro spectra_time, pnum, mjs0, time, dseq, wl, spt 
    ; len(time)
    a=size(time)
    time_len=a[1]
    ; len(wl)
    get_w, mjs0, pnum, wl
    b = size(wl)
    wl_len = b[1]
    ; create array with dimensions, time, wl
    spt=fltarr(time_len, wl_len)
    ; iterate over every time, produce spectrum and save
    for i = 0,time_len-1 do begin
        spectra, pnum, mjs0, time[i], dseq[i,*,*], sp
        spt[i,*] = sp
    endfor
end

;save, filename='frog_pnum.sav', wl, spt