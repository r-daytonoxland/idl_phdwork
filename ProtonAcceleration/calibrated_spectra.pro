pro calibrated_spectra, start, duration, interval
;
;
;

rows = duration/interval

extract_datetime, start, datetime
intervals, datetime, interval, duration, start_times

data = fltarr(403, rows)

for i=0, rows -1 do begin
    read_tim, start_times[i], interval, mjs0, time, dseq, icount, /nophot
    spectra, 3, mjs0, time, dseq, spectrum
    data[*, i] = [starttimes[i], spectrum]
endfor

get_w, mjs0, 3, wl
header = ['timestamp', string(wl)]

fname_maker, start, 'spectrum', 'csv', fname, lib=lib
write_csv, fname, data, header
end