pro calibrated_spectra, start, duration, interval, fpath

rows = duration/interval

extract_datetime, start, datetime
intervals, datetime, interval, duration, start_times

fname_maker, start, 'spectrum', 'csv', fname, lib=lib

read_tim, start, 0.1/60., mjs0, time, dseq, icount, /nophot
get_w, mjs0, 3, wl

header = ['timestamp', wl]

data = fltarr[402, rows]

for i=0, rows -1 do begin
    read_tim, start_times[i], interval, mjs0, time, dseq, icount, /nophot
    spectra, 3, mjs0, time, dseq, spectrum
    data[*, i] = [starttimes[i], spectrum]
endfor

write_csv, fname, data, header
end