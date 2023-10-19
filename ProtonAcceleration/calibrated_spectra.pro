pro calibrated_spectra, start, duration, interval
; Creates CSV files with the calibrated spectra (proton panel is chosen) in the desired intervals
; File header is 'timestamp' followed by each wavelength, first column is the start time of the interval and the following columns are the spectrum intensity values
; Inputs
;   start (string) : The start time of the desired time period in form 'dd/mm/yyyy hh:mm:ss'
;   duration (int) : The duration of the entire time period in s
;   interval (int) : The length of each interval (time the spectrum will be integrated over) in s
; Keywords
;   lib=lib : Saves the output file in ~/lib instead of the working directory
; Outputs
;   CSV file of form yyyymmdd_h_m_spectrum.csv 

rows = duration/interval

extract_datetime, start, datetime
intervals, datetime, interval, duration, start_times

data = fltarr(403, rows)

for i=0, rows -1 do begin
    read_tim, start_times[i], interval/3600., mjs0, time, dseq, icount, /nophot
    spectra, 3, mjs0, time, dseq, spectrum
    data[*, i] = [string(start_times[i]), string(spectrum)]
    print, 'Interval completed'
endfor

get_w, mjs0, 3, wl
header = ['timestamp', string(wl)]

fname_maker, start, 'spectrum', 'csv', fname, lib=lib
write_csv, fname, data, header
print, 'File saved at' + fname
end