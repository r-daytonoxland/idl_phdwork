pro calibrated_spectra, start, duration, interval, lib=lib
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

data = fltarr(rows, 402)

for i = 0, rows - 1 do begin
    read_tim, start_times[i], interval/3600., mjs0, time, dseq, icount, /nophot
    spectra, 3, mjs0, time, dseq, spectrum
    data[i,*] = [spectrum]
    print, 'Interval completed'
endfor

get_w, mjs0, 3, wl
columns = ['Wavelength', start_times]
datas =[[wl], [transpose(data)]]

fname_maker, start, 'spectrum', 'csv', fname, lib=lib
write_csv, fname, transpose(datas), header=columns
print, 'File saved at' + fname
end

pro runall
; Run all times of interest
; With Pc1
calibrated_spectra, '13/12/2021 06:10:00', 50*60, 10*60, /lib
calibrated_spectra, '15/12/2021 08:30:00', 30*60, 10*60, /lib
; WithoutPc1
calibrated_spectra, '06/12/2021 04:40:00', 60*60, 10*60, /lib
calibrated_spectra, '13/12/2021 14:10:00', 80*60, 10*60, /lib
calibrated_spectra, '11/12/2021 02:10:00', 110*60, 10*60, /lib
calibrated_spectra, '12/12/2021 05:00:00', 60*60, 10*60, /lib
calibrated_spectra, '14/12/2021 04:20:00', 40*60, 10*60, /lib

end