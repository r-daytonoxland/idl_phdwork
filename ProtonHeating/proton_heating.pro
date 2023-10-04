pro proton_intensity, mjs0, dseq, result
; Calculates the intensity of the peak of the proton spectrum (assume pnum = 3) from index 210:270 (works for 2021 season)
; Inputs
;       mjs0 (int) : The mjs value from read_tim
;       dseq (array) : The dseq from read_tim
; Outputs
;       result (float): The integrated peak intensity

get_p, mjs0, dseq, 3, panel, /percentff
slice, panel, space_integrated, scount, slice=[0.5, 0.3]
result = total(space_integrated[*, 210:270], 2)

end

pro proton_blueshift, mjs0, dseq, blueshift
;
pnum = 3
get_p, mjs0, dseq, pnum, panel, /percentff
get_w, mjs0, pnum, wl
slice, panel, space_integrated, scount
peak = space_integrated[timeindices[0]:timeindices[1],210:270]
; a = size(peak)
; b = a[1]
; pp holds ithe wavelength [0,*] and intensity [1,*] of the spectrum at 'time' [i]
pp = fltarr(2, 61)
pp[0,*] = wl[210:270]

pp[1,*] = peak[j,*]
centre = total(pp[1,*] * pp[0,*])/total(pp[1,*])
blueshift = 6562.81 - centre

end

pro oh_temperature
;

end

pro extract_datetime, input_string, datetime

; Input string in the 'dd/mm/yyyy hh:mm:ss' format

; Split the string into date and time components
date_time_components = input_string.Split(' ')

; Extract date components
date_components = date_time_components[0].Split('/')

day = Long(date_components[0])
month = Long(date_components[1])
year = Long(date_components[2])

; Extract time components
time_components = date_time_components[1].split(':')

hour = Long(time_components[0])
minute = Long(time_components[1])
second = Long(time_components[2])

datetime = [year, month, day, hour, minute, second]

end

pro intervals, datetime, interval_length, total_length, start_times

;Splits a longer time into intervals and returns the start of each interval in 'dd/mm/yyyy hh:mm:ss' format
; Inputs-
;       Datetime (array(integers)) : Output from extract_datetime, array with [year, month, day, hour, minute, second]
;       interval_length (integer)  : Length of each interval in seconds
;       total_length (integer)     : Length of total period of interest
; Outputs-
;       start_times (array(string)): start of each interval in 'dd/mm/yyyy hh:mm:ss' format

tt_mjs, datetime[0], datetime[1], datetime[2], datetime[3], datetime[4], datetime[5], 00, mjs

number_intervals = total_length / interval_length

start_times = strarr(number_intervals)

for i = 0 ,number_intervals - 1 do begin
	startimemjs = mjs + (i * interval_length)
	dat2str, startimemjs, startimestr, /other
	start_times[i] = startimestr
endfor

end

pro fname_maker, time_string, result_type, file_suffix, fname
; Generates a filename from the starting input string (time) and the result eg 'intensity' or 'spectrum'
; Inputs
;       time_string (string) : The start time of the data in 'dd/mm/yyyy hh:mm:ss' form
;       result_type (string) : A string describing the contents of the file e.g. 'spectrum'
;       file_suffix (string) : The file type suffix e.g. 'txt'
; Outputs
;       fname (string) : Hopefully a suitable, standard filename for your result

extract_datetime, time_string, datetime

dt = string(datetime)

fname = dt[0] + dt[1] + dt[2] + '_' + dt[3] + '_' + dt[4] + '_' + result_type + '.' + file_suffix
fname = fname.compress()

end

pro create_timeseries, panel, start, interval_length, total_length, intensity=intensity, blueshift=blueshift, temperature=temperature

; Creates a timeseries doing one of a few processes I hope??
; Input:
;       panel (integer) = 2 for OH panel and 3 for proton panel
;       time_int (integer) = Length of time in (s) for 

extract_datetime, start, datetime
intervals, datetime, interval_length, total_length, start_times

a = size(start_times)
b = a[-1]  ; Gets the number of startimes i.e. intervals from start_times

;Reading in data...
for i = 0, b-1 do begin
        read_tim, start_times[i], interval_length / 3600., mjs0, time, dseq, icount, /nophot, tadd = interval_length
        if keyword_set(intensity) then begin
                fname_maker, start, 'intensity', 'txt', fname
                proton_intensity, mjs0, dseq, intensity
                openw, 1, fname, /append
                printf, 1, intensity
                close, 1
        endif
        if keyword_set(blueshift) then begin
                fname_maker, start, 'blueshift', 'txt', fname
                proton_blueshift, blueshift
                openw, 1, fname, /append
                printf, 1, blueshift
                close, 1
        endif
        if keyword_set(temperature) then begin
                fname_maker, start, 'temperature', 'txt', fname
                oh_temperature, temperature
                openw, 1, fname, /append
                printf, 1, temperature
                close, 1
        endif

endfor

print, ' Files saved at ' + fname

;Getting panel...
; a = size(time)
; dim = a[1]          ; The length of time in 1/2 second steps
; nosteps = dim / 120   ; Number of minutes-1 

;Integrating in spatial direction...
;peakint = total(space_integrated[*,210:270], 2)
;peakint = peakint/120. ; mean


; ;Error handling...
; blank = fltarr(60) + !Values.F_NaN

; ;Loop
; for i=0,nosteps do begin
;         catch, error_status
;         if error_status ne 0 then begin
;                 openw, 1, fname, /append
;                 printf, 1, blank, form='(60f10.5)'
;                 close, 1
;                 continue
;                 catch, /cancel
;         endif
;         ;timeindices = [where(time eq i*60.),where(time eq ((i+1)*60.)-0.5)] ; Start,end of minute
;         ;minutedata = peakint[timeindices[0]:timeindices[1]] ; Data for that minute	

; 	;Background remove
; 	;x = findgen(120)
; 	;background = linfit(x, minutedata)
; 	;timeseries = minutedata - (background[0] + x*background[1])

; 	;Make and save autocorrelation for minute of data
;         a = size(space_integrated)
;         openw, 1, fname, /append
;         if a[1] eq 120 then begin                           ; Check there is 60s of data
;                 acorr = a_correlate(timeseries, lags)
;                 printf, 1, acorr, form='(60f10.5)'
;         endif else begin
;                 printf, 1, blank, form='(60f10.5)'          ; If not print row of NaNs
;         endelse
;         close, 1
; endfor

end