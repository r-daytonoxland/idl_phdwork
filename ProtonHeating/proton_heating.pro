pro fname_maker, time_string, result_type, interval_length, file_suffix, lib=lib, fname
; Generates a filename from the starting input string (time) and the result eg 'intensity' or 'spectrum'
; Inputs
;       time_string (string) : The start time of the data in 'dd/mm/yyyy hh:mm:ss' form
;       result_type (string) : A string describing the contents of the file e.g. 'spectrum'
;       interval_length (int): Length of each time interval in (s)
;       file_suffix (string) : The file type suffix e.g. 'txt'
; Keywords
;       lib : Saves the file to the lib directory in home 
; Outputs
;       fname (string) : Hopefully a suitable, standard filename for your result

extract_datetime, time_string, datetime

dt = string(datetime)

if keyword_set(lib) then begin
        fname = '~/lib/' + dt[0] + dt[1] + dt[2] + '_' + dt[3] + '_' + dt[4] + '_' + result_type + string(interval_length) + '.' + file_suffix
endif else begin
                fname = dt[0] + dt[1] + dt[2] + '_' + dt[3] + '_' + dt[4] + '_' + result_type + string(interval_length) + '.' + file_suffix
endelse

fname = fname.compress()

end

pro sp_func, X, A, F
; Input function for fitting the combined spectrum in the OH panel
; Inputs
;       X is wl (the wavelengths of the panel)
;       A is the guess free parameters A[0] = Temperature1, A[1] = Temperature2, A[2] = Intensity1, A[3] = Intensity2, A[4] = Background
; Outputs
;       A is the updated parameters
;       F is the function y values

pnum = 1
get_wlrange, pnum, wls

synth_oh, wls, A[0], ohwl, ohint, width, upperv=9
convolve_sp, ohwl, ohint, 0.6d, X, sp1

synth_oh, wls, A[1], ohwl, ohint, width, upperv=5
convolve_sp, ohwl, ohint, 0.6d, X, sp2

F = (A[2] * sp1) + (A[3] * sp2) + A[4]

end

pro sp_func_o, X, A, F
; Input function for fitting the combined spectrum in the OH panel
; Inputs
;       X is wl (the wavelengths of the panel)
;       A is the guess free parameters A[0] = Temperature, A[1] = Intensity, A[2] = Background
; Outputs
;       A is the updated parameters
;       F is the function y values

pnum = 2
get_wlrange, pnum, wls

synth_oh, wls, A[0], ohwl, ohint, width, upperv=9
convolve_sp, ohwl, ohint, 0.6d, X, sp

F = (A[1] * sp) + A[2]

end

pro ohpanel_fit, wl, sp, A, result
; Fits the two spectra for the oh panel separately
; Inputs
;       wl vector (402) : the wavelengths
;       sp vector (402) : the spectrum to be fit
; Outputs
;       A vector (5) : the output free parameters A[0] = Temperature1, A[1] = Temperature2, A[2] = Intensity1, A[3] = Intensity2, A[4] = Background
;       result vector (402) : the output best fit model spectrum

A = [200, 200, 0.05, 0.05, 0.0015]  ; Guess inputs
fita = [1, 1, 1, 1, 1]  ; Fit everything
weights = fltarr(402) + 1  ; Don't weight

result = curvefit(wl, sp, weights, A, sigma, function_name='sp_func', /noderivative, itmax=100, /double, fita=fita)

end

pro opanel_fit, wl, sp, A, result

A = [200d, 0.05d, 0.0015d]
fita = [1, 1, 1]
weights = fltarr(402) + 1

result = curvefit(wl, sp, weights, A, sigma, function_name='sp_func_o', /noderivative, itmax=100, /double, fita=fita)

end

pro proton_intensity, mjs0, dseq, result
; Calculates the intensity of the peak of the proton spectrum (assume pnum = 3) from index 210:270 (works for 2021 season)
; Inputs
;       mjs0 (double) : The mjs value from read_tim
;       dseq (float array) : The dseq from read_tim
; Outputs
;       result (double): The integrated peak intensity

get_p, mjs0, dseq, 3, panel, /percentff
slice, panel, space_integrated, scount, slice=[0.5, 0.3]
result = total(space_integrated[*, 210:270], 2)

end

pro proton_blueshift, mjs0, dseq, blueshift
; Calculates the blueshift as the wavelength difference from the proton profile COM from the rest Halpha wavelength
; Inputs
;       mjs0 (double) : Ths mjs value from read_tim
;       dseq (float array) : The dseq from read_tim
; Outputs
;       blueshift (double): The integrated peak intensity

pnum = 3
get_p, mjs0, dseq, pnum, panel, /percentff
get_w, mjs0, pnum, wl
slice, panel, space_integrated, scount
peak = space_integrated[210:270]
wls = wl[210:270]

; Compensates for calibration fuck up... but badly
peak -= mean(peak)
peak = -peak
for i=0,60 do if peak[i] le 0 then peak[i] = 0
; UPDATE this and make it not shite

centre = total(peak * wls)/total(peak)
blueshift = 6562.81 - centre

end

pro oh_fitting, mjs0, time, dseq, A, opanel=opanel, ohpanel=ohpanel
; Gets the temperatures and intensities of the OH band in the O panel or the combined bands in the OH panel.
; Inputs
;       mjs0 (double) : Ths mjs value from read_tim
;       time (float array) : The time from read_tim
;       dseq (float array) : The dseq from read_tim
; Outputs
;       A is either A[0] = Temperature, A[1] = Intensity, A[2] = Background for /opanel
;                or A[0] = Temperature1, A[1] = Temperature2, A[2] = Intensity1, A[3] = Intensity2, A[4] = Background for /ohpanel

av = reform(total(dseq, 1) / double(n_elements(time)), [1, 512, 512])
spectra, 2, mjs0, time[0], av, sp
get_w, mjs0, 2, wl

if keyword_set(ohpanel) then begin
        ohpanel_fit, wl, sp, A, result
endif

if keyword_set(opanel) then begin 
        opanel_fit, wl, sp, A, result
endif

end

pro extract_datetime, input_string, datetime
; Gets a datetime array out of a read_tim style datetime string
; Inputs
;       input_string (string) : datetime in the 'dd/mm/yyyy hh:mm:ss' format
; Outputs
;       datetime (array) : datetime in array form [yyyy, mm, dd, hh, mm, ss]

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

pro create_timeseries, start, interval_length, total_length, lib=lib
; Creates a timeseries of intensity, blueshift, and/or OH temperature for a chosen time period and interval length
; Inputs
;       start (string) : The start time of the period of interest in form 'dd/mm/yyyy hh:mm:ss'
;       interval_length (double) : The length of each interval in (s)
;       total_lentgh (double) : The total duration of the period of interst in (s) - must be a multiple of interval_length or it will break and be annoying
; Keywords
;       lib : saves to lib directory
; Outputs
;       Creates file

print, 'Getting datetimes...'
extract_datetime, start, datetime
intervals, datetime, interval_length, total_length, start_times

a = n_elements(start_times)

fname_maker, start, 'allparams', interval_length, 'txt', fname_gen, lib=lib

for i = 0, b-1 do begin
        print, 'Reading tim...'
        read_tim, start_times[i], interval_length / 3600., mjs0, time, dseq, icount, /nophot, tadd = interval_length

        ; Need to reduce space_integrated to just one averaged spectrum
        dseq = reform(total(dseq, 1), [1, 512, 512])

        proton_intensity, mjs0, dseq, intensity
        proton_blueshift, mjs0, dseq, blueshift
        oh_fitting, mjs0, time, dseq, A, /opanel
        oh_fitting, mjs0, time, dseq, B, /ohpanel

endfor

end

; pro create_timeseries, start, interval_length, total_length, intensity=intensity, blueshift=blueshift, temperature=temperature, lib=lib
; ; Creates a timeseries of intensity, blueshift, and/or OH temperature for a chosen time period and interval length
; ; Inputs
; ;       start (string) : The start time of the period of interest in form 'dd/mm/yyyy hh:mm:ss'
; ;       interval_length (double) : The length of each interval in (s)
; ;       total_lentgh (double) : The total duration of the period of interst in (s) - must be a multiple of interval_length or it will break and be annoying
; ; Keywords
; ;       intensity : Proton peak intensities will be calculated
; ;       blueshift : The blueshift of the COM of the proton peak will be calculated as the difference from the rest Ha wavelength
; ;       temperature : This will calculate the rotational temperature of the OH profiles from the OH panel (2)
; ;       lib : Saves outputs to lib directory
; ; Outputs
; ;       Creates separate file for each of the keywords used with filename of form (below), in the working directory or in lib/ depending on keywords
; ;               'yyyymmdd_hh_mm_intensity.txt' 

; print, 'Getting datetimes...'
; extract_datetime, start, datetime
; intervals, datetime, interval_length, total_length, start_times

; a = size(start_times)
; b = a[-1]  ; Gets the number of startimes i.e. intervals from start_times

; fname_maker, start, 'datatype', interval_length, 'txt', fname_gen, lib=lib

; ;Reading in data...
; for i = 0, b-1 do begin
;         print, 'Reading tim...'
;         read_tim, start_times[i], interval_length / 3600., mjs0, time, dseq, icount, /nophot, tadd = interval_length

;         ; Need to reduce space_integrated to just one averaged spectrum
;         dseq = reform(total(dseq, 1), [1, 512, 512])

;         if keyword_set(intensity) then begin
;                 fname_maker, start, 'intensity', interval_length, 'txt', lib=lib, fname
;                 print, 'Calculating intensity...'
;                 proton_intensity, mjs0, dseq, intensity
;                 openw, 1, fname, /append
;                 printf, 1, intensity
;                 close, 1
;                 print, 'Intensity saved'
;         endif
;         if keyword_set(blueshift) then begin
;                 fname_maker, start, 'blueshift', interval_length, 'txt', lib=lib, fname
;                 print, 'Calculating blueshift...'
;                 proton_blueshift, mjs0, dseq, blueshift
;                 openw, 1, fname, /append
;                 printf, 1, blueshift
;                 close, 1
;                 print, 'Blueshift saved'
;         endif
;         if keyword_set(temperature) then begin
;                 fname_maker, start, 'temperature', interval_length, 'txt', lib=lib, fname
;                 print, 'Calculating temperature...'
;                 oh_83_temperature, mjs0, time, dseq, temperature
;                 openw, 1, fname, /append
;                 printf, 1, temperature
;                 close, 1
;                 print, 'Temperature saved'
;         endif

; endfor

; print, ' Files saved at ' + fname_gen

; ; ;Error handling...
; ; blank = fltarr(60) + !Values.F_NaN

; ; ;Loop
; ; for i=0,nosteps do begin
; ;         catch, error_status
; ;         if error_status ne 0 then begin
; ;                 openw, 1, fname, /append
; ;                 printf, 1, blank, form='(60f10.5)'
; ;                 close, 1
; ;                 continue
; ;                 catch, /cancel
; ;         endif
; ;
; ;         openw, 1, fname, /append
; ;         if a[1] eq 120 then begin                           ; Check there is 60s of data
; ;                 acorr = a_correlate(timeseries, lags)
; ;                 printf, 1, acorr, form='(60f10.5)'
; ;         endif else begin
; ;                 printf, 1, blank, form='(60f10.5)'          ; If not print row of NaNs
; ;         endelse
; ;         close, 1
; ; endfor

; end