pro readacfs, fname, len, acorr
; Reads in the acfs we just made
; Updated 28/04/2022
;
; Inputs
; fname is the filename, use filepath '/stp/raid1/workdir/rado1g21/data/'
; len in length of data in mins
;
; Outputs
; acorr, a [60,norows] fltarray containing the autocorrelation data


acorr = fltarr(60, len)

openr, 1, fname
readf, 1, acorr
close, 1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;start='15/12/2021 05:00:00'
;year=2021 & month=12 & day=15 & hour=05 & min=00 & sec=00 & millisec=00
;tt_mjs, year, month, day, hour, min, sec, millisec, mjs
;len = 20 * 60 ;20 minutes in seconds
;totalen = 4 * 60 ;minutes
;totalensecs = totalen * 60
;nochonks = totalensecs / len

pro getchonks, year, month, day, hour, min, sec, millisec, minsinchonk, totlength, startimes

tt_mjs, year, month, day, hour, min, sec, millisec, mjs
lensecs = totlength * 60
secsinchonk = minsinchonk * 60

nochonks = totlength/minsinchonk

startimes = strarr(nochonks)

for i=0,nochonks-1 do begin
	startimemjs = mjs + (i*secsinchonk)
	dat2str, startimemjs, startimestr, /other
	startimes[i] = startimestr
endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro makeffts, start, length, fname 

panel = 3 ; Choose H  panel
print, 'Reading in data...'
read_tim, start, (length)/60., mjs0, time, dseq, icount, /nophot

a = size(time)
dim = a[1]          ; The length of time in 1/2 second steps
nosteps = dim/120   ; Number of minutes-1 
print, 'Getting panel...'
get_p, mjs0, dseq, 3, panel, /percentff

;Doing a bunch of initialisation
space_integrated = total(panel, 2)/144
peakint = total(space_integrated[*,210:270], 2)
peakint = peakint/120. ; mean

blank = fltarr(120) + !Values.F_NaN

for i=0,nosteps do begin

	catch, error_status
        if error_status ne 0 then begin
                openw, 1, fname, /append
                printf, 1, blank, form='(120f10.5)'
                close, 1
                continue
                catch, /cancel
        endif

	timeindices = [where(time eq i*60.),where(time eq ((i+1)*60.)-0.5)] ; Start,end of minute
        minutedata = peakint[timeindices[0]:timeindices[1]] ; Data for that minute

        ;Background remove
        x = findgen(120)
        background = linfit(x, minutedata)
        timeseries = minutedata - (background[0] + x*background[1])

	a=size(timeseries)
	hwindow = hanning(120, /double)	
	openw, 1, fname, /append
	if a[1] eq 120 then begin
		minfft = fft(hwindow*timeseries)
		printf, 1, real_part(minfft), form='(120f10.5)'
	endif else begin
		printf, 1, blank, form='(120f10.5)'
	endelse
	close, 1 
endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro readffts, fname, len, fftdat
; Reads in the acfs we just made
; Updated 28/04/2022
;
; Inputs
; fname is the filename, use filepath '/stp/raid1/workdir/rado1g21/data/'
; norows the number of rowa (length-1)
;
; Outputs
; acorr, a [60,norows] fltarray containing the autocorrelation data


fftdat = fltarr(120, len)

openr, 1, fname
readf, 1, fftdat
close, 1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro make_intensity_ACFs, start, length, fname

; Procedure to make ACFs for every minute in a certain period for lags consistent with Pc1 pulsations i.e. 0-15s
; Updated 19/05/2022
;
; Inputs
; Start - Start time in form 'dd/mm/yyyy hh:mm:ss'
; Length - in minutes
; fname is the filename, use filepath '/stp/raid1/workdir/rado1g21/data/'
;
; Outputs
; peakint - timeseries of integrated proton lines


panel = 3 ; Choose H  panel
print, 'Reading in data...'
read_tim, start, (length)/60., mjs0, time, dseq, icount, /nophot
print, 'Getting panel...'
a = size(time)
dim = a[1]          ; The length of time in 1/2 second steps
nosteps = dim/120   ; Number of minutes-1 
get_p, mjs0, dseq, 3, panel, /percentff

;Integrate to get time series
;space_integrated = total(panel, 2)/144
slice, panel, space_integrated, scount, slice=[0.5, 0.3]
peakint = total(space_integrated[*,210:270], 2)
peakint = peakint/120. ; mean

lags = indgen(60)
acorr = fltarr(nosteps+1, 60)
print, 'Making  ACFs'

;Error handling
blank = fltarr(60) + !Values.F_NaN

;Autocorrelation loop
for i=0,nosteps do begin
        catch, error_status
        if error_status ne 0 then begin
                openw, 1, fname, /append
                printf, 1, blank, form='(60f10.5)'
                close, 1
                continue
                catch, /cancel
        endif
        timeindices = [where(time eq i*60.),where(time eq ((i+1)*60.)-0.5)] ; Start,end of minute
        minutedata = peakint[timeindices[0]:timeindices[1]] ; Data for that minute	

	;Background remove
	x = findgen(120)
	background = linfit(x, minutedata)
	timeseries = minutedata - (background[0] + x*background[1])

	;Make and save autocorrelation for minute of data
        a = size(timeseries)
        openw, 1, fname, /append
        if a[1] eq 120 then begin                           ; Check there is 60s of data
                acorr = a_correlate(timeseries, lags)
                printf, 1, acorr, form='(60f10.5)'
        endif else begin
                printf, 1, blank, form='(60f10.5)'          ; If not print row of NaNs
        endelse
        close, 1
endfor
end

;;;;;;;;;;;;;;;;;

pro make_blueshift_ACFs, start, length, fname

pnum = 3 ; Choose H  panel
print, 'Reading in data...'
read_tim, start, (length)/60., mjs0, time, dseq, icount, /nophot
print, 'Getting panel...'
a = size(time)
dim = a[1]          ; The length of time in 1/2 second steps
nosteps = dim/120 
; Get panel and wavelength arrays
get_p, mjs0, dseq, pnum, panel, /percentff
get_w, mjs0, pnum, wl
; Integrate panel to get 1D spectrum
slice, panel, space_integrated, scount
;Lags
lags = indgen(60)
;Error handling
blank = fltarr(60) + !Values.F_NaN

;Autocorrelation loop
for i=0,nosteps do begin
        catch, error_status
        if error_status ne 0 then begin
                openw, 1, fname, /append
                printf, 1, blank, form='(60f10.5)'
                close, 1
                continue
                catch, /cancel
        endif
        timeindices = [where(time eq i*60.),where(time eq ((i+1)*60.)-0.5)] ; Start,end of minute
	peak = space_integrated[timeindices[0]:timeindices[1],210:270]
	a = size(peak)
	b = a[1]
	; pp holds ithe wavelength [0,*] and intensity [1,*] of the spectrum at 'time' [i]
	pp = fltarr(2, 61)
	pp[0,*] = wl[210:270]
	blueshifts = fltarr(b)

	for j=0,b-1 do begin
        	pp[1,*] = peak[j,*]
        	centre = total(pp[1,*] * pp[0,*])/total(pp[1,*])
        	blueshifts[j] = 6562.81 - centre
	endfor

        ;Background remove
        x = findgen(120)
        background = linfit(x, blueshifts)
        timeseries = blueshifts - (background[0] + x*background[1])

        ;Make and save autocorrelation for minute of data
        a = size(timeseries)
        openw, 1, fname, /append
        if a[1] eq 120 then begin                           ; Check there is 60s of data
                acorr = a_correlate(timeseries, lags)
                printf, 1, acorr, form='(60f10.5)'
        endif else begin
                printf, 1, blank, form='(60f10.5)'          ; If not print row of NaNs
        endelse
        close, 1
endfor
end


;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
