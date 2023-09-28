pro total_box, im, x, y, dim, tot_box 

; Procedure integrates the auroral brightness in the box defined by the input.
; Inputs;
;        im - the ASK image
;        x, y - The pixel location of the magnetic zenith
;        dim - the pizel dimensions of the box, 10px default????

r = dim/2
; round x and y to whole pixel
xx = fix(x) & yy = fix(y)
box = im[*,xx-r:xx+r, yy-r:yy+r]
tot_box = total(box)

end

;;;;;;;;;;;;

pro get_box_timeseries, fn, length, x, y, timeseries

dim=10
read_mv, fn, length, ims
timeseries = make_array(length)

for i=0,length-1 do begin
	total_box, ims[i,*,*], x, y, dim, tot_box
	timeseries[i] = tot_box
endfor
end

;;;;;;;;;;;;;

pro read_askdat, fname, duration, length, dat

nosteps = duration/length

dat = fltarr(length, nosteps)

openr, 1, fname
readf, 1, dat
close, 1

end

;;;;;;;;;;;;;

pro ASK_ACF, start, duration, length, x, y, fname

nosteps = duration/length
lags = indgen(length)

for i=0, nosteps-1 do begin
        get_box_timeseries, start, length, x, y, timeseries
        start += length

	range=findgen(length)
	background = linfit(range, timeseries)
	timeseries = timeseries - (background[0] + range*background[1])

        acf = a_correlate(timeseries, lags)

        openw, 1, fname, /append
        printf, 1, acf, form='(400f10.5)'
        close, 1
endfor
end

;;;;;;;;;;;;;

pro ASK_bootstrap, start, length, x, y, max

bootstrap = 1000

lags = indgen(length)
max_amplitudes=fltarr(bootstrap,length)

get_box_timeseries, start, length, x, y, timeseries
range=findgen(length)
background = linfit(range, timeseries)
timeseries = timeseries - (background[0] + range*background[1])

for i=0, bootstrap-1 do begin
	b = randomu(dseed,length)
	a = timeseries[sort(b)]
	acf = a_correlate(a, lags)
	max_amplitudes[i,*] = acf
endfor

max = max(abs(max_amplitudes), dimension=1)

plot, lags/20., max

end

;;;;;;;;;;;;;;

pro ASK_FFT, start, duration, length, x, y, fname

nosteps = duration/length
lags = indgen(length)
fn = start
hwindow = hanning(length, /double)

for i=0, nosteps-1 do begin
        get_box_timeseries, fn, length, x, y, timeseries
        fn += length

        range=findgen(length)
        background = linfit(range, timeseries)
        timeseries = timeseries - (background[0] + range*background[1])

        fft = fft(timeseries*hwindow)

        openw, 1, fname, /append
        printf, 1, abs(fft), form='(' + strtrim(length,1) + 'f10.5)'
        close, 1
endfor
end

;;;;;;;;;;;;;;;

pro quickgo, x, y
; Quick setup while I'm working on certain files and times
; Read_vs, copy in megablock of interest
; Irgf_senith, get zenith pixel numbers
; times - length of ACF or FFT time-step in seconds

read_vs,file='20160102060659r1.txt'
igrf_zenith,2016.01d,78.153d,16.029d,0.4d,200d,az,el
get_cnv, cnv
conv_xy_ae,x,y,az,el,cnv,/back


end

;;;;;;;;;;;;;;;

pro ASK_spectrogram, start, duration, length, x, y, fname

fn = start
hwindow = hanning(length, /double, alpha=0.54) ;Hann window a=0.54 apparently

for i=0, duration-1 do begin
        get_box_timeseries, fn, length, x, y, timeseries
        fn += 1

        range=findgen(length)
        background = linfit(range, timeseries)
        timeseries = timeseries - (background[0] + range*background[1])

        fft = fft(timeseries*hwindow)

        openw, 1, fname, /append
        printf, 1, real_part(fft), form='(' + strtrim(length,1) + 'f10.5)'
        close, 1
endfor

end

