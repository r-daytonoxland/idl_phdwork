pro readfile, file, data

; Select a text file and open for reading
OPENR, lun, file, /GET_LUN
; Read one line at a time, saving the result into array
array = ''
line = ''
WHILE NOT EOF(lun) DO BEGIN & $
  READF, lun, line & $
  array = [array, line] & $
ENDWHILE
; Close the file and free the file unit
FREE_LUN, lun

data = float(array)

end

pro readsimple, filename, array

n_lines = file_lines(filename)
array = fltarr(n_lines)

openr, lun, filename, /get_lun
readf, lun, array
free_lun, lun

end

pro readall, file_start, interval_length, intensity, blueshift, temperature, lib=lib

path = '~/lib/'

if keyword_set(lib) then begin
  fnames = [path + file_start + 'intensity' + strtrim(interval_length, 2) + '.txt', path + file_start + 'blueshift' + strtrim(interval_length, 2) + '.txt', path + file_start + 'temperature' + strtrim(interval_length, 2) + '.txt']
endif else begin
  fnames = [file_start + 'intensity' + strtrim(interval_length, 2) + '.txt', file_start + 'blueshift' + strtrim(interval_length, 2) + '.txt', file_start + 'temperature' + strtrim(interval_length, 2) + '.txt']
endelse

readsimple, fnames[0], intensity
readsimple, fnames[1], blueshift
readsimple, fnames[2], temperature

end

pro readone, file_start, interval_length, index, data, lib=lib

path = '~/lib/'

if keyword_set(lib) then begin
  fnames = [path + file_start + 'intensity' + strtrim(interval_length, 2) + '.txt', path + file_start + 'blueshift' + strtrim(interval_length, 2) + '.txt', path + file_start + 'temperature' + strtrim(interval_length, 2) + '.txt']
endif else begin
  fnames = [file_start + 'intensity' + strtrim(interval_length, 2) + '.txt', file_start + 'blueshift' + strtrim(interval_length, 2) + '.txt', file_start + 'temperature' + strtrim(interval_length, 2) + '.txt']
endelse

if index eq 0 then begin
  readsimple, fnames[0], intensity
endif
if index eq 1 then begin
  readsimple, fnames[1], blueshift
endif
if index eq 2 then begin
  readsimple, fnames[2], temperature
endif

end