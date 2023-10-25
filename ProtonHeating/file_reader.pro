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

pro readall, file_start, intensity, blueshift, temperature, lib=lib

path = '~/lib/'

if keyword_set(lib) then begin
  fnames = [path + file_start + 'intensity.txt', path + file_start + 'blueshift.txt', path + file_start + 'temperature.txt']
endif else begin
  fnames = [file_start + 'intensity.txt', file_start + 'blueshift.txt', file_start + 'temperature.txt']
endelse

readsimple, fnames[0], intensity
readsimple, fnames[1], blueshift
readsimple, fnames[2], temperature

end