FUNCTION read_specfit_file, filename, cols=cols
;
; Routine to read a file output by specfit, with parameters like OH temperature
; etc.
;
; Inputs:
;  filename - The file to be read.
;
; Outputs:
;  dat - The data read from the file.
;
; Keywords:
;  cols - Set to a named variable to receive the column headings (string array).
;

errortags=['-','p','n','N','v','t','f','pv','Nv']

openr,lun,filename,/get_lun

; Read (and throw away) the header
a=' '
while ((strmid(a,0,3) ne 'mjs') && ~eof(lun)) do readf,lun,a
; Extract column headings and find the errortag (which is non-numeric)
cols=strsplit(a,/extract)
errloc=where(strcmp(cols,'errortag',8) eq 1,c)
if (c ne 1) then message,'Unable to find errortag'
errloc=errloc[0]

; Read the rest of the data
while ~eof(lun) do begin
 readf,lun,a ; a is a whole line of the file as a string
 b=strsplit(a,/extract)
 if n_elements(b) ne 23 then begin
  message,'Line has '+string(n_elements(b),form='(i0)')+' columns, expected 23.'
 endif else begin

  ; Convert error tag to number
  b[errloc]=where(strcmp(errortags,b[errloc]))

  ; Find bad data (column width not big enough) to set to NaN later
  missing=where((strcmp(b,'***',3) eq 1) or (strcmp(b,'Inf',3) eq 1),cmiss)
  if (cmiss gt 0) then b[missing]='0'
  
  ; Convert data to doubles and set missing values to NaN
  db=double(b)
  if (cmiss gt 0) then db[missing]=!values.F_NaN
  
  if (n_elements(dat) eq 0) then dat=db else dat=[[dat],[db]]
  
 endelse
endwhile
close,lun
free_lun,lun

return,dat

end
