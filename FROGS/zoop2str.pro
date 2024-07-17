PRO zoop_str2frame, str, megablock, frame, print=print
;
; Does the reverse of zoop_frame2str, and extracts the megablock name and frame
; number from a short string, e.g. a zoo image filename.
;
; Input:
;  str - A string encoding megablock and frame, as returned by zoop_frame2str.
;       If the string is longer than normal, e.g. because it has a filename
;       extension, then this routine will still work fine.
;
; Outputs:
;  megablock - The megablock name, as a string. Only the date-time part.
;  frame - The frame number, as a long.
;
; Currently only works for scalar inputs.
;

alpha=strjoin(['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'])
numalpha=strjoin(['0','1','2','3','4','5','6','7','8','9',alpha])

yr=strpos(numalpha,strmid(str,0,1))+2000
mo=strpos(alpha,strmid(str,1,1))+1
da=strpos(numalpha,strmid(str,2,1))
hr=strpos(alpha,strmid(str,3,1))
mise=strmid(str,4,4)
megablock=string(yr,mo,da,hr,mise,form='(i04,i02,i02,i02,a4)')

frame=0L
;reads,strmid(str,8,4),frame,form='(Z04)'
frame=long(strmid(str,8,5))

if keyword_set(print) then print,megablock,frame

end