FUNCTION hitran_intensity, v, E, S, T, Qfile=Qfile
;
; Calculates the temperature-dependant line intensity, see equation 4
; of https://hitran.org/docs/definitions-and-units/
;
; Inputs:
;  v - Line wavenumber, cm^-1
;  E - Lower-state energy of the transition (cm^-1)
;  S - Spectral line intensity (cm^-1/(molecule cm^-2)) at T=296K
;  All of the above are as read in from HITRAN, in cgs units
;  T - Desired temperature in K
; Output:
;  S(T), i.e. line intensity at T
;
; Keywords:
;  Qfile - The location of the file containing data on Q, the total
;      internal partition sum, probably generated from TIPS. This should
;      have 2 columns, first temperature in K, second Q. If the file
;      has already been read by this routine, it is not reread.
;
; All inputs can be vector arrays. v, E, S should all have the same number
; of elements, m. T can have a different number of elements, n.
; The output is an array of dimensions (m,n), i.e. (lines, temperatures).
;

common _hitran_intensity, ht_Q_file,ht_QT,ht_Q

if n_elements(Qfile) lt 1 then Qfile='~/lib/HITRAN_H2O_Q.dat'
; Read in the file of total internal partition sum data
if n_elements(ht_Q_file) le 0 || ht_Q_file ne Qfile then begin
 openr,lun,Qfile,/get_lun
 ht_QT=dblarr(10000)
 ht_Q=dblarr(10000)
 i=0
 while not(eof(lun)) do begin
  thisT=0d
  thisQ=0d
  readf,lun,thisT,thisQ
  ht_QT[i]=thisT
  ht_Q[i]=thisQ
  i+=1
 endwhile
 ht_QT=ht_QT[0:i-1]
 ht_Q=ht_Q[0:i-1]
 close,lun
 ht_Q_file=Qfile
endif

c2=1.4387769d0 ; hc/k in cgs units, cm K
Tref=296d0 ; Reference temperature in K
Qref=interpol(ht_Q,ht_QT,Tref) ; total internal partition function at Tref
Q=interpol(ht_Q,ht_QT,T) ; Q at all of our T

n=n_elements(v)
m=n_elements(T)

; part of the equation with dimensions (n,m)
part1=exp((-c2*E) # (1d0/T))*(1-exp((-c2*v) # (1d0/T)))/((dblarr(n)+1) # Q)

; part of the equation with dimensions (n) only
part2=S*Qref/(exp(-c2*E/Tref)*(1-exp(-c2*v/Tref)))

return,part1*(part2 # (dblarr(m)+1))

end
