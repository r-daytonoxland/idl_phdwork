PRO read_hitran, parfile=parfile, force=force
;
; Reads a HITRAN ".par" file in 160-byte fixed width format into the common
; block, for use with other HITRAN-related routines.
;
; Keyword parfile is the file to read, default ~/lib/HITRAN_H2O.par
;
; If parfile matches the name of the file already read in to the common block
; then this routine returns without doing anything, unless the force keyword
; is set, in which case it will reread the file.
;
; Makes use of information at https://hitran.org/docs/definitions-and-units/
;

if n_elements(parfile) lt 1 then parfile='~/lib/HITRAN_H2O.par'
n=40000L ; Maximum 40000 lines in the file
mol=1 ; 1 for H2O, 2 for CO2

common _hitran,ht_file,ht_v,ht_S,ht_A,ht_gair,ht_gself,ht_E,ht_nair,ht_dair

if (n_elements(force) eq 0 and n_elements(ht_file) gt 0) && (ht_file eq parfile) then return ; Don't bother reading again
ht_file=parfile

; Things read from file:
mid=0 ; Molecule ID
iid=0 ; Isotopologue ID
v=0d  ; Wavenumber in vacuum (cm^-1)
S=0d  ; Spectral line intensity (cm^-1/(molecule cm^-2)) at T=296K
A=0d  ; Einstein coefficient (s^-1)
gair=0d ; Air-broadened HWHM (cm^-1/atm) at T=296K and P=1atm
gself=0d ; Self-broadened HWHM (cm^-1/atm) at T=296K and P=1atm
E=0d  ; Lower-state energy of the transition (cm^-1)
nair=0d ; Coefficient of temperature dependence of air-broadened half width
dair=0d ; Pressure shift (cm^-1/atm) at T=296K and P=1atm of the line position

ht_v=dblarr(n)
ht_S=dblarr(n)
ht_A=dblarr(n)
ht_gair=dblarr(n)
ht_gself=dblarr(n)
ht_E=dblarr(n)
ht_nair=dblarr(n)
ht_dair=dblarr(n)

openr,1,parfile
i=0L
while not(eof(1)) do begin
 aa=''
 readf,1,aa
 reads,aa,mid,iid,v,S,A,gair,gself,E,nair,dair,form='(i2,i1,f12.6,e10.3,e10.3,f5.4,f5.3,f10.4,f4.2,f8.6)'
 if (mid ne mol) or (iid ne 1) then continue ; Only care about the first isotopologue of the molecule we want
 ht_v[i]=v
 ht_S[i]=S
 ht_A[i]=A
 ht_gair[i]=gair
 ht_gself[i]=gself
 ht_E[i]=E
 ht_nair[i]=nair
 ht_dair[i]=dair
 i+=1
endwhile
close,1

ht_v=ht_v[0:i-1]
ht_S=ht_S[0:i-1]
ht_A=ht_A[0:i-1]
ht_gair=ht_gair[0:i-1]
ht_gself=ht_gself[0:i-1]
ht_E=ht_E[0:i-1]
ht_nair=ht_nair[0:i-1]
ht_dair=ht_dair[0:i-1]

end
