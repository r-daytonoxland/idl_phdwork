PRO pregen_n2, Ts, minwl, maxwl, N2_i, N2_w, vibband
;
; Pre-calculates the N2 synthetic spectra for temperatures given in Ts, and
; puts these in a common block. Then returns N2_i, N2_w, and vibband.
; If the pre-calculation has already been done, just returns the values from
; the common block.
;  N2_i - Intensities of N2 1P lines
;  N2_w - Wavelengths of N2 1P lines (in A)
;  vibband - String array of vibrational bands of the N2 1P lines
; Note vibband is a 1d array. The others are 2d arrays, first dimension temperature.
;

common pregen_n2_cmn, n2s_Ts, n2s_minwl, n2s_maxwl, n2s_N2_i, n2s_N2_w, n2s_vibband

if (n_elements(n2s_Ts) gt 0) then if (array_equal(Ts,n2s_Ts) and minwl eq n2s_minwl and maxwl eq n2s_maxwl) then begin
 N2_i=n2s_N2_i
 N2_w=n2s_N2_w
 vibband=n2s_vibband
 return
endif

nT=n_elements(Ts)
N2_i=dblarr(nT,40000)*!values.F_NaN
N2_w=N2_i
mc=0
for i=0,nT-1 do begin
 print,string(13b)+"Generating N2 spectra: "+string(i,form='(i0)')+'/'+string(nT,form='(i0)')+" temperatures...",form='(a,$)'
 prog_n2_1p_v3_exp,90d,Ts[i],_N2_i,_N2_w,minwl,maxwl,/newc
 s=where(_N2_w gt minwl and _N2_w lt maxwl,c)
 N2_i[i,0:c-1]=_N2_i[s]
 N2_w[i,0:c-1]=_N2_w[s]
 ; Note for determining the vibrational band we assume the wavelengths are the same for every temperature...
 ivib=array_indices(_N2_w,s)
 vibband=strjoin(strtrim(transpose([[string(ivib[1,*],form='(i2)')],[string(ivib[2,*],form='(i2)')]]),2),'-')
 if c gt mc then mc=c
endfor
print,string(13b)+"Generating N2 spectra: "+string(i,form='(i0)')+'/'+string(nT,form='(i0)')+" temperatures...",form='(a,$)'
N2_i=N2_i[*,0:mc-1]
N2_w=N2_w[*,0:mc-1]
n2s_Ts=Ts
n2s_minwl=minwl
n2s_maxwl=maxwl
n2s_N2_i=N2_i
n2s_N2_w=N2_w
n2s_vibband=vibband
print,"...done."

end
