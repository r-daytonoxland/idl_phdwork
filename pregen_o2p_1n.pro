PRO pregen_o2p_1n, Ts, minwl, maxwl, O2p_i, O2p_w
;
; Pre-calculates the O2+ 1N synthetic spectra for temperatures given in Ts, and
; puts these in a common block. Then returns O2p_i and O2p_w.
; If the pre-calculation has already been done, just returns O2p_i and O2p_w
; from the common block.
;

common pregen_o2p_1n_cmn, o2p_Ts, o2p_minwl, o2p_maxwl, o2p_O2p_i, o2p_O2p_w

if (n_elements(o2p_Ts) gt 0) then if (array_equal(Ts,o2p_Ts) and minwl eq o2p_minwl and maxwl eq o2p_maxwl) then begin
 O2p_i=o2p_O2p_i
 O2p_w=o2p_O2p_w
 return
endif

nT=n_elements(Ts)
O2p_i=dblarr(nT,40000)*!values.F_NaN
O2p_w=O2p_i
mc=0
for i=0,nT-1 do begin
 print,string(13b)+"Generating O2+ 1N spectra: "+string(i,form='(i0)')+'/'+string(nT,form='(i0)')+" temperatures...",form='(a,$)'
 prog_o2_1n_v3_exp,90d,Ts[i],_O2p_i,_O2p_w,/newc,/nocorr;,minwl,maxwl
 s=where(_O2p_w gt minwl and _O2p_w lt maxwl,c)
 O2p_i[i,0:c-1]=_O2p_i[s]
 O2p_w[i,0:c-1]=_O2p_w[s]
 if c gt mc then mc=c
endfor
print,string(13b)+"Generating O2+ 1N spectra: "+string(i,form='(i0)')+'/'+string(nT,form='(i0)')+" temperatures...",form='(a,$)'
O2p_i=O2p_i[*,0:mc]
O2p_w=O2p_w[*,0:mc]
o2p_Ts=Ts
o2p_minwl=minwl
o2p_maxwl=maxwl
o2p_O2p_i=O2p_i
o2p_O2p_w=O2p_w
print,"...done."

end
