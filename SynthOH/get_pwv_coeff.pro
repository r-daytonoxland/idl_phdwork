FUNCTION get_pwv_coeff, wl, plot=plot
;
; Get's a PWV coefficient as defined in Chadney et al 2017, by fitting a
; straight line to the H2O optical depth for different values of PWV.
;
; Inputs:
;  wl - The wavelength of the line we want the factor for, in A
;       Can be an array. Should be wavelength *in air* here.
; Returns:
;  Coefficient to quickly calculate transmission for an arbitrary
;  PWV value. See Chadney 2017, doi:10.5194/angeo-35-481-2017.
;
; Keywords:
;  plot - Make a plot of optical depth vs PWV, with the fit,
;       in the current plotting location. Each wavelength will
;       be shown in a different colour. The plot keyword is ignored
;       if you've asked for more than 40 wavelengths.
;

; Convert wavelengths in air to wavenumber
v=1/(air2vac(wl)*1d-8)
nwl=n_elements(v)

;if nwl gt 3 then print,'Calculating PWV coefficients for '+string(nwl,form='(i0)')+' wavelengths...',form='(a,$)'

; Get temperature (T, in K), pressure (P, in atm), and H2O number density
; (nH2O, in m^-3) altitude profiles. Also dz, altitude step in m.
; All corresponding to 1mm of PWV.
restore,file='~/lib/NYA_water_vapour.idl'

pwvs=[0.1,1,2,5,10,20,40]

odepth=dblarr(n_elements(pwvs),nwl)

for i=0,n_elements(pwvs)-1 do begin

 if nwl gt 3 then print,string(13b)+'Calculating PWV coefficients for '+string(nwl,form='(i0)')+' wavelengths...'+string(i*100/n_elements(pwvs),form='(i0)')+'%',form='(a,$)'

 ; Scale H2O number density by the PWV
 this_nH2O=nH2O*pwvs[i]

 ; Partial pressure of H2O, in atm
 pp=this_nH2O*constant(/kb)*T/1.01325d5

 ; Absorption cross section, dimensions (P, wl), converted to m^2
 sigma=hitran_absorption(v, T, P, pp)*1d-4

 ; Optical depth - here we do the integration in altitude.
 ; Output from total will have dimensions (nwl).
 odepth[i,*]=total(sigma*((this_nH2O*dz) # (dblarr(nwl)+1)),1)

endfor
if nwl gt 3 then print,string(13b)+'Calculating PWV coefficients for '+string(nwl,form='(i0)')+' wavelengths...'+string(i*100/n_elements(pwvs),form='(i0)')+'%',form='(a,$)'

; Now we have the optical depth, fit a straight line with PWV.
; We force it to go through zero, so the fit is trivial.
coeff=mean(odepth / (pwvs # (dblarr(nwl)+1)), dim=1)

if n_elements(plot) gt 0 and nwl le 40 then begin
 loadct,0,/silent
 plot,fltarr(3),/nodata,xran=[min(pwvs),max(pwvs)],/xsty,xtitle='PWV, mm',yran=[0,max(odepth)],/ysty,ytitle='Optical depth'
 loadct,39,/silent
 cols=round((findgen(nwl)+1)*254/nwl)
 for i=0L,nwl-1 do begin
  oplot,pwvs,odepth[*,i],psym=1,color=cols[i]
  oplot,!x.crange,!x.crange*coeff[i],color=cols[i]
 endfor
endif

if nwl gt 3 then print,string(13b)+'Calculating PWV coefficients for '+string(nwl,form='(i0)')+' wavelengths...done.'
return,coeff

end
