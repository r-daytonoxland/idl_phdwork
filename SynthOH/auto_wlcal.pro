pro auto_wlcal,mjs0,time,dseq,pnum,min_wl,max_wl,p0,wout,$
               mid=mid,hydroxyl=hydroxyl,width=width,test=test,$
	       solar=solar,kurucz2005=kurucz2005,ctype=ctype, $
	       nooffset=nooffset, slice=slice, $
	       correct_banding=correct_banding
;
; Interactive routine to wavelength calibrate the spectrograph
;
; IMPORTANT: Modified by DKW May 2015, to generate a new style of parameters using a sinusoidal
; component (type 2) instead of a 3rd order polynomial (type 1). For type 2 parameters
; wout contains 5 elements.
;
; Type 1: wl = w0 + w1(x - p0) + w2(x - p0)^2 + w3(x - p0)^3
; Type 2: wl = w0 + w1(x - p0) + w2 sin(w3(x - p0) + w4)
; (These are definied in wl_function.pro)
;
; Inputs:
;  mjs0,time,dseq  variables as produced by a read_tim for the time to use as calibration
;                   If calibrating using the solar spectrum you should obviously choose a period
;                   in daylight. This is usually easier than using the hydroxyl lines, which can
;                   be difficult to see particularly at the shorter wavelengths. Hydroxyl lines
;                   are a result of airglow so a clear period at night should be chosen. It is
;                   probably a good idea to attempt the calibration using both methods and compare
;                   the results.
;  pnum            panel number
;  min_wl, max_wl  supposed minimum and maximum wavelengths of filter (in angstroms)
;  p0              central pixel to use in the calibration. If keyword mid is set this is an output,
;                   and equals half the pixel width of the panel.
;
; Outputs:
;  wout            calibration coefficients produced (4 element array)
;                   wout[0] = central wavelength
;
; Keywords:
;  width           width to use when convolving theoretical spectrum with sif spectrum (default 0.7)
;  hydroxyl        set this keyword to calibrate with a hydroxyl air glow spectrum instead of
;                   a solar spectrum. As of May 2015 Hydroxyl is the default.
;                   Set to 2 to use "synth_oh" instead of "ohspectra".
;                   If set to something larger than 4, then synth_oh is used
;                   and hydroxyl is passed through to synth_oh in the upperv
;                   keyword, to select a specific upper vibrational level.
;  solar           set this keyword to calibrate with the solar Fraunhofer spectrum instead
;                   of the OH spectrum.
;  kurucz2005      set this instead of solar, to use the Kurucz 2005 reference solar spectrum.
;  test            if this keyword is set then wout and p0 become inputs, and the calibration is
;                   used to plot the theoretical spectrum and the calibrated sif spectrum. This is
;                   for checking coefficients determined at a different time or using a different
;                   method etc.
;  ctype           Specify the calibration formula type. Original 3rd order
;                   polynomial is type 1. May 2015 formula with a sinusoidal
;                   component is type 2. See wl_function.
;  slice           Set to a 2-element array to specify the slicing. If not
;                   set then get_def_slc is used to choose the slicing instead.
;  correct_banding If set then the spectrum is cleaned with the remove_banding routine,
;                  correct_banding is passed through to remove_banding as the type keyword.
;
; DKW added nooffset keyword 22/02/2022 (Twosday!), passed through to spectra.
; DKW added slice keyword, 15/11/2023.
;
if not(keyword_set(ctype)) then ctype=2
get_size, mjs0,nx,ny
get_panels,mjs0,np,xx,yy
p_width=float(yy[pnum-1,1]+1-yy[pnum-1,0])
;x_width=xx[pnum-1,1]+1-xx[pnum-1,0]
;left=xx[pnum-1,0]+fix(x_width*0.2)
;right=xx[pnum-1,0]+fix(x_width*0.8)
;if ((nx/2) le left) then slice=[0.2,0.2] else begin
; if ((nx/2) ge right) then slice=[0.8,0.2] else slice=[((nx/2)-xx[pnum-1,0])/float(x_width),0.2]
;endelse
if (n_elements(slice) ne 2) then begin
 get_def_slc,mjs0,slices
 slice=reform(slices[pnum-1,*])
endif
if keyword_set(mid) then p0=p_width/2.0
if not(keyword_set(width)) then width=0.7
if keyword_set(solar) then begin
 type='solar'
 read_sol
 common solarsp, solar_wl, solar_in, solar_cnt, solar_flux
 wl=solar_wl
 in=solar_in
endif else if keyword_set(kurucz2005) then begin
 type='solar'
 read_sol,/kurucz2005
 common solarsp, solar_wl, solar_in
 wl = solar_wl
 in = solar_in
endif else begin
 type='hydroxyl'
 if keyword_set(hydroxyl) and hydroxyl ge 2 then begin
  if (hydroxyl gt 4) then $
   synth_oh, [min_wl,max_wl], 190d, wl,in,dummy,ohlabels,upperv=hydroxyl $
  else $
   synth_oh, [min_wl,max_wl], 190d, wl,in,dummy,ohlabels
  ; Select the wavelengths for only bright enough lines, and merge e and f components
  ohsel=where(in gt 1d-3)
  ef=strmid(ohlabels[ohsel],0,1,/reverse)
  eonly=where(ef eq 'e',comp=fonly)
  sorte=sort(ohlabels[ohsel[eonly]])
  sortf=sort(ohlabels[ohsel[fonly]])
  ohwl=((wl[ohsel[eonly[sorte]]]*in[ohsel[eonly[sorte]]])+(wl[ohsel[fonly[sortf]]]*in[ohsel[fonly[sortf]]]))/(in[ohsel[eonly[sorte]]]+in[ohsel[fonly[sortf]]])
 endif else begin
  ohspectra,min_wl,max_wl, wl,in
  common oh, ohband, ohk, ohp1, ohp2
  ohwl=[ohp1,ohp2]
  ;wl=findgen(max_wl-min_wl)+min_wl
 endelse
endelse
sifpx=findgen(p_width)
sifwl=(findgen(p_width)*(max_wl-min_wl)/p_width)+min_wl
if not(keyword_set(test)) then convolve_sp,wl,in,width,sifwl,solarspec
mjs_tt,mjs0,yr,mo,da,hr,mi,se
tstring=string(da,mo,yr,hr,mi,se,form='(2(i2.2,"/"),i4.4, i3.2,2(":",i2.2))')
tlength=(((time[n_elements(time)-1]-time[0])/n_elements(time))+time[n_elements(time)-1]-time[0])/60.0
print,'slice=',slice
print,'timeint=',tstring
print,'length=',tlength
spectra,pnum,mjs0,time,dseq,sifspectrum,str,slice=slice,/percentff,nooffset=nooffset;,timeint=tstring,length=tlength
if n_elements(correct_banding) eq 1 then begin
 remove_banding, sifspectrum, spectrummasked, type=correct_banding, /force
 sifspectrum=spectrummasked
endif
;sifspectrum=smooth(sifspectrum,5)
pos=[0.07,0.1,0.95,0.95]
x_ran=[min_wl,max_wl]
if not(keyword_set(test)) then begin
 y_ran=[min([min(solarspec/max(solarspec)),min(sifspectrum/max(sifspectrum))])-0.02,1.02]
 y_ransol=[min(solarspec/max(solarspec))-0.02,1.02]
 y_ransif=[min(sifspectrum/max(sifspectrum))-0.02,1.02]
 window,2
 !p.position=pos
 plot,sifwl,solarspec/max(solarspec),xrange=x_ran,xstyle=1,yrange=y_ransol,ystyle=1,title=type+' spectrum'
 print, 'Click on distinctive points on the sif spectrum.'
 print, 'You should be able to see them in the '+type+' spectrum as well.'
 print, 'Click in the left of the window to stop.'
 window,1
 !p.position=pos
 plot,sifpx,sifspectrum/max(sifspectrum),xrange=[0,p_width],xstyle=1,yrange=y_ransif,ystyle=1,title='sif spectrum'
 xin=p_width/2
 sifpos=fltarr(2)
 while xin gt 0 do begin
  cursor,xin,yin,/data,/down
  if xin gt 0 then begin
   if (strcmp(type,'hydroxyl') eq 1) then begin
    ; If we're using a hydroxyl spectrum then fit a Gaussian
    ; to find a more accurate value for the line centre
    dummy=min(abs(sifpx-xin),s)
    ran=[max([0,s[0]-10]),min([s[0]+10,n_elements(sifpx)-1])]
    dummy=gaussfit(sifpx[ran[0]:ran[1]],sifspectrum[ran[0]:ran[1]]/max(sifspectrum),Gfit,nterms=4,estimates=[yin,xin,1.7,median(sifspectrum/max(sifspectrum))])
    loadct,39,/silent
    oplot,sifpx[ran[0]:ran[1]],dummy,color=250
    loadct,0,/silent
    xin=Gfit[1]
   endif
   print,xin
   sifpos=[[sifpos],[xin,yin]]
  endif
  plots,sifpos(0,*),sifpos(1,*),/data,psym=4
 endwhile
 sifpos=sifpos(*,1:*)
 print, 'Now match these points tn the '+type+' spectrum.'
 window,2
 !p.position=pos
 plot,sifwl,solarspec/max(solarspec),xrange=x_ran,xstyle=1,yrange=y_ransol,ystyle=1,title=type+' spectrum'
 xin=((max_wl-min_wl)/2)+min_wl
 solarpos=fltarr(2)
 x_norm=(pos[2]-pos[0])/p_width
 x_norm2=pos[0]
 y_norm=(pos[3]-pos[1])/(y_ransol[1]-y_ransol[0])
 y_norm2=pos[1]-(y_norm*y_ransol[0])
 for i=1,n_elements(sifpos(0,*)) do begin
  wset,1
  if i ge 2 then plots,(sifpos(0,i-2)*x_norm)+x_norm2,(sifpos(1,i-2)*y_norm)+y_norm2,psym=4,/norm
  loadct,1,/silent
  plots,(sifpos(0,i-1)*x_norm)+x_norm2,(sifpos(1,i-1)*y_norm)+y_norm2,color=200,psym=4,/norm
  loadct,0,/silent
  wset,2
  cursor,xin,yin,/data,/down
  if (strcmp(type,'hydroxyl') eq 1) then begin
   ; Match to actual OH wavelength
   dummy=min(abs(ohwl-xin),s)
   xin=ohwl[s]
  endif
  print,xin
  solarpos=[[solarpos],[xin,yin]]
  plots,solarpos(0,*),solarpos(1,*),/data,psym=4
 endfor
 print, 'Calculating...'
 solarpos=solarpos(*,1:*)
 sifpix2=sifpos[0,*];-p0
 solarpos2=solarpos[0,*]
 ;wout=poly_fit(sifpix2,solarpos2,3)
 case ctype of
  1: begin
      wout=[double(ctype),double(p0),((max_wl-min_wl)/2.0)+min_wl,dblarr(3)]
      fita=[0,0,1,1,1,1]
      noder=1
     end
  2: begin
      wout=[double(ctype),double(p0),((max_wl-min_wl)/2.0)+min_wl,dblarr(4)]
      fita=[0,0,1,1,0,0,0]
      noder=0
      dummy=curvefit(reform(sifpix2),reform(solarpos2),dblarr(n_elements(sifpix2))+1,wout,fita=fita,function_name='wl_function',/double)
      wout[4:5]=0.1
      fita=[0,0,1,1,1,1,1]
     end
  else: begin
      wl_function,ctype,np,/numparams
      wout=[double(ctype),double(p0),dblarr(np-2)]
      fita=[0,0,intarr(np-2)+1]
      noder=0
     end
 endcase
 dummy=curvefit(reform(sifpix2),reform(solarpos2),dblarr(n_elements(sifpix2))+1,wout,fita=fita,noder=noder,function_name='wl_function',/double)
 ;wout=wout[2:*]
endif
;w0=wout[0]
;w1=wout[1]
;w2=wout[2]
;w3=wout[3]
;sifwl2=w0+((findgen(p_width)-p0)*w1)+(((findgen(p_width)-p0)^2)*w2)+(((findgen(p_width)-p0)^3)*w3)
wl_function,findgen(p_width),wout,sifwl2;[double(ctype),double(p0),wout],sifwl2
convolve_sp,wl,in,width,sifwl2,solarspec2
if not(keyword_set(test)) then begin
 window,3
 plot,sifpos[0,*],solarpos2,psym=4,title='Fit Done',xtitle='Pixel number',ytitle='Wavelength',yrange=x_ran,xrange=[0,p_width],xstyle=1
 oplot,sifpx,sifwl2
endif else y_ran=[min([min(solarspec2/max(solarspec2)),min(sifspectrum/max(sifspectrum))])-0.02,1.02]
loadct,1,/silent
window,4
plot,sifwl2,sifspectrum/max(sifspectrum),xrange=x_ran,yrange=y_ran,xstyle=1,ystyle=1,color=250,title='Blue: '+type+' spectrum, White: calibrated sif spectrum' 
oplot,sifwl2,(solarspec2/max(solarspec2)),color=200
loadct,0,/silent
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;