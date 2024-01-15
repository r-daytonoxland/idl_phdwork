Forward_Function Sigma4F1,Pi4F1,Sigma4F2,Pi4F2,Sigma4F3,Pi4F3,Sigma4F4,Pi4F4,nuva2leai,njp

pro prog_O2_1N_v3_exp,lines,temp,all_O2_i,all_O2_w,newc=newc,nocorr=nocorr

;Outputs are: 
;       all_O2_w    - Wavelengths of line spectra      - [Lambda,v1,v2]
;       all_O2_i    - Intensities of line spectra      - [   Int,v1,v2]
;       instru_grid - Wavelengths of instrument grid   - {Lambda]
;       int_colv    - Convolved intensities
;                     with instrument function         - [   Int]

;MAIN INPUTS

                 brcs = 48.     ;Number of branches
                 v1_1 = 0.      ;V' First vibrational level [ v' from 0 to 3]
                 v1_2 = 3.      ;V'  Last vibrational level [ v' from 0 to 3]
                 v2_1 = 0.      ;V" First vibrational level [ v" from 0 to 10]
                 v2_2 = 10.     ;V"  Last vibrational level [ v" from 0 to 10]
                 plt1 = 1.      ;plot spectra   [1-on 0-off]
                 plt2 = 1.      ;plot reference [1-on 0-off] *FIX or REMOVE
                 plt3 = 1.      ;plot HL        [1-on 0-off] *FIX or REMOVE

;READ DATA OF VIBRATIONAL LEVELS
O2_vib1, O2_1N_vib_ints 

;SET UNCONVOLVED ARRAYS
all_O2_w=dblarr(lines*brcs,4,11)
all_O2_i=dblarr(lines*brcs,4,11)

;if (v1_2-v1_1+1. gt 2) or (v2_2-v2_1+1. gt 2) then $
; print, string((v1_2-v1_1+1.)*(v2_2-v2_1+1.), form='("   Calculating ",i3.0," O2+ 1N bands...")')
;print,round(v1_2-v1_1+1.)*round(v2_2-v2_1+1.)," O2+ 1N bands to calculate..."
;SET V' AND V'' LOOP
for v1= v1_1,v1_2 do begin  ;Upper b-State vibrational band
for v2= v2_1,v2_2 do begin  ;Lower a-State vibrational band

;CREATE SPECTRA
O2_1N,v1,v2,lines,temp, newc=newc, nocorr=nocorr, $
    R11,Q11,R21,Q21,R31,Q31,R41,Q41,P11,R12,P21,R22,P31,R32,P41,R42,Q12,$
    P12,Q22,P22,Q32,P32,Q42,P42,R13,R23,R33,R43,Q13,Q23,Q33,Q43,P13,P23,$
    P33,P43,R14,R24,R34,R44,Q14,Q24,Q34,Q44,P14,P24,P34,P44,R11_i,Q11_i,$
    R21_i,Q21_i,R31_i,Q31_i,R41_i,Q41_i,P11_i,R12_i,P21_i,R22_i,P31_i,$
    R32_i,P41_i,R42_i,Q12_i,P12_i,Q22_i,P22_i,Q32_i,P32_i,Q42_i,P42_i,$
    R13_i,R23_i,R33_i,R43_i,Q13_i,Q23_i,Q33_i,Q43_i,P13_i,P23_i,P33_i,$
    P43_i,R14_i,R24_i,R34_i,R44_i,Q14_i,Q24_i,Q34_i,Q44_i,P14_i,P24_i,$
    P34_i,P44_i   ;,al

;MODEL vs THEORY PLOT
;if (plt2 eq 1) and (v1 eq 0) and (v2 eq 0) then Dieke_2P,P1,R1,P2,R2,Q2,P3,R3,Q3,al

;GROUP BAND SPECTRA TOGETHER
group_all1,lines,brcs,R11,Q11,R21,Q21,R31,Q31,R41,Q41,P11,R12,P21,R22,P31,R32,P41,R42,Q12,$
    P12,Q22,P22,Q32,P32,Q42,P42,R13,R23,R33,R43,Q13,Q23,Q33,Q43,P13,P23,$
    P33,P43,R14,R24,R34,R44,Q14,Q24,Q34,Q44,P14,P24,P34,P44,R11_i,Q11_i,$
    R21_i,Q21_i,R31_i,Q31_i,R41_i,Q41_i,P11_i,R12_i,P21_i,R22_i,P31_i,$
    R32_i,P41_i,R42_i,Q12_i,P12_i,Q22_i,P22_i,Q32_i,P32_i,Q42_i,P42_i,$
    R13_i,R23_i,R33_i,R43_i,Q13_i,Q23_i,Q33_i,Q43_i,P13_i,P23_i,P33_i,$
    P43_i,R14_i,R24_i,R34_i,R44_i,Q14_i,Q24_i,Q34_i,Q44_i,P14_i,P24_i,$
    P34_i,P44_i,group_w,group_i

all_O2_w(*,v1,v2)=group_w                            ;Add wavelength to array
all_O2_i(*,v1,v2)=group_i*O2_1N_vib_ints(v2,v1)     ;Add normalised profile * IBC3 intensity

endfor     ;end of lower state v''
endfor     ;end of higher state v'



;;CREATE WAVELENGTH GRID
;instru_grid=minw+lindgen((maxw-minw)/grid)*grid
;
;;INSTRUMENT FUNCTION
;if (v1_2-v1_1+1. gt 2) or (v2_2-v2_1+1. gt 2) then print, " Convolving O2+ 1N...", round(maxw-minw), " Angstrom wavelength range"
;convolve_sp, all_O2_w, all_O2_i, instrument_function, instru_grid, int_colv
;
;;Outputs are: 
;;
;;       instru_grid - Wavelengths
;;       int_colv    - Convolved intensities
;
;
;;IF YOU WANT PLOTS
;if plt1 eq 1 then hk1,instru_grid,int_colv
;;if plt3 eq 1 then plot_HL
;

;print, v1_2, v1_1, v2_2
;TRYING THIS OUT ###############################################################
;if (v1_2 and v1_1 eq 1) and (v2_2 eq 0) then begin

;investigate_10,R11,Q11,R21,Q21,R31,Q31,R41,Q41,P11,R12,P21,R22,P31,R32,P41,$
;    R42,Q12,P12,Q22,P22,Q32,P32,Q42,P42,R13,R23,R33,R43,Q13,Q23,Q33,Q43,$
;    P13,P23,P33,P43,R14,R24,R34,R44,Q14,Q24,Q34,Q44,P14,P24,P34,P44
;endif
;
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro O2_vib1, O2_1N_vib_ints ;THESE ARE THE VIBRATIONAL BAND INTENSITIES
                            ;Vallance Jones 1974, can be modelled with T_vib
; v''0      1      2      3      4      5     6     7     8     9    10
O2_1N_vib_ints=$                                                             ;v'
[[  3.98,  3.47,  1.90,   .85,   .33,   .12,  .04,  .01,  .00,  .00,  .00],$ ;0
 [  5.90,   .23,   .38,   .86,   .74,   .43,  .21,  .09,  .04,  .01,  .01],$ ;1
 [  1.84,  1.08,   .65,   .02,   .08,   .19,  .18,  .12,  .07,  .03,  .01],$ ;2
 [  0.20,   .90,   .06,   .23,   .07,   .00,  .02,  .04,  .04,  .03,  .02]]  ;3
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro group_all1,lines,brcs,R11,Q11,R21,Q21,R31,Q31,R41,Q41,P11,R12,P21,R22,P31,R32,P41,R42,Q12,$
    P12,Q22,P22,Q32,P32,Q42,P42,R13,R23,R33,R43,Q13,Q23,Q33,Q43,P13,P23,$
    P33,P43,R14,R24,R34,R44,Q14,Q24,Q34,Q44,P14,P24,P34,P44,R11_i,Q11_i,$
    R21_i,Q21_i,R31_i,Q31_i,R41_i,Q41_i,P11_i,R12_i,P21_i,R22_i,P31_i,$
    R32_i,P41_i,R42_i,Q12_i,P12_i,Q22_i,P22_i,Q32_i,P32_i,Q42_i,P42_i,$
    R13_i,R23_i,R33_i,R43_i,Q13_i,Q23_i,Q33_i,Q43_i,P13_i,P23_i,P33_i,$
    P43_i,R14_i,R24_i,R34_i,R44_i,Q14_i,Q24_i,Q34_i,Q44_i,P14_i,P24_i,$
    P34_i,P44_i,group_w,group_i

group_w=dblarr(brcs*lines)
group_w(        0 : (1*lines)-1) =R11
group_w((1*lines) : (2*lines)-1) =Q11
group_w((2*lines) : (3*lines)-1) =R21
group_w((3*lines) : (4*lines)-1) =Q21
group_w((4*lines) : (5*lines)-1) =R31
group_w((5*lines) : (6*lines)-1) =Q31
group_w((6*lines) : (7*lines)-1) =R41
group_w((7*lines) : (8*lines)-1) =Q41
group_w((8*lines) : (9*lines)-1) =P11
group_w((9*lines) :(10*lines)-1) =R12
group_w((10*lines):(11*lines)-1) =P21
group_w((11*lines):(12*lines)-1) =R22
group_w((12*lines):(13*lines)-1) =P31
group_w((13*lines):(14*lines)-1) =R32
group_w((14*lines):(15*lines)-1) =P41
group_w((15*lines):(16*lines)-1) =R42
group_w((16*lines):(17*lines)-1) =Q12
group_w((17*lines):(18*lines)-1) =P12
group_w((18*lines):(19*lines)-1) =Q22
group_w((19*lines):(20*lines)-1) =P22
group_w((20*lines):(21*lines)-1) =Q32
group_w((21*lines):(22*lines)-1) =P32
group_w((22*lines):(23*lines)-1) =Q42
group_w((23*lines):(24*lines)-1) =P42
group_w((24*lines):(25*lines)-1) =R13
group_w((25*lines):(26*lines)-1) =R23
group_w((26*lines):(27*lines)-1) =R33
group_w((27*lines):(28*lines)-1) =R43
group_w((28*lines):(29*lines)-1) =Q13
group_w((29*lines):(30*lines)-1) =Q23
group_w((30*lines):(31*lines)-1) =Q33
group_w((31*lines):(32*lines)-1) =Q43
group_w((32*lines):(33*lines)-1) =P13
group_w((33*lines):(34*lines)-1) =P23
group_w((34*lines):(35*lines)-1) =P33
group_w((35*lines):(36*lines)-1) =P43
group_w((36*lines):(37*lines)-1) =R14
group_w((37*lines):(38*lines)-1) =R24
group_w((38*lines):(39*lines)-1) =R34
group_w((39*lines):(40*lines)-1) =R44
group_w((40*lines):(41*lines)-1) =Q14
group_w((41*lines):(42*lines)-1) =Q24
group_w((42*lines):(43*lines)-1) =Q34
group_w((43*lines):(44*lines)-1) =Q44
group_w((44*lines):(45*lines)-1) =P14
group_w((45*lines):(46*lines)-1) =P24
group_w((46*lines):(47*lines)-1) =P34
group_w((47*lines):(48*lines)-1) =P44

group_i=dblarr(brcs*lines)
group_i(        0 : (1*lines)-1) =R11_i
group_i((1*lines) : (2*lines)-1) =Q11_i
group_i((2*lines) : (3*lines)-1) =R21_i
group_i((3*lines) : (4*lines)-1) =Q21_i
group_i((4*lines) : (5*lines)-1) =R31_i
group_i((5*lines) : (6*lines)-1) =Q31_i
group_i((6*lines) : (7*lines)-1) =R41_i
group_i((7*lines) : (8*lines)-1) =Q41_i
group_i((8*lines) : (9*lines)-1) =P11_i
group_i((9*lines) :(10*lines)-1) =R12_i
group_i((10*lines):(11*lines)-1) =P21_i
group_i((11*lines):(12*lines)-1) =R22_i
group_i((12*lines):(13*lines)-1) =P31_i
group_i((13*lines):(14*lines)-1) =R32_i
group_i((14*lines):(15*lines)-1) =P41_i
group_i((15*lines):(16*lines)-1) =R42_i
group_i((16*lines):(17*lines)-1) =Q12_i
group_i((17*lines):(18*lines)-1) =P12_i
group_i((18*lines):(19*lines)-1) =Q22_i
group_i((19*lines):(20*lines)-1) =P22_i
group_i((20*lines):(21*lines)-1) =Q32_i
group_i((21*lines):(22*lines)-1) =P32_i
group_i((22*lines):(23*lines)-1) =Q42_i
group_i((23*lines):(24*lines)-1) =P42_i
group_i((24*lines):(25*lines)-1) =R13_i
group_i((25*lines):(26*lines)-1) =R23_i
group_i((26*lines):(27*lines)-1) =R33_i
group_i((27*lines):(28*lines)-1) =R43_i
group_i((28*lines):(29*lines)-1) =Q13_i
group_i((29*lines):(30*lines)-1) =Q23_i
group_i((30*lines):(31*lines)-1) =Q33_i
group_i((31*lines):(32*lines)-1) =Q43_i
group_i((32*lines):(33*lines)-1) =P13_i
group_i((33*lines):(34*lines)-1) =P23_i
group_i((34*lines):(35*lines)-1) =P33_i
group_i((35*lines):(36*lines)-1) =P43_i
group_i((36*lines):(37*lines)-1) =R14_i
group_i((37*lines):(38*lines)-1) =R24_i
group_i((38*lines):(39*lines)-1) =R34_i
group_i((39*lines):(40*lines)-1) =R44_i
group_i((40*lines):(41*lines)-1) =Q14_i
group_i((41*lines):(42*lines)-1) =Q24_i
group_i((42*lines):(43*lines)-1) =Q34_i
group_i((43*lines):(44*lines)-1) =Q44_i
group_i((44*lines):(45*lines)-1) =P14_i
group_i((45*lines):(46*lines)-1) =P24_i
group_i((46*lines):(47*lines)-1) =P34_i
group_i((47*lines):(48*lines)-1) =P44_i
     sr=sort(group_w)
group_w=group_w(sr)
group_i=group_i(sr)

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;pro plot_HL
;window,11
; plot, indgen(10),HL_P(indgen(10)),yrange=[0,20],xrange=[0,5]
;oplot, indgen(10),HL_R(indgen(10)),linestyle=1
;oplot, indgen(10),HL_Q(indgen(10))

;print, HL_P(indgen(10))
;print, HL_R(indgen(10))
;print, HL_Q(indgen(10))

;end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro O2_1N,v1,v2,lines,temp,newc=newc,nocorr=nocorr, $
    R11,Q11,R21,Q21,R31,Q31,R41,Q41,P11,R12,P21,R22,P31,R32,P41,R42,Q12,$
    P12,Q22,P22,Q32,P32,Q42,P42,R13,R23,R33,R43,Q13,Q23,Q33,Q43,P13,P23,$
    P33,P43,R14,R24,R34,R44,Q14,Q24,Q34,Q44,P14,P24,P34,P44,R11_i,Q11_i,$
    R21_i,Q21_i,R31_i,Q31_i,R41_i,Q41_i,P11_i,R12_i,P21_i,R22_i,P31_i,$
    R32_i,P41_i,R42_i,Q12_i,P12_i,Q22_i,P22_i,Q32_i,P32_i,Q42_i,P42_i,$
    R13_i,R23_i,R33_i,R43_i,Q13_i,Q23_i,Q33_i,Q43_i,P13_i,P23_i,P33_i,$
    P43_i,R14_i,R24_i,R34_i,R44_i,Q14_i,Q24_i,Q34_i,Q44_i,P14_i,P24_i,$
    P34_i,P44_i


if keyword_set(newc) then begin
;#################################NEW CONSTANTS######################################
;using new constants - put in new constants here david
; Hui Liu 2014

    	  ;L&G      	    		    ;Krupener 1987	    ;L&G/Hansen (1983)	    ;Somewhere?   	    		    	
b_T_e    =16587.9478d  	    	    	    ;146759.0d	    	    ;16587.9478d 	    ;145961.7816d    	    	
b_w_e    =1197.0173d  	    	    	    ;1196.913d	    	    ;1197.0173d 	    ;1197.02d      	    	
b_wexe   =17.1720d  	    	    	    ;17.13546d	    	    ;17.1720d 	    	    ;17.172d       	    	
b_weye   =1.177d-2  	    	    	    ;0.0d	    	    ;1.177d-2 	    	    ;1.18d-2  	    	    
b_weze   =-9.92d-4  	    	    	    ;0.0d	    	    ;-9.92d-4  	    	    ;-1.0d-3  	    	    	

b_B_e    =1.2876556d  	    	    	    ;1.287297d	    	    ;1.2876556d 	    ;1.28766d       	    	
b_alpha_e=2.19209d-2  	    	    	    ;0.02206747d	    ;2.19209d-2 	    ;2.192d-2     	    	  
b_gamma_e=-1.2762d-4  	    	    	    ;0.0d	    	    ;-1.2762d-4  	    ;-1.28d-4 	  	    	

b_D_e    =5.9992d-6  	    	    	    ;5.9992d-6	    	    ;5.9992d-6 	    	    ;5.81d-6      	    	
b_SS_e   =-1.127d-7  	    	    	    ;-1.127d-7	     	    ;-1.127d-7  	    ;0.185d-6     	    	


a_T_e    =0.434d  	    	    	    ;130161d	    	    ;0.434d 	    	    ;129374.2679d   	    	
a_w_e    =1035.13d  	    	    	    ;1035.534d	    	    ;1035.13d 	    	    ;1035.13d    	    	
a_wexe   =10.115d  	    	    	    ;10.32194d	    	    ;10.115d 	    	    ;10.115d       	    
a_weye   =-3.31d-2  	    	    	    ;0.0d	    	    ;-3.31d-2  	    	    ;-3.31d-2 	    	    	
a_weze   =2.1d-4  	    	    	    ;0.0d	    	    ;2.1d-4 	    	    ;2.1d-4

a_B_e    =1.1047580d  	    	    	    ;1.104320d	    	    ;1.1047580d 	    ;1.10476d             	
a_alpha_e=1.54762d-2  	    	    	    ;0.0154546d	    	    ;1.54762d-2 	    ;1.548d-2     	    	
a_gamma_e=1.2d-5  	    	    	    ;0.0d	    	    ;1.2d-5 	    	    ;1.2d-5   	    	    	
a_delta_e=-5.006d-6  	    	    	    ;0.0d	    	    ;-5.006d-6  	    t;-5.0d-6
	    	    
a_D_e    =5.0392d-6  	    	    	    ;5.0392d-6 	    	    ;5.0392d-6 	    	    ;4.88d-6        	    	
a_SS_e   =-1.20d-8  	    	    	    ;-1.20d-8 		    ;-1.20d-8 	    	    ;-0.095d-6   ;O2  	    	

;^source of originals is h&h footnotes 

Bu=1.4378  ;O2 ground state rotational constant  ;Lofthus 1972

;Electronic & vibronic energy------------------------------------------------
deltaE=  b_T_e-a_T_e ;  - 8.0
deltaG=( b_w_e*(v1+0.5) - b_wexe*(v1+0.5)^2 + b_weye*(v1+0.5)^3 + b_weze*(v1+0.5)^4 )  - $
       ( a_w_e*(v2+0.5) - a_wexe*(v2+0.5)^2 + a_weye*(v2+0.5)^3 + a_weze*(v2+0.5)^4 )



;print, 'newc'
;print, 'deltaE :', deltaE, ' deltaG :', deltaG
;Upper 4Sigma State vibrational settings ---------------------------------------

b_B_v = b_B_e - b_alpha_e *(v1+0.5) + b_gamma_e*(v1+0.5)^2
b_D_v = b_D_e-b_SS_e    *(v1+0.5)


eps = 0.1487
ga =  -0.00033

;Lower 4Pi State vibrational settings-------------------------------------

a_B_v = a_B_e - a_alpha_e *(v2+0.5) + a_gamma_e*(v2+0.5)^2 + a_delta_e*(v2+0.5)^3

a_D_v = a_D_e-a_SS_e*(v2+0.5) 
   

;HANSEN 1983 - current best
A_v = (-47.78421d) + (-2.348d-2*(v2+0.5))  + (9.095d-3*(v2+0.5)^2)

;Albritton 1977
;A_v = (-47.78561d) + (-1.868d-2*(v2+0.5))  + (-7.967d-3*(v2+0.5)^2)


Y = A_v/B_B_v

 
;print, 'newc'
;print, 'deltaE :', deltaE, ' deltaG :', deltaG


endif else begin
;#################################OLD CONSTANTS######################################

b_T_e    = 49552.0    ;O2 from Herzberg 1979
b_w_e    = 1196.77    ;O2
b_wexe   = 17.09      ;O2
b_B_e    = 1.28729    ;O2 from Nevin 1940
b_alpha_e= 0.02206    ;O2
b_D_e    = 5.81d-6    ;O2
b_SS_e   = 0.185d-6   ;O2

a_T_e    = 32964.0   ;O2 from Herzberg 1979
a_w_e    = 1035.69   ;O2
a_wexe   = 10.39     ;O2
a_B_e    = 1.10466   ;O2 from Nevin 1940    
a_alpha_e= 0.01575   ;O2
a_D_e    = 4.88d-6   ;O2
a_SS_e   = -0.095d-6 ;O2

Bu=1.4378  ;O2 ground state rotational constant  ;Lofthus 1972

;Electronic & vibronic energy------------------------------------------------
 deltaE=  b_T_e-a_T_e
 deltaG=  b_w_e*(v1+0.5)-b_wexe*(v1+0.5)^2.-a_w_e*(v2+0.5)+a_wexe*(v2+0.5)^2.
 
;print, 'oldc'
;print, 'deltaE :', deltaE, ' deltaG :', deltaG

;Upper 4Sigma State vibrational settings ---------------------------------------
 b_B_v = b_B_e-b_alpha_e *(v1+0.5)
 b_D_v = b_D_e-b_SS_e    *(v1+0.5)
   eps = 0.1487
    ga = -0.00033
;Lower 4Pi State vibrational settings-------------------------------------
 a_B_v = a_B_e-a_alpha_e *(v2+0.5)
 a_D_v = a_D_e-a_SS_e    *(v2+0.5)
   A_v = -47.79d; -  0.068d*(v2+0.5)     $ 47.79
                ;- (2.5d-3)*(v2+0.5)^2. $
                ;- (2.6d-4)*(v2+0.5)^3.
     Y = A_v/B_B_v

endelse


if lines gt 100. then begin
print, 'too many lines to reference'
stop
endif

;LOAD correcting factors
if not(keyword_set(nocorr)) then restore,'O2_E_correction.idl'
;restore, 'O2_E_correction_test.idl'

if keyword_set(nocorr) then begin
    COR_11B = fltarr(100)
    COR_22B = fltarr(100)
    COR_33B = fltarr(100)
    COR_44B = fltarr(100)
endif

;print, min(COR_11B), max(COR_11B), mean(COR_11B), '--'

;states,b_B_v,b_D_v,eps,ga,a_B_v,a_D_v,Y,lines,cor_11b,cor_22b,cor_33b,cor_44b

;Load HL factors
O2_HL, lines, $
R11_J, R11_new,Q11_J, Q11_new,R21_J, R21_new,Q21_J, Q21_new,$
R31_J, R31_new,Q31_J, Q31_new,R41_J, R41_new,Q41_J, Q41_new,$
P11_J, P11_new,R12_J, R12_new,P21_J, P21_new,R22_J, R22_new,$
P31_J, P31_new,R32_J, R32_new,P41_J, P41_new,R42_J, R42_new,$
Q12_J, Q12_new,P12_J, P12_new,Q22_J, Q22_new,P22_J, P22_new,$
Q32_J, Q32_new,P32_J, P32_new,Q42_J, Q42_new,P42_J, P42_new,$
R13_J, R13_new,R23_J, R23_new,R33_J, R33_new,R43_J, R43_new,$
Q13_J, Q13_new,Q23_J, Q23_new,Q33_J, Q33_new,Q43_J, Q43_new,$
P13_J, P13_new,P23_J, P23_new,P33_J, P33_new,P43_J, P43_new,$
R14_J, R14_new,R24_J, R24_new,R34_J, R34_new,R44_J, R44_new,$
Q14_J, Q14_new,Q24_J, Q24_new,Q34_J, Q34_new,Q44_J, Q44_new,$
P14_J, P14_new,P24_J, P24_new,P34_J, P34_new,P44_J, P44_new 

;F1-F1-------------------------------------------------------------
;P-branch F1-F1 levels (deltaJ = +1)
P11  =fltarr(lines)                                                   ;#
P11_i=fltarr(lines)                                                   ;#;J is J''
for J=P11_J,lines-1. do begin                                         ;#
P11  (J)=Sigma4F1(b_B_v,b_D_v,eps,ga,J-.5)-Pi4F1(a_B_v,a_D_v,Y,J+.5,cor_11b(J)) ;#
P11  (J)=nuva2leai(deltaE+deltaG+P11(J))                              ;#
P11_i(J)=P11_new(J)*njp(J+.5,temp)                                    ;#
endfor
;R-branch F1-F1 levels (deltaJ = -1)
R11  =fltarr(lines)                                                   ;#
R11_i=fltarr(lines)                                                   ;#;J is J''
for J=R11_J,lines-1. do begin                                         ;#
R11  (J)=Sigma4F1(b_B_v,b_D_v,eps,ga,J+1.5)-Pi4F1(a_B_v,a_D_v,Y,J+.5,cor_11b(J)) ;#
R11  (J)=nuva2leai(deltaE+deltaG+R11(J))                              ;#
R11_i(J)=R11_new(J)*njp(J+.5,temp)                                    ;#
endfor
;Q-branch F1-F1 levels (deltaJ =  0)
Q11  =fltarr(lines)                                                   ;#
Q11_i=fltarr(lines)                                                   ;#;J is J''
for J=Q11_J,lines-1. do begin                                         ;#
Q11  (J)=Sigma4F1(b_B_v,b_D_v,eps,ga,J+0.5)-Pi4F1(a_B_v,a_D_v,Y,J+.5,cor_11b(J)) ;#
Q11  (J)=nuva2leai(deltaE+deltaG+Q11(J))                              ;#
Q11_i(J)=Q11_new(J)*njp(J+.5,temp)                                    ;#
endfor

;F2-F2-------------------------------------------------------------
;P-branch F2-F2 levels (deltaJ = +1)
P22  =fltarr(lines)                                                   ;#
P22_i=fltarr(lines)                                                   ;#;J is J''
for J=P22_J,lines-1. do begin                                         ;#
P22  (J)=Sigma4F2(b_B_v,b_D_v,eps,ga,J-.5)-Pi4F2(a_B_v,a_D_v,Y,J+.5,cor_22b(J)) ;#
P22  (J)=nuva2leai(deltaE+deltaG+P22(J))                              ;#
P22_i(J)=P22_new(J)*njp(J+.5,temp)                                    ;#
endfor
;R-branch F2-F2 levels (deltaJ = -1)
R22  =fltarr(lines)                                                   ;#
R22_i=fltarr(lines)                                                   ;#;J is J''
for J=R22_J,lines-1. do begin                                         ;#
R22  (J)=Sigma4F2(b_B_v,b_D_v,eps,ga,J+1.5)-Pi4F2(a_B_v,a_D_v,Y,J+.5,cor_22b(J)) ;#
R22  (J)=nuva2leai(deltaE+deltaG+R22(J))                              ;#
R22_i(J)=R22_new(J)*njp(J+.5,temp)                                    ;#
endfor
;Q-branch F2-F2 levels (deltaJ =  0)
Q22  =fltarr(lines)                                                   ;#
Q22_i=fltarr(lines)                                                   ;#;J is J''
for J=Q22_J,lines-1. do begin                                         ;#
Q22  (J)=Sigma4F2(b_B_v,b_D_v,eps,ga,J+0.5)-Pi4F2(a_B_v,a_D_v,Y,J+.5,cor_22b(J)) ;#
Q22  (J)=nuva2leai(deltaE+deltaG+Q22(J))                              ;#
Q22_i(J)=Q22_new(J)*njp(J+.5,temp)                                    ;#
endfor

;F3-F3-------------------------------------------------------------
;P-branch F3-F3 levels (deltaJ = +1)
P33  =fltarr(lines)                                                   ;#
P33_i=fltarr(lines)                                                   ;#;J is J''
for J=P33_J,lines-1. do begin                                         ;#
P33  (J)=Sigma4F3(b_B_v,b_D_v,eps,ga,J-.5)-Pi4F3(a_B_v,a_D_v,Y,J+.5,cor_33b(J)) ;#
P33  (J)=nuva2leai(deltaE+deltaG+P33(J))                              ;#
P33_i(J)=P33_new(J)*njp(J+.5,temp)                                    ;#
endfor
;R-branch F3-F3 levels (deltaJ = -1)
R33  =fltarr(lines)                                                   ;#
R33_i=fltarr(lines)                                                   ;#;J is J''
for J=R33_J,lines-1. do begin                                         ;#
R33  (J)=Sigma4F3(b_B_v,b_D_v,eps,ga,J+1.5)-Pi4F3(a_B_v,a_D_v,Y,J+.5,cor_33b(J)) ;#
R33  (J)=nuva2leai(deltaE+deltaG+R33(J))                              ;#
R33_i(J)=R33_new(J)*njp(J+.5,temp)                                    ;#
endfor
;Q-branch F3-F3 levels (deltaJ =  0)
Q33  =fltarr(lines)                                                   ;#
Q33_i=fltarr(lines)                                                   ;#;J is J''
for J=Q33_J,lines-1. do begin                                         ;#
Q33  (J)=Sigma4F3(b_B_v,b_D_v,eps,ga,J+0.5)-Pi4F3(a_B_v,a_D_v,Y,J+.5,cor_33b(J)) ;#
Q33  (J)=nuva2leai(deltaE+deltaG+Q33(J))                              ;#
Q33_i(J)=Q33_new(J)*njp(J+.5,temp)                                    ;#
endfor

;F4-F4-------------------------------------------------------------
;P-branch F4-F4 levels (deltaJ = +1)
P44  =fltarr(lines)                                                   ;#
P44_i=fltarr(lines)                                                   ;#;J is J''
for J=P44_J,lines-1. do begin                                         ;#
P44  (J)=Sigma4F4(b_B_v,b_D_v,eps,ga,J-.5)-Pi4F4(a_B_v,a_D_v,Y,J+.5,cor_44b(J)) ;#
P44  (J)=nuva2leai(deltaE+deltaG+P44(J))                              ;#
P44_i(J)=P44_new(J)*njp(J+.5,temp)                                    ;#
endfor
;R-branch F4-F4 levels (deltaJ = -1)
R44  =fltarr(lines)                                                   ;#
R44_i=fltarr(lines)                                                   ;#;J is J''
for J=R44_J,lines-1. do begin                                         ;#
R44  (J)=Sigma4F4(b_B_v,b_D_v,eps,ga,J+1.5)-Pi4F4(a_B_v,a_D_v,Y,J+.5,cor_44b(J)) ;#
R44  (J)=nuva2leai(deltaE+deltaG+R44(J))                              ;#
R44_i(J)=R44_new(J)*njp(J+.5,temp)                                    ;#
endfor
;Q-branch F4-F4 levels (deltaJ =  0)
Q44  =fltarr(lines)                                                   ;#
Q44_i=fltarr(lines)                                                   ;#;J is J''
for J=Q44_J,lines-1. do begin                                         ;#
Q44  (J)=Sigma4F4(b_B_v,b_D_v,eps,ga,J+0.5)-Pi4F4(a_B_v,a_D_v,Y,J+.5,cor_44b(J)) ;#
Q44  (J)=nuva2leai(deltaE+deltaG+Q44(J))                              ;#
Q44_i(J)=Q44_new(J)*njp(J+.5,temp)                                    ;#
endfor

;F2-F1-------------------------------------------------------------
;P-branch F2-F1 levels (deltaJ = +1)
P21  =fltarr(lines)                                                   ;#
P21_i=fltarr(lines)                                                   ;#;J is J''
for J=P21_J,lines-1. do begin                                         ;#
P21  (J)=Sigma4F2(b_B_v,b_D_v,eps,ga,J-.5)-Pi4F1(a_B_v,a_D_v,Y,J+.5,cor_11b(J))   +0.0026*J*J - 0.0348*J + 0.0344;+0.06*J ;#
P21  (J)=nuva2leai(deltaE+deltaG+P21(J))                              ;#
P21_i(J)=P21_new(J)*njp(J+.5,temp)                                    ;#
endfor
;R-branch F2-F1 levels (deltaJ = -1)
R21  =fltarr(lines)                                                   ;#
R21_i=fltarr(lines)                                                   ;#;J is J''
for J=R21_J,lines-1. do begin                                         ;#
R21  (J)=Sigma4F2(b_B_v,b_D_v,eps,ga,J+1.5)-Pi4F1(a_B_v,a_D_v,Y,J+.5,cor_11b(J))   +0.0026*J*J - 0.0348*J + 0.0344;+0.06*J ;#
R21  (J)=nuva2leai(deltaE+deltaG+R21(J))                              ;#
R21_i(J)=R21_new(J)*njp(J+.5,temp)                                    ;#
endfor
;Q-branch F2-F1 levels (deltaJ =  0)
Q21  =fltarr(lines)                                                   ;#
Q21_i=fltarr(lines)                                                   ;#;J is J''
for J=Q21_J,lines-1. do begin                                         ;#
Q21  (J)=Sigma4F2(b_B_v,b_D_v,eps,ga,J+0.5)-Pi4F1(a_B_v,a_D_v,Y,J+.5,cor_11b(J))   +0.0026*J*J - 0.0348*J + 0.0344;+0.06*J ;#
Q21  (J)=nuva2leai(deltaE+deltaG+Q21(J))                              ;#
Q21_i(J)=Q21_new(J)*njp(J+.5,temp)                                    ;#
endfor

;F3-F1-------------------------------------------------------------
;P-branch F3-F1 levels (deltaJ = +1)
P31  =fltarr(lines)                                                   ;#
P31_i=fltarr(lines)                                                   ;#;J is J''
for J=P31_J,lines-1. do begin                                         ;#
P31  (J)=Sigma4F3(b_B_v,b_D_v,eps,ga,J-.5)-Pi4F1(a_B_v,a_D_v,Y,J+.5,cor_11b(J)) +0.0001*J*J*J - 0.0009*J*J + 0.0094*J - 0.031;#
P31  (J)=nuva2leai(deltaE+deltaG+P31(J))                              ;#
P31_i(J)=P31_new(J)*njp(J+.5,temp)                                    ;#
endfor
;R-branch F3-F1 levels (deltaJ = -1)
R31  =fltarr(lines)                                                   ;#
R31_i=fltarr(lines)                                                   ;#;J is J''
for J=R31_J,lines-1. do begin                                         ;#
R31  (J)=Sigma4F3(b_B_v,b_D_v,eps,ga,J+1.5)-Pi4F1(a_B_v,a_D_v,Y,J+.5,cor_11b(J)) +0.0001*J*J*J - 0.0009*J*J + 0.0094*J - 0.031;#
R31  (J)=nuva2leai(deltaE+deltaG+R31(J))                              ;#
R31_i(J)=R31_new(J)*njp(J+.5,temp)                                    ;#
endfor
;Q-branch F3-F1 levels (deltaJ =  0)
Q31  =fltarr(lines)                                                   ;#
Q31_i=fltarr(lines)                                                   ;#;J is J''
for J=Q31_J,lines-1. do begin                                         ;#
Q31  (J)=Sigma4F3(b_B_v,b_D_v,eps,ga,J+0.5)-Pi4F1(a_B_v,a_D_v,Y,J+.5,cor_11b(J)) +0.0001*J*J*J - 0.0009*J*J + 0.0094*J - 0.031;#
Q31  (J)=nuva2leai(deltaE+deltaG+Q31(J))                              ;#
Q31_i(J)=Q31_new(J)*njp(J+.5,temp)                                    ;#
endfor

;F4-F1-------------------------------------------------------------
;P-branch F4-F1 levels (deltaJ = +1)
P41  =fltarr(lines)                                                   ;#
P41_i=fltarr(lines)                                                   ;#;J is J''
for J=P41_J,lines-1. do begin                                         ;#
P41  (J)=Sigma4F4(b_B_v,b_D_v,eps,ga,J-.5)-Pi4F1(a_B_v,a_D_v,Y,J+.5,cor_11b(J))   +0.0002*J*J*J - 0.001*J*J + 0.016*J - 0.0745;#
P41  (J)=nuva2leai(deltaE+deltaG+P41(J))                              ;#
P41_i(J)=P41_new(J)*njp(J+.5,temp)                                    ;#
endfor
;R-branch F4-F1 levels (deltaJ = -1)
R41  =fltarr(lines)                                                   ;#
R41_i=fltarr(lines)                                                   ;#;J is J''
for J=R41_J,lines-1. do begin                                         ;#
R41  (J)=Sigma4F4(b_B_v,b_D_v,eps,ga,J+1.5)-Pi4F1(a_B_v,a_D_v,Y,J+.5,cor_11b(J))   +0.0002*J*J*J - 0.001*J*J + 0.016*J - 0.0745;#
R41  (J)=nuva2leai(deltaE+deltaG+R41(J))                              ;#
R41_i(J)=R41_new(J)*njp(J+.5,temp)                                    ;#
endfor
;Q-branch F4-F1 levels (deltaJ =  0)
Q41  =fltarr(lines)                                                   ;#
Q41_i=fltarr(lines)                                                   ;#;J is J''
for J=Q41_J,lines-1. do begin                                         ;#
Q41  (J)=Sigma4F4(b_B_v,b_D_v,eps,ga,J+0.5)-Pi4F1(a_B_v,a_D_v,Y,J+.5,cor_11b(J))   +0.0002*J*J*J - 0.001*J*J + 0.016*J - 0.0745;#
Q41  (J)=nuva2leai(deltaE+deltaG+Q41(J))                              ;#
Q41_i(J)=Q41_new(J)*njp(J+.5,temp)                                    ;#
endfor

;F1-F2-------------------------------------------------------------
;P-branch F1-F2 levels (deltaJ = +1)
P12  =fltarr(lines)                                                   ;#
P12_i=fltarr(lines)                                                   ;#;J is J''
for J=P12_J,lines-1. do begin                                         ;#
P12  (J)=Sigma4F1(b_B_v,b_D_v,eps,ga,J-.5)-Pi4F2(a_B_v,a_D_v,Y,J+.5,cor_22b(J))    -0.0025*J*J + 0.0348*J - 0.0966;-0.0066*J*J + 0.0052*J + 0.1093 ;#
P12  (J)=nuva2leai(deltaE+deltaG+P12(J))                              ;#
P12_i(J)=P12_new(J)*njp(J+.5,temp)                                    ;#
endfor
;R-branch F1-F2 levels (deltaJ = -1)
R12  =fltarr(lines)                                                   ;#
R12_i=fltarr(lines)                                                   ;#;J is J''
for J=R22_J,lines-1. do begin                                         ;#
R12  (J)=Sigma4F1(b_B_v,b_D_v,eps,ga,J+1.5)-Pi4F2(a_B_v,a_D_v,Y,J+.5,cor_22b(J))   -0.0025*J*J + 0.0348*J - 0.0966;-0.011*J*J + 0.049*J - 0.0795 ;#
R12  (J)=nuva2leai(deltaE+deltaG+R12(J))                              ;#
R12_i(J)=R12_new(J)*njp(J+.5,temp)                                    ;#
endfor
;Q-branch F1-F2 levels (deltaJ =  0)
Q12  =fltarr(lines)                                                   ;#
Q12_i=fltarr(lines)                                                   ;#;J is J''
for J=Q12_J,lines-1. do begin                                         ;#
Q12  (J)=Sigma4F1(b_B_v,b_D_v,eps,ga,J+0.5)-Pi4F2(a_B_v,a_D_v,Y,J+.5,cor_22b(J))   -0.0025*J*J + 0.0348*J - 0.0966;-0.0122*J*J + 0.0983*J - 0.1879 ;#
Q12  (J)=nuva2leai(deltaE+deltaG+Q12(J))                              ;#
Q12_i(J)=Q12_new(J)*njp(J+.5,temp)                                    ;#
endfor

;F3-F2-------------------------------------------------------------
;P-branch F2-F2 levels (deltaJ = +1)
P32  =fltarr(lines)                                                   ;#
P32_i=fltarr(lines)                                                   ;#;J is J''
for J=P32_J,lines-1. do begin                                         ;#
P32  (J)=Sigma4F3(b_B_v,b_D_v,eps,ga,J-.5)-Pi4F2(a_B_v,a_D_v,Y,J+.5,cor_22b(J))   +0.00007*J*J*J - 0.0014*J*J + 0.0194*J - 0.0451;#
P32  (J)=nuva2leai(deltaE+deltaG+P32(J))                              ;#
P32_i(J)=P32_new(J)*njp(J+.5,temp)                                    ;#
endfor
;R-branch F3-F2 levels (deltaJ = -1)
R32  =fltarr(lines)                                                   ;#
R32_i=fltarr(lines)                                                   ;#;J is J''
for J=R32_J,lines-1. do begin                                         ;#
R32  (J)=Sigma4F3(b_B_v,b_D_v,eps,ga,J+1.5)-Pi4F2(a_B_v,a_D_v,Y,J+.5,cor_22b(J))   +0.00007*J*J*J - 0.0014*J*J + 0.0194*J - 0.0451 ;#
R32  (J)=nuva2leai(deltaE+deltaG+R32(J))                              ;#
R32_i(J)=R32_new(J)*njp(J+.5,temp)                                    ;#
endfor
;Q-branch F3-F2 levels (deltaJ =  0)
Q32  =fltarr(lines)                                                   ;#
Q32_i=fltarr(lines)                                                   ;#;J is J''
for J=Q32_J,lines-1. do begin                                         ;#
Q32  (J)=Sigma4F3(b_B_v,b_D_v,eps,ga,J+0.5)-Pi4F2(a_B_v,a_D_v,Y,J+.5,cor_22b(J))   +0.00007*J*J*J - 0.0014*J*J + 0.0194*J - 0.0451 ;#
Q32  (J)=nuva2leai(deltaE+deltaG+Q32(J))                              ;#
Q32_i(J)=Q32_new(J)*njp(J+.5,temp)                                    ;#
endfor

;F4-F2-------------------------------------------------------------
;P-branch F4-F2 levels (deltaJ = +1)
P42  =fltarr(lines)                                                   ;#
P42_i=fltarr(lines)                                                   ;#;J is J''
for J=P42_J,lines-1. do begin                                         ;#
P42  (J)=Sigma4F4(b_B_v,b_D_v,eps,ga,J-.5)-Pi4F2(a_B_v,a_D_v,Y,J+.5,cor_22b(J))    +0.0001*J*J*J - 0.0018*J*J + 0.0265*J - 0.0581;#
P42  (J)=nuva2leai(deltaE+deltaG+P42(J))                              ;#
P42_i(J)=P42_new(J)*njp(J+.5,temp)                                    ;#
endfor
;R-branch F4-F2 levels (deltaJ = -1)
R42  =fltarr(lines)                                                   ;#
R42_i=fltarr(lines)                                                   ;#;J is J''
for J=R42_J,lines-1. do begin                                         ;#
R42  (J)=Sigma4F4(b_B_v,b_D_v,eps,ga,J+1.5)-Pi4F2(a_B_v,a_D_v,Y,J+.5,cor_22b(J))    +0.0001*J*J*J - 0.0018*J*J + 0.0265*J - 0.0581 ;#
R42  (J)=nuva2leai(deltaE+deltaG+R42(J))                              ;#
R42_i(J)=R42_new(J)*njp(J+.5,temp)                                    ;#
endfor
;Q-branch F4-F2 levels (deltaJ =  0)
Q42  =fltarr(lines)                                                   ;#
Q42_i=fltarr(lines)                                                   ;#;J is J''
for J=Q42_J,lines-1. do begin                                         ;#
Q42  (J)=Sigma4F4(b_B_v,b_D_v,eps,ga,J+0.5)-Pi4F2(a_B_v,a_D_v,Y,J+.5,cor_22b(J))    +0.0001*J*J*J - 0.0018*J*J + 0.0265*J - 0.0581 ;#
Q42  (J)=nuva2leai(deltaE+deltaG+Q42(J))                              ;#
Q42_i(J)=Q42_new(J)*njp(J+.5,temp)                                    ;#
endfor

;F1-F3-------------------------------------------------------------
;P-branch F1-F3 levels (deltaJ = +1)
P13  =fltarr(lines)                                                   ;#
P13_i=fltarr(lines)                                                   ;#;J is J''
for J=P13_J,lines-1. do begin                                         ;#
P13  (J)=Sigma4F1(b_B_v,b_D_v,eps,ga,J-.5)-Pi4F3(a_B_v,a_D_v,Y,J+.5,cor_33b(J))    -0.0054*J*J + 0.0888*J - 0.3239  ;-0.0222*J*J + 0.1665*J - 0.2817 ;#
P13  (J)=nuva2leai(deltaE+deltaG+P13(J))                              ;#
P13_i(J)=P13_new(J)*njp(J+.5,temp)                                    ;#
endfor
;R-branch F1-F3 levels (deltaJ = -1)
R13  =fltarr(lines)                                                   ;#
R13_i=fltarr(lines)                                                   ;#;J is J''
for J=R13_J,lines-1. do begin                                         ;#
R13  (J)=Sigma4F1(b_B_v,b_D_v,eps,ga,J+1.5)-Pi4F3(a_B_v,a_D_v,Y,J+.5,cor_33b(J))   -0.0054*J*J + 0.0888*J - 0.3239  ;-0.0225*J*J + 0.1722*J - 0.2563 ;#
R13  (J)=nuva2leai(deltaE+deltaG+R13(J))                              ;#
R13_i(J)=R13_new(J)*njp(J+.5,temp)                                    ;#
endfor
;Q-branch F1-F3 levels (deltaJ =  0)
Q13  =fltarr(lines)                                                   ;#
Q13_i=fltarr(lines)                                                   ;#;J is J''
for J=Q13_J,lines-1. do begin                                         ;#
Q13  (J)=Sigma4F1(b_B_v,b_D_v,eps,ga,J+0.5)-Pi4F3(a_B_v,a_D_v,Y,J+.5,cor_33b(J))   -0.0054*J*J + 0.0888*J - 0.3239  ;-0.0234*J*J + 0.1955*J - 0.4428 ;#
Q13  (J)=nuva2leai(deltaE+deltaG+Q13(J))                              ;#
Q13_i(J)=Q13_new(J)*njp(J+.5,temp)                                    ;#
endfor

;F2-F3-------------------------------------------------------------
;P-branch F2-F3 levels (deltaJ = +1)
P23  =fltarr(lines)                                                   ;#
P23_i=fltarr(lines)                                                   ;#;J is J''
for J=P23_J,lines-1. do begin                                         ;#
P23  (J)=Sigma4F2(b_B_v,b_D_v,eps,ga,J-.5)-Pi4F3(a_B_v,a_D_v,Y,J+.5,cor_33b(J)) -0.00004*J*J*J - 0.0002*J*J + 0.0046*J - 0.0186;#
P23  (J)=nuva2leai(deltaE+deltaG+P23(J))                              ;#
P23_i(J)=P22_new(J)*njp(J+.5,temp)                                    ;#
endfor
;R-branch F2-F3 levels (deltaJ = -1)
R23  =fltarr(lines)                                                   ;#
R23_i=fltarr(lines)                                                   ;#;J is J''
for J=R23_J,lines-1. do begin                                         ;#
R23  (J)=Sigma4F2(b_B_v,b_D_v,eps,ga,J+1.5)-Pi4F3(a_B_v,a_D_v,Y,J+.5,cor_33b(J)) -0.00004*J*J*J - 0.0002*J*J + 0.0046*J - 0.0186;#
R23  (J)=nuva2leai(deltaE+deltaG+R23(J))                              ;#
R23_i(J)=R23_new(J)*njp(J+.5,temp)                                    ;#
endfor
;Q-branch F2-F3 levels (deltaJ =  0)
Q23  =fltarr(lines)                                                   ;#
Q23_i=fltarr(lines)                                                   ;#;J is J''
for J=Q23_J,lines-1. do begin                                         ;#
Q23  (J)=Sigma4F2(b_B_v,b_D_v,eps,ga,J+0.5)-Pi4F3(a_B_v,a_D_v,Y,J+.5,cor_33b(J)) -0.00004*J*J*J - 0.0002*J*J + 0.0046*J - 0.0186;#
Q23  (J)=nuva2leai(deltaE+deltaG+Q23(J))                              ;#
Q23_i(J)=Q23_new(J)*njp(J+.5,temp)                                    ;#
endfor

;F4-F3-------------------------------------------------------------
;P-branch F4-F3 levels (deltaJ = +1)
P43  =fltarr(lines)                                                   ;#
P43_i=fltarr(lines)                                                   ;#;J is J''
for J=P43_J,lines-1. do begin                                         ;#
P43  (J)=Sigma4F4(b_B_v,b_D_v,eps,ga,J-.5)-Pi4F3(a_B_v,a_D_v,Y,J+.5,cor_33b(J))   +0.00005*J*J*J - 0.0002*J*J + 0.0018*J + 0.009;#
P43  (J)=nuva2leai(deltaE+deltaG+P43(J))                              ;#
P43_i(J)=P43_new(J)*njp(J+.5,temp)                                    ;#
endfor
;R-branch F4-F3 levels (deltaJ = -1)
R43  =fltarr(lines)                                                   ;#
R43_i=fltarr(lines)                                                   ;#;J is J''
for J=R43_J,lines-1. do begin                                         ;#
R43  (J)=Sigma4F4(b_B_v,b_D_v,eps,ga,J+1.5)-Pi4F3(a_B_v,a_D_v,Y,J+.5,cor_33b(J))   +0.00005*J*J*J - 0.0002*J*J + 0.0018*J + 0.009 ;#
R43  (J)=nuva2leai(deltaE+deltaG+R43(J))                              ;#
R43_i(J)=R43_new(J)*njp(J+.5,temp)                                    ;#
endfor
;Q-branch F4-F3 levels (deltaJ =  0)
Q43  =fltarr(lines)                                                   ;#
Q43_i=fltarr(lines)                                                   ;#;J is J''
for J=Q43_J,lines-1. do begin                                         ;#
Q43  (J)=Sigma4F4(b_B_v,b_D_v,eps,ga,J+0.5)-Pi4F3(a_B_v,a_D_v,Y,J+.5,cor_33b(J))   +0.00005*J*J*J - 0.0002*J*J + 0.0018*J + 0.009 ;#
Q43  (J)=nuva2leai(deltaE+deltaG+Q43(J))                              ;#
Q43_i(J)=Q43_new(J)*njp(J+.5,temp)                                    ;#
endfor

;F1-F4-------------------------------------------------------------
;P-branch F1-F4 levels (deltaJ = +1)
P14  =fltarr(lines)                                                   ;#
P14_i=fltarr(lines)                                                   ;#;J is J''
for J=P14_J,lines-1. do begin                                         ;#
P14  (J)=Sigma4F1(b_B_v,b_D_v,eps,ga,J-.5)-Pi4F4(a_B_v,a_D_v,Y,J+.5,cor_44b(J)) -0.0002*J*J*J + 0.0009*J*J - 0.0058*J - 0.0445;#
P14  (J)=nuva2leai(deltaE+deltaG+P14(J))                              ;#
P14_i(J)=P14_new(J)*njp(J+.5,temp)                                    ;#
endfor
;R-branch F1-F4 levels (deltaJ = -1)
R14  =fltarr(lines)                                                   ;#
R14_i=fltarr(lines)                                                   ;#;J is J''
for J=R14_J,lines-1. do begin                                         ;#
R14  (J)=Sigma4F1(b_B_v,b_D_v,eps,ga,J+1.5)-Pi4F4(a_B_v,a_D_v,Y,J+.5,cor_44b(J)) -0.0002*J*J*J + 0.0009*J*J - 0.0058*J - 0.0445;#
R14  (J)=nuva2leai(deltaE+deltaG+R14(J))                              ;#
R14_i(J)=R14_new(J)*njp(J+.5,temp)                                    ;#
endfor
;Q-branch F1-F4 levels (deltaJ =  0)
Q14  =fltarr(lines)                                                   ;#
Q14_i=fltarr(lines)                                                   ;#;J is J''
for J=Q14_J,lines-1. do begin                                         ;#
Q14  (J)=Sigma4F1(b_B_v,b_D_v,eps,ga,J+0.5)-Pi4F4(a_B_v,a_D_v,Y,J+.5,cor_44b(J)) -0.0002*J*J*J + 0.0009*J*J - 0.0058*J - 0.0445;#
Q14  (J)=nuva2leai(deltaE+deltaG+Q14(J))                              ;#
Q14_i(J)=Q14_new(J)*njp(J+.5,temp)                                    ;#
endfor

;F2-F4-------------------------------------------------------------
;P-branch F2-F4 levels (deltaJ = +1)
P24  =fltarr(lines)                                                   ;#
P24_i=fltarr(lines)                                                   ;#;J is J''
for J=P24_J,lines-1. do begin                                         ;#
P24  (J)=Sigma4F2(b_B_v,b_D_v,eps,ga,J-.5)-Pi4F4(a_B_v,a_D_v,Y,J+.5,cor_44b(J)) -0.0001*J*J*J + 0.0027*J*J - 0.0377*J + 0.0678;#
P24  (J)=nuva2leai(deltaE+deltaG+P24(J))                              ;#
P24_i(J)=P24_new(J)*njp(J+.5,temp)                                    ;#
endfor
;R-branch F2-F4 levels (deltaJ = -1)
R24  =fltarr(lines)                                                   ;#
R24_i=fltarr(lines)                                                   ;#;J is J''
for J=R24_J,lines-1. do begin                                         ;#
R24  (J)=Sigma4F2(b_B_v,b_D_v,eps,ga,J+1.5)-Pi4F4(a_B_v,a_D_v,Y,J+.5,cor_44b(J)) -0.0001*J*J*J + 0.0027*J*J - 0.0377*J + 0.0678;#
R24  (J)=nuva2leai(deltaE+deltaG+R24(J))                              ;#
R24_i(J)=R24_new(J)*njp(J+.5,temp)                                    ;#
endfor
;Q-branch F2-F4 levels (deltaJ =  0)
Q24  =fltarr(lines)                                                   ;#
Q24_i=fltarr(lines)                                                   ;#;J is J''
for J=Q24_J,lines-1. do begin                                         ;#
Q24  (J)=Sigma4F2(b_B_v,b_D_v,eps,ga,J+0.5)-Pi4F4(a_B_v,a_D_v,Y,J+.5,cor_44b(J)) -0.0001*J*J*J + 0.0027*J*J - 0.0377*J + 0.0678;#
Q24  (J)=nuva2leai(deltaE+deltaG+Q24(J))                              ;#
Q24_i(J)=Q24_new(J)*njp(J+.5,temp)                                    ;#
endfor

;F3-F4-------------------------------------------------------------
;P-branch F3-F4 levels (deltaJ = +1)
P34  =fltarr(lines)                                                   ;#
P34_i=fltarr(lines)                                                   ;#;J is J''
for J=P34_J,lines-1. do begin                                         ;#
P34  (J)=Sigma4F3(b_B_v,b_D_v,eps,ga,J-.5)-Pi4F4(a_B_v,a_D_v,Y,J+.5,cor_44b(J)) -0.00007*J*J*J + 0.0009*J*J - 0.0065*J - 0.0406;#
P34  (J)=nuva2leai(deltaE+deltaG+P34(J))                              ;#
P34_i(J)=P34_new(J)*njp(J+.5,temp)                                    ;#
endfor
;R-branch F3-F4 levels (deltaJ = -1)
R34  =fltarr(lines)                                                   ;#
R34_i=fltarr(lines)                                                   ;#;J is J''
for J=R34_J,lines-1. do begin                                         ;#
R34  (J)=Sigma4F3(b_B_v,b_D_v,eps,ga,J+1.5)-Pi4F4(a_B_v,a_D_v,Y,J+.5,cor_44b(J)) -0.00007*J*J*J + 0.0009*J*J - 0.0065*J - 0.0406;#
R34  (J)=nuva2leai(deltaE+deltaG+R34(J))                              ;#
R34_i(J)=R34_new(J)*njp(J+.5,temp)                                    ;#
endfor
;Q-branch F3-F4 levels (deltaJ =  0)
Q34  =fltarr(lines)                                                   ;#
Q34_i=fltarr(lines)                                                   ;#;J is J''
for J=Q34_J,lines-1. do begin                                         ;#
Q34  (J)=Sigma4F3(b_B_v,b_D_v,eps,ga,J+0.5)-Pi4F4(a_B_v,a_D_v,Y,J+.5,cor_44b(J)) -0.00007*J*J*J + 0.0009*J*J - 0.0065*J - 0.0406;#
Q34  (J)=nuva2leai(deltaE+deltaG+Q34(J))                              ;#
Q34_i(J)=Q34_new(J)*njp(J+.5,temp)                                    ;#
endfor

;------------------------------------------------------------------

R11(where(R11 eq 0))=!values.F_NaN  ;Set zero wavelengths to nothing
Q11(where(Q11 eq 0))=!values.F_NaN
R21(where(R21 eq 0))=!values.F_NaN
Q21(where(Q21 eq 0))=!values.F_NaN
R31(where(R31 eq 0))=!values.F_NaN
Q31(where(Q31 eq 0))=!values.F_NaN
R41(where(R41 eq 0))=!values.F_NaN
Q41(where(Q41 eq 0))=!values.F_NaN
P11(where(P11 eq 0))=!values.F_NaN
R12(where(R12 eq 0))=!values.F_NaN
P21(where(P21 eq 0))=!values.F_NaN
R22(where(R22 eq 0))=!values.F_NaN
P31(where(P31 eq 0))=!values.F_NaN
R32(where(R32 eq 0))=!values.F_NaN
P41(where(P41 eq 0))=!values.F_NaN
R42(where(R42 eq 0))=!values.F_NaN
Q12(where(Q12 eq 0))=!values.F_NaN
P12(where(P12 eq 0))=!values.F_NaN
Q22(where(Q22 eq 0))=!values.F_NaN
P22(where(P22 eq 0))=!values.F_NaN
Q32(where(Q32 eq 0))=!values.F_NaN
P32(where(P32 eq 0))=!values.F_NaN
Q42(where(Q42 eq 0))=!values.F_NaN
P42(where(P42 eq 0))=!values.F_NaN
R13(where(R13 eq 0))=!values.F_NaN
R23(where(R23 eq 0))=!values.F_NaN
R33(where(R33 eq 0))=!values.F_NaN
R43(where(R43 eq 0))=!values.F_NaN
Q13(where(Q13 eq 0))=!values.F_NaN
Q23(where(Q23 eq 0))=!values.F_NaN
Q33(where(Q33 eq 0))=!values.F_NaN
Q43(where(Q43 eq 0))=!values.F_NaN
P13(where(P13 eq 0))=!values.F_NaN
P23(where(P23 eq 0))=!values.F_NaN
P33(where(P33 eq 0))=!values.F_NaN
P43(where(P43 eq 0))=!values.F_NaN
R14(where(R14 eq 0))=!values.F_NaN
R24(where(R24 eq 0))=!values.F_NaN
R34(where(R34 eq 0))=!values.F_NaN
R44(where(R44 eq 0))=!values.F_NaN
Q14(where(Q14 eq 0))=!values.F_NaN
Q24(where(Q24 eq 0))=!values.F_NaN
Q34(where(Q34 eq 0))=!values.F_NaN
Q44(where(Q44 eq 0))=!values.F_NaN
P14(where(P14 eq 0))=!values.F_NaN
P24(where(P24 eq 0))=!values.F_NaN
P34(where(P34 eq 0))=!values.F_NaN
P44(where(P44 eq 0))=!values.F_NaN

;Spectral linear correction factor [A]
; al=0.8

;P1=P1-al
;R1=R1-al
;P2=P2-al
;Q2=Q2-al
;R2=R2-al
;P3=P3-al
;Q3=Q3-al
;R3=R3-al

;NORMALISE
tt=total(R11_i+Q11_i+R21_i+Q21_i+R31_i+Q31_i+R41_i+Q41_i+ $
         P11_i+R12_i+P21_i+R22_i+P31_i+R32_i+P41_i+R42_i+ $
         Q12_i+P12_i+Q22_i+P22_i+Q32_i+P32_i+Q42_i+P42_i+ $
         R13_i+R23_i+R33_i+R43_i+Q13_i+Q23_i+Q33_i+Q43_i+ $
         P13_i+P23_i+P33_i+P43_i+R14_i+R24_i+R34_i+R44_i+ $
         Q14_i+Q24_i+Q34_i+Q44_i+P14_i+P24_i+P34_i+P44_i)

for J=0,lines-1. do R11_i(J)=R11_i(J)/tt
for J=0,lines-1. do Q11_i(J)=Q11_i(J)/tt
for J=0,lines-1. do R21_i(J)=R21_i(J)/tt
for J=0,lines-1. do Q21_i(J)=Q21_i(J)/tt
for J=0,lines-1. do R31_i(J)=R31_i(J)/tt
for J=0,lines-1. do Q31_i(J)=Q31_i(J)/tt
for J=0,lines-1. do R41_i(J)=R41_i(J)/tt
for J=0,lines-1. do Q41_i(J)=Q41_i(J)/tt
for J=0,lines-1. do P11_i(J)=P11_i(J)/tt
for J=0,lines-1. do R12_i(J)=R12_i(J)/tt
for J=0,lines-1. do P21_i(J)=P21_i(J)/tt
for J=0,lines-1. do R22_i(J)=R22_i(J)/tt
for J=0,lines-1. do P31_i(J)=P31_i(J)/tt
for J=0,lines-1. do R32_i(J)=R32_i(J)/tt
for J=0,lines-1. do P41_i(J)=P41_i(J)/tt
for J=0,lines-1. do R42_i(J)=R42_i(J)/tt
for J=0,lines-1. do Q12_i(J)=Q12_i(J)/tt
for J=0,lines-1. do P12_i(J)=P12_i(J)/tt
for J=0,lines-1. do Q22_i(J)=Q22_i(J)/tt
for J=0,lines-1. do P22_i(J)=P22_i(J)/tt
for J=0,lines-1. do Q32_i(J)=Q32_i(J)/tt
for J=0,lines-1. do P32_i(J)=P32_i(J)/tt
for J=0,lines-1. do Q42_i(J)=Q42_i(J)/tt
for J=0,lines-1. do P42_i(J)=P42_i(J)/tt
for J=0,lines-1. do R13_i(J)=R13_i(J)/tt
for J=0,lines-1. do R23_i(J)=R23_i(J)/tt
for J=0,lines-1. do R33_i(J)=R33_i(J)/tt
for J=0,lines-1. do R43_i(J)=R43_i(J)/tt
for J=0,lines-1. do Q13_i(J)=Q13_i(J)/tt
for J=0,lines-1. do Q23_i(J)=Q23_i(J)/tt
for J=0,lines-1. do Q33_i(J)=Q33_i(J)/tt
for J=0,lines-1. do Q43_i(J)=Q43_i(J)/tt
for J=0,lines-1. do P13_i(J)=P13_i(J)/tt
for J=0,lines-1. do P23_i(J)=P23_i(J)/tt
for J=0,lines-1. do P33_i(J)=P33_i(J)/tt
for J=0,lines-1. do P43_i(J)=P43_i(J)/tt
for J=0,lines-1. do R14_i(J)=R14_i(J)/tt
for J=0,lines-1. do R24_i(J)=R24_i(J)/tt
for J=0,lines-1. do R34_i(J)=R34_i(J)/tt
for J=0,lines-1. do R44_i(J)=R44_i(J)/tt
for J=0,lines-1. do Q14_i(J)=Q14_i(J)/tt
for J=0,lines-1. do Q24_i(J)=Q24_i(J)/tt
for J=0,lines-1. do Q34_i(J)=Q34_i(J)/tt
for J=0,lines-1. do Q44_i(J)=Q44_i(J)/tt
for J=0,lines-1. do P14_i(J)=P14_i(J)/tt
for J=0,lines-1. do P24_i(J)=P24_i(J)/tt
for J=0,lines-1. do P34_i(J)=P34_i(J)/tt
for J=0,lines-1. do P44_i(J)=P44_i(J)/tt

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function njp,J,temp

    h   = 6.6606876e-34       ;[J][s]
    c   = (2.99792458e8)*1e2  ;[cm][s-2]
    k   = 1.3806503e-23       ;[J][K-1]
    Bu  = 1.98954             ;+/- .00003 N2 X ground state rotational constant
    val = exp(-( (Bu*J*(J+1.))*h*c)/(k*temp))
return,val
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function nuva2leai, wavenumber_vacuum ;Wavenumber in vacuum to wavelength in air
;
; This is function to replace Olli's version in prog_N2_1P_v2_exp.pro.
; In this version, coefficients come from Ciddor (1996), referenced in
; Chadney et al., Angeo, 2017 (water vapour paper).
;
; Input:
;  wavenumber in vacuum, in cm^-1
; Output:
;  wavelength in Angstroms

k0=238.0185d
k1=5792105d
k2=57.362d
k3=167917.0d

wl_vac = (1d0/wavenumber_vacuum)*1d8

s2 = 1/(wl_vac*1d-4)^2
n = (((k1/(k0-s2)) + (k3/(k2-s2)))/1d8) + 1
wl_air = wl_vac/n

;print, 'Using new Wl2air.'
return, wl_air

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; NEVIN 1938 part I
function Pi4F1, B_v, D_v, Y, J,cor
delta = (6.*Y*(Y+4.))/(2.*Y*(Y-4.)+8.*J*(J+1.)+56.)
 Z1 = (Y*(Y-4.))+4.*J*(J+1.)+23/9.+(2*delta)/9.
;val =  B_v*  (J*(J+1)-1.5*sqrt(Z1))  -D_v*(J-1.5)^2.*(J-0.5)^2.    -2.*((Y*(Y-1.)-2.*J*(J+1.))/(Y*(Y-4.)+4.*J*(J+1.)))-10.
val =  B_v*  (J*(J+1)-1.5*sqrt(Z1))  -D_v*(J-1.5)^2.*(J-0.5)^2.    -cor;-12.65 +J/10. ;-10./J - 10.
return, val
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function Pi4F2, B_v, D_v, Y, J,cor
delta = (6.*Y*(Y+4.))/(2.*Y*(Y-4.)+8.*J*(J+1.)+56.)
 Z1 = (Y*(Y-4.))+4.*J*(J+1.)-5.-2.*delta
;val =  B_v*  (J*(J+1)-0.5*sqrt(Z1))  -D_v*(J-0.5)^2.*(J+0.5)^2.    +2.*((Y*(Y-1.)-2.*J*(J+1.))/(Y*(Y-4.)+4.*J*(J+1.)))
val =  B_v*  (J*(J+1)-0.5*sqrt(Z1))  -D_v*(J-0.5)^2.*(J+0.5)^2.    -cor;-5.9
return, val
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function Pi4F3, B_v, D_v, Y, J,cor
delta = (6.*Y*(Y+4.))/(2.*Y*(Y-4.)+8.*J*(J+1.)+56.)
 Z1 = (Y*(Y-4.))+4.*J*(J+1.)-5.-2.*delta
;val =  B_v*  (J*(J+1)+0.5*sqrt(Z1))  -D_v*(J+0.5)^2.*(J+1.5)^2.    +2.*((Y*(Y-1.)-2.*J*(J+1.))/(Y*(Y-4.)+4.*J*(J+1.)))
val =  B_v*  (J*(J+1)+0.5*sqrt(Z1))  -D_v*(J+0.5)^2.*(J+1.5)^2.    -cor;  -10./J
return, val
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function Pi4F4, B_v, D_v, Y, J,cor
 alpha=-0.99  ;v"=0  .66  v"=1   .75 v"=2
 beta = 0.038 ;v"=0  .027 v"=1  .034 v"=2
delta = (6.*Y*(Y+4.))/(2.*Y*(Y-4.)+8.*J*(J+1.)+56.)
 Z1 = (Y*(Y-4.))+4.*J*(J+1.)+23/9.+(2*delta)/9.
;val =  B_v*  (J*(J+1)+1.5*sqrt(Z1))  -D_v*(J+1.5)^2.*(J+2.5)^2.+alpha+beta*J  -2.*((Y*(Y-1.)-2.*J*(J+1.))/(Y*(Y-4.)+4.*J*(J+1.)))
val =  B_v*  (J*(J+1)+1.5*sqrt(Z1))  -D_v*(J+1.5)^2.*(J+2.5)^2.+alpha+beta*J  -cor;+7.-J/5.
return, val
end

;############################################################################################################
;NOTE I CHANGED all (- D_v) terms to (+ D_v) for 4 equations below to match the form found in Nevin 1938 pg. 474.
;############################################################################################################
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function Sigma4F1, B_v, D_v, epsilon, gamma,J   ; K=J-1.5
val = B_v*(J-1.5)*(J-0.5)  +  D_v*(J-1.5)^2.*(J-0.5)^2.  -1.5*epsilon*((2.*J-3.)/(2.*J))+3.*gamma*(J-1.5)
return, val                            ;Nevin 1938 pg. 474 expanded!!!
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function Sigma4F2, B_v, D_v, epsilon, gamma,J   ; K=J-0.5
val = B_v*(J-0.5)*(J+0.5)  +  D_v*(J-0.5)^2.*(J+0.5)^2.  +1.5*epsilon*((2.*J+5.)/(2.*(J+1.)))+gamma*(J-3.5)
return, val
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function Sigma4F3, B_v, D_v, epsilon, gamma,J   ; K=J+0.5
val = B_v*(J+0.5)*(J+1.5)  +  D_v*(J+0.5)^2.*(J+1.5)^2.  +1.5*epsilon*((2.*J-3.)/(2.*J))-gamma*(J+4.5)
return, val
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function Sigma4F4, B_v, D_v, epsilon, gamma,J   ; K=J+1.5
val = B_v*(J+1.5)*(J+2.5)  +  D_v*(J+1.5)^2.*(J+2.5)^2.  -1.5*epsilon*((2.*J+5.)/(2.*(J+1.)))-3.*gamma*(J+2.5)
return, val
end
;##############################################################################
function do_interpol1, lines, x_in        ; can use spline, quadratic
a=x_in                                   ; but linear works best
x_out=interpol(a,indgen(size(a,/dimensions)), indgen(lines))
return,x_out
end
;###################################################################
function do_interpol, lines, x_in        ; can use spline, quadratic
a=x_in(2,where(x_in(2,*) gt 0.))         ; but linear works best
x_out=interpol(a,indgen(size(a,/dimensions)), indgen(lines))
return,x_out
end
;###################################################################
function apol, lines, first_J, interpolated
value=fltarr(lines)
for k=0,(lines-1-first_J)/2. do value(k*2.+first_J)=interpolated(k)
value(where(value le 0.))=0.0  ;le used as lt can give -1 as value
return, value                  ;meaning last index set to zero
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro states,B1_v,D1_v,epsilon,gamma,B2_v,D2_v,Y,lines,cor_11b,cor_22b,cor_33b,cor_44b
window
tc=5
 plot, fltarr(1),xstyle=1,ystyle=1,xrange=[0,10],yrange=[-100,2e2],title='Quad state energy level splitting',charsize=1.5,ytitle='[cm-1]',xtitle='J'
oplot,indgen(lines)+0.5,sigma4F1(B1_v, D1_v, epsilon, gamma, indgen(lines)+0.5),color=250,thick=tc
oplot,indgen(lines)+0.5,sigma4F2(B1_v, D1_v, epsilon, gamma, indgen(lines)+0.5),color=250,thick=tc
oplot,indgen(lines)+0.5,sigma4F3(B1_v, D1_v, epsilon, gamma, indgen(lines)+0.5),color=250,thick=tc
oplot,indgen(lines)+0.5,sigma4F4(B1_v, D1_v, epsilon, gamma, indgen(lines)+0.5),color=250,thick=tc
oplot,indgen(lines)+0.5,Pi4F1(B2_v, D2_v, Y, indgen(lines)+0.5,cor_11b(indgen(lines))),color=150,thick=tc
oplot,indgen(lines)+0.5,Pi4F2(B2_v, D2_v, Y, indgen(lines)+0.5,cor_22b(indgen(lines))),color=150,thick=tc
oplot,indgen(lines)+0.5,Pi4F3(B2_v, D2_v, Y, indgen(lines)+0.5,cor_33b(indgen(lines))),color=150,thick=tc
oplot,indgen(lines)+0.5,Pi4F4(B2_v, D2_v, Y, indgen(lines)+0.5,cor_44b(indgen(lines))),color=150,thick=tc
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hk1,instru_grid,int_colv

   window,4,xsize=1400,ysize=400
   plot, instru_grid,int_colv,xstyle=1,ystyle=1,yrange=[0.,max(int_colv)]
;Hities spectrum
   restore, 'Sif_O2_ref.idl'
   ;oplot, wl_ref,spectrum2/500.+0.1,linestyle=1
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro investigate_10,R11,Q11,R21,Q21,R31,Q31,R41,Q41,P11,R12,P21,R22,P31,R32,P41,$
    R42,Q12,P12,Q22,P22,Q32,P32,Q42,P42,R13,R23,R33,R43,Q13,Q23,Q33,Q43,$
    P13,P23,P33,P43,R14,R24,R34,R44,Q14,Q24,Q34,Q44,P14,P24,P34,P44

restore, filename='HL_data.idl'
print,  'INVESTIGATE'
ap11=fltarr(36)
ar11=fltarr(36)
aq11=fltarr(36)
ap22=fltarr(36)
ar22=fltarr(36)
aq22=fltarr(36)
ap33=fltarr(36)
ar33=fltarr(36)
aq33=fltarr(36)
ap44=fltarr(36)
ar44=fltarr(36)
aq44=fltarr(36)

;for i=0,n_elements(ap11)-1. do ap11(i-1)=P11__(1,i-1)  -P11(i-1)
;for i=0,n_elements(ap11)-1. do ar11(i-1)=R11__(1,i-1)  -R11(i-1)
;for i=0,n_elements(ap11)-1. do aq11(i-1)=Q11__(1,i-1)  -Q11(i-1)
;for i=1,18 do ap22(i-1)=P22__(1,2*(i-1))-P22(2*(i-1))
;for i=1,18 do ar22(i-1)=R22__(1,2*(i-1))-R22(2*(i-1))
;for i=1,18 do aq22(i-1)=Q22__(1,2*i-1)  -Q22(2*i-1)
;for i=1,18 do ap33(i-1)=P33__(1,2*i-1)  -P33(2*i-1)
;for i=1,18 do ar33(i-1)=R33__(1,2*i-1)  -R33(2*i-1)
;for i=1,18 do aq33(i-1)=Q33__(1,2*(i-1))-Q33(2*(i-1))
;for i=1,18 do ap44(i-1)=P44__(1,2*(i-1))-P44(2*(i-1))
;for i=1,18 do ar44(i-1)=R44__(1,2*(i-1))-R44(2*(i-1))
;for i=1,18 do aq44(i-1)=Q44__(1,2*i-1)  -Q44(2*i-1)

;for i=1,18 do ap11(i-1)=1/(P11__(1,2*i-1)  *1e-8) -1/(P11(2*i-1)  *1e-8)
;for i=1,18 do ar11(i-1)=1/(R11__(1,2*i-1)  *1e-8) -1/(R11(2*i-1)  *1e-8)
;for i=1,18 do aq11(i-1)=1/(Q11__(1,2*(i-1))*1e-8) -1/(Q11(2*(i-1))*1e-8)
;for i=1,18 do ap22(i-1)=1/(P22__(1,2*(i-1))*1e-8) -1/(P22(2*(i-1))*1e-8)
;for i=1,18 do ar22(i-1)=1/(R22__(1,2*(i-1))*1e-8) -1/(R22(2*(i-1))*1e-8)
;for i=1,18 do aq22(i-1)=1/(Q22__(1,2*i-1)  *1e-8) -1/(Q22(2*i-1)  *1e-8)
;for i=1,18 do ap33(i-1)=1/(P33__(1,2*i-1)  *1e-8) -1/(P33(2*i-1)  *1e-8)
;for i=1,18 do ar33(i-1)=1/(R33__(1,2*i-1)  *1e-8) -1/(R33(2*i-1)  *1e-8)
;for i=1,18 do aq33(i-1)=1/(Q33__(1,2*(i-1))*1e-8) -1/(Q33(2*(i-1))*1e-8)
;for i=1,18 do ap44(i-1)=1/(P44__(1,2*(i-1))*1e-8) -1/(P44(2*(i-1))*1e-8)
;for i=1,18 do ar44(i-1)=1/(R44__(1,2*(i-1))*1e-8) -1/(R44(2*(i-1))*1e-8)
;for i=1,18 do aq44(i-1)=1/(Q44__(1,2*i-1)  *1e-8) -1/(Q44(2*i-1)  *1e-8)

for i=0,35 do ap11(i)=1/(P11__(1,i)*1e-8)-1/(P11(i)*1e-8)
for i=0,35 do ar11(i)=1/(R11__(1,i)*1e-8)-1/(R11(i)*1e-8)
for i=0,35 do aq11(i)=1/(Q11__(1,i)*1e-8)-1/(Q11(i)*1e-8)
for i=0,35 do ap22(i)=1/(P22__(1,i)*1e-8)-1/(P22(i)*1e-8)
for i=0,35 do ar22(i)=1/(R22__(1,i)*1e-8)-1/(R22(i)*1e-8)
for i=0,35 do aq22(i)=1/(Q22__(1,i)*1e-8)-1/(Q22(i)*1e-8)
for i=0,35 do ap33(i)=1/(P33__(1,i)*1e-8)-1/(P33(i)*1e-8)
for i=0,35 do ar33(i)=1/(R33__(1,i)*1e-8)-1/(R33(i)*1e-8)
for i=0,35 do aq33(i)=1/(Q33__(1,i)*1e-8)-1/(Q33(i)*1e-8)
for i=0,35 do ap44(i)=1/(P44__(1,i)*1e-8)-1/(P44(i)*1e-8)
for i=0,35 do ar44(i)=1/(R44__(1,i)*1e-8)-1/(R44(i)*1e-8)
for i=0,35 do aq44(i)=1/(Q44__(1,i)*1e-8)-1/(Q44(i)*1e-8)

;ap11(where(~finite(ap11), /null))=!values.f_nan
;ar11(where(~finite(ar11), /null))=!values.f_nan
;aq11(where(~finite(aq11), /null))=!values.f_nan
;ap22(where(~finite(ap22), /null))=!values.f_nan
;ar22(where(~finite(ar22), /null))=!values.f_nan
;aq22(where(~finite(aq22), /null))=!values.f_nan
;ap33(where(~finite(ap33), /null))=!values.f_nan
;ar33(where(~finite(ar33), /null))=!values.f_nan
;aq33(where(~finite(aq33), /null))=!values.f_nan
;ap44(where(~finite(ap44), /null))=!values.f_nan
;ar44(where(~finite(ar44), /null))=!values.f_nan
;aq44(where(~finite(aq44), /null))=!values.f_nan

cor_11=fltarr(36)
cor_22=fltarr(36)
cor_33=fltarr(36)
cor_44=fltarr(36)
for i=0,35 do cor_11(i)=mean([ap11(i),ar11(i),aq11(i)], /nan)
for i=0,35 do cor_22(i)=mean([ap22(i),ar22(i),aq22(i)], /nan)
for i=0,35 do cor_33(i)=mean([ap33(i),ar33(i),aq33(i)], /nan)
for i=0,35 do cor_44(i)=mean([ap44(i),ar44(i),aq44(i)], /nan)

cor_11(where(~finite(cor_11), /null))=0.
cor_22(where(~finite(cor_22), /null))=0.

;Correction factors done up till here

window,2
!p.multi=[0,2,2]
;!p.multi=0
lines=100.

cor_11b=do_interpol1(lines,cor_11)
plot,fltarr(1),xrange=[0,lines],yrange=[min(cor_11b),max(cor_11b)],ystyle=1
;oplot,indgen(18),1/(ap11*1e-8)
;oplot,indgen(18),1/(ar11*1e-8)
;oplot,indgen(18),1/(aq11*1e-8)
  oplot, indgen(lines),cor_11b(0:lines-1.),linestyle=1

cor_22b=do_interpol1(lines,cor_22)
plot,fltarr(1),xrange=[0,lines],yrange=[min(cor_22b),max(cor_22b)],ystyle=1
;oplot,indgen(18),ap22
;oplot,indgen(18),ar22
;oplot,indgen(18),aq22
  oplot, indgen(lines),cor_22b(0:lines-1.),linestyle=1

cor_33b=do_interpol1(lines,cor_33)
plot,fltarr(1),xrange=[0,lines],yrange=[min(cor_33b),max(cor_33b)],ystyle=1
;oplot,indgen(18),ap33
;oplot,indgen(18),ar33
;oplot,indgen(18),aq33
  oplot, indgen(lines),cor_33b(0:lines-1.),linestyle=1

cor_44b=do_interpol1(lines,cor_44)
plot,fltarr(1),xrange=[0,lines],yrange=[min(cor_44b),max(cor_44b)],ystyle=1
;oplot,indgen(18),ap44,color=150
;oplot,indgen(18),ar44,color=150
;oplot,indgen(18),aq44,color=150
  oplot, indgen(lines),cor_44b(0:lines-1.),linestyle=1

;THIS IS USED TO CREATE AN INPUT FILE TO BE USED EARLIER in code
; used for 1-1 2-2 3-3 4-4 PQR transitions when no correction of any kind
; is applied to the 1-1 2-2 3-3 4-4 energy calculations.
;save, cor_11b,cor_22b,cor_33b,cor_44b, filename='O2_E_correction_test.idl'


;CORRECTIONS for other bands than main branches
cor_12=fltarr(36)
cor_13=fltarr(36)
cor_14=fltarr(36)
cor_21=fltarr(36)
cor_23=fltarr(36)
cor_24=fltarr(36)
cor_31=fltarr(36)
cor_32=fltarr(36)
cor_34=fltarr(36)
cor_41=fltarr(36)
cor_42=fltarr(36)
cor_43=fltarr(36)
for i=0,35 do cor_12(i)=mean([1/(P12__(1,i)*1e-8)-1/(P12(i)*1e-8),1/(R12__(1,i)*1e-8)-1/(R12(i)*1e-8),1/(Q12__(1,i)*1e-8)-1/(Q12(i)*1e-8)], /nan)
for i=0,35 do cor_13(i)=mean([1/(P13__(1,i)*1e-8)-1/(P13(i)*1e-8),1/(R13__(1,i)*1e-8)-1/(R13(i)*1e-8),1/(Q13__(1,i)*1e-8)-1/(Q13(i)*1e-8)], /nan)
for i=0,35 do cor_14(i)=mean([1/(P14__(1,i)*1e-8)-1/(P14(i)*1e-8),1/(R14__(1,i)*1e-8)-1/(R14(i)*1e-8),1/(Q14__(1,i)*1e-8)-1/(Q14(i)*1e-8)], /nan)
for i=0,35 do cor_21(i)=mean([1/(P21__(1,i)*1e-8)-1/(P21(i)*1e-8),1/(R21__(1,i)*1e-8)-1/(R21(i)*1e-8),1/(Q21__(1,i)*1e-8)-1/(Q21(i)*1e-8)], /nan)
for i=0,35 do cor_23(i)=mean([1/(P23__(1,i)*1e-8)-1/(P23(i)*1e-8),1/(R23__(1,i)*1e-8)-1/(R23(i)*1e-8),1/(Q23__(1,i)*1e-8)-1/(Q23(i)*1e-8)], /nan)
for i=0,35 do cor_24(i)=mean([1/(P24__(1,i)*1e-8)-1/(P24(i)*1e-8),1/(R24__(1,i)*1e-8)-1/(R24(i)*1e-8),1/(Q24__(1,i)*1e-8)-1/(Q24(i)*1e-8)], /nan)
for i=0,35 do cor_31(i)=mean([1/(P31__(1,i)*1e-8)-1/(P31(i)*1e-8),1/(R31__(1,i)*1e-8)-1/(R31(i)*1e-8),1/(Q31__(1,i)*1e-8)-1/(Q31(i)*1e-8)], /nan)
for i=0,35 do cor_32(i)=mean([1/(P32__(1,i)*1e-8)-1/(P32(i)*1e-8),1/(R32__(1,i)*1e-8)-1/(R32(i)*1e-8),1/(Q32__(1,i)*1e-8)-1/(Q32(i)*1e-8)], /nan)
for i=0,35 do cor_34(i)=mean([1/(P34__(1,i)*1e-8)-1/(P34(i)*1e-8),1/(R34__(1,i)*1e-8)-1/(R34(i)*1e-8),1/(Q34__(1,i)*1e-8)-1/(Q34(i)*1e-8)], /nan)
for i=0,35 do cor_41(i)=mean([1/(P41__(1,i)*1e-8)-1/(P41(i)*1e-8),1/(R41__(1,i)*1e-8)-1/(R41(i)*1e-8),1/(Q41__(1,i)*1e-8)-1/(Q41(i)*1e-8)], /nan)
for i=0,35 do cor_42(i)=mean([1/(P42__(1,i)*1e-8)-1/(P42(i)*1e-8),1/(R42__(1,i)*1e-8)-1/(R42(i)*1e-8),1/(Q42__(1,i)*1e-8)-1/(Q42(i)*1e-8)], /nan)
for i=0,35 do cor_43(i)=mean([1/(P43__(1,i)*1e-8)-1/(P43(i)*1e-8),1/(R43__(1,i)*1e-8)-1/(R43(i)*1e-8),1/(Q43__(1,i)*1e-8)-1/(Q43(i)*1e-8)], /nan)

; HOW THIS IS DONE. Quickest way: copy printed values to Excel, get polynomial fit, add/subtract function from the energy levels.

!p.multi=0

;ANGSTROMS:
;print, '           11                                     22'
;for i=1,18 do print,P11__(1,2*i-1)  -P11(2*i-1)  ,R11__(1,2*i-1)  -R11(2*i-1)  ,Q11__(1,2*(i-1))-Q11(2*(i-1)),$
;                    P22__(1,2*(i-1))-P22(2*(i-1)),R22__(1,2*(i-1))-R22(2*(i-1)),Q22__(1,2*i-1)  -Q22(2*i-1)
;print, '           33                                     44'
;for i=1,18 do print,P33__(1,2*i-1)  -P33(2*i-1)  ,R33__(1,2*i-1)  -R33(2*i-1)  ,Q33__(1,2*(i-1))-Q33(2*(i-1)),$
;                    P44__(1,2*(i-1))-P44(2*(i-1)),R44__(1,2*(i-1))-R44(2*(i-1)),Q44__(1,2*i-1)  -Q44(2*i-1)
;print, '           12                                     13'
;for i=1,18 do print,P12__(1,2*i-1)  -P12(2*i-1)  ,R12__(1,2*i-1)  -R12(2*i-1)  ,Q12__(1,2*(i-1))-Q12(2*(i-1)),$
;                    P13__(1,2*i-1)  -P13(2*i-1)  ,R13__(1,2*i-1)  -R13(2*i-1)  ,Q13__(1,2*(i-1))-Q13(2*(i-1))
;print, '           14                                     21'
;for i=1,18 do print,P14__(1,2*i-1)  -P14(2*i-1)  ,R14__(1,2*i-1)  -R14(2*i-1)  ,Q14__(1,2*(i-1))-Q14(2*(i-1)),$
;                    P21__(1,2*(i-1))-P21(2*(i-1)),R21__(1,2*(i-1))-R21(2*(i-1)),Q21__(1,2*i-1)  -Q21(2*i-1)
;print, '           23                                     24'
;for i=1,18 do print,P23__(1,2*(i-1))-P23(2*(i-1)),R23__(1,2*(i-1))-R23(2*(i-1)),Q23__(1,2*i-1)  -Q23(2*i-1),$
;                    P24__(1,2*(i-1))-P24(2*(i-1)),R24__(1,2*(i-1))-R24(2*(i-1)),Q24__(1,2*i-1)  -Q24(2*i-1)
;print, '           31                                     32'
;for i=1,18 do print,P31__(1,2*i-1)  -P31(2*i-1)  ,R31__(1,2*i-1)  -R31(2*i-1)  ,Q31__(1,2*(i-1))-Q31(2*(i-1)),$
;                    P32__(1,2*i-1)  -P32(2*i-1)  ,R32__(1,2*i-1)  -R32(2*i-1)  ,Q32__(1,2*(i-1))-Q32(2*(i-1))
;print, '           34                                     41'
;for i=1,18 do print,P34__(1,2*i-1)  -P34(2*i-1)  ,R34__(1,2*i-1)  -R34(2*i-1)  ,Q34__(1,2*(i-1))-Q34(2*(i-1)),$
;                    P41__(1,2*(i-1))-P41(2*(i-1)),R41__(1,2*(i-1))-R41(2*(i-1)),Q41__(1,2*i-1)  -Q41(2*i-1)
;print, '           42                                     43'
;for i=1,18 do print,P42__(1,2*(i-1))-P42(2*(i-1)),R42__(1,2*(i-1))-R42(2*(i-1)),Q42__(1,2*i-1)  -Q42(2*i-1),$
;                    P43__(1,2*(i-1))-P43(2*(i-1)),R43__(1,2*(i-1))-R43(2*(i-1)),Q43__(1,2*i-1)  -Q43(2*i-1)

;WAVENUMBER:
print, '           11                                     22'
for i=1,18 do print,1/(P11__(1,2*i-1)  *1e-8)-1/(P11(2*i-1)  *1e-8),1/(R11__(1,2*i-1)  *1e-8)-1/(R11(2*i-1)  *1e-8),1/(Q11__(1,2*(i-1))*1e-8)-1/(Q11(2*(i-1))*1e-8),$
                    1/(P22__(1,2*(i-1))*1e-8)-1/(P22(2*(i-1))*1e-8),1/(R22__(1,2*(i-1))*1e-8)-1/(R22(2*(i-1))*1e-8),1/(Q22__(1,2*i-1)  *1e-8)-1/(Q22(2*i-1)  *1e-8)
print, '           33                                     44'
for i=1,18 do print,1/(P33__(1,2*i-1)  *1e-8)-1/(P33(2*i-1)  *1e-8),1/(R33__(1,2*i-1)  *1e-8)-1/(R33(2*i-1)  *1e-8),1/(Q33__(1,2*(i-1))*1e-8)-1/(Q33(2*(i-1))*1e-8),$
                    1/(P44__(1,2*(i-1))*1e-8)-1/(P44(2*(i-1))*1e-8),1/(R44__(1,2*(i-1))*1e-8)-1/(R44(2*(i-1))*1e-8),1/(Q44__(1,2*i-1)  *1e-8)-1/(Q44(2*i-1)  *1e-8)
print, '           12                                     13'
for i=1,18 do print,1/(P12__(1,2*i-1)  *1e-8)-1/(P12(2*i-1)  *1e-8),1/(R12__(1,2*i-1)  *1e-8)-1/(R12(2*i-1)  *1e-8),1/(Q12__(1,2*(i-1))*1e-8)-1/(Q12(2*(i-1))*1e-8),$
                    1/(P13__(1,2*i-1)  *1e-8)-1/(P13(2*i-1)  *1e-8),1/(R13__(1,2*i-1)  *1e-8)-1/(R13(2*i-1)  *1e-8),1/(Q13__(1,2*(i-1))*1e-8)-1/(Q13(2*(i-1))*1e-8)
print, '           14                                     21'
for i=1,18 do print,1/(P14__(1,2*i-1)  *1e-8)-1/(P14(2*i-1)  *1e-8),1/(R14__(1,2*i-1)  *1e-8)-1/(R14(2*i-1)  *1e-8),1/(Q14__(1,2*(i-1))*1e-8)-1/(Q14(2*(i-1))*1e-8),$
                    1/(P21__(1,2*(i-1))*1e-8)-1/(P21(2*(i-1))*1e-8),1/(R21__(1,2*(i-1))*1e-8)-1/(R21(2*(i-1))*1e-8),1/(Q21__(1,2*i-1)  *1e-8)-1/(Q21(2*i-1)  *1e-8)
print, '           23                                     24'
for i=1,18 do print,1/(P23__(1,2*(i-1))*1e-8)-1/(P23(2*(i-1))*1e-8),1/(R23__(1,2*(i-1))*1e-8)-1/(R23(2*(i-1))*1e-8),1/(Q23__(1,2*i-1)  *1e-8)-1/(Q23(2*i-1)  *1e-8),$
                    1/(P24__(1,2*(i-1))*1e-8)-1/(P24(2*(i-1))*1e-8),1/(R24__(1,2*(i-1))*1e-8)-1/(R24(2*(i-1))*1e-8),1/(Q24__(1,2*i-1)  *1e-8)-1/(Q24(2*i-1)  *1e-8)
print, '           31                                     32'
for i=1,18 do print,1/(P31__(1,2*i-1)  *1e-8)-1/(P31(2*i-1)  *1e-8),1/(R31__(1,2*i-1)  *1e-8)-1/(R31(2*i-1)  *1e-8),1/(Q31__(1,2*(i-1))*1e-8)-1/(Q31(2*(i-1))*1e-8),$
                    1/(P32__(1,2*i-1)  *1e-8)-1/(P32(2*i-1)  *1e-8),1/(R32__(1,2*i-1)  *1e-8)-1/(R32(2*i-1)  *1e-8),1/(Q32__(1,2*(i-1))*1e-8)-1/(Q32(2*(i-1))*1e-8)
print, '           34                                     41'
for i=1,18 do print,1/(P34__(1,2*i-1)  *1e-8)-1/(P34(2*i-1)  *1e-8),1/(R34__(1,2*i-1)  *1e-8)-1/(R34(2*i-1)  *1e-8),1/(Q34__(1,2*(i-1))*1e-8)-1/(Q34(2*(i-1))*1e-8),$
                    1/(P41__(1,2*(i-1))*1e-8)-1/(P41(2*(i-1))*1e-8),1/(R41__(1,2*(i-1))*1e-8)-1/(R41(2*(i-1))*1e-8),1/(Q41__(1,2*i-1)  *1e-8)-1/(Q41(2*i-1)  *1e-8)
print, '           42                                     43'
for i=1,18 do print,1/(P42__(1,2*(i-1))*1e-8)-1/(P42(2*(i-1))*1e-8),1/(R42__(1,2*(i-1))*1e-8)-1/(R42(2*(i-1))*1e-8),1/(Q42__(1,2*i-1)  *1e-8)-1/(Q42(2*i-1)  *1e-8),$
                    1/(P43__(1,2*(i-1))*1e-8)-1/(P43(2*(i-1))*1e-8),1/(R43__(1,2*(i-1))*1e-8)-1/(R43(2*(i-1))*1e-8),1/(Q43__(1,2*i-1)  *1e-8)-1/(Q43(2*i-1)  *1e-8)



print, ' ---------- '
print,P11__(1,2*18-1)-P11(2*18-1),R11__(1,2*18-1)-R11(2*18-1),Q11__(1,2*17)  -Q11(2*17),$
      P22__(1,2*17)  -P22(2*17)  ,R22__(1,2*17)  -R22(2*17)  ,Q22__(1,2*18-1)-Q22(2*18-1)
print,P33__(1,2*18-1)-P33(2*18-1),R33__(1,2*18-1)-R33(2*18-1),Q33__(1,2*17)  -Q33(2*17),$
      P44__(1,2*17)  -P44(2*17)  ,R44__(1,2*17)  -R44(2*17)  ,Q44__(1,2*18-1)-Q44(2*18-1)
print,P12__(1,2*18-1)-P12(2*18-1),R12__(1,2*18-1)-R12(2*18-1),Q12__(1,2*17)  -Q12(2*17),$
      P13__(1,2*18-1)-P13(2*18-1),R13__(1,2*18-1)-R13(2*18-1),Q13__(1,2*17)  -Q13(2*17)
print,P14__(1,2*18-1)-P14(2*18-1),R14__(1,2*18-1)-R14(2*18-1),Q14__(1,2*17)  -Q14(2*17),$
      P21__(1,2*17)  -P21(2*17)  ,R21__(1,2*17)  -R21(2*17)  ,Q21__(1,2*18-1)-Q21(2*18-1)
print,P23__(1,2*17)  -P23(2*17)  ,R23__(1,2*17)  -R23(2*17)  ,Q23__(1,2*18-1)-Q23(2*18-1),$
      P24__(1,2*17)  -P24(2*17)  ,R24__(1,2*17)  -R24(2*17)  ,Q24__(1,2*18-1)-Q24(2*18-1)
print,P31__(1,2*18-1)-P31(2*18-1),R31__(1,2*18-1)-R31(2*18-1),Q31__(1,2*17)  -Q31(2*17),$
      P32__(1,2*18-1)-P32(2*18-1),R32__(1,2*18-1)-R32(2*18-1),Q32__(1,2*17)  -Q32(2*17)
print,P34__(1,2*18-1)-P34(2*18-1),R34__(1,2*18-1)-R34(2*18-1),Q34__(1,2*17)  -Q34(2*17),$
      P41__(1,2*17)  -P41(2*17)  ,R41__(1,2*17)  -R41(2*17)  ,Q41__(1,2*18-1)-Q41(2*18-1)
print,P42__(1,2*17)  -P42(2*17)  ,R42__(1,2*17)  -R42(2*17)  ,Q42__(1,2*18-1)-Q42(2*18-1),$
      P43__(1,2*17)  -P43(2*17)  ,R43__(1,2*17)  -R43(2*17)  ,Q43__(1,2*18-1)-Q43(2*18-1)


end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;###################################################################
pro plot_br, index,x2,y2,li,pol_J,pol_I
plot,index(0,*),index(2,*), xrange=[0,x2], yrange=[0,y2]
if li gt 34 then oplot,indgen(li)*2. +pol_J,pol_I,color=100
end

;##############################################################################
pro O2_HL, lines, $
R11_J, R11_new,Q11_J, Q11_new,R21_J, R21_new,Q21_J, Q21_new,$
R31_J, R31_new,Q31_J, Q31_new,R41_J, R41_new,Q41_J, Q41_new,$
P11_J, P11_new,R12_J, R12_new,P21_J, P21_new,R22_J, R22_new,$
P31_J, P31_new,R32_J, R32_new,P41_J, P41_new,R42_J, R42_new,$
Q12_J, Q12_new,P12_J, P12_new,Q22_J, Q22_new,P22_J, P22_new,$
Q32_J, Q32_new,P32_J, P32_new,Q42_J, Q42_new,P42_J, P42_new,$
R13_J, R13_new,R23_J, R23_new,R33_J, R33_new,R43_J, R43_new,$
Q13_J, Q13_new,Q23_J, Q23_new,Q33_J, Q33_new,Q43_J, Q43_new,$
P13_J, P13_new,P23_J, P23_new,P33_J, P33_new,P43_J, P43_new,$
R14_J, R14_new,R24_J, R24_new,R34_J, R34_new,R44_J, R44_new,$
Q14_J, Q14_new,Q24_J, Q24_new,Q34_J, Q34_new,Q44_J, Q44_new,$
P14_J, P14_new,P24_J, P24_new,P34_J, P34_new,P44_J, P44_new 

;Digitised from Henriksen 1987
;O2+ 1N Honl-London factors for **** -> J" <- ****
;format = Array[3, 35] 
;
;R11 Q11 R21 Q21 R31 Q31 R41 Q41 ;All branches
;P11 R12 P21 R22 P31 R32 P41 R42
;Q12 P12 Q22 P22 Q32 P32 Q42 P42
;R13 R23 R33 R43 Q13 Q23 Q33 Q43
;P13 P23 P33 P43 R14 R24 R34 R44
;Q14 Q24 Q34 Q44 P14 P24 P34 P44

;Manipulates HL factors

;lines =200.  ;input
if (lines lt 35.) then print, "less than 35 lines!"

R11 = [ $
[ .5 ,    0.0,  0.0  ],$
[1.5 ,    0.0,  0.0  ],$
[2.5 ,    0.0,  0.0  ],$
[3.5 , 5584.1,   .350],$
[4.5 ,    0.0,  0.0  ],$
[5.5 , 5583.5,   .928],$
[6.5 ,    0.0,  0.0  ],$
[7.5 , 5582.3,  1.718],$
[8.5 ,    0.0,  0.0  ],$
[9.5 , 5580.6,  2.722],$
[10.5,    0.0,  0.0  ],$
[11.5, 5578.4,  3.932],$
[12.5,    0.0,  0.0  ],$
[13.5, 5575.6,  5.334],$
[14.5,    0.0,  0.0  ],$
[15.5, 5572.4,  6.907],$
[16.5,    0.0,  0.0  ],$
[17.5, 5568.6,  8.625],$
[18.5,    0.0,  0.0  ],$
[19.5, 5564.4, 10.466],$
[20.5,    0.0,  0.0  ],$
[21.5, 5559.8, 12.406],$
[22.5,    0.0,  0.0  ],$
[23.5, 5554.6, 14.424],$
[24.5,    0.0,  0.0  ],$
[25.5, 5549.1, 16.503],$
[26.5,    0.0,  0.0  ],$
[27.5, 5543.1, 18.628],$
[28.5,    0.0,  0.0  ],$
[29.5, 5536.7, 20.787],$
[30.5,    0.0,  0.0  ],$
[31.5, 5530.0, 22.970],$
[32.5,    0.0,  0.0  ],$
[33.5, 5522.8, 25.170],$
[34.5,    0.0,  0.0  ],$
[35.5, 5515.2, 27.380]]

Q11 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5,    0.0, 0.0   ],$
[ 2.5, 5585.8,  .825 ],$
[ 3.5,    0.0, 0.0   ],$
[ 4.5, 5587.0, 2.198 ],$
[ 5.5,    0.0, 0.0   ],$
[ 6.5, 5587.7, 3.818 ],$
[ 7.5,    0.0, 0.0   ],$
[ 8.5, 5587.8, 5.809 ],$
[ 9.5,    0.0, 0.0   ],$
[10.5, 5587.4, 8.189 ],$
[11.5,    0.0, 0.0   ],$
[12.5, 5586.4, 10.942],$
[13.5,    0.0, 0.0   ],$
[14.5, 5584.9, 14.033],$
[15.5,    0.0, 0.0   ],$
[16.5, 5583.0, 17.419],$
[17.5,    0.0, 0.0   ],$
[18.5, 5580.5, 21.052],$
[19.5,    0.0, 0.0   ],$
[20.5, 5577.6, 24.889],$
[21.5,    0.0, 0.0   ],$
[22.5, 5574.2, 28.887],$
[23.5,    0.0, 0.0   ],$
[24.5, 5570.4, 33.011],$
[25.5,    0.0, 0.0   ],$
[26.5, 5566.1, 37.230],$
[27.5,    0.0, 0.0   ],$
[28.5, 5561.4, 41.521],$
[29.5,    0.0, 0.0   ],$
[30.5, 5556.3, 45.863],$
[31.5,    0.0, 0.0   ],$
[32.5, 5550.7, 50.240],$
[33.5,    0.0, 0.0   ],$
[34.5, 5544.8, 54.640],$
[35.5,    0.0, 0.0   ]]

R21 = [ $
[ .5,     0.0, 0.0   ],$
[1.5,     0.0, 0.0   ],$
[2.5,  5581.8, 0.263 ],$
[3.5,     0.0, 0.0   ],$
[4.5,  5579.9, 1.014 ],$
[5.5,     0.0, 0.0   ],$
[6.5,  5577.4, 1.863 ],$
[7.5,     0.0, 0.0   ],$
[8.5,  5574.4, 2.714 ],$
[9.5,     0.0, 0.0   ],$
[10.5, 5570.8, 3.511 ],$
[11.5,    0.0, 0.0   ],$
[12.5, 5566.8, 4.222 ],$
[13.5,    0.0, 0.0   ],$
[14.5, 5562.2, 4.821 ],$
[15.5,    0.0, 0.0   ],$
[16.5, 5557.2, 5.320 ],$
[17.5,    0.0, 0.0   ],$
[18.5, 5551.6, 5.705 ],$
[19.5,    0.0, 0.0   ],$
[20.5, 5545.7, 5.990 ],$
[21.5,    0.0, 0.0   ],$
[22.5, 5539.2, 6.189 ],$
[23.5,    0.0, 0.0   ],$
[24.5, 5532.4, 6.313 ],$
[25.5,    0.0, 0.0   ],$
[26.5, 5525.1, 6.375 ],$
[27.5,    0.0, 0.0   ],$
[28.5, 5517.4, 6.387 ],$
[29.5,    0.0, 0.0   ],$
[30.5, 5509.3, 6.359 ],$
[31.5,    0.0, 0.0   ],$
[32.5, 5500.8, 6.299 ],$
[33.5,    0.0, 0.0   ],$
[34.5, 5491.9, 6.215 ],$
[35.5,    0.0, 0.0   ]]

Q21 = [ $
[ 0.5,    0.0, 0.0    ],$
[ 1.5,    0.0, 0.0    ],$
[ 2.5,    0.0, 0.0    ],$
[ 3.5, 5584.0, 2.773  ],$
[ 4.5,    0.0, 0.0    ],$
[ 5.5, 5583.4, 4.914  ],$
[ 6.5,    0.0, 0.0    ],$
[ 7.5, 5582.2, 6.874  ],$
[ 8.5,    0.0, 0.0    ],$
[ 9.5, 5580.5, 8.655  ],$
[10.5,    0.0, 0.0    ],$
[11.5, 5518.2, 10.218 ],$
[12.5,    0.0, 0.0    ],$
[13.5, 5575.4, 11.536 ],$
[14.5,    0.0, 0.0    ],$
[15.5, 5572.1, 12.602 ],$
[16.5,    0.0, 0.0    ],$
[17.5, 5568.4, 13.424 ],$
[18.5,    0.0, 0.0    ],$
[19.5, 5564.1, 14.023 ],$
[20.5,    0.0, 0.0    ],$
[21.5, 5559.4, 14.427 ],$
[22.5,    0.0, 0.0    ],$
[23.5, 5554.3, 14.664 ],$
[24.5,    0.0, 0.0    ],$
[25.5, 5548.6, 14.765 ],$
[26.5,    0.0, 0.0    ],$
[27.5, 5542.6, 14.754 ],$
[28.5,    0.0, 0.0    ],$
[29.5, 5536.1, 14.654 ],$
[30.5,    0.0, 0.0    ],$
[31.5, 5529.3, 14.486 ],$
[32.5,    0.0, 0.0    ],$
[33.5, 5522.0, 14.266 ],$
[34.5,    0.0, 0.0    ],$
[35.5, 5514.3, 14.006 ]]

R31 = [ $
[  .5,    0.0, 0.0    ],$
[ 1.5,    0.0, 0.0    ],$
[ 2.5,    0.0, 0.0    ],$
[ 3.5, 5577.0,  .334  ],$
[ 4.5,    0.0, 0.0    ],$
[ 5.5, 5573.3,  .683  ],$
[ 6.5,    0.0, 0.0    ],$
[ 7.5, 5569.0,  .942  ],$
[ 8.5,    0.0, 0.0    ],$
[ 9.5, 5564.1, 1.102  ],$
[10.5,    0.0, 0.0    ],$
[11.5, 5558.8, 1.175  ],$
[12.5,    0.0, 0.0    ],$
[13.5, 5552.9, 1.183  ],$
[14.5,    0.0, 0.0    ],$
[15.5, 5546.6, 1.144  ],$
[16.5,    0.0, 0.0    ],$
[17.5, 5539.7, 1.077  ],$
[18.5,    0.0, 0.0    ],$
[19.5, 5532.5,  .995  ],$
[20.5,    0.0, 0.0    ],$
[21.5, 5524.7,  .906  ],$
[22.5,    0.0, 0.0    ],$
[23.5, 5516.6,  .818  ],$
[24.5,    0.0, 0.0    ],$
[25.5, 5508.0,  .733  ],$
[26.5,    0.0, 0.0    ],$
[27.5, 5499.0,  .654  ],$
[28.5,    0.0, 0.0    ],$
[29.5, 5489.6,  .581  ],$
[30.5,    0.0, 0.0    ],$
[31.5, 5479.8,  .516  ],$
[32.5,    0.0, 0.0    ],$
[33.5, 5469.7,  .458  ],$
[34.5,    0.0, 0.0    ],$
[35.5, 5459.1,  .406  ]]

Q31 = [ $
[  .5,    0.0, 0.0    ],$
[ 1.5,    0.0, 0.0    ],$
[ 2.5, 5581.8,  .900  ],$
[ 3.5,    0.0, 0.0    ],$
[ 4.5, 5519.9, 2.102  ],$
[ 5.5,    0.0, 0.0    ],$
[ 6.5, 5577.4, 2.827  ],$
[ 7.5,    0.0, 0.0    ],$
[ 8.5, 5574.4, 3.221  ],$
[ 9.5,    0.0, 0.0    ],$
[10.5, 5570.8, 3.370  ],$
[11.5,    0.0, 0.0    ],$
[12.5, 5566.7, 3.343  ],$
[13.5,    0.0, 0.0    ],$
[14.5, 5562.2, 3.200  ],$
[15.5,    0.0, 0.0    ],$
[16.5, 5557.1, 2.987  ],$
[17.5,    0.0, 0.0    ],$
[18.5, 5551.5, 2.738  ],$
[19.5,    0.0, 0.0    ],$
[20.5, 5545.5, 2.479  ],$
[21.5,    0.0, 0.0    ],$
[22.5, 5539.0, 2.225  ],$
[23.5,    0.0, 0.0    ],$
[24.5, 5532.1, 1.984  ],$
[25.5,    0.0, 0.0    ],$
[26.5, 5524.8, 1.762  ],$
[27.5,    0.0, 0.0    ],$
[28.5, 5517.0, 1.561  ],$
[29.5,    0.0, 0.0    ],$
[30.5, 5508.8, 1.382  ],$
[31.5,    0.0, 0.0    ],$
[32.5, 5500.2, 1.222  ],$
[33.5,    0.0, 0.0    ],$
[34.5, 5491.2, 1.082  ],$
[35.5,    0.0, 0.0    ]]

R41 =[ $
[  .5,    0.0, 0.0    ],$
[ 1.5,    0.0, 0.0    ],$
[ 2.5, 5574.9,  .026  ],$
[ 3.5,    0.0, 0.0    ],$
[ 4.5, 4569.9,  .081  ],$
[ 5.5,    0.0, 0.0    ],$
[ 6.5, 5564.3,  .112  ],$
[ 7.5,    0.0, 0.0    ],$
[ 8.5, 5558.2,  .121  ],$
[ 9.5,    0.0, 0.0    ],$
[10.5, 5551.6,  .115  ],$
[11.5,    0.0, 0.0    ],$
[12.5, 5544.4,  .101  ],$
[13.5,    0.0, 0.0    ],$
[14.5, 5536.8,  .086  ],$
[15.5,    0.0, 0.0    ],$
[16.5, 5528.7,  .071  ],$
[17.5,    0.0, 0.0    ],$
[18.5, 5520.1,  .057  ],$
[19.5,    0.0, 0.0    ],$
[20.5, 5511.1,  .045  ],$
[21.5,    0.0, 0.0    ],$
[22.5, 5501.7,  .036  ],$
[23.5,    0.0, 0.0    ],$
[24.5, 5491.9,  .028  ],$
[25.5,    0.0, 0.0    ],$
[26.5, 5481.6,  .022  ],$
[27.5,    0.0, 0.0    ],$
[28.5, 5470.9,  .017  ],$
[29.5,    0.0, 0.0    ],$
[30.5, 5459.9,  .014  ],$
[31.5,    0.0, 0.0    ],$
[32.5, 5448.5,  .011  ],$
[33.5,    0.0, 0.0    ],$
[34.5, 5436.7,  .009  ],$
[35.5,    0.0, 0.0   ]]

Q41 =[ $
[   .5,    0.0, 0.0   ],$
[  1.5,    0.0, 0.0   ],$
[  2.5,    0.0, 0.0   ],$
[  3.5, 5577.2,  .273 ],$
[  4.5,    0.0, 0.0   ],$
[  5.5, 5573.4,  .393 ],$
[  6.5,    0.0, 0.0   ],$
[  7.5, 5569.1,  .417 ],$
[  8.5,    0.0, 0.0   ],$
[  9.5, 5564.2,  .390 ],$
[ 10.5,    0.0, 0.0   ],$
[ 11.5, 5558.9,  .341 ],$
[ 12.5,    0.0, 0.0   ],$
[ 13.5, 5553.0,  .285 ],$
[ 14.5,    0.0, 0.0   ],$
[ 15.5, 5546.6,  .232 ],$
[ 16.5,    0.0, 0.0   ],$
[ 17.5, 5539.8,  .185 ],$
[ 18.5,    0.0, 0.0   ],$
[ 19.5, 5532.5,  .147 ],$
[ 20.5,    0.0, 0.0   ],$
[ 21.5, 5524.7,  .115 ],$
[ 22.5,    0.0, 0.0   ],$
[ 23.5, 5516.5,  .090 ],$
[ 24.5,    0.0, 0.0   ],$
[ 25.5, 5507.9,  .071 ],$
[ 26.5,    0.0, 0.0   ],$
[ 27.5, 5498.8,  .055 ],$
[ 28.5,    0.0, 0.0   ],$
[ 29.5, 5489.4,  .043 ],$
[ 30.5,    0.0, 0.0   ],$
[ 31.5, 5479.5,  .034 ],$
[ 32.5,    0.0, 0.0   ],$
[ 33.5, 5469.3,  .027 ],$
[ 34.5,    0.0, 0.0   ],$
[ 35.5, 5458.6,  .021 ]]

P11 =[ $
[   .5,    0.0,  0.0   ],$
[  1.5,    0.0,  0.0   ],$
[  2.5,    0.0,  0.0   ],$
[  3.5, 5588.0,  2.040 ],$
[  4.5,    0.0,  0.0   ],$
[  5.5, 5590.6,  2.593 ],$
[  6.5,    0.0,  0.0   ],$
[  7.5, 5592.5,  3.463 ],$
[  8.5,    0.0,  0.0   ],$
[  9.5, 5593.9,  4.567 ],$
[ 10.5,    0.0,  0.0   ],$
[ 11.5, 5594.8,  5.874 ],$
[ 12.5,    0.0,  0.0   ],$
[ 13.5, 5595.1,  7.361 ],$
[ 14.5,    0.0,  0.0   ],$
[ 15.5, 5595.0,  9.003 ],$
[ 16.5,    0.0,  0.0   ],$
[ 17.5, 5594.3, 10.775 ],$
[ 18.5,    0.0,  0.0   ],$
[ 19.5, 5593.2, 12.655 ],$
[ 20.5,    0.0,  0.0   ],$
[ 21.5, 5591.5, 14.620 ],$
[ 22.5,    0.0,  0.0   ],$
[ 23.5, 5589.4, 16.652 ],$
[ 24.5,    0.0,  0.0   ],$
[ 25.5, 5586.9, 18.736 ],$
[ 26.5,    0.0,  0.0   ],$
[ 27.5, 5583.9, 20.858 ],$
[ 28.5,    0.0,  0.0   ],$
[ 29.5, 5580.5, 23.008 ],$
[ 30.5,    0.0,  0.0   ],$
[ 31.5, 5576.6, 25.178 ],$
[ 32.5,    0.0,  0.0   ],$
[ 33.5, 5572.3, 27.360 ],$
[ 34.5,    0.0,  0.0   ],$
[ 35.5, 5567.6, 29.551 ]]

R12 =[ $
[  .5,    0.0, 0.0   ],$
[ 1.5, 5600.0,  .325 ],$
[ 2.5,    0.0, 0.0   ],$
[ 3.5, 5600.1, 1.193 ],$
[ 4.5,    0.0, 0.0   ],$
[ 5.5, 5599.8, 2.167 ],$
[ 6.5,    0.0, 0.0   ],$
[ 7.5, 5599.0, 3.156 ],$
[ 8.5,    0.0, 0.0   ],$
[ 9.5, 5597.7, 4.098 ],$
[10.5,    0.0, 0.0   ],$
[11.5, 5596.1, 4.950 ],$
[12.5,    0.0, 0.0   ],$
[13.5, 5594.0, 5.687 ],$
[14.5,    0.0, 0.0   ],$
[15.5, 5591.4, 6.296 ],$
[16.5,    0.0, 0.0   ],$
[17.5, 5588.5, 6.778 ],$
[18.5,   0.0,  0.0   ],$
[19.5, 5585.1, 7.141 ],$
[20.5,    0.0, 0.0   ],$
[21.5, 5581.4, 7.398 ],$
[22.5,    0.0, 0.0   ],$
[23.5, 5577.2, 7.564 ],$
[24.5,    0.0, 0.0   ],$
[25.5, 5572.6, 7.653 ],$
[26.5,    0.0, 0.0   ],$
[27.5, 5567.7, 7.679 ],$
[28.5,    0.0, 0.0   ],$
[29.5, 5562.3, 7.655 ],$
[30.5,    0.0, 0.0   ],$
[31.5, 5556.6, 7.591 ],$
[32.5,    0.0, 0.0   ],$
[33.5, 5550.5, 7.497 ],$
[34.5,    0.0, 0.0   ],$
[35.5, 5544.0, 7.379 ]]

P21 =[ $
[  .5,    0.0, 0.0   ],$
[ 1.5,    0.0, 0.0   ],$
[ 2.5, 5585.7, 3.678 ],$
[ 3.5,    0.0, 0.0   ],$
[ 4.5, 5586.9, 4.193 ],$
[ 5.5,    0.0, 0.0   ],$
[ 6.5, 5587.5, 4.979 ],$
[ 7.5,    0.0, 0.0   ],$
[ 8.5, 5587.6, 5.780 ],$
[ 9.5,    0.0, 0.0   ],$
[10.5, 5587.2, 6.516 ],$
[11.5,    0.0, 0.0   ],$
[12.5, 5586.2, 7.148 ],$
[13.5,    0.0, 0.0   ],$
[14.5, 5584.8, 7.660 ],$
[15.5,    0.0, 0.0   ],$
[16.5, 5582.8, 8.051 ],$
[17.5,    0.0, 0.0   ],$
[18.5, 5580.3, 8.328 ],$
[19.5,    0.0, 0.0   ],$
[20.5, 5577.3, 8.503 ],$
[21.5,    0.0, 0.0   ],$
[22.5, 5573.9, 8.592 ],$
[23.5,    0.0, 0.0   ],$
[24.5, 5570.0, 8.608 ],$
[25.5,    0.0, 0.0   ],$
[26.5, 5565.7, 8.566 ],$
[27.5,    0.0, 0.0   ],$
[28.5, 5560.9, 8.479 ],$
[29.5,    0.0, 0.0   ],$
[30.5, 5555.7, 8.356 ],$
[31.5,    0.0, 0.0   ],$
[32.5, 5550.0, 8.207 ],$
[33.5,    0.0, 0.0   ],$
[34.5, 5543.9, 8.038 ],$
[35.5,    0.0, 0.0   ]]

R22 =[ $
[  .5,    0.0,  0.0   ],$
[ 1.5,    0.0,  0.0   ],$
[ 2.5, 5597.6,   .071 ],$
[ 3.5,    0.0,  0.0   ],$
[ 4.5, 5595.9,   .056 ],$
[ 5.5,    0.0,  0.0   ],$
[ 6.5, 5593.8,   .008 ],$
[ 7.5,    0.0,  0.0   ],$
[ 8.5, 5591.2,   .020 ],$
[ 9.5,    0.0,  0.0   ],$
[10.5, 5588.2,   .179 ],$
[11.5,    0.0,  0.0   ],$
[12.5, 5584.7,   .551 ],$
[13.5,    0.0,  0.0   ],$
[14.5, 5580.8,  1.172 ],$
[15.5,    0.0,  0.0   ],$
[16.5, 5576.5,  2.055 ],$
[17.5,    0.0,  0.0   ],$
[18.5, 5571.8,  3.191 ],$
[19.5,    0.0,  0.0   ],$
[20.5, 5566.7,  4.560 ],$
[21.5,    0.0,  0.0   ],$
[22.5, 5561.2,  6.134 ],$
[23.5,    0.0,  0.0   ],$
[24.5, 5555.2,  7.883 ],$
[25.5,    0.0,  0.0   ],$
[26.5, 5548.9,  9.777 ],$
[27.5,    0.0,  0.0   ],$
[28.5, 5542.2, 11.790 ],$ 
[29.5,    0.0,  0.0   ],$
[30.5, 5535.2, 13.897 ],$
[31.5,    0.0,  0.0   ],$
[32.5, 5527.7, 16.078 ],$
[33.5,    0.0,  0.0   ],$
[34.5, 5519.9, 18.316 ],$
[35.5,    0.0,  0.0   ]]
 
P31 =[ $
[   .5,    0.0, 0.0   ],$
[  1.5,    0.0, 0.0   ],$
[  2.5,    0.0, 0.0   ],$
[  3.5, 5584.1, 2.230 ],$
[  4.5,    0.0, 0.0   ],$
[  5.5, 5583.4, 2.489 ],$
[  6.5,    0.0, 0.0   ],$
[  7.5, 5582.2, 2.584 ],$
[  8.5,    0.0, 0.0   ],$
[  9.5, 5580.5, 2.564 ],$
[ 10.5,    0.0, 0.0   ],$
[ 11.5, 5578.2, 2.460 ],$
[ 12.5,    0.0, 0.0   ],$
[ 13.5, 5575.4, 2.301 ],$
[ 14.5,    0.0, 0.0   ],$
[ 15.5, 5572.1, 2.113 ],$
[ 16.5,    0.0, 0.0   ],$
[ 17.5, 5568.3, 1.913 ],$
[ 18.5,    0.0, 0.0   ],$
[ 19.5, 5564.0, 1.714 ],$
[ 20.5,    0.0, 0.0   ],$
[ 21.5, 5559.3, 1.526 ],$
[ 22.5,    0.0, 0.0   ],$
[ 23.5, 5554.1, 1.351 ],$
[ 24.5,    0.0, 0.0   ],$
[ 25.5, 5548.4, 1.193 ],$
[ 26.5,    0.0, 0.0   ],$
[ 27.5, 5542.3, 1.051 ],$
[ 28.5,    0.0, 0.0   ],$
[ 29.5, 5535.8,  .926 ],$
[ 30.5,    0.0, 0.0   ],$
[ 31.5, 5528.8,  .816 ],$
[ 32.5,    0.0, 0.0   ],$
[ 33.5, 5521.4,  .719 ],$
[ 34.5,    0.0, 0.0   ],$
[ 35.5, 5513.6,  .635 ]]

R32 =[ $
[   .5,    0.0,  0.0   ],$
[  1.5, 5596.0,   .184 ],$
[  2.5,     0.0, 0.0   ],$
[  3.5, 5592.9,   .709 ],$
[  4.5,     0.0, 0.0   ],$
[  5.5, 5589.4,  1.414 ],$
[  6.5,     0.0, 0.0   ],$
[  7.5, 5585.5,  2.247 ],$
[  8.5,     0.0, 0.0   ],$
[  9.5, 5581.1,  3.141 ],$
[ 10.5,     0.0, 0.0   ],$
[ 11.5, 5576.3,  4.032 ],$
[ 12.5,     0.0, 0.0   ],$
[ 13.5, 5571.1,  4.866 ],$
[ 14.5,     0.0, 0.0   ],$
[ 15.5, 5565.4,  5.610 ],$
[ 16.5,     0.0, 0.0   ],$
[ 17.5, 5559.4,  6.242 ],$
[ 18.5,     0.0, 0.0   ],$
[ 19.5, 5552.9,  6.759 ],$
[ 20.5,     0.0, 0.0   ],$
[ 21.5, 5546.1,  7.162 ],$
[ 22.5,     0.0, 0.0   ],$
[ 23.5, 5538.8,  7.461 ],$
[ 24.5,     0.0, 0.0   ],$
[ 25.5, 5531.2,  7.669 ],$
[ 26.5,     0.0, 0.0   ],$
[ 27.5, 5523.2,  7.799 ],$
[ 28.5,     0.0, 0.0   ],$
[ 29.5, 5514.8,  7.863 ],$
[ 30.5,     0.0, 0.0   ],$
[ 31.5, 5506.0,  7.873 ],$
[ 32.5,     0.0, 0.0   ],$
[ 33.5, 5496.9,  7.839 ],$
[ 34.5,     0.0, 0.0   ],$
[ 35.5, 5487.4,  7.770 ]]

P41 =[ $
[  .5,    0.0, 0.0   ],$
[ 1.5,    0.0, 0.0   ],$
[ 2.5, 5582.0,  .308 ],$
[ 3.5,    0.0, 0.0   ],$
[ 4.5, 5580.0,  .414 ],$
[ 5.5,    0.0, 0.0   ],$
[ 6.5, 5577.5,  .401 ],$
[ 7.5,    0.0, 0.0   ],$
[ 8.5, 5574.5,  .355 ],$
[ 9.5,    0.0, 0.0   ],$
[10.5, 5570.9,  .299 ],$
[11.5,    0.0, 0.0   ],$
[12.5, 5566.9,  .244 ],$
[13.5,    0.0, 0.0   ],$
[14.5, 5562.3,  .195 ],$
[15.5,    0.0, 0.0   ],$
[16.5, 5557.2,  .153 ],$
[17.5,    0.0, 0.0   ],$
[18.5, 5551.6,  .120 ],$
[19.5,    0.0, 0.0   ],$
[20.5, 5545.5,  .093 ],$
[21.5,    0.0, 0.0   ],$
[22.5, 5539.0,  .072 ],$
[23.5,    0.0, 0.0   ],$
[24.5, 5532.1,  .056 ],$
[25.5,    0.0, 0.0   ],$
[26.5, 5524.7,  .044 ],$
[27.5,    0.0, 0.0   ],$
[28.5, 5516.8,  .034 ],$
[29.5,    0.0, 0.0   ],$
[30.5, 5508.6,  .027 ],$
[31.5,    0.0, 0.0   ],$
[32.5, 5499.9,  .021 ],$
[33.5,    0.0, 0.0   ],$
[34.5, 5490.8,  .017 ],$
[35.5,    0.0, 0.0  ]]

R42 =[ $
[   .5,    0.0, 0.0   ],$
[  1.5,    0.0, 0.0   ],$
[  2.5, 5590.7,  .502 ],$
[  3.5,    0.0, 0.0   ],$
[  4.5, 5585.9,  .904 ],$
[  5.5,    0.0, 0.0   ],$
[  6.5, 5580.6, 1.175 ],$
[  7.5,    0.0, 0.0   ],$
[  8.5, 5574.9, 1.328 ],$
[  9.5,    0.0, 0.0   ],$
[ 10.5, 5568.8, 1.386 ],$
[ 11.5,    0.0, 0.0   ],$
[ 12.5, 5562.2, 1.375 ],$
[ 13.5,    0.0, 0.0   ],$
[ 14.5, 5555.2, 1.316 ],$
[ 15.5,    0.0, 0.0   ],$
[ 16.5, 5547.9, 1.228 ],$
[ 17.5,    0.0, 0.0   ],$
[ 18.5, 5540.1, 1.126 ],$
[ 19.5,    0.0, 0.0   ],$
[ 20.5, 5531.9, 1.020 ],$
[ 21.5,    0.0, 0.0   ],$
[ 22.5, 5523.3,  .915 ],$
[ 23.5,    0.0, 0.0   ],$
[ 24.5, 5514.4,  .816 ],$
[ 25.5,    0.0, 0.0   ],$
[ 26.5, 5505.1,  .725 ],$
[ 27.5,    0.0, 0.0   ],$
[ 28.5, 5495.4,  .643 ],$
[ 29.5,    0.0, 0.0   ],$
[ 30.5, 5485.3,  .569 ],$
[ 31.5,    0.0, 0.0   ],$
[ 32.5, 5474.9,  .503 ],$
[ 33.5,    0.0, 0.0   ],$
[ 34.5, 5464.1,  .445 ],$
[ 35.5,    0.0, 0.0   ]]

Q12 = [$
[ 0.5,    0.0,  0.0   ],$
[ 1.5,    0.0,  0.0   ],$
[ 2.5, 5601.7,  2.001 ],$
[ 3.5,    0.0,  0.0   ],$
[ 4.5, 5603.1,  3.879 ],$
[ 5.5,    0.0,  0.0   ],$
[ 6.5, 5604.1,  5.747 ],$
[ 7.5,    0.0,  0.0   ],$
[ 8.5, 5604.7,  7.532 ],$
[ 9.5,    0.0,  0.0   ],$
[10.5, 5604.8,  9.161 ],$
[11.5,    0.0,  0.0   ],$
[12.5, 5604.5, 10.581 ],$
[13.5,    0.0,  0.0   ],$
[14.5, 5603.7, 11.76  ],$
[15.5,    0.0,  0.0   ],$
[16.5, 5602.5, 12.710 ],$
[17.5,    0.0,  0.0   ],$
[18.5, 5600.9, 13.426 ],$
[19.5,    0.0,  0.0   ],$
[20.5, 5598.9, 13.937 ],$
[21.5,    0.0,  0.0   ],$
[22.5, 5596.4, 14.269 ],$
[23.5,    0.0,  0.0   ],$
[24.5, 5593.6, 14.451 ],$
[25.5,    0.0,  0.0   ],$
[26.5, 5590.3, 14.509 ],$
[27.5,    0.0,  0.0   ],$
[28.5, 5586.7, 14.467 ],$
[29.5,    0.0,  0.0   ],$
[30.5, 5582.6, 14.347 ],$
[31.5,    0.0,  0.0   ],$
[32.5, 5578.2, 14.166 ],$
[33.5,    0.0,  0.0   ],$
[34.5, 5573.3, 13.939 ],$
[35.5,    0.0,  0.0   ]]

P12 =[ $
[  .5,    0.0, 0.0   ],$
[ 1.5,    0.0, 0.0   ],$
[ 2.5,    0.0, 0.0   ],$
[ 3.5, 5604.0, 1.778 ],$
[ 4.5,    0.0, 0.0   ],$
[ 5.5, 5606.8, 2.650 ],$
[ 6.5,    0.0, 0.0   ],$
[ 7.5, 5609.2, 3.474 ],$
[ 8.5,    0.0, 0.0   ],$
[ 9.5, 5611.1, 4.234 ],$
[10.5,    0.0, 0.0   ],$
[11.5, 5612.6, 4.905 ],$
[12.5,    0.0, 0.0   ],$
[13.5, 5613.6, 5.471 ],$
[14.5,    0.0, 0.0   ],$
[15.5, 5614.2, 5.927 ],$
[16.5,    0.0, 0.0   ],$
[17.5, 5614.3, 6.275 ],$
[18.5,    0.0, 0.0   ],$
[19.5, 5614.1, 6.526 ],$
[20.5,    0.0, 0.0   ],$
[21.5, 5613.4, 6.691 ],$
[22.5,    0.0, 0.0   ],$
[23.5, 5612.3, 6.783 ],$
[24.5,    0.0, 0.0   ],$
[25.5, 5610.8, 6.814 ],$
[26.5,    0.0, 0.0   ],$
[27.5, 5608.8, 6.797 ],$
[28.5,    0.0, 0.0   ],$
[29.5, 5606.5, 6.741 ],$
[30.5,    0.0, 0.0   ],$
[31.5, 5603.7, 6.655 ],$
[32.5,    0.0, 0.0   ],$
[33.5, 5600.5, 6.547 ],$
[34.5,    0.0, 0.0   ],$
[35.5, 5596.9, 6.422]]

Q22 = [ $
[  .5,    0.0,  0.0   ],$
[ 1.5, 5599.8,   .148 ],$
[ 2.5,    0.0,  0.0   ],$
[ 3.5, 5600.0,   .192 ],$
[ 4.5,    0.0,  0.0   ],$
[ 5.5, 5599.6,   .073 ],$
[ 6.5,    0.0,  0.0   ],$
[ 7.5, 5598.8,   .001 ],$
[ 8.5,    0.0,  0.0   ],$
[ 9.5, 5597.6,   .168 ],$
[10.5,    0.0,  0.0   ],$
[11.5, 5595.9,   .727 ],$
[12.5,    0.0,  0.0   ],$
[13.5, 5593.8,  1.774 ],$
[14.5,    0.0,  0.0   ],$
[15.5, 5591.2,  3.349 ],$
[16.5,    0.0,  0.0   ],$
[17.5, 5588.2,  5.449 ],$
[18.5,    0.0,  0.0   ],$
[19.5, 5584.8,  8.039 ],$
[20.5,    0.0,  0.0   ],$
[21.5, 5581.0, 11.066 ],$
[22.5,    0.0,  0.0   ],$
[23.5, 5576.8, 14.468 ],$
[24.5,    0.0,  0.0   ],$
[25.5, 5572.2, 18.186 ],$
[26.5,    0.0,  0.0   ],$
[27.5, 5567.1, 22.161 ],$
[28.5,    0.0,  0.0   ],$
[29.5, 5561.7, 26.341 ],$
[30.5,    0.0,  0.0   ],$
[31.5, 5555.9, 30.684 ],$
[32.5,    0.0,  0.0   ],$
[33.5, 5549.7, 35.151 ],$
[34.5,    0.0,  0.0   ],$
[35.5, 5543.1, 39.713]]


P22 = [ $
[  .5,    0.0,  0.0   ],$
[ 1.5,    0.0,  0.0   ],$
[ 2.5, 5601.5,   .096 ],$
[ 3.5,    0.0,  0.0   ],$
[ 4.5, 5603.0,   .068 ],$
[ 5.5,    0.0,  0.0   ],$
[ 6.5, 5604.0,   .006 ],$
[ 7.5,    0.0,  0.0   ],$
[ 8.5, 5604.5,   .032 ],$
[ 9.5,    0.0,  0.0   ],$
[10.5, 5604.6,   .232 ],$
[11.5,    0.0,  0.0   ],$
[12.5, 5604.3,   .666 ],$
[13.5,    0.0,  0.0   ],$
[14.5, 5603.5,  1.363 ],$
[15.5,    0.0,  0.0   ],$
[16.5, 5602.3,  2.329 ],$
[17.5,    0.0,  0.0   ],$
[18.5, 5600.7,  3.550 ],$
[19.5,    0.0,  0.0   ],$
[20.5, 5598.6,  5.002 ],$
[21.5,    0.0,  0.0   ],$
[22.5, 5596.1,  6.654 ],$
[23.5,    0.0,  0.0   ],$
[24.5, 5593.2,  8.476 ],$
[25.5,    0.0,  0.0   ],$
[26.5, 5589.9, 10.436 ],$
[27.5,    0.0,  0.0   ],$
[28.5, 5586.1, 12.508 ],$
[29.5,    0.0,  0.0   ],$
[30.5, 5582.0, 14.668 ],$
[31.5,    0.0,  0.0   ],$
[32.5, 5577.4, 16.896 ],$
[33.5,    0.0,  0.0   ],$
[34.5, 5572.5, 19.175 ],$
[35.5,    0.0,  0.0   ]]

Q32 = [ $
[  .5,    0.0,  0.0   ],$
[ 1.5,    0.0,  0.0   ],$
[ 2.5, 5597.7,  1.442 ],$
[ 3.5,    0.0,  0.0   ],$
[ 4.5, 5595.9,  2.911 ],$
[ 5.5,    0.0,  0.0   ],$
[ 6.5, 5593.8,  4.686 ],$
[ 7.5,    0.0,  0.0   ],$
[ 8.5, 5591.2,  6.639 ],$
[ 9.5,    0.0,  0.0   ],$
[10.5, 5588.2,  8.626 ],$
[11.5,    0.0,  0.0   ],$
[12.5, 5584.7, 10.521 ],$
[13.5,    0.0,  0.0   ],$
[14.5, 5580.8, 12.234 ],$
[15.5,    0.0,  0.0   ],$
[16.5, 5576.4, 13.710 ],$
[17.5,    0.0,  0.0   ],$
[18.5, 5571.7, 14.927 ],$
[19.5,    0.0,  0.0   ],$
[20.5, 5566.5, 15.888 ],$
[21.5,    0.0,  0.0   ],$
[22.5, 5561.0, 16.610 ],$
[23.5,    0.0,  0.0   ],$
[24.5, 5555.0, 17.119 ],$
[25.5,    0.0,  0.0   ],$
[26.5, 5548.6, 17.444 ],$
[27.5,    0.0,  0.0   ],$
[28.5, 5541.9, 17.614 ],$
[29.5,    0.0,  0.0   ],$
[30.5, 5534.7, 17.657 ],$
[31.5,    0.0,  0.0   ],$
[32.5, 5527.2, 17.596 ],$
[33.5,    0.0,  0.0   ],$
[34.5, 5519.2, 17.453 ],$
[35.5,    0.0,  0.0  ]]

P32 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5, 5600.1, 2.166 ],$
[ 2.5,    0.0, 0.0   ],$
[ 3.5, 5600.0, 1.784 ],$
[ 4.5,    0.0, 0.0   ],$
[ 5.5, 5599.6, 2.628 ],$
[ 6.5,    0.0, 0.0   ],$
[ 7.5, 5598.8, 3.646 ],$
[ 8.5,    0.0, 0.0   ],$
[ 9.5, 5597.6, 4.721 ],$
[10.5,    0.0, 0.0   ],$
[11.5, 5595.9, 5.771 ],$
[12.5,    0.0, 0.0   ],$
[13.5, 5593.7, 6.737 ],$
[14.5,    0.0, 0.0   ],$
[15.5, 5591.2, 7.581 ],$
[16.5,    0.0, 0.0   ],$
[17.5, 5588.2, 8.285 ],$
[18.5,    0.0, 0.0   ],$
[19.5, 5584.7, 8.846 ],$
[20.5,    0.0, 0.0   ],$
[21.5, 5580.9, 9.271 ],$
[22.5,    0.0, 0.0   ],$
[23.5, 5576.6, 9.574 ],$
[24.5,    0.0, 0.0   ],$
[25.5, 5571.9, 9.770 ],$
[26.5,    0.0, 0.0   ],$
[27.5, 5566.8, 9.875 ],$
[28.5,    0.0, 0.0   ],$
[29.5, 5561.3, 9.907 ],$
[30.5,    0.0, 0.0   ],$
[31.5, 5555.4, 9.878 ],$
[32.5,    0.0, 0.0   ],$
[33.5, 5549.1, 9.800 ],$
[34.5,    0.0, 0.0   ],$
[35.5, 5542.4, 9.685]]


Q42 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5, 5596.1, 1.177 ],$
[ 2.5,    0.0, 0.0   ],$
[ 3.5, 5593.1, 2.345 ],$
[ 4.5,    0.0, 0.0   ],$
[ 5.5, 5589.6, 3.068 ],$
[ 6.5,    0.0, 0.0   ],$
[ 7.5, 5585.6, 3.477 ],$
[ 8.5,    0.0, 0.0   ],$
[ 9.5, 5581.2, 3.638 ],$
[10.5,    0.0, 0.0   ],$
[11.5, 5576.4, 3.614 ],$
[12.5,    0.0, 0.0   ],$
[13.5, 5571.2, 3.465 ],$
[14.5,    0.0, 0.0   ],$
[15.5, 5565.5, 3.238 ],$
[16.5,    0.0, 0.0   ],$
[17.5, 5559.4, 2.970 ],$
[18.5,    0.0, 0.0   ],$
[19.5, 5553.0, 2.689 ],$
[20.5,    0.0, 0.0   ],$
[21.5, 5546.1, 2.412 ],$
[22.5,    0.0, 0.0   ],$
[23.5, 5538.8, 2.150 ],$
[24.5,    0.0, 0.0   ],$
[25.5, 5531.1, 1.908 ],$
[26.5,    0.0, 0.0   ],$
[27.5, 5523.0, 1.689 ],$
[28.5,    0.0, 0.0   ],$
[29.5, 5514.5, 1.493 ],$
[30.5,    0.0, 0.0   ],$
[31.5, 5505.7, 1.319 ],$
[32.5,    0.0, 0.0   ],$
[33.5, 5496.5, 1.166 ],$
[34.5,    0.0, 0.0   ],$
[35.5, 5486.9, 1.031]]


P42 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5,    0.0, 0.0   ],$
[ 2.5, 5597.8, 1.887 ],$
[ 3.5,    0.0, 0.0   ],$
[ 4.5, 5596.1, 2.182 ],$
[ 5.5,    0.0, 0.0   ],$
[ 6.5, 5593.9, 2.379 ],$
[ 7.5,    0.0, 0.0   ],$
[ 8.5, 5591.3, 2.450 ],$
[ 9.5,    0.0, 0.0   ],$
[10.5, 5588.3, 2.416 ],$
[11.5,    0.0, 0.0   ],$
[12.5, 5584.8, 2.307 ],$
[13.5,    0.0, 0.0   ],$
[14.5, 5580.9, 2.150 ],$
[15.5,    0.0, 0.0   ],$
[16.5, 5576.5, 1.968 ],$
[17.5,    0.0, 0.0   ],$
[18.5, 5571.8, 1.779 ],$
[19.5,    0.0, 0.0   ],$
[20.5, 5566.6, 1.593 ],$
[21.5,    0.0, 0.0   ],$
[22.5, 5561.0, 1.417 ],$
[23.5,    0.0, 0.0   ],$
[24.5, 5554.9, 1.255 ],$
[25.5,    0.0, 0.0   ],$
[26.5, 5548.5, 1.109 ],$
[27.5,    0.0, 0.0   ],$
[28.5, 5541.7, 0.979 ],$
[29.5,    0.0, 0.0   ],$
[30.5, 5534.5, 0.863 ],$
[31.5,    0.0, 0.0   ],$
[32.5, 5526.8, 0.761 ],$
[33.5,    0.0, 0.0   ],$
[34.5, 5518.8, 0.672 ],$
[35.5,    0.0, 0.0  ]]

R13 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5, 5615.8,  .810 ],$
[ 2.5,    0.0, 0.0   ],$
[ 3.5, 5616.1, 1.446 ],$
[ 4.5,    0.0, 0.0   ],$
[ 5.5, 5616.0, 1.891 ],$
[ 6.5,    0.0, 0.0   ],$
[ 7.5, 5615.6, 2.160 ],$
[ 8.5,    0.0, 0.0   ],$
[ 9.5, 5614.9, 2.283 ],$
[10.5,    0.0, 0.0   ],$
[11.5, 5613.8, 2.291 ],$
[12.5,    0.0, 0.0   ],$
[13.5, 5612.4, 2.218 ],$
[14.5,    0.0, 0.0   ],$
[15.5, 5610.6, 2.093 ],$
[16.5,    0.0, 0.0   ],$
[17.5, 5608.5, 1.939 ],$
[18.5,    0.0, 0.0   ],$
[19.5, 5606.0, 1.773 ],$
[20.5,    0.0, 0.0   ],$
[21.5, 5603.1, 1.605 ],$
[22.5,    0.0, 0.0   ],$
[23.5, 5599.9, 1.445 ],$
[24.5,    0.0, 0.0   ],$
[25.5, 5596.4, 1.295 ],$
[26.5,    0.0, 0.0   ],$
[27.5, 5592.5, 1.157 ],$
[28.5,    0.0, 0.0   ],$
[29.5, 5588.2, 1.033 ],$
[30.5,    0.0, 0.0   ],$
[31.5, 5583.5,  .922 ],$
[32.5,    0.0, 0.0   ],$
[33.5, 5578.5,  .823 ],$
[34.5,    0.0, 0.0   ],$
[35.5, 5573.2,  .735 ]]

Q13 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5,    0.0, 0.0   ],$
[ 2.5, 5617.6, 1.782 ],$
[ 3.5,    0.0, 0.0   ],$
[ 4.5, 5619.3, 2.749 ],$
[ 5.5,    0.0, 0.0   ],$
[ 6.5, 5620.6, 3.394 ],$
[ 7.5,    0.0, 0.0   ],$
[ 8.5, 5621.6, 3.753 ],$
[ 9.5,    0.0, 0.0   ],$
[10.5, 5622.3, 3.881 ],$
[11.5,    0.0, 0.0   ],$
[12.5, 5622.6, 3.835 ],$
[13.5,    0.0, 0.0   ],$
[14.5, 5622.6, 3.673 ],$
[15.5,    0.0, 0.0   ],$
[16.5, 5622.2, 3.438 ],$
[17.5,    0.0, 0.0   ],$
[18.5, 5621.4, 3.166 ],$
[19.5,    0.0, 0.0   ],$
[20.5, 5620.3, 2.883 ],$
[21.5,    0.0, 0.0   ],$
[22.5, 5618.8, 2.604 ],$
[23.5,    0.0, 0.0   ],$
[24.5, 5617.0, 2.339 ],$
[25.5,    0.0, 0.0   ],$
[26.5, 5614.8, 2.094 ],$
[27.5,    0.0, 0.0   ],$
[28.5, 5612.2, 1.871 ],$
[29.5,    0.0, 0.0   ],$
[30.5, 5609.2, 1.670 ],$
[31.5,    0.0, 0.0   ],$
[32.5, 5605.9, 1.491 ],$
[33.5,    0.0, 0.0   ],$
[34.5, 5602.2, 1.331 ],$
[35.5,    0.0, 0.0  ]]

R23 = [ $
[  .5, 5614.6,  .106 ],$
[ 1.5,    0.0, 0.0   ],$
[ 2.5, 5613.5,  .585 ],$
[ 3.5,    0.0, 0.0   ],$
[ 4.5, 5612.0, 1.313 ],$
[ 5.5,    0.0, 0.0   ],$
[ 6.5, 5610.2, 2.229 ],$
[ 7.5,    0.0, 0.0   ],$
[ 8.5, 5608.1, 3.255 ],$
[ 9.5,    0.0, 0.0   ],$
[10.5, 5605.6, 4.310 ],$
[11.5,    0.0, 0.0   ],$
[12.5, 5602.7, 5.326 ],$
[13.5,    0.0, 0.0   ],$
[14.5, 5599.5, 6.252 ],$
[15.5,    0.0, 0.0   ],$
[16.5, 5596.0, 7.055 ],$
[17.5,    0.0, 0.0   ],$
[18.5, 5592.1, 7.724 ],$
[19.5,    0.0, 0.0   ],$
[20.5, 5587.9, 8.257 ],$
[21.5,    0.0, 0.0   ],$
[22.5, 5583.3, 8.663 ],$
[23.5,    0.0, 0.0   ],$
[24.5, 5578.3, 8.956 ],$
[25.5,    0.0, 0.0   ],$
[26.5, 5573.0, 9.150 ],$
[27.5,    0.0, 0.0   ],$
[28.5, 5567.4, 9.260 ],$
[29.5,    0.0, 0.0   ],$
[30.5, 5561.3, 9.301 ],$
[31.5,    0.0, 0.0   ],$
[32.5, 5555.0, 9.286 ],$
[33.5,    0.0, 0.0   ],$
[34.5, 5548.2, 9.226 ],$
[35.5,    0.0, 0.0  ]]

Q23 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5, 5615.6,   .376],$
[ 2.5,    0.0, 0.0   ],$
[ 3.5, 5615.9,  1.477],$
[ 4.5,    0.0, 0.0   ],$
[ 5.5, 5615.9,  3.003],$
[ 6.5,    0.0, 0.0   ],$
[ 7.5, 5615.5,  4.820],$
[ 8.5,    0.0, 0.0   ],$
[ 9.5, 5614.7,  6.777],$
[10.5,    0.0, 0.0   ],$
[11.5, 5613.6,  8.726],$
[12.5,    0.0, 0.0   ],$
[13.5, 5612.2, 10.553],$
[14.5,    0.0, 0.0   ],$
[15.5, 5610.4, 12.179],$
[16.5,    0.0, 0.0   ],$
[17.5, 5608.2, 13.562],$
[18.5,    0.0, 0.0   ],$
[19.5, 5605.7, 14.691],$
[20.5,    0.0, 0.0   ],$
[21.5, 5602.8, 15.572],$
[22.5,    0.0, 0.0   ],$
[23.5, 5599.6, 16.227],$
[24.5,    0.0, 0.0   ],$
[25.5, 5595.9, 16.683],$
[26.5,    0.0, 0.0   ],$
[27.5, 5591.9, 16.968],$
[28.5,    0.0, 0.0   ],$
[29.5, 5587.6, 17.110],$
[30.5,    0.0, 0.0   ],$
[31.5, 5582.8, 17.135],$
[32.5,    0.0, 0.0   ],$
[33.5, 5577.7, 17.065],$
[34.5,    0.0, 0.0   ],$
[35.5, 5572.3, 16.920]]

R33 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5, 5611.8,  .288 ],$
[ 2.5,    0.0, 0.0   ],$
[ 3.5, 5608.9,  .224 ],$
[ 4.5,    0.0, 0.0   ],$
[ 5.5, 5605.7,  .093 ],$
[ 6.5,    0.0, 0.0   ],$
[ 7.5, 5602.1,  .004 ],$
[ 8.5,    0.0, 0.0   ],$
[ 9.5, 5598.2,  .052 ],$
[10.5,    0.0, 0.0   ],$
[11.5, 5594.0,  .311 ],$
[12.5,    0.0, 0.0   ],$
[13.5, 5589.4,  .823 ],$
[14.5,    0.0, 0.0   ],$
[15.5, 5584.5, 1.606 ],$
[16.5,    0.0, 0.0   ],$
[17.5, 5579.2, 2.655 ],$
[18.5,    0.0, 0.0   ],$
[19.5, 5573.6, 3.949 ],$
[20.5,    0.0, 0.0   ],$
[21.5, 5567.6, 5.460 ],$
[22.5,    0.0, 0.0   ],$
[23.5, 5561.3, 7.158 ],$
[24.5,    0.0, 0.0   ],$
[25.5, 5554.6, 9.011 ],$
[26.5,    0.0, 0.0   ],$
[27.5, 5547.6,10.991 ],$
[28.5,    0.0, 0.0   ],$
[29.5, 5540.2,13.072 ],$
[30.5,    0.0, 0.0   ],$
[31.5, 5532.5,15.233 ],$
[32.5,    0.0, 0.0   ],$
[33.5, 5524.4,17.455 ],$
[34.5,    0.0, 0.0   ],$
[35.5, 5516.0,19.723 ]]

Q33 = [ $
[  .5, 5614.8, 1.331 ],$
[ 1.5,    0.0, 0.0   ],$
[ 2.5, 5613.5,  .594 ],$
[ 3.5,    0.0, 0.0   ],$
[ 4.5, 5612.0,  .343 ],$
[ 5.5,    0.0, 0.0   ],$
[ 6.5, 5610.2,  .084 ],$
[ 7.5,    0.0, 0.0   ],$
[ 8.5, 5608.1,  .008 ],$
[ 9.5,    0.0, 0.0   ],$
[10.5, 5605.6,  .286 ],$
[11.5,    0.0, 0.0   ],$
[12.5, 5602.7, 1.039 ],$
[13.5,    0.0, 0.0   ],$
[14.5, 5599.5, 2.327 ],$
[15.5,    0.0, 0.0   ],$
[16.5, 5595.9, 4.162 ],$
[17.5,    0.0, 0.0   ],$
[18.5, 5592.0, 6.518 ],$
[19.5,    0.0, 0.0   ],$
[20.5, 5587.7, 9.344 ],$
[21.5,    0.0, 0.0   ],$
[22.5, 5583.1,12.580 ],$
[23.5,    0.0, 0.0   ],$
[24.5, 5578.1,16.161 ],$
[25.5,    0.0, 0.0   ],$
[26.5, 5572.7,20.026 ],$
[27.5,    0.0, 0.0   ],$
[28.5, 5567.0,24.121 ],$
[29.5,    0.0, 0.0   ],$
[30.5, 5560.9,28.397 ],$
[31.5,    0.0, 0.0   ],$
[32.5, 5554.4,32.814 ],$
[33.5,    0.0, 0.0   ],$
[34.5, 5547.6,37.3385],$
[35.5,    0.0, 0.0   ]]

R43 = [ $
[  .5, 5610.8,  .562 ],$
[ 1.5,    0.0, 0.0   ],$
[ 2.5, 5606.5, 1.423 ],$
[ 3.5,    0.0, 0.0   ],$
[ 4.5, 5601.9, 2.287 ],$
[ 5.5,    0.0, 0.0   ],$
[ 6.5, 5597.0, 3.138 ],$
[ 7.5,    0.0, 0.0   ],$
[ 8.5, 5591.7, 3.939 ],$
[ 9.5,    0.0, 0.0   ],$
[10.5, 5586.1, 4.656 ],$
[11.5,    0.0, 0.0   ],$
[12.5, 5580.1, 5.270 ],$
[13.5,    0.0, 0.0   ],$
[14.5, 5573.8, 5.772 ],$
[15.5,    0.0, 0.0   ],$
[16.5, 5567.1, 6.163 ],$
[17.5,    0.0, 0.0   ],$
[18.5, 5560.2, 6.451 ],$
[19.5,    0.0, 0.0   ],$
[20.5, 5552.8, 6.648 ],$
[21.5,    0.0, 0.0   ],$
[22.5, 5545.2, 6.767 ],$
[23.5,    0.0, 0.0   ],$
[24.5, 5537.1, 6.820 ],$
[25.5,    0.0, 0.0   ],$
[26.5, 5528.8, 6.821 ],$
[27.5,    0.0, 0.0   ],$
[28.5, 5520.1, 6.780 ],$
[29.5,    0.0, 0.0   ],$
[30.5, 5511.0, 6.706 ],$
[31.5,    0.0, 0.0   ],$
[32.5, 5501.6, 6.607 ],$
[33.5,    0.0, 0.0   ],$
[34.5, 5491.9, 6.489 ],$
[35.5,    0.0, 0.0   ]]

Q43 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5, 5611.9, 2.029 ],$
[ 2.5,    0.0, 0.0   ],$
[ 3.5, 5609.0, 3.836 ],$
[ 4.5,    0.0, 0.0   ],$
[ 5.5, 5605.8, 5.731 ],$
[ 6.5,    0.0, 0.0   ],$
[ 7.5, 5602.2, 7.577 ],$
[ 8.5,    0.0, 0.0   ],$
[ 9.5, 5598.3, 9.280 ],$
[10.5,    0.0, 0.0   ],$
[11.5, 5594.1,10.777 ],$
[12.5,    0.0, 0.0   ],$
[13.5, 5589.5,12.033 ],$
[14.5,    0.0, 0.0   ],$
[15.5, 5584.5,13.039 ],$
[16.5,    0.0, 0.0   ],$
[17.5, 5579.2,13.806 ],$
[18.5,    0.0, 0.0   ],$
[19.5, 5573.6,14.355 ],$
[20.5,    0.0, 0.0   ],$
[21.5, 5567.6,14.715 ],$
[22.5,    0.0, 0.0   ],$
[23.5, 5561.2,14.915 ],$
[24.5,    0.0, 0.0   ],$
[25.5, 5554.5,14.982 ],$
[26.5,    0.0, 0.0   ],$
[27.5, 5547.4,14.944 ],$
[28.5,    0.0, 0.0   ],$
[29.5, 5539.9,14.823 ],$
[30.5,    0.0, 0.0   ],$
[31.5, 5532.1,14.637 ],$
[32.5,    0.0, 0.0   ],$
[33.5, 5524.0,14.402 ],$
[34.5,    0.0, 0.0   ],$
[35.5, 5515.4,14.131]]

P13 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5,    0.0, 0.0   ],$
[ 2.5,    0.0, 0.0   ],$
[ 3.5, 5620.0,	.792 ],$
[ 4.5,    0.0, 0.0   ],$
[ 5.5, 5623.2, 1.188 ],$
[ 6.5,    0.0, 0.0   ],$
[ 7.5, 5625.9, 1.433 ],$
[ 8.5,    0.0, 0.0   ],$
[ 9.5, 5628.4, 1.559 ],$
[10.5,    0.0, 0.0   ],$
[11.5, 5630.4, 1.592 ],$
[12.5,    0.0, 0.0   ],$
[13.5, 5632.2, 1.559 ],$
[14.5,    0.0, 0.0   ],$
[15.5, 5633.5, 1.483 ],$
[16.5,    0.0, 0.0   ],$
[17.5, 5634.5, 1.381 ],$
[18.5,    0.0, 0.0   ],$
[19.5, 5635.2, 1.268 ],$
[20.5,    0.0, 0.0   ],$
[21.5, 5635.4, 1.152 ],$
[22.5,    0.0, 0.0   ],$
[23.5, 5635.3, 1.039 ],$
[24.5,    0.0, 0.0   ],$
[25.5, 5634.8,	.933 ],$
[26.5,    0.0, 0.0   ],$
[27.5, 5634.0,	.835 ],$
[28.5,    0.0, 0.0   ],$
[29.5, 5632.7,	.746 ],$
[30.5,    0.0, 0.0   ],$
[31.5, 5631.1,	.667 ],$
[32.5,    0.0, 0.0   ],$
[33.5, 5629.1,	.596 ],$
[34.5,    0.0, 0.0   ],$
[35.5, 5626.7,	.532 ]]


P23 = [ $
[  .5,    0.0,  0.0   ],$
[ 1.5,    0.0,  0.0   ],$
[ 2.5, 5617.4,	 .237 ],$
[ 3.5,    0.0,  0.0   ],$
[ 4.5, 5019.1,	 .854 ],$
[ 5.5,    0.0,  0.0   ],$
[ 6.5, 5620.5,	1.637 ],$
[ 7.5,    0.0,  0.0   ],$
[ 8.5, 5621.5,	2.522 ],$
[ 9.5,    0.0,  0.0   ],$
[10.5, 5622.2,	3.439 ],$
[11.5,    0.0,  0.0   ],$
[12.5, 5622.5,	4.326 ],$
[13.5,    0.0,  0.0   ],$
[14.5, 5622.4,	5.136 ],$
[15.5,    0.0,  0.0   ],$
[16.5, 5622.0,	5.842 ],$
[17.5,    0.0,  0.0   ],$
[18.5, 5621.2,	6.431 ],$
[19.5,    0.0,  0.0   ],$
[20.5, 5620.0,	6.902 ],$
[21.5,    0.0,  0.0   ],$
[22.5, 5618.5,	7.263 ],$
[23.5,    0.0,  0.0   ],$
[24.5, 5616.6,	7.523 ],$
[25.5,    0.0,  0.0   ],$
[26.5, 5614.3,	7.698 ],$
[27.5,    0.0,  0.0   ],$
[28.5, 5611.6,	7.799 ],$
[29.5,    0.0,  0.0   ],$
[30.5, 5608.6,	7.839 ],$
[31.5,    0.0,  0.0   ],$
[32.5, 5605.2,	7.830 ],$
[33.5,    0.0,  0.0   ],$
[34.5, 5601.4,	7.781 ],$
[35.5,    0.0,	0.0  ]]

P33 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5, 5615.9,	.497 ],$
[ 2.5,    0.0, 0.0   ],$
[ 3.5, 5616.0,	.225 ],$
[ 4.5,    0.0, 0.0   ],$
[ 5.5, 5615.9,	.095 ],$
[ 6.5,    0.0, 0.0   ],$
[ 7.5, 5615.5,	.005 ],$
[ 8.5,    0.0, 0.0   ],$
[ 9.5, 5614.8,  .049 ],$
[10.5,    0.0, 0.0   ],$
[11.5, 5613.6,	.302 ],$
[12.5,    0.0, 0.0   ],$
[13.5, 5612.2,	.813 ],$
[14.5,    0.0, 0.0   ],$
[15.5, 5610.3, 1.599 ],$
[16.5,    0.0, 0.0   ],$
[17.5, 5608.2, 2.656 ],$
[18.5,    0.0, 0.0   ],$
[19.5, 5605.6, 3.965 ],$
[20.5,    0.0, 0.0   ],$
[21.5, 5602.7, 5.496 ],$
[22.5,    0.0, 0.0   ],$
[23.5, 5599.4, 7.217 ],$
[24.5,    0.0, 0.0   ],$
[25.5, 5595.7, 9.097 ],$
[26.5,    0.0, 0.0   ],$
[27.5, 5591.6,11.105 ],$
[28.5,    0.0, 0.0   ],$
[29.5, 5587.2,13.216 ],$
[30.5,    0.0, 0.0   ],$
[31.5, 5582.4,15.407 ],$
[32.5,    0.0, 0.0   ],$
[33.5, 5577.2,17.660 ],$
[34.5,    0.0, 0.0   ],$
[35.5, 5571.6,19.959]]


P43 = [ $
[  .5,    0.0,  0.0   ],$
[ 1.5,    0.0,  0.0   ],$
[ 2.5, 5613.7,	1.380 ],$
[ 3.5,    0.0,  0.0   ],$
[ 4.5, 5612.2,	2.454 ],$
[ 5.5,    0.0,  0.0   ],$
[ 6.5, 5610.4,	3.519 ],$
[ 7.5,    0.0,  0.0   ],$
[ 8.5, 5608.2,	4.524 ],$
[ 9.5,    0.0,  0.0   ],$
[10.5, 5605.7,	5.428 ],$
[11.5,    0.0,  0.0   ],$
[12.5, 5602.8,	6.204 ],$
[13.5,    0.0,  0.0   ],$
[14.5, 5599.6,	6.841 ],$
[15.5,    0.0,  0.0   ],$
[16.5, 5596.0, 	7.340 ],$
[17.5,    0.0,  0.0   ],$
[18.5, 5592.1,	7.710 ],$
[19.5,    0.0,  0.0   ],$
[20.5, 5587.8,	7.966 ],$
[21.5,    0.0,  0.0   ],$
[22.5, 5583.1,	8.124 ],$
[23.5,    0.0,  0.0   ],$
[24.5, 5578.0,	8.201 ],$
[25.5,    0.0,  0.0   ],$
[26.5, 5572.6,	8.212 ],$
[27.5,    0.0,  0.0   ],$
[28.5, 5566.8,	8.170 ],$
[29.5,    0.0,  0.0   ],$
[30.5, 5560.6,	8.086 ],$
[31.5,    0.0,  0.0   ],$
[32.5, 5554.1,	7.972 ],$
[33.5,    0.0,  0.0   ],$
[34.5, 5547.1,  7.833 ],$
[35.5,    0.0,  0.0  ]]

R14 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5, 5631.5,	.758 ],$
[ 2.5,    0.0, 0.0   ],$
[ 3.5, 5632.0,	.707 ],$
[ 4.5,    0.0, 0.0   ],$
[ 5.5, 5632.2,  .644 ],$
[ 6.5,    0.0, 0.0   ],$
[ 7.5, 5632.2,  .562 ],$
[ 8.5,    0.0, 0.0   ],$
[ 9.5, 5632.0,	.474 ],$
[10.5,    0.0, 0.0   ],$
[11.5, 5631.6,	.390 ],$
[12.5,    0.0, 0.0   ],$
[13.5, 5630.8,	.316 ],$
[14.5,    0.0, 0.0   ],$
[15.5, 5629.8,	.252 ],$
[16.5,    0.0, 0.0   ],$
[17.5, 5628.6,	.200 ],$
[18.5,    0.0, 0.0   ],$
[19.5, 5627.0,	.158 ],$
[20.5,    0.0, 0.0   ],$
[21.5, 5625.1,	.125 ],$
[22.5,    0.0, 0.0   ],$
[23.5, 5622.9,	.099 ],$
[24.5,    0.0, 0.0   ],$
[25.5, 5620.3,	.079 ],$
[26.5,    0.0, 0.0   ],$
[27.5, 5617.5,	.063 ],$
[28.5,    0.0, 0.0   ],$
[29.5, 5614.3,	.050 ],$
[30.5,    0.0, 0.0   ],$
[31.5, 5610.8,  .040 ],$
[32.5,    0.0, 0.0   ],$
[33.5, 5606.9,	.033 ],$
[34.5,    0.0, 0.0   ],$
[35.5, 5602.7,	.027]]

R24 = [ $
[  .5, 5630.2,	1.742 ],$
[ 1.5,    0.0,  0.0   ],$
[ 2.5, 5629.2,	2.122 ],$
[ 3.5,    0.0,  0.0   ],$
[ 4.5, 5628.0,  2.453 ],$
[ 5.5,    0.0,  0.0   ],$
[ 6.5, 5626.6,	2.644 ],$
[ 7.5,    0.0,  0.0   ],$
[ 8.5, 5624.9,	2.703 ],$
[ 9.5,    0.0,  0.0   ],$
[10.5, 5622.9,	2.656 ],$
[11.5,    0.0,  0.0   ],$
[12.5, 5620.8,	2.534 ],$
[13.5,    0.0,  0.0   ],$
[14.5, 5618.3,	2.365 ],$
[15.5,    0.0,  0.0   ],$
[16.5, 5615.6,	2.172 ],$
[17.5,    0.0,  0.0   ],$
[18.5, 5612.5,	1.972 ],$
[19.5,    0.0,  0.0   ],$
[20.5, 5609.2,	1.775 ],$
[21.5,    0.0,  0.0   ],$
[22.5, 5605.5,	1.590 ],$
[23.5,    0.0,  0.0   ],$
[24.5, 5601.6,	1.418 ],$
[25.5,    0.0,  0.0   ],$
[26.5, 5597.3,	1.263 ],$
[27.5,    0.0,  0.0   ],$
[28.5, 5592.7,	1.123 ],$
[29.5,    0.0,  0.0   ],$
[30.5, 5587.8,	 .999 ],$
[31.5,    0.0,  0.0   ],$
[32.5, 5582.5,	 .889 ],$
[33.5,    0.0,  0.0   ],$
[34.5, 5576.9,	 .792 ],$
[35.5,    0.0,	0.0  ]]

R34 = [ $
[  .5,    0.0,  0.0   ],$
[ 1.5, 5627.4,	1.634 ],$
[ 2.5,    0.0,  0.0   ],$
[ 3.5, 5624.7,	2.592 ],$
[ 4.5,    0.0,  0.0   ],$
[ 5.5, 5621.8,	3.565 ],$
[ 6.5,    0.0,  0.0   ],$
[ 7.5, 5618.6,	4.504 ],$
[ 8.5,    0.0,  0.0   ],$
[ 9.5, 5615.2,	5.366 ],$
[10.5,    0.0,  0.0   ],$
[11.5, 5611.6,	6.118 ],$
[12.5,    0.0,  0.0   ],$
[13.5, 5607.7,	6.745 ],$
[14.5,    0.0,  0.0   ],$
[15.5, 5603.5,	7.243 ],$
[16.5,    0.0,  0.0   ],$
[17.5, 5599.0,	7.618 ],$
[18.5,    0.0,  0.0   ],$
[19.5, 5594.3,	7.882 ],$
[20.5,    0.0,  0.0   ],$
[21.5, 5589.2,	8.048 ],$
[22.5,    0.0,  0.0   ],$
[23.5, 5583.9,	8.134 ],$
[24.5,    0.0,  0.0   ],$
[25.5, 5578.2,	8.152 ],$
[26.5,    0.0,  0.0   ],$
[27.5, 5572.2,	8.118 ],$
[28.5,    0.0,  0.0   ],$
[29.5, 5565.8,	8.041 ],$
[30.5,    0.0,  0.0   ],$
[31.5, 5559.2,	7.932 ],$
[32.5,    0.0,  0.0   ],$
[33.5, 5552.2,	7.799 ],$
[34.5,    0.0,  0.0   ],$
[35.5, 5544.9,	7.647]]

R44 = [ $
[  .5, 5626.5,  .256 ],$
[ 1.5,    0.0, 0.0   ],$
[ 2.5, 5622.2,  .722 ],$
[ 3.5,    0.0, 0.0   ],$
[ 4.5, 5617.8, 1.347 ],$
[ 5.5,    0.0, 0.0   ],$
[ 6.5, 5613.2, 2.165 ],$
[ 7.5,    0.0, 0.0   ],$
[ 8.5, 5608.4, 3.186 ],$
[ 9.5,    0.0, 0.0   ],$
[10.5, 5603.3, 4.404 ],$
[11.5,    0.0, 0.0   ],$
[12.5, 5598.0, 5.807 ],$
[13.5,    0.0, 0.0   ],$
[14.5, 5592.4, 7.373 ],$
[15.5,    0.0, 0.0   ],$
[16.5, 5586.5, 9.079 ],$
[17.5,    0.0, 0.0   ],$
[18.5, 5580.3,10.903 ],$
[19.5,    0.0, 0.0   ],$
[20.5, 5573.9,12.820 ],$
[21.5,    0.0, 0.0   ],$
[22.5, 5567.1,14.813 ],$
[23.5,    0.0, 0.0   ],$
[24.5, 5560.1,16.863 ],$
[25.5,    0.0, 0.0   ],$
[26.5, 5552.7,18.958 ],$
[27.5,    0.0, 0.0   ],$
[28.5, 5545.0,21.085 ],$
[29.5,    0.0, 0.0   ],$
[30.5, 5537.0,23.235 ],$
[31.5,    0.0, 0.0   ],$
[32.5, 5528.6,25.401 ],$
[33.5,    0.0, 0.0   ],$
[34.5, 5519.9,27.578 ],$
[35.5,    0.0, 0.0  ]]

Q14 = [ $
[  .5,     0.0, 0.0   ],$
[ 1.5,     0.0, 0.0   ],$
[ 2.5,  5633.3,  .717 ],$
[ 3.5,     0.0, 0.0   ],$
[ 4.5,  5635.3,  .818 ],$
[ 5.5,     0.0, 0.0   ],$
[ 6.5,  5637.0,  .799 ],$
[ 7.5,     0.0, 0.0   ],$
[ 8.5,  5638.5,  .722 ],$
[ 9.5,     0.0, 0.0   ],$
[10.5,  5639.8,  .621 ],$
[11.5,     0.0, 0.0   ],$
[12.5,  5640.8,  .518 ],$
[13.5,     0.0, 0.0   ],$
[14.5,  5641.5,  .423 ],$
[15.5,     0.0, 0.0   ],$
[16.5,  5641.9,  .341 ],$
[17.5,     0.0, 0.0   ],$
[18.5,  5642.1,  .273 ],$
[19.5,     0.0, 0.0   ],$
[20.5,  5641.9,  .217 ],$
[21.5,     0.0, 0.0   ],$
[22.5,  5641.4,  .173 ],$
[23.5,     0.0, 0.0   ],$
[24.5,  5640.6,  .138 ],$
[25.5,     0.0, 0.0   ],$
[26.5,  5639.4,  .110 ],$
[27.5,     0.0, 0.0   ],$
[28.5,  5637.9,  .088 ],$
[29.5,     0.0, 0.0   ],$
[30.5,  5636.1,  .071 ],$
[31.5,     0.0, 0.0   ],$
[32.5,  5633.9,  .057 ],$
[33.5,     0.0, 0.0   ],$
[34.5,  5631.4,  .047 ],$
[35.5,     0.0, 0.0  ]]


Q24 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5, 5631.3, 1.319 ],$
[ 2.5,    0.0, 0.0   ],$
[ 3.5, 5631.8, 2.596 ],$
[ 4.5,    0.0, 0.0   ],$
[ 5.5, 5632.1, 3.390 ],$
[ 6.5,    0.0, 0.0   ],$
[ 7.5, 5632.1, 3.847 ],$
[ 8.5,    0.0, 0.0   ],$
[ 9.5, 5631.9, 4.037 ],$
[10.5,    0.0, 0.0   ],$
[11.5, 5631.4, 4.028 ],$
[12.5,    0.0, 0.0   ],$
[13.5, 5630.6, 3.880 ],$
[14.5,    0.0, 0.0   ],$
[15.5, 5629.6, 3.647 ],$
[16.5,    0.0, 0.0   ],$
[17.5, 5628.3, 3.367 ],$
[18.5,    0.0, 0.0   ],$
[19.5, 5626.7, 3.069 ],$
[20.5,    0.0, 0.0   ],$
[21.5, 5624.7, 2.774 ],$
[22.5,    0.0, 0.0   ],$
[23.5, 5611.5, 2.492 ],$
[24.5,    0.0, 0.0   ],$
[25.5, 5619.9, 2.230 ],$
[26.5,    0.0, 0.0   ],$
[27.5, 5616.9, 1.991 ],$
[28.5,    0.0, 0.0   ],$
[29.5, 5613.7, 1.776 ],$
[30.5,    0.0, 0.0   ],$
[31.5, 5610.0, 1.584 ],$
[32.5,    0.0, 0.0   ],$
[33.5, 5606.1, 1.414 ],$
[34.5,    0.0, 0.0   ],$
[35.5, 5601.8, 1.263 ]]


Q34 = [ $
[  .5, 5630.5,   .002 ],$
[ 1.5,    0.0,  0.0   ],$
[ 2.5, 5629.3,  2.025 ],$
[ 3.5,    0.0,  0.0   ],$
[ 4.5, 5628.0,  3.991 ],$
[ 5.5,    0.0,  0.0   ],$
[ 6.5, 5626.6,  5.926 ],$
[ 7.5,    0.0,  0.0   ],$
[ 8.5, 5624.9,  7.757 ],$
[ 9.5,    0.0,  0.0   ],$
[10.5, 5622.9,  9.410 ],$
[11.5,    0.0,  0.0   ],$
[12.5, 5620.7, 10.835 ],$
[13.5,    0.0,  0.0   ],$
[14.5, 5618.2, 12.012 ],$
[15.5,    0.0,  0.0   ],$
[16.5, 5615.5, 12.940 ],$
[17.5,    0.0,  0.0   ],$
[18.5, 5612.4, 13.635 ],$
[19.5,    0.0,  0.0   ],$
[20.5, 5609.0, 14.124 ],$
[21.5,    0.0,  0.0   ],$
[22.5, 5605.4, 14.436 ],$
[23.5,    0.0,  0.0   ],$
[24.5, 5601.3, 14.598 ],$
[25.5,    0.0,  0.0   ],$
[26.5, 5597.0, 14.640 ],$
[27.5,    0.0,  0.0   ],$
[28.5, 5592.3, 14.584 ],$
[29.5,    0.0,  0.0   ],$
[30.5, 5587.3, 14.453 ],$
[31.5,    0.0,  0.0   ],$
[32.5, 5581.9, 14.262 ],$
[33.5,    0.0,  0.0   ],$
[34.5, 5576.2, 14.028 ],$
[35.5,    0.0,  0.0  ]]

Q44 = [ $
[   5,    0.0,  0.0   ],$
[ 1.5, 5627.6,   .285 ],$
[ 2.5,    0.0,  0.0   ],$
[ 3.5, 5624.8,  1.239 ],$
[ 4.5,    0.0,  0.0   ],$
[ 5.5, 5621.9,  2.589 ],$
[ 6.5,    0.0,  0.0   ],$
[ 7.5, 5618.7,  4.360 ],$
[ 8.5,    0.0,  0.0   ],$
[ 9.5, 5615.3,  6.553 ],$
[10.5,    0.0,  0.0   ],$
[11.5, 5611.7,  9.151 ],$
[12.5,    0.0,  0.0   ],$
[13.5, 5607.8, 12.116 ],$
[14.5,    0.0,  0.0   ],$
[15.5, 5603.6, 15.402 ],$
[16.5,    0.0,  0.0   ],$
[17.5, 5599.1, 18.958 ],$
[18.5,    0.0,  0.0   ],$
[19.5, 5594.3, 22.737 ],$
[20.5,    0.0,  0.0   ],$
[21.5, 5589.2, 26.692 ],$
[22.5,    0.0,  0.0   ],$
[23.5, 5583.8, 30.785 ],$
[24.5,    0.0,  0.0   ],$
[25.5, 5578.1, 34.983 ],$
[26.5,    0.0,  0.0   ],$
[27.5, 5572.0, 39.259 ],$
[28.5,    0.0,  0.0   ],$
[29.5, 5565.6, 43.592 ],$
[30.5,    0.0,  0.0   ],$
[31.5, 5558.8, 47.964 ],$
[32.5,    0.0,  0.0   ],$
[33.5, 5551.8, 52.362 ],$
[34.5,    0.0,  0.0   ],$
[35.5, 5544.3, 56.776]]


P14 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5,    0.0, 0.0   ],$
[ 2.5,    0.0, 0.0   ],$
[ 3.5, 5635.9,	.171 ],$
[ 4.5,    0.0, 0.0   ],$
[ 5.5, 5639.4,	.231 ],$
[ 6.5,    0.0, 0.0   ],$
[ 7.5, 5642.6,	.241 ],$
[ 8.5,    0.0, 0.0   ],$
[ 9.5, 5645.6,	.226 ],$
[10.5,    0.0, 0.0   ],$
[11.5, 5648.3,	.199 ],$
[12.5,    0.0, 0.0   ],$
[13.5, 5650.7,	.168 ],$
[14.5,    0.0, 0.0   ],$
[15.5, 5652.9,	.139 ],$
[16.5,    0.0, 0.0   ],$
[17.5, 5654.8,	.113 ],$
[18.5,    0.0, 0.0   ],$
[19.5, 5656.1,	.091 ],$
[20.5,    0.0, 0.0   ],$
[21.5, 5657.6,	.073 ],$
[22.5,    0.0, 0.0   ],$
[23.5, 5658.5,	.059 ],$
[24.5,    0.0, 0.0   ],$
[25.5, 5659.1,	.047 ],$
[26.5,    0.0, 0.0   ],$
[27.5, 5659.3,	.038 ],$
[28.5,    0.0, 0.0   ],$
[29.5, 5659.2,	.030 ],$
[30.5,    0.0, 0.0   ],$
[31.5, 5658.8,	.025 ],$
[32.5,    0.0, 0.0   ],$
[33.5, 5658.0,	.020 ],$
[34.5,    0.0, 0.0   ],$
[35.5, 5656.8,	.016]]

P24 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5,    0.0, 0.0   ],$
[ 2.5, 5633.2,  .298 ],$
[ 3.5,    0.0, 0.0   ],$
[ 4.5, 5635.2,  .806 ],$
[ 5.5,    0.0, 0.0   ],$
[ 6.5, 5636.9, 1.163 ],$
[ 7.5,    0.0, 0.0   ],$
[ 8.5, 5638.4, 1.381 ],$
[ 9.5,    0.0, 0.0   ],$
[10.5, 5639.6, 1.486 ],$
[11.5,    0.0, 0.0   ],$
[12.5, 5640.6, 1.505 ],$
[13.5,    0.0, 0.0   ],$
[14.5, 5641.3, 1.464 ],$
[15.5,    0.0, 0.0   ],$
[16.5, 5641.7, 1.386 ],$
[17.5,    0.0, 0.0   ],$
[18.5, 5641.8, 1.287 ],$
[19.5,    0.0, 0.0   ],$
[20.5, 5641.6, 1.179 ],$
[21.5,    0.0, 0.0   ],$
[22.5, 5641.1, 1.070 ],$
[23.5,    0.0, 0.0   ],$
[24.5, 5640.2,	.965 ],$
[25.5,    0.0, 0.0   ],$
[26.5, 5639.0,	.867 ],$
[27.5,    0.0, 0.0   ],$
[28.5, 5637.4,	.777 ],$
[29.5,    0.0, 0.0   ],$
[30.5, 5635.5,  .695 ],$
[31.5,    0.0, 0.0   ],$
[32.5, 5633.2,	.622 ],$
[33.5,    0.0, 0.0   ],$
[34.5, 5630.6,	.554 ],$
[35.5,    0.0, 0.0  ]]

P34 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5, 5631.6,  .003 ],$
[ 2.5,    0.0, 0.0   ],$
[ 3.5, 5631.9,  .695 ],$
[ 4.5,    0.0, 0.0   ],$
[ 5.5, 5632.1, 1.581 ],$
[ 6.5,    0.0, 0.0   ],$
[ 7.5, 5632.1, 2.486 ],$
[ 8.5,    0.0, 0.0   ],$
[ 9.5, 5631.9, 3.343 ],$
[10.5,    0.0, 0.0   ],$
[11.5, 5631.4, 4.114 ],$
[12.5,    0.0, 0.0   ],$
[13.5, 5630.6, 4.775 ],$
[14.5,    0.0, 0.0   ],$
[15.5, 5629.6, 5.317 ],$
[16.5,    0.0, 0.0   ],$
[17.5, 5628.2, 5.743 ],$
[18.5,    0.0, 0.0   ],$
[19.5, 5626.6, 6.063 ],$
[20.5,    0.0, 0.0   ],$
[21.5, 5624.6, 6.287 ],$
[22.5,    0.0, 0.0   ],$
[23.5, 5622.3, 6.432 ],$
[24.5,    0.0, 0.0   ],$
[25.5, 5619.6, 6.509 ],$
[26.5,    0.0, 0.0   ],$
[27.5, 5616.6, 6.531 ],$
[28.5,    0.0, 0.0   ],$
[29.5, 5613.3, 6.510 ],$
[30.5,    0.0, 0.0   ],$
[31.5, 5609.6, 6.454 ],$
[32.5,    0.0, 0.0   ],$
[33.5, 5605.5, 6.373 ],$
[34.5,    0.0, 0.0   ],$
[35.5, 5601.1, 6.271 ]]

P44 = [ $
[  .5,    0.0, 0.0   ],$
[ 1.5,    0.0, 0.0   ],$
[ 2.5, 5629.4,  .116 ],$
[ 3.5,    0.0, 0.0   ],$
[ 4.5, 5628.2,  .585 ],$
[ 5.5,    0.0, 0.0   ],$
[ 6.5, 5626.7, 1.302 ],$
[ 7.5,    0.0, 0.0   ],$
[ 8.5, 5625.0, 2.251 ],$
[ 9.5,    0.0, 0.0   ],$
[10.5, 5623.1, 3.423 ],$
[11.5,    0.0, 0.0   ],$
[12.5, 5620.8, 4.802 ],$
[13.5,    0.0, 0.0   ],$
[14.5, 5618.3, 6.363 ],$
[15.5,    0.0, 0.0   ],$
[16.5, 5615.6, 8.082 ],$
[17.5,    0.0, 0.0   ],$
[18.5, 5612.5, 9.930 ],$
[19.5,    0.0, 0.0   ],$
[20.5, 5609.1,11.883 ],$
[21.5,    0.0, 0.0   ],$
[22.5, 5605.3,13.919 ],$
[23.5,    0.0, 0.0   ],$
[24.5, 5601.3,16.017 ],$
[25.5,    0.0, 0.0   ],$
[26.5, 5596.9,18.163 ],$
[27.5,    0.0, 0.0   ],$
[28.5, 5592.1,20.343 ],$
[29.5,    0.0, 0.0   ],$
[30.5, 5587.0,22.548 ],$
[31.5,    0.0, 0.0   ],$
[32.5, 5581.6,24.768 ],$
[33.5,    0.0, 0.0   ],$
[34.5, 5575.8,26.999 ],$
[35.5,    0.0, 0.0   ]]

;SAVE
;R11__=R11
;Q11__=Q11
;R21__=R21
;Q21__=Q21
;R31__=R31
;Q31__=Q31
;R41__=R41
;Q41__=Q41
;P11__=P11
;R12__=R12
;P21__=P21
;R22__=R22
;P31__=P31
;R32__=R32
;P41__=P41
;R42__=R42
;Q12__=Q12
;P12__=P12
;Q22__=Q22
;P22__=P22
;Q32__=Q32
;P32__=P32
;Q42__=Q42
;P42__=P42
;R13__=R13
;R23__=R23
;R33__=R33
;R43__=R43
;Q13__=Q13
;Q23__=Q23
;Q33__=Q33
;Q43__=Q43
;P13__=P13
;P23__=P23
;P33__=P33
;P43__=P43
;R14__=R14
;R24__=R24
;R34__=R34
;R44__=R44
;Q14__=Q14
;Q24__=Q24
;Q34__=Q34
;Q44__=Q44
;P14__=P14
;P24__=P24
;P34__=P34
;P44__=P44

;save,R11__,Q11__,R21__,Q21__,R31__,Q31__,R41__,Q41__,P11__,R12__,P21__,R22__,P31__,R32__,P41__,R42__,$
;     Q12__,P12__,Q22__,P22__,Q32__,P32__,Q42__,P42__,R13__,R23__,R33__,R43__,Q13__,Q23__,Q33__,Q43__,$
;     P13__,P23__,P33__,P43__,R14__,R24__,R34__,R44__,Q14__,Q24__,Q34__,Q44__,P14__,P24__,P34__,P44__,$
;     filename='HL_data.idl' ;save digitised data to file


;EXTRAPOLATE TO J>34
if lines gt 34 then begin
R11_=do_interpol(lines,R11)
Q11_=do_interpol(lines,Q11)
R21_=do_interpol(lines,R21)
Q21_=do_interpol(lines,Q21)
R31_=do_interpol(lines,R31)
Q31_=do_interpol(lines,Q31)
R41_=do_interpol(lines,R41)
Q41_=do_interpol(lines,Q41)

P11_=do_interpol(lines,P11)
R12_=do_interpol(lines,R12)
P21_=do_interpol(lines,P21)
R22_=do_interpol(lines,R22)
P31_=do_interpol(lines,P31)
R32_=do_interpol(lines,R32)
P41_=do_interpol(lines,P41)
R42_=do_interpol(lines,R42)

Q12_=do_interpol(lines,Q12)
P12_=do_interpol(lines,P12)
Q22_=do_interpol(lines,Q22)
P22_=do_interpol(lines,P22)
Q32_=do_interpol(lines,Q32)
P32_=do_interpol(lines,P32)
Q42_=do_interpol(lines,Q42)
P42_=do_interpol(lines,P42)

R13_=do_interpol(lines,R13)
R23_=do_interpol(lines,R23)
R33_=do_interpol(lines,R33)
R43_=do_interpol(lines,R43)
Q13_=do_interpol(lines,Q13)
Q23_=do_interpol(lines,Q23)
Q33_=do_interpol(lines,Q33)
Q43_=do_interpol(lines,Q43)

P13_=do_interpol(lines,P13)
P23_=do_interpol(lines,P23)
P33_=do_interpol(lines,P33)
P43_=do_interpol(lines,P43)
R14_=do_interpol(lines,R14)
R24_=do_interpol(lines,R24)
R34_=do_interpol(lines,R34)
R44_=do_interpol(lines,R44)

Q14_=do_interpol(lines,Q14)
Q24_=do_interpol(lines,Q24)
Q34_=do_interpol(lines,Q34)
Q44_=do_interpol(lines,Q44)
P14_=do_interpol(lines,P14)
P24_=do_interpol(lines,P24)
P34_=do_interpol(lines,P34)
P44_=do_interpol(lines,P44)
endif

;FIRST non-zero J
R11_J=min(where(R11(2,*) ne 0.))
Q11_J=min(where(Q11(2,*) ne 0.))
R21_J=min(where(R21(2,*) ne 0.))
Q21_J=min(where(Q21(2,*) ne 0.))
R31_J=min(where(R31(2,*) ne 0.))
Q31_J=min(where(Q31(2,*) ne 0.))
R41_J=min(where(R41(2,*) ne 0.))
Q41_J=min(where(Q41(2,*) ne 0.))
P11_J=min(where(P11(2,*) ne 0.))
R12_J=min(where(R12(2,*) ne 0.))
P21_J=min(where(P21(2,*) ne 0.))
R22_J=min(where(R22(2,*) ne 0.))
P31_J=min(where(P31(2,*) ne 0.))
R32_J=min(where(R32(2,*) ne 0.))
P41_J=min(where(P41(2,*) ne 0.))
R42_J=min(where(R42(2,*) ne 0.))
Q12_J=min(where(Q12(2,*) ne 0.))
P12_J=min(where(P12(2,*) ne 0.))
Q22_J=min(where(Q22(2,*) ne 0.))
P22_J=min(where(P22(2,*) ne 0.))
Q32_J=min(where(Q32(2,*) ne 0.))
P32_J=min(where(P32(2,*) ne 0.))
Q42_J=min(where(Q42(2,*) ne 0.))
P42_J=min(where(P42(2,*) ne 0.))
R13_J=min(where(R13(2,*) ne 0.))
R23_J=min(where(R23(2,*) ne 0.))
R33_J=min(where(R33(2,*) ne 0.))
R43_J=min(where(R43(2,*) ne 0.))
Q13_J=min(where(Q13(2,*) ne 0.))
Q23_J=min(where(Q23(2,*) ne 0.))
Q33_J=min(where(Q33(2,*) ne 0.))
Q43_J=min(where(Q43(2,*) ne 0.))
P13_J=min(where(P13(2,*) ne 0.))
P23_J=min(where(P23(2,*) ne 0.))
P33_J=min(where(P33(2,*) ne 0.))
P43_J=min(where(P43(2,*) ne 0.))
R14_J=min(where(R14(2,*) ne 0.))
R24_J=min(where(R24(2,*) ne 0.))
R34_J=min(where(R34(2,*) ne 0.))
R44_J=min(where(R44(2,*) ne 0.))
Q14_J=min(where(Q14(2,*) ne 0.))
Q24_J=min(where(Q24(2,*) ne 0.))
Q34_J=min(where(Q34(2,*) ne 0.))
Q44_J=min(where(Q44(2,*) ne 0.))
P14_J=min(where(P14(2,*) ne 0.))
P24_J=min(where(P24(2,*) ne 0.))
P34_J=min(where(P34(2,*) ne 0.))
P44_J=min(where(P44(2,*) ne 0.))

as=1
if as eq 0 then begin
;DO A PLOT #########################################
window,xsize=1300, ysize=800
loadct,39
!p.multi=[0,8,6]
hlt=lines ;x-axis
hgt=lines ;y-axis

plot_br, R11, hlt, hgt, lines, R11_J, R11_
plot_br, Q11, hlt, hgt, lines, Q11_J, Q11_
plot_br, R21, hlt, hgt, lines, R21_J, R21_
plot_br, Q21, hlt, hgt, lines, Q21_J, Q21_
plot_br, R31, hlt, hgt, lines, R31_J, R31_
plot_br, Q31, hlt, hgt, lines, Q31_J, Q31_
plot_br, R41, hlt, hgt, lines, R41_J, R41_
plot_br, Q41, hlt, hgt, lines, Q41_J, Q41_

plot_br, P11, hlt, hgt, lines, P11_J, P11_
plot_br, R12, hlt, hgt, lines, R12_J, R12_
plot_br, P21, hlt, hgt, lines, P21_J, P21_
plot_br, R22, hlt, hgt, lines, R22_J, R22_
plot_br, P31, hlt, hgt, lines, P31_J, P31_
plot_br, R32, hlt, hgt, lines, R32_J, R32_
plot_br, P41, hlt, hgt, lines, P41_J, P41_
plot_br, R42, hlt, hgt, lines, R42_J, R42_

plot_br, Q12, hlt, hgt, lines, Q12_J, Q12_
plot_br, P12, hlt, hgt, lines, P12_J, P12_
plot_br, Q22, hlt, hgt, lines, Q22_J, Q22_
plot_br, P22, hlt, hgt, lines, P22_J, P22_
plot_br, Q32, hlt, hgt, lines, Q32_J, Q32_
plot_br, P32, hlt, hgt, lines, P32_J, P32_
plot_br, Q42, hlt, hgt, lines, Q42_J, Q42_
plot_br, P42, hlt, hgt, lines, P42_J, P42_

plot_br, R13, hlt, hgt, lines, R13_J, R13_
plot_br, R23, hlt, hgt, lines, R23_J, R23_
plot_br, R33, hlt, hgt, lines, R33_J, R33_
plot_br, R43, hlt, hgt, lines, R43_J, R43_
plot_br, Q13, hlt, hgt, lines, Q13_J, Q13_
plot_br, Q23, hlt, hgt, lines, Q23_J, Q23_
plot_br, Q33, hlt, hgt, lines, Q33_J, Q33_
plot_br, Q43, hlt, hgt, lines, Q43_J, Q43_

plot_br, P13, hlt, hgt, lines, P13_J, P13_
plot_br, P23, hlt, hgt, lines, P23_J, P23_
plot_br, P33, hlt, hgt, lines, P33_J, P33_
plot_br, P43, hlt, hgt, lines, P43_J, P43_
plot_br, R14, hlt, hgt, lines, R14_J, R14_
plot_br, R24, hlt, hgt, lines, R24_J, R24_
plot_br, R34, hlt, hgt, lines, R34_J, R34_
plot_br, R44, hlt, hgt, lines, R44_J, R44_

plot_br, Q14, hlt, hgt, lines, Q14_J, Q14_
plot_br, Q24, hlt, hgt, lines, Q24_J, Q24_
plot_br, Q34, hlt, hgt, lines, Q34_J, Q34_
plot_br, Q44, hlt, hgt, lines, Q44_J, Q44_
plot_br, P14, hlt, hgt, lines, P14_J, P14_
plot_br, P24, hlt, hgt, lines, P24_J, P24_
plot_br, P34, hlt, hgt, lines, P34_J, P34_
plot_br, P44, hlt, hgt, lines, P44_J, P44_
endif

;NEW ARRAY WITH EXTRAPOLATED POINTS
if lines gt 34. then begin
R11_new=apol(lines, R11_J, R11_)
Q11_new=apol(lines, Q11_J, Q11_)
R21_new=apol(lines, R21_J, R21_)
Q21_new=apol(lines, Q21_J, Q21_)
R31_new=apol(lines, R31_J, R31_)
Q31_new=apol(lines, Q31_J, Q31_)
R41_new=apol(lines, R41_J, R41_)
Q41_new=apol(lines, Q41_J, Q41_)

P11_new=apol(lines, P11_J, P11_)
R12_new=apol(lines, R12_J, R12_)
P21_new=apol(lines, P21_J, P21_)
R22_new=apol(lines, R22_J, R22_)
P31_new=apol(lines, P31_J, P31_)
R32_new=apol(lines, R32_J, R32_)
P41_new=apol(lines, P41_J, P41_)
R42_new=apol(lines, R42_J, R42_)

Q12_new=apol(lines, Q12_J, Q12_)
P12_new=apol(lines, P12_J, P12_)
Q22_new=apol(lines, Q22_J, Q22_)
P22_new=apol(lines, P22_J, P22_)
Q32_new=apol(lines, Q32_J, Q32_)
P32_new=apol(lines, P32_J, P32_)
Q42_new=apol(lines, Q42_J, Q42_)
P42_new=apol(lines, P42_J, P42_)

R13_new=apol(lines, R13_J, R13_)
R23_new=apol(lines, R23_J, R23_)
R33_new=apol(lines, R33_J, R33_)
R43_new=apol(lines, R43_J, R43_)
Q13_new=apol(lines, Q13_J, Q13_)
Q23_new=apol(lines, Q23_J, Q23_)
Q33_new=apol(lines, Q33_J, Q33_)
Q43_new=apol(lines, Q43_J, Q43_)

P13_new=apol(lines, P13_J, P13_)
P23_new=apol(lines, P23_J, P23_)
P33_new=apol(lines, P33_J, P33_)
P43_new=apol(lines, P43_J, P43_)
R14_new=apol(lines, R14_J, R14_)
R24_new=apol(lines, R24_J, R24_)
R34_new=apol(lines, R34_J, R34_)
R44_new=apol(lines, R44_J, R44_)

Q14_new=apol(lines, Q14_J, Q14_)
Q24_new=apol(lines, Q24_J, Q24_)
Q34_new=apol(lines, Q34_J, Q34_)
Q44_new=apol(lines, Q44_J, Q44_)
P14_new=apol(lines, P14_J, P14_)
P24_new=apol(lines, P24_J, P24_)
P34_new=apol(lines, P34_J, P34_)
P44_new=apol(lines, P44_J, P44_)
endif

ar=1
if ar eq 0 then begin ; #########################################
fi=lines-5.
se=lines-1. ;199
print, R11_J, R11_new(fi:se)
print, Q11_J, Q11_new(fi:se)
print, R21_J, R21_new(fi:se)
print, Q21_J, Q21_new(fi:se)
print, R31_J, R31_new(fi:se)
print, Q31_J, Q31_new(fi:se)
print, R41_J, R41_new(fi:se)
print, Q41_J, Q41_new(fi:se)

print, P11_J, P11_new(fi:se)
print, R12_J, R12_new(fi:se)
print, P21_J, P21_new(fi:se)
print, R22_J, R22_new(fi:se)
print, P31_J, P31_new(fi:se)
print, R32_J, R32_new(fi:se)
print, P41_J, P41_new(fi:se)
print, R42_J, R42_new(fi:se)

print, Q12_J, Q12_new(fi:se)
print, P12_J, P12_new(fi:se)
print, Q22_J, Q22_new(fi:se)
print, P22_J, P22_new(fi:se)
print, Q32_J, Q32_new(fi:se)
print, P32_J, P32_new(fi:se)
print, Q42_J, Q42_new(fi:se)
print, P42_J, P42_new(fi:se)

print, R13_J, R13_new(fi:se)
print, R23_J, R23_new(fi:se)
print, R33_J, R33_new(fi:se)
print, R43_J, R43_new(fi:se)
print, Q13_J, Q13_new(fi:se)
print, Q23_J, Q23_new(fi:se)
print, Q33_J, Q33_new(fi:se)
print, Q43_J, Q43_new(fi:se)

print, P13_J, P13_new(fi:se)
print, P23_J, P23_new(fi:se)
print, P33_J, P33_new(fi:se)
print, P43_J, P43_new(fi:se)
print, R14_J, R14_new(fi:se)
print, R24_J, R24_new(fi:se)
print, R34_J, R34_new(fi:se)
print, R44_J, R44_new(fi:se)

print, Q14_J, Q14_new(fi:se)
print, Q24_J, Q24_new(fi:se)
print, Q34_J, Q34_new(fi:se)
print, Q44_J, Q44_new(fi:se)
print, P14_J, P14_new(fi:se)
print, P24_J, P24_new(fi:se)
print, P34_J, P34_new(fi:se)
print, P44_J, P44_new(fi:se)
endif

end

;###################################################################
