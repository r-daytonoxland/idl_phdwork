

@N2_balance_4VKband


;-----------------------------VIBRATIONAL AND ROTATIONAL FUNCTIONS -----------------------

;Vibrational Function

function func_Gv,we,wexe,weye,weze,weae,v
Gv=(we*(v+0.5))-(wexe*(v+0.5)^2)+(weye*(v+0.5)^3)+(weze*(v+0.5)^4)+(weae*(v+0.5)^5)
return, Gv
end


;Rotational Functions


;#1
function func_Dv,De,betae,v

    Dv=De-(betae*(v+0.5))
    return, Dv

end

;#2

function func_Bv, Be,ae,ye,dele,v

    Bv=Be-(ae*(v+0.5))-(ye*(v+0.5)^2)+(dele*(v+0.5)^3)
    return, Bv

end


;#3

function lambda_v,A_lambda,v

   if (v ge 0) and (v le 13) then return, A_lambda[v] $
   else return, (-1.335 +0.007*v)

end


;#4

function func_FvJ, Bv,Dv,J,v

    ;Bv=func_Bv(Be,ae,ye,dele,v)
    ;Dv=_func_Dv(De,betae,v)
   
    Fv=(Bv*J*(J+1))-(Dv*(J^2)*((J+1)^2))
    return, Fv

end



;#5

function func_F1, Bv,Dv,A_gamma,A_lambda,N,v

    ;Bv=Be-(ae*(v+0.5))-(ye*(v+0.5)^2)+(dele*(v+0.5)^3)
    ;Dv=De-(betae*(v+0.5))
    
;    Bv=func_Bv(Be,ae,ye,dele,v)
 ;   Dv=_func_Dv(De,betae,v)
    
    Fv=(Bv*N*(N+1))+(Dv*(N^2)*((N+1)^2))  + ((2*N+3)*Bv) - lambda_v(A_lambda,v)- $
       sqrt((((2*N+3)*Bv)^2 )  +  (lambda_v(A_lambda,v)^2)  - (2*lambda_v(A_lambda,v)*Bv)) + (A_gamma*(N+1))

    return, Fv

end


;#6

function func_F2, Bv,Dv,N,v

    
    ;Bv=func_Bv(Be,ae,ye,dele,v)
    ;Dv=_func_Dv(De,betae,v)
    
    Fv=(Bv*N*(N+1))+(Dv*(N^2)*((N+1)^2)) 

    return, Fv

end



;#7

function func_F3, Bv,Dv,A_gamma,A_lambda,N,v

        
    ;Bv=func_Bv(Be,ae,ye,dele,v)
    ;Dv=_func_Dv(De,betae,v)
    
    Fv=(Bv*N*(N+1))+(Dv*(N^2)*((N+1)^2))  - ((2*N-1)*Bv) - lambda_v(A_lambda,v)+ $
       sqrt((((2*N-1)*Bv)^2 )  +  (lambda_v(A_lambda,v)^2)  - (2*lambda_v(A_lambda,v)*Bv)) - (A_gamma*N)

    return, Fv

end




;________________________________________________________________________________________

pro prog_n2_VKband_v1,Ts,wv_VK,int_VK

;Program to get the VK bands of N2 molecule -Srimoyee 2/4/25


;Inputs: Ts:Temperature (K)
;Outputs: wv_VK[n_elements(Ts),vL,vU,J] : (nm)wavelengths as a function of Ts-Temperature, vL-Lower vibrational band, vU-Upper vibrational band, J-Number of transitional lines
;	: int_VK[n_elements(Ts),vL,vU,J]:intensities NOTE:intensities are not normalized







;Franck-Condon and Transition Moment Files

restore, 'ElecTransMoment_VKband.sav'
;file has Re_VK[v'',v'] ; v''->vL, v'->vU

Re_VK=Re_VK/0.393456  ; convert from a.u. to debye, 1 debye =0.393456 a.u 

nvibL_AX= (size(Re_VK,/dimensions))[0] ;number of lower vibrational levels

nvibU_AX= (size(Re_VK,/dimensions))[1] ;number of upper vibrational levels



restore,'FCfactors_VKband.sav'
;file has FC_VK[vL,vU] ; v''->vL, v'->vU


;Rotational branches
n_lines=101 


N2_balance_4VKband,Ts,N0

;N0[temperature,vU], vU=22



;-----------------------------------------------------------

;Constants
c  = 2.99792458e10              ;cm s-1
h  = 6.6606876e-27           ;[ergs][s]
k  = 1.3806503e-16           ;[ergs][K-1]
;-----------------------------------------------------------

;From NIST
;Te	minimum electronic energy (cm-1)
;we	vibrational constant - first term (cm-1)
;wexe	vibrational constant - second term (cm-1)
;weye	vibrational constant - third term (cm-1)

;Be	rotational constant in equilibrium position (cm-1)
;ae	rotational constant - first term (cm-1)
;ye	rotation-vibration interaction constant (cm-1)
;De	centrifugal distortion  (cm-1)
;beta_e rotational constant - first term, centrifugal force (cm-1)
;re internuclear distance (angstrom)

;------------------For A3Sigma_u_pls-------------------------
;NIST data
Te_A=50203.6
we_A=1460.64
wexe_A=13.87
weye_A=0.0103	
;Be_A=1.4546
;ae_A=0.0180
;ye_A=-8.8e-5
;De_A=6.15e-6  
;beta_A=	 	
;re_A=1.2866	

;Gilmore data

;Te_A=49754.8
;Vibrational constants
;we_A=1460.48
;wexe_A=13.775
;weye_A=-1.175d-2
weze_A= 1.41d-4
weae_A=-7.29d-5	

;Rotational Constants
Be_A=1.45499
ae_A=1.8385d-2
ye_A=1.24d-5
dele_A=-6.7d-6

De_A=6.15e-6  ;NIST data
betae_A= 3.9d-8	;Justin Yonker Master's Thesis 

A_gamma=-0.003  ;cm-1
A_lambda=[-1.326 ,-1.320,-1.314,-1.307,-1.300,-1.293,-1.285,-1.276,-1.268,-1.259,-1.249,-1.239,-1.229,-1.219] ;in cm-1
 	



;-------------------For X1Sigma_u_pls------------------------
Te_X=0.

;Vibrationsl constants -NIST
we_X=2358.57
wexe_X=14.324
weye_X=	-2.26e-3
weze_X=-2.4d-4
weae_X=0d

;Rotational constants
Be_X= 1.998241
ae_X= 0.017318
ye_X=-3.3d-5
dele_X=0d

De_X= 5.76d-6  ;NIST data
betae_X=4.6d-9	;Justin Yonker Master's Thesis  		 	

re_X= 1.097685	
;_______________________________________________________________________________________





N=indgen(n_lines)      ; N = 0,1,2,3,....100
J1=N+1                 ; J1= 1,2,3,4,....101 (N+1)
J2=indgen(n_lines-1)+1 ; J2= 1,2,3,....100 (=N without the 0 line)
J3=indgen(n_lines-1)   ; J3= 0,1,2,....99  (=N-1)

PQ=dblarr(nvibL_AX,nvibU_AX,n_elements(J1))
PP=dblarr(nvibL_AX,nvibU_AX,n_elements(J2))
RR=dblarr(nvibL_AX,nvibU_AX,n_elements(J2))
RQ=dblarr(nvibL_AX,nvibU_AX,n_elements(J3))


I_PQ=dblarr(nvibL_AX,nvibU_AX,n_elements(J1))
I_PP=dblarr(nvibL_AX,nvibU_AX,n_elements(J2))
I_RR=dblarr(nvibL_AX,nvibU_AX,n_elements(J2))
I_RQ=dblarr(nvibL_AX,nvibU_AX,n_elements(J3))

wv_VK=  dblarr(n_elements(Ts),nvibL_AX,nvibU_AX,n_elements(J1)+n_elements(J2)+n_elements(J2)+n_elements(J3))
int_VK= dblarr(n_elements(Ts),nvibL_AX,nvibU_AX,n_elements(J1)+n_elements(J2)+n_elements(J2)+n_elements(J3))


for tt=0, n_elements(Ts)-1 do begin 
  
  for vU= 0,nvibU_AX-1 do begin
   
    Gv_A=func_Gv(we_A,wexe_A,weye_A,weze_A,weae_A,vU)
    Bv_A=func_Bv(Be_A,ae_A,ye_A,dele_A,vU)
    Dv_A=func_Dv(De_A,betae_A,vU)
    
    ;Qr= k*Ts(tt)/(h*c*Bv_A)   ;MODIFY THIS LATER TO 
    

   for vL= 0,nvibL_AX-1 do begin
       
        Gv_X=func_Gv(we_X,wexe_X,weye_X,weze_X,weae_X,vL)
        Bv_X=func_Bv(Be_X,ae_X,ye_X,dele_X,vL)
        Dv_X=func_Dv(De_X,betae_X,vL)
  
     
        ;Rotational Branches
	
	;#1-------------------------------PQ branch-------------------------------------------------------------------
	
	F1=func_F1(Bv_A,Dv_A,A_gamma,A_lambda,J1,vU)  ;J1(nn)
	phi=(2-(J1 mod 2))/3.
	FvJ=func_FvJ(Bv_X,Dv_X,J1,vL)                 ;J1(nn)

	Qr=total( (2*J1+1)*exp(-(F1*h*c)/(k*Ts(tt))))


        PQ(vL,vU,*)=Te_A+Gv_A+F1-Te_X-Gv_X-FvJ

	    
        ;Intensity
	I_PQ(vL,vU,*)=(64* ((!pi*PQ[vL,vU,*])^4) *N0[tt,vU]*c/(3*Qr)) * ((Re_VK(vL,vU)^2)*FC_VK(vL,vU)*0.5*(J1+0.5)*phi)   *exp(- ((F1+Gv_A)*h*c)/(k*Ts(tt)) )  


       
       
       
       ;#2--------------------------------- PP branch-----------------------------------------------------------------
            
	F2=func_F2(Bv_A,Dv_A,J2,vU)  ;J2(nn)
	phi=(1+(J2 mod 2))/3.
        FvJ=func_FvJ(Bv_X,Dv_X,J2+1,vL)  ;J2(nn)+1
        Qr=total( (2*J2+1)*exp(-(F2*h*c)/(k*Ts(tt))))
	    
        PP(vL,vU,*)=Te_A+Gv_A+F2- Te_X-Gv_X-FvJ
			
			
       ;Intensity
	  
        I_PP(vL,vU,*)=(64* ((!pi*PP[vL,vU,*])^4) *N0[tt,vU]*c/(3*Qr)) * ((Re_VK(vL,vU)^2)*FC_VK(vL,vU)*0.5*(J2+1)*phi)   *exp(- ((F2+Gv_A)*h*c)/(k*Ts(tt)) )  


        
      ;#3 ------------------------------RR branch--------------------------------------------------------------------
	  
	  FvJ=func_FvJ(Bv_X,Dv_X,J2-1,vL)  ;J2(nn)-1

          RR(vL,vU,*)=Te_A+Gv_A+F2- Te_X-Gv_X-FvJ
	
	
	 ;Intensity
	  
	  I_RR(vL,vU,*)=(64* ((!pi*RR[vL,vU,*])^4) *N0[tt,vU]*c/(3*Qr)) * ((Re_VK(vL,vU)^2)*FC_VK(vL,vU)*0.5*J2*phi)   *exp(- ((F2+Gv_A)*h*c)/(k*Ts(tt)) )  
		
           

 

     ;#4 ------------------------------RQ branch-----------------------------------------------------------------------


   
       F3=func_F3(Bv_A,Dv_A,A_gamma,A_lambda,J3,vU); J3(nn)
       phi=(2-(J3 mod 2))/3.
       FvJ=func_FvJ(Bv_X,Dv_X,J3,vL)        ;   J3(nn)
       Qr=total( (2*J3+1)*exp(-(F3*h*c)/(k*Ts(tt))))
      
       RQ(vL,vU,*)=Te_A+Gv_A+F3- Te_X-Gv_X-FvJ  

      ;Intensity
        
       I_RQ(vL,vU,*)=(64* ((!pi*RQ[vL,vU,*])^4) *N0[tt,vU]*c/(3*Qr)) * ((Re_VK(vL,vU)^2)*FC_VK(vL,vU)*0.5*(J3+0.5)*phi)   *exp(- ((F3+Gv_A)*h*c)/(k*Ts(tt)) )  
	
	
	
     ;-----------------------------------------------------------------------------------------------------------------------------------
     
     
     ;Combine all 4 branches together and sort according to wavelengths
     
     
     wvl=[reform(PQ(vL,vU,*)), reform(PP(vL,vU,*)), reform(RR(vL,vU,*)),reform(RQ(vL,vU,*))]   
     int=[reform(I_PQ(vL,vU,*)), reform(I_PP(vL,vU,*)), reform(I_RR(vL,vU,*)),reform(I_RQ(vL,vU,*))]   

     wv_VK[tt,vL,vU,*]=1.e7/wvl[sort(wvl)]
     int_VK[tt,vL,vU,*]=int[sort(wvl)]

     
     


   endfor
   
 endfor



endfor


    
    






end
