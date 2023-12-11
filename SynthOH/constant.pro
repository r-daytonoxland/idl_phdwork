FUNCTION constant, kB=kB,amu=amu,mu0=mu0,q=q,c=c,h=h,Re=Re,me=me,Na=Na,e0=e0
;
; Returns the value of a physical constant as a double, mostly in SI units.
; You should set ONE keyword only, to determine the returned constant:
;  kB - Boltzmann constant, in units of J/K
;  amu - Atomic mass unit, in kg
;  mu0 - Magnetic permeability of free space, in H/m (= N/A^2)
;  e0 - Electric permittivity of free space, in F/m
;  q - Electron charge, in C
;  c - Speed of light, in m/s
;  h - Planck's constant, in Js
;  Re - Earth's mean radius, in m
;  me - Electron mass, in kg
;  Na - Avogadro's number

if keyword_set(kB) then return, 1.38064852d-23 ; J/K
if keyword_set(amu) then return, 1.660539040d-27 ; kg
if keyword_set(mu0) then return, !dpi*4d-7 ; H/m
if keyword_set(q) then return, 1.6021766208d-19 ; C
if keyword_set(c) then return, 2.99792458d8 ; m/s
if keyword_set(h) then return, 6.626069934d-34 ; Js
if keyword_set(Re) then return,6371.0d3 ; m
if keyword_set(me) then return,9.10938356d-31 ; kg
if keyword_set(Na) then return,6.02214076d23
if keyword_set(e0) then return, 8.8541878128d-12 ; F/m
return,!values.F_NaN

end
