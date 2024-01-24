FUNCTION oh_label, qnum
;
; A routine to make an OH line label from the quantum numbers of the upper
; and lower states as given in the data from Brooke et al. 2016
; Inputs:
;  qnum - The quantum numbers, a 6-elements vector with these elements:
;   0: v1 - upper state vibrational level
;   1: v2 - lower state vibrational level
;   2: J1 - upper state J
;   3: J2 - lower state J
;   4: F1 - upper state F level
;   5: F2 - lower state F level
; These are the first 6 columns in the Brooke et al data file.
;

if qnum[4] eq qnum[5] then s=qnum[4] else s=(qnum[4]*10)+qnum[5]
letters=['P','Q','R']
return, string(qnum[0],qnum[1],letters[1+qnum[2]-qnum[3]],s,qnum[3]+qnum[4]-1.5,form='("OH ",i0,"-",i0," ",a1,"!D",i0,"!N(",i0,")")')

end
