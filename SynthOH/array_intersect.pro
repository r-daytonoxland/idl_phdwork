FUNCTION array_intersect, array1, array2, index=index
;
; Returns all the elements of array1 which are also present in array2. Note
; that the function returns a vector array, as array1 and array2 are converted
; to single dimension arrays internally. If none of the elements of array1
; are present in array2 then NaN is returned.
;
; If the index keyword is set then indices into array1 are returned instead
; of the actual values. NOTE: If the values of array1 or array2 are not unique,
; then only the last occurence of the search value is returned.
;

all=[reform(array1,n_elements(array1)),reform(array2,n_elements(array2))]
allf=[replicate(0b,n_elements(array1)),replicate(1b,n_elements(array2))]
s=sort(all)
all=all[s]
allf=allf[s]
s2=where((all eq shift(all,-1)) and (allf ne shift(allf,-1)),c)
if (c le 0) then return,!values.F_NaN else if keyword_set(index) then return,s[s2] else return,all[s2]

end
