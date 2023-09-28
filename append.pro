function append, a, b

isa = n_elements(a) ne 0
isb = n_elements(b) ne 0
if ~isa && ~isb then $
message, 'One of the input vectors must exist!'

if ~isa && isb then return, b
if ~isb && isa then return, a
  
if n_elements(a) eq 0 then return, b
if n_elements(b) eq 0 then return, a
return, [a, b]
end

