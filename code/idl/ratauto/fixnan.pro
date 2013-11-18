pro fixnan, array, mask

if n_elements(flag) gt 0 then array[where(mask gt 0)] = !values.f_nan
if n_elements(fixdefect) gt 0 then begin
  array[where(mask gt 0)] = !values.f_nan
  ct = 1
  while ct gt 0 do begin 
    bad = where(finite(array) eq 0, ct)
    for m = 0, ct-1 do begin
      mm = bad[m]
      array[mm] = median(array[mm-2:mm+2])
    endfor
  endwhile
endif

end

