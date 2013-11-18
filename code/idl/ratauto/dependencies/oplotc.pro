pro oplotc, x, y, color=color, _EXTRA=ex

  if n_elements(color) eq 1 then color = intarr(n_elements(x))+color

  for i = 0, n_elements(x)-1 do begin
     oplot, [x[i]], [y[i]], color=color[i], _EXTRA=ex
  endfor

end

