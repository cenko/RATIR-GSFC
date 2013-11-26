function unique, arr, count=count

  uniquearr = arr
  n = n_elements(arr)
  nunique = 0L

  for i = 0L, n-1 do begin
    known = 0

    for j = 0L, nunique - 1 do begin
      if arr[i] eq uniquearr[j] then known = 1
    endfor

    if known eq 0 then begin
      uniquearr[nunique] = arr[i]
      nunique = nunique + 1
    endif
  endfor
  uniquearr = uniquearr[0:nunique-1]

  count = intarr(nunique)
  for u = 0L, nunique-1 do begin
    dum = where(arr eq uniquearr[u], ctu)
    count[u] = ctu
  endfor

  return, uniquearr

end
