function fpr, val, fmt
  if n_elements(val) eq 0 then return, ''
  if n_elements(fmt) eq 0 then return, clip(float(val))
  nleft = fix(fmt) > 1
  ndec = fix((fmt+0.001 - fix(fmt)) * 10)

  return, string(val, $
                 format='(f'+clip(nleft+ndec+(ndec gt 0))+'.'+clip(ndec)+')')
end

