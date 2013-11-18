function clipzerostr, strarrin
  strarrout = strarr(n_elements(strarrin))
  for p = 0, n_elements(strarrin)-1 do begin
    strarrout[p] = clipzero(strarrin[p])
  endfor
  return, strarrout
end

function clipzero, strin
  if n_elements(strin) gt 1 then return, clipzerostr(strin)
  if n_elements(strin) eq 0 then fail
  str = string(strin)
  str = strtrim(str,2)
  if strpos(str,'.') eq -1 then return, str
  curlen = strlen(str)
  epos = strpos(str,'e')
  i = (curlen-1)
  if epos gt 0 then begin
     i = epos-1
     ee = strmid(str,epos)
  endif
  while strmid(str, i, 1) eq '0' and i gt 0 do i = i - 1
  if strmid(str, i, 1) eq '.' and i gt 0 then i = i - 1
  str = strmid(str, 0, i+1)
  if n_elements(ee) gt 0 then str = str + ee
  return, str
end
