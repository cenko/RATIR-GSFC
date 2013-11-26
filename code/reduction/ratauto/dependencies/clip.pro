function clip, strin, maxlen
  str = string(strin)
  if n_params() eq 1 then return, strtrim(str,2)

  str = strmid(strtrim(str,2),0,maxlen)
  curlen = strlen(str)
  if curlen lt maxlen then $
    for i = 0, (maxlen-curlen-1) do $
      str = str + ' '

  return, str
end