function getsexpath
  common sexblock, sexpath
  if n_elements(sexpath) eq 0 then sexpath = '' ; default
  return, sexpath
end

