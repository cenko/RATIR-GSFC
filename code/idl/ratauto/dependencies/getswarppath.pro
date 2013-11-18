function getswarppath
  common swarpblock, swarppath
  if n_elements(swarppath) eq 0 then swarppath = '' ; default
  return, swarppath
end

