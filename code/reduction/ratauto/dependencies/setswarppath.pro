pro setswarppath, swarppathin
  common swarpblock, swarppath
  swarppath = swarppathin
  if strlen(swarppath) gt 0 and strmid(swarppath,strlen(swarppath)-1) ne '/' then swarppath += '/'
end

