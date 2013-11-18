pro setsexpath, sexpathin
  common sexblock, sexpath
  sexpath = sexpathin
  if strlen(sexpath) gt 0 and strmid(sexpath,strlen(sexpath)-1) ne '/' then sexpath += '/'
end


