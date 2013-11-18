function flux2mag, flux, filter

  ;units of microJy

  mag = 2.5*(zeropt(filter)-alog10(flux))
  negflux = where(flux le 0, ct)
  if ct gt 0 then mag[negflux] = 99
  return, mag

end
