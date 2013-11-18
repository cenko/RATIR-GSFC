function mag2flux, mag, filter

  ;units of microJy

  return, 10.0^(zeropt(filter)-mag/2.5)

end
