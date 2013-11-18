pro writespec, filename, spec, h

  if n_elements(spec) eq 0 then begin
    print, 'No spectral data; not writing file to disk.'
  endif

  openw, 1, filename
  inline = ''
  for i = 0, n_elements(h)-1 do begin
    printf, 1, '# '+h[i]
  endfor

  skyyes = total(spec.sky) ne 0
  uncyes = total(spec.unc) ne 0
  pixyes = total(spec.pix) ne 0
  respyes = total(spec.resp) ne 0
  if respyes eq 0 then printf, 1, '# wavelength     flux         sky_flux     flux_unc     pixel'
  if respyes eq 1 then printf, 1, '# wavelength     flux         sky_flux     flux_unc     pixel    response'
  for i = 0, n_elements(spec)-1 do begin
    if respyes eq 0 then printf, 1, fpr(spec[i].wav,6.4), spec[i].flux, spec[i].sky, spec[i].unc, '  ', clip(long(spec[i].pix),7)
    if respyes eq 1 then printf, 1, fpr(spec[i].wav,6.4), spec[i].flux, spec[i].sky, spec[i].unc, '  ', clip(long(spec[i].pix),7), spec[i].resp

    ;if skyyes eq 0 and pixyes eq 0 then printf, 1, spec[i].wav, spec[i].flux
    ;if skyyes eq 1 and pixyes eq 0 then printf, 1, spec[i].wav, spec[i].flux, spec[i].sky
    ;if skyyes eq 1 and pixyes eq 1 then printf, 1, spec[i].wav, spec[i].flux, spec[i].sky, spec[i].unc
    ;if skyyes eq 1 and pixyes eq 1 then printf, 1, spec[i].wav, spec[i].flux, spec[i].sky, spec[i].pix
  endfor


  close, 1
end
