function readspec, filename, h, header=header, timenorm=timenorm, cutoff=cutoff

  ; reads an ASCII spectrum in Perley's standard format, analogous to a fits file.

  close, 1
  if file_test(filename) eq 0 then begin
     print, 'File '+filename+' not found!'
     return, ['']
  endif
  if strpos(filename,'.fit') ge 0 then begin
     print, 'Cannot process fits files.'
     return, ['']
  endif

  nl = countlines(filename)  
  spec = replicate({wav:0.0d, flux:0.0, sky:0.0, unc:0.0, pix:0L, resp:0.0},nl)
  header = ['']

  openr, 1, filename
  inline = ''
  i = 0
  for l = 0l, nl-1 do begin
     readf, 1, inline
     inline = clip(inline)
     if strlen(inline) eq 0 then continue
     if strmid(inline,0,1) eq '#' then begin
        inline = clip(strmid(inline,1))
        hline = processheaderline(inline)
        if hline ne '' then begin
          header = [header, clip(hline,78)]  ; might not be strictly valid
        endif
     endif else begin
        splitline = strsplit(inline, /extract)
        ni = n_elements(splitline)
        if ni eq 0 then continue
        if ni lt 2 then begin
           print, 'Problem reading this line:', inline
        endif else begin
           spec[i].wav = float(splitline[0])
           spec[i].flux = float(splitline[1])
           if ni ge 3 then spec[i].sky = float(splitline[2])
           if ni ge 4 then spec[i].unc = float(splitline[3])
           if ni ge 5 then spec[i].pix = float(splitline[4])
           if ni ge 6 then spec[i].resp = float(splitline[5])
        endelse
        i += 1
     endelse
  endfor

  spec = spec[0:i-1]

  if keyword_set(timenorm) then begin
       tn = sxpar(header,'ELAPTIME')
       spec.flux /= tn
       spec.unc /= tn
       spec.sky /= tn
  endif

  dwav = spec.wav - shift(spec.wav,+1)
  dwav[0] = dwav[1]
  negdwav = where(dwav lt 0, ctnegdwav)
  if ctnegdwav gt 0 and ctnegdwav lt n_elements(dwav) then begin
      print, 'Warning: Wavelength solution for this spectrum is double-valued.'
      print, '     Removing points below the inflection. '
   ; This is a temporary solution; you won't always want this behavior.  Work on the wavelength solver to
   ;   recognize and stop such bad solutions earlier in the process.
      spec = spec[where(dwav gt 0.)]
  endif

  if n_elements(cutoff) ge 1 then begin
     spec = spec[where(spec.wav gt cutoff[0])]
  endif

  close, 1
  h = header
  return, spec

end

