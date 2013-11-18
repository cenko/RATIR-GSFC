function headspec, filename

  ; reads the header of an ASCII spectrum in Perley's standard format, analogous to a fits file.

  close, 1
  nl = countlines(filename)  

  openr, 1, filename
  inline = ''
  header = ['']

  i = 0
  for l = 0, nl-1 do begin
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
        if l gt 1 then break
     endelse
  endfor

  close, 1
  h = header
  return, h

end

