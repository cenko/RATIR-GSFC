function readcalibration, calfile, filter=filter, magcolnum=magcolnum, unccolnum=unccolnum, f2=f2

  ; Read in calibration star file

  ; f2 is a second filter.
  imag = -1
  imag2 = -1
  iunc2 = -1

  if file_test(calfile) then begin
     ns = countlines(calfile)-1
     cat = replicate({ra:0.D, dec:0.D, mag:0., unc:0., mag2:0., unc2:0.}, ns)

     inline = ''
     openr, 1, calfile
     readf, 1, inline
     headings = strsplit(strmid(inline,1),/extract)
     lowheadings = strlowcase(headings)
     nh = n_elements(headings)

     ira = (where(lowheadings eq 'ra' or lowheadings eq '#ra' or lowheadings eq 'raj2000' or lowheadings eq '#raj2000', ct)) [0]
     if ct eq 0 then ira = 0

     idec = (where(lowheadings eq 'dec' or lowheadings eq '#dec' or lowheadings eq 'dej2000', ct)) [0]
     if ct eq 0 then idec = 1

     if n_elements(magcolnum) eq 0 then begin
        if n_elements(filter) eq 0 then begin
          imag = idec + 1
          print, 'No filter specified, using first column after declination ('+headings[imag]+')'
        endif else begin
          imag = where(headings eq filter or headings eq filter+'mag', ct)
          if ct eq 0 then begin
             imag = idec+1
             if n_elements(headings) le imag then begin
                print, 'No magnitude information in header.'
                imag = -1
             endif else begin
                print, 'Cannot match filter name in header ('+filter+').  Using first column after declination. ('+headings[imag]+')'
             endelse
          endif
        endelse
     endif else begin
        imag = magcolnum
     endelse

     if n_elements(f2) gt 0 then begin
        imag2 = where(headings eq f2 or headings eq f2+'mag', ct)
        iunc2 = where(headings eq f2+'unc' or headings eq f2+'err' or headings eq f2+'_u' or headings eq 'e'+f2 or headings eq 'e_'+f2, ct)
     endif

   
     if n_elements(unccolnum) eq 0 and n_elements(filter) gt 0 then begin
        iunc = where(headings eq filter+'unc' or headings eq filter+'err' or headings eq filter+'_u' or headings eq 'e'+filter or headings eq 'e_'+filter or headings eq 'e'+filter+'mag' or headings eq 'e_'+filter+'mag', ct)
        if ct eq 0 then iunc = -1
     endif else begin
        if n_elements(unccolnum) gt 0 then iunc = unccolnum else iunc = -1
     endelse

     for s = 0, ns-1 do begin
        readf, 1, inline
        inline = clip(inline)
        if strlen(inline) eq 0 then continue
        if strmid(inline,0,1) eq '#' then continue
        inarr = strsplit(inline,/extract)
        inra = inarr[ira]
        indec = inarr[idec]
        if strpos(inra,':')  gt 0 then cat[s].ra  = 15.*ten(inra) else cat[s].ra = double(inra)
        if strpos(indec,':') gt 0 then cat[s].dec =   ten(indec)  else cat[s].dec = double(indec)
        if imag gt -1 then if clip(inarr[imag]) ne '-' then cat[s].mag = float(inarr[imag])
        if iunc gt 0 then if clip(inarr[iunc]) ne '-' then cat[s].unc = float(inarr[iunc])
        if imag2 gt -1 then begin
          if clip(inarr[imag2]) ne '-' then cat[s].mag2 = float(inarr[imag2])
          if iunc gt 0 then if clip(inarr[iunc]) ne '-' then cat[s].unc2 = float(inarr[iunc2])
        endif
     endfor
     close, 1

     good = where(cat.mag ne 0.0, ct)
     if ct eq 0 then print, 'Found no magnitudes in calibration file: wrong column?' else cat = cat[good]

     return, cat

  endif else begin
     print, 'Cannot find calibration file.'
     return, -1
  endelse


end

