pro fringecorrect, filename, fringefile=fringefile, fringeimg=fringeimg, median=median, mode=mode, prefix=prefix, scale=scale, iterate=iterate, itersec=itersec

  if n_elements(fringeimg) eq 0 then begin
    if n_elements(fringefile) eq 0 then fringefile = 'fringe.fits'
    if strpos(fringefile,'.fits') eq -1 then begin
      print, 'Fringe file must be of fits (.fits) format'
      return
    endif

    fringeimg = float(readfits(fringefile,fringehead,/silent))
  endif
  if n_elements(filename) gt 1 then begin
    fringeimg = float(readfits(fringefile,fringehead,/silent))
    if n_elements(scale) eq n_elements(filename) then begin
      for i = 0, n_elements(filename)-1 do $
         fringecorrect, filename[i], fringefile=fringefile, fringeimg=fringeimg, $
                            median=median, mode=mode, prefix=prefix, scale=scale[i], iterate=iterate, itersec=itersec
    endif else begin
      for i = 0, n_elements(filename)-1 do $
         fringecorrect, filename[i], fringefile=fringefile, fringeimg=fringeimg, $
                              median=median, mode=mode, prefix=prefix, scale=scale, iterate=iterate, itersec=itersec
    endelse
    return
  endif else begin
  if strpos(filename,'.fit') eq -1 then begin
    if strpos(filename, '.txt') eq -1 and $
       strpos(filename,'.lis') eq -1 then begin
      print, 'File must be fits (.fits) or a list (.txt, .lis)'
      return
    endif

    iline = ''
    openr, 1, filename
    while not eof(1) do begin
      readf, 1, iline
      fringecorrect, iline, fringefile=fringefile, fringeimg=fringeimg, $
                            median=median, mode=mode, prefix=prefix, iterate=iterate, itersec=itersec
    endwhile

    close, 1
    return
  endif
  endelse


  ;print, filename

  tstart = systime(/seconds)

  if n_elements(prefix) eq 0 then prefix = 'f'

  if n_elements(scale) eq 0 then begin
  if keyword_set(median) eq 0 and keyword_set(mode) eq 0 then begin
    median = 1
    ;return
  endif

  file = float(readfits(filename,fhead,/silent))

  med = median(file)
  sclip = file[where(file lt 2 * med)]
  
  if keyword_set(median) then begin
    med = median(sclip)
    ;print, filename, ' : median = ', strtrim(med,2)
    avg = med
  endif
  if keyword_set(mode) then begin
    h = histogram(sclip,binsize=5)
    mode = min(sclip) + 5*where(h eq max(h)) + 2.5
    ;print, filename, ' :   mode = ', strtrim(mode,2)
    avg = mode[0]
  endif
  endif 
 
  ;corr = file - avg * fringeimg  + avg


  if n_elements(scale) eq 0 then imscale = avg else imscale = scale

  print, filename, ' :   scale = ', strtrim(imscale,2)

  rescale = 1.0
  if n_elements(iterate) gt 0 then begin
     if n_elements(itersec) eq 0 then isec = [0,((size(corr)) [1])-1,0,((size(corr)) [2])-1] $
                                 else isec = itersec
     if isec[1] lt 0 then isec[1] = (size(file)) [1]+isec[1]
     if isec[3] lt 0 then isec[3] = (size(file)) [2]+isec[3]

     sec       =      file[isec[0]:isec[1],isec[2]:isec[3]]
     fringesec = fringeimg[isec[0]:isec[1],isec[2]:isec[3]]
     maskobjects, sec, /nan
     
     stepfactor = 0.5^(1+findgen(11)) ;[0.5,0.25,0.125,0.0625,0.03125,0.015625...]
     niter = n_elements(stepfactor)
     skydev = sigclipabsdev(sec - imscale*fringesec)
     for i = 0, niter-1 do begin
        step = stepfactor[i]
        skydevp = sigclipabsdev(sec - imscale*(rescale+step)*fringesec)
        skydevm = sigclipabsdev(sec - imscale*(rescale-step)*fringesec)
        print, clip(rescale), ':=', fpr(skydev,4.5),'    ', clip(rescale+step), ':=', fpr(skydevp,4.5), '    ', clip(rescale-step), ':=', fpr(skydevm,4.5)
        if      skydevp lt skydev then begin
            rescale = rescale + step
            skydev = skydevp
        endif else begin
        if skydevm lt skydev then begin
            rescale = rescale - step
            skydev = skydevm
        endif
        endelse
     endfor
     print, ' x ', rescale
  endif

  corr = file - imscale*rescale*fringeimg


  sxaddpar, fhead, 'FRINGEFL', fringefile, 'Fringe file'
  sxaddpar, fhead, 'FRINGESC', imscale*rescale, 'Fringe scale'

  writefits, prefix+filename, corr, fhead

  ;print, string(systime(/seconds)-tstart,format='(f6.2)'),' seconds to normalize ',filename

end
