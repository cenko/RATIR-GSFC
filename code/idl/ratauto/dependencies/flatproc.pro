pro flatproc, filename, flatname, crop=crop, flatminval=flatminval, flatmaxval=flatmaxval

; now corrects for different windowings (LTV); soon correct for different binnings (LTM)

if n_elements(filename) eq 1 then begin
  filearr = strsplit(filename,'.', /extract)
  fileroot = filearr[0]
  fileext = filearr[1]

  if fileext eq 'cat' or fileext eq 'lis' or fileext eq 'list' or fileext eq 'txt' then begin
    filenames = grabcolumn(filename,0,/str)
  endif else  begin
    filenames = filename
  endelse
endif else begin
  filenames = filename
endelse

nfile = n_elements(filenames)
flat = mrdfits(flatname, 0, hflat, /silent)

flattype = sxpar(hflat,'FLATTYPE')
med = median(flat)
if med lt 0.5 or med gt 2.0 then begin
  print, 'Warning: flat is not normalized to zero'
endif

for f = 0, nfile-1 do begin
   data = mrdfits(filenames[f],0,h, /silent)
   if n_elements(data) ne n_elements(flat) then begin
      print, 'WARNING - Data and flat field are different dimensions!!'
         ; runs, but doesn't work; check this later (is this still true????????)
      dataoriginx = -sxpar(h,'LTV1')
      dataoriginy = -sxpar(h,'LTV2')
      flatoriginx = -sxpar(hflat,'LTV1')
      flatoriginy = -sxpar(hflat,'LTV2')
      nxdata = (size(data)) [1]
      nxflat = (size(flat)) [1]
      nyflat = (size(flat)) [2]
      pflat = data*0. + !values.f_nan
      xdata = lindgen(n_elements(data)) mod nxdata  ; 1D array
      ydata = lindgen(n_elements(data)) / nxdata
      xflat = xdata + flatoriginx-dataoriginx ; x pixels on the flat corresponding to (xdata,ydata)
      yflat = ydata + flatoriginy-dataoriginy ; y pixels on the flat corresponding to (xdata,ydata)
      noflatpix = where(xflat lt 0 or xflat ge nxflat or yflat gt 0 or yflat ge nyflat, ct, complement=flatpix)
      ; stop   ; cautiously  not stopping
      if n_elements(flatpix) eq 0 then begin
        print, 'WARNING - Flat processing failed - no overlap!'
      endif
      pflat[xdata[flatpix],ydata[flatpix]] = flat[xflat[flatpix],yflat[flatpix]]
      if ct gt 0 then begin
        print, 'WARNING - Flat does not fully overlap science image!  Cannot flatfield all pixels.'
        pflat[xdata[noflatpix],ydata[noflatpix]] =  1
      endif
      if total(finite(pflat)) eq 0 then begin
        print, 'WARNING - Flat processing failed - no overlap!'
        return
      endif
   endif else begin
      pflat = flat
   endelse

   ; determine crop zones
;
;   if n_elements(crop) eq 1 and size(crop, /tname) eq 'STRING' then begin
;      if crop eq 'auto' or crop eq 'autoy' then findcroplines, pflat, xcrop, ycrop, thresh=0.2
;      if      crop eq 'auto' then crop = [xcrop[0], ycrop[0], xcrop[1], ycrop[1]] $
;      else if crop eq 'autoy' then crop = [0, ycrop[0], ((size(fdata)) [1])-1, ycrop[1]]
;   endif

   ; set values too low/high to NaNs
   if n_elements(flatminval) gt 0 then begin
      w = where(pflat lt flatminval, ct, complement=goodsignal)
      if ct gt 0 then pflat[w] = !values.f_nan
   endif else begin
      w = where(pflat lt 0.1, complement=goodsignal) ; for sky count determination only
   endelse
   if n_elements(flatmaxval) gt 0 then begin
      w = where(pflat gt flatmaxval, ct)
       if ct gt 0 then pflat[w] = !values.f_nan   ; Nans can complicate life, though useful as fix markers.
   endif
 
   fdata = data / pflat 

;   if n_elements(crop) gt 0 then begin
;      if total(crop eq 0) ne n_elements(crop) then begin
;        s = size(fdata)
;        if crop[2] le 0 then crop[2] = s[1] + crop[2]-1  ; can use negatives to count from right/top edges
;        if crop[3] le 0 then crop[3] = s[2] + crop[3]-1
;        fdata = fdata[crop[0]:crop[2],crop[1]:crop[3]]
;        sxaddpar, h, 'LTV1', sxpar(h,'LTV1')-crop[0]
;        sxaddpar, h, 'LTV2', sxpar(h,'LTV2')-crop[1]
;        sxaddpar, h, 'CROPX0', crop[0]
;        sxaddpar, h, 'CROPX1', crop[2]
;        sxaddpar, h, 'CROPY0', crop[1]
;        sxaddpar, h, 'CROPY1', crop[3]
;        endif
;   endif
   sxaddpar, h, 'FLATFLD', flatname
   sxaddpar, h, 'FLATTYPE', flattype
   skycts = median(fdata[goodsignal])
   sxaddpar, h, 'COUNTS', skycts
   sxaddpar, h, 'SKYCTS', skycts, 'Sky counts'
   if sxpar(h,'EXPTIME') gt 0 then $
      sxaddpar, h, 'CTRATE', skycts/sxpar(h,'EXPTIME'), after='SKYCTS', 'Sky counts per second'

   get_date, now
   hist = 'Processed by flatproc '+now
   sxaddhist, hist, h
   
   lastslashpos = strpos(filenames[f],'/',/reverse_search)
   if lastslashpos gt 0 then $
     outname = strmid(filenames[f],0,lastslashpos) + '/' + 'f' + strmid(filenames[f],lastslashpos+1) $
   else  $
     outname = 'f'+filenames[f]

   mwrfits, fdata, outname, h, /create
endfor

end
