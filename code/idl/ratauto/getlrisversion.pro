function getlrisversion, filename, datadir=datadir, quiet=quiet
   ; Look for an example file and find the date to infer the camera and file format
   ; (MJD upgrade dates are approximate)

   if n_elements(filename) eq 0 then begin
      if n_elements(datadir) eq 0 then datadir = './'
      files = (findfile(datadir+'r??????_????'+'.fits'))
      if files[0] eq '' then files = (findfile(datadir+'b??????_????'+'.fits'))
      if files[0] eq '' then files = (findfile(datadir+'*r??????_????'+'.fits')) 
      if files[0] eq '' then files = (findfile(datadir+'*b??????_????'+'.fits')) 
      if files[0] eq '' then files = (findfile(datadir+'lred????'+'.fits')) 
      if files[0] eq '' then files = (findfile(datadir+'lblue????'+'.fits')) 
      if files[0] eq '' then files = (findfile(datadir+'*lred????'+'.fits')) 
      if files[0] eq '' then files = (findfile(datadir+'*lblue????'+'.fits')) 
   endif else begin
      files = (findfile(filename))
   endelse

   if n_elements(files) eq 0 then return, 0   ; 0 = no files found

   vers = intarr(3)
   for i = 0, 2 do begin
     if i eq 0 then try = files[n_elements(files)/2]
     if i eq 1 then try = files[n_elements(files)-1]
     if i eq 2 then try = files[0]
     if try ne '' then begin
        h = headfits(try[0])
        mjd = float(sxpar(h,'MJD-OBS'))
        if mjd gt 0 and mjd lt 54952 then vers[i] = 1
        if mjd ge 54952 and mjd lt 55562 then vers[i] = 2
        if mjd ge 55562 then vers[i] = 3
     endif 
   endfor

   if max(vers) le 0 then begin
     if keyword_set(quiet) eq 0 then print, 'Cannot determine LRIS version.'  ; -1 = unknown
     return, -1
   endif
   vers = vers[where(vers gt 0)]
   if max(vers)-min(vers) ne 0 then begin
     print, 'Conflicting LRIS version numbers!'
     print, 'This directory seems to include data from multiple runs with different red CCDs.'
     return, -2                               ; -2 = conflict
   endif

   return, vers[0]

end

