pro readstandard, name, inwav, flux, unit=unit, path=path, check=check, result=result
   ; Reads in a standard star spectrum from absolute calibration database.
   ; If wav is specified, interpolate flux onto wav.
   ; Returns in microJy.

    ; should probably enable loading just by RA/DEC after all.

   ; check just makes sure the standard is real; it doesn't actually load it.

   ;common specstdpath, path

   if n_elements(path) eq 0 then path = '~/progs/idl/lris/refdata/' ; slash included this time...

   if n_elements(unit) eq 0 then unit = 'flambda/ang'
  
   standardtable = readstandardtable(path+'specstandards.dat')

   if size(name,/tn) ne 'STRING' then print, 'Warning: Standard name is not a string.'
   name = clip(name) ; non-string input should be used with caution

   flux = -1
   result=0

   nstd = n_elements(standardtable)
   for n = 0, nstd-1 do begin
      allnames = [standardtable[n].name]
      if strlen(standardtable[n].altnames) gt 1 then $
         allnames = [allnames, strsplit(standardtable[n].altnames,',',/extract)]
      m = where(name eq allnames, nmatch)
      ;print, clip(nmatch,3), allnames
      if nmatch eq 0 then continue

      file = standardtable[n].file
      indir = (strsplit(file, '/', /extract)) [0] ; the top-level directory 

      result = file_test(path+file) ; just check if the file exists.
      if result eq 0 then begin
         print, 'Standard reference spectrum ', path+file, ' not found.'
         return
      endif
      if keyword_set(check) then return
 
      if indir eq 'ctiostan' or indir eq 'hststan' or indir eq 'okestan' or indir eq 'wdstan' then begin
         readcol, path+file, wav, flam_16, flux_mJy, format='f,f,f', /silent
         flux = flux_mJy * 1000. ; we want uJy
      endif
      if indir eq 'calspec' then begin
         readcol, path+file, wav, flam, format='f,f', /silent
         flux = flam / 1d-8 / 1d-29 / 2.998e10 * (wav*1d-8)^2  ; we want uJy
      endif
      if indir eq 'ing' then begin
         nl = countlines(path+file)
         openr, 10, path+file
         iline = ''
         inunit = ''
         wav = fltarr(nl)
         f = fltarr(nl)
         i = 0L
         while not eof(10) do begin
            readf, 10, iline
            if strlen(iline) lt 3 then continue
            if strmid(iline,0,1) eq '*' then continue
            if strmid(iline,0,3) eq 'SET' then begin
               if strpos(iline,'UNIT') gt 0 then begin
                 p = strpos(iline,'"')
                 inunit = strmid(iline,p+1)
                 p = strpos(inunit,'"')
                 inunit = strlowcase(strmid(inunit,0,p))
               endif
               continue
            endif
            if inunit eq '' then print, 'Unit not recognized?'
            idata = strsplit(iline,/extract)
            wav[i] = idata[0]
            f[i] = idata[1]
            i += 1
         endwhile
         close, 10
         wav = wav[0:i-1] 
         f = f[0:i-1]
         if inunit eq 'ujy' or inunit eq 'micro-janskys' then flux = f
         if inunit eq 'mjy' or inunit eq 'milli-janskys' then flux = f*1000.   ; we want uJy
         if inunit eq 'ab magnitudes' then flux = 10.0^(-(f-23.9)/2.5)
      endif

      break
   endfor
   if keyword_set(check) then return

   if n_elements(flux) eq 1 then begin
      print, 'Standard name '+name+' not recognized.'
      return
   endif

   ; units don't get changed if not recognized which could confuse user
   if unit ne 'uJy' and unit ne 'fnu' and unit ne 'flambda/cgs' and unit ne 'flambda/ang' then begin
     print, '  WARNING:  Unit not recognized; assuming uJy'
     return
   endif
   if unit eq 'uJy'         then flux *= 1
   if unit eq 'fnu'         then flux *= 1d-29
   if unit eq 'flambda/cgs' then flux *= 1d-29 * 2.998e10 / (wav*1d-8)^2
   if unit eq 'flambda/ang' then flux *= 1d-29 * 2.998e10 / (wav*1d-8)^2 * 1d-8

   ;readcol, 'fltt7987.oke', wav, flux, fluxunc, fluxbin skip=13

   ;readcol, path+'vma2.tab', wav, flux, skip=15

   if n_elements(inwav) gt 0 then begin
     oldflux = flux
     flux = interpol(flux, wav, inwav)   ; interpolate it onto wav

     ; Extrapolate values above the maximum assuming RJ law
     wabovemax = where(inwav gt max(wav), ct)
     weqmax = (where(wav eq max(wav))) [0]
     if ct gt 0 then flux[wabovemax] = oldflux[weqmax] * (inwav[wabovemax]/max(wav))^(-4.)

     ; the extrapolation for below max is less clear.  implement this when you actually need to.

   endif else begin
     inwav = wav
   endelse

end

