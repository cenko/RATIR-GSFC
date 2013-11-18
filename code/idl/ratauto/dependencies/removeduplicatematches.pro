pro removeduplicatematches, xmatch, wavmatch, residual

   uwavmatch = unique(wavmatch)       
   if n_elements(uwavmatch) lt n_elements(wavmatch) then begin
      outmatchi = intarr(n_elements(uwavmatch))
      for iu = 0, n_elements(uwavmatch)-1 do begin
         wm = where(wavmatch eq uwavmatch[iu], ct)
         if ct ge 2 then begin
            outwm = (where(residual[wm] eq min(residual[wm]))) [0]
         endif else begin
            outwm = wm
         endelse
         outmatchi[iu] = wm[outwm]
      endfor
      xmatch = xmatch[outmatchi]
      wavmatch = wavmatch[outmatchi]
   endif

   uxmatch = unique(xmatch)
   if n_elements(uxmatch) lt n_elements(xmatch) then begin
      outmatchi = intarr(n_elements(uxmatch))
      for iu = 0, n_elements(uxmatch)-1 do begin
         wm = where(xmatch eq uxmatch[iu], ct)
         if ct ge 2 then begin
            outwm = (where(residual[wm] eq min(residual[wm]))) [0]
         endif else begin
            outwm = wm
         endelse
         outmatchi[iu] = wm[outwm]
      endfor
      xmatch = xmatch[outmatchi]
      wavmatch = wavmatch[outmatchi]
   endif

end

