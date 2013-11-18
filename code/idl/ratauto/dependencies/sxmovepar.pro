pro sxmovepar, header, key, before=before, after=after

   val = sxpar(header, key, count=count, /silent)
   if count eq 0 then return
   sxdelpar, header, key
   if n_elements(after) gt 0 then begin 
      sxaddpar, header, key, val, after=after
   endif else if n_elements(before) gt 0 then begin
      sxaddpar, header, key, val, before=before
   endif else begin
      sxaddpar, header, key, val ; put at end
   endelse

end

