pro sxrenamepar, header, key, newname
   if n_elements(newname) eq 0 then begin
     print, 'New name?'
     return
   endif

   if key eq newname then return
   hold = sxpar(header, key, count=count, /silent)
   if count eq 0 then return
 
   sxaddpar, header, newname, hold, after=key
   sxdelpar, header, key

   
end

