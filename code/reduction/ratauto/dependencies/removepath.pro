function removepath, filename
   n = n_elements(filename) 
   if n gt 1 then begin
      out = strarr(n)
      for i = 0, n-1 do out[i] = removepath(filename[i])
      return, out
   endif else begin
     slashpos = strpos(filename,'/',/reverse_search)
     return, strmid(filename,slashpos+1)
   endelse
end

