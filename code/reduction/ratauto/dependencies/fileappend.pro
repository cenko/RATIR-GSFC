function fileappend, filename, prefix, suffix
   ; Adds a prefix or suffix onto a file, even if it has a directory.

   n = n_elements(filename)
   if n gt 1 then begin
      outfilenames = strarr(n)
      for i = 0, n-1 do begin
         outfilenames[i] = fileappend(filename[i],prefix,suffix)
      endfor
      return, outfilenames
   endif else begin

      if n_elements(prefix) eq 0 then prefix = ''
      if n_elements(suffix) eq 0 then suffix = ''

      slashpos = strpos(filename,'/',/reverse_search)
      dotpos = strpos(filename,'.',/reverse_search)
      if dotpos eq -1 then dotpos = strlen(filename) ; unusual to not have a dot
      fileroot = strmid(filename,slashpos+1,dotpos-slashpos-1)

      dir  = strmid(filename,0,slashpos+1)
      file = strmid(filename,slashpos+1,dotpos-slashpos-1)
      ext  = strmid(filename,dotpos)

      return, dir + prefix + file + suffix + ext
   endelse

end

