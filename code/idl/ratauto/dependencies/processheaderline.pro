function processheaderline, inline

   ; process a single line in a FITS-like ASCII header

   eqpos = strpos(inline, '=')
   if eqpos le 0 then return, ''

   hline = ''
   hline += clip(strmid(inline,0,eqpos),8)
   hline += '= '
   slpos = strpos(inline, ' /', /reverse_search) ; comment separator
   if slpos lt eqpos then slpos = strlen(inline)

   val = clip(strmid(inline,eqpos+1,slpos-(eqpos+1)))
   if valid_num(val) eq 0 then begin
       ; enclose string values in quotes for full compatibility with sxpar
       if val eq '' then begin
          val = "''"
       endif else begin
         if strmid(clip(val),0,1) ne "'" then val = "'" + val
         if strmid(clip(val),strlen(clip(val))-1,1) ne "'" then val = val + "'"
       endelse
   endif 
   hline += val

   if slpos lt strlen(inline) then begin
      comment = clip(strmid(inline,slpos+2))
      if strlen(clip(comment)) gt 1 then hline += ' / ' + comment
   endif

   return, hline

end

