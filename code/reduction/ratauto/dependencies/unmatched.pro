function unmatched, filelist, prefix
   outfilelist = filelist
   for f = 0, n_elements(outfilelist)-1 do begin
      slashpos = strpos(filelist[f],'/',/reverse_search)
      if slashpos ge 0 then outfilelist[f] = strmid(filelist[f],0,slashpos+1) + prefix + strmid(filelist[f],slashpos+1)
   endfor
   for c = 0, n_elements(outfilelist)-1 do begin
      if outfilelist[c] eq '' then check = '' $
                              else check = findfile(outfilelist[c])
      if check ne '' then outfilelist[c] = ''
   endfor
   unm = where(outfilelist ne '', ctunmatched)
   if ctunmatched eq 0 then return, ''
   return, filelist[unm]
end

