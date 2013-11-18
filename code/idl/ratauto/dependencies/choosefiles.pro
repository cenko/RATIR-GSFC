function repstr1, in, subout, subin
   p = strpos(in,subout)
   l = strlen(subout)
   frag1 = strmid(in, 0,p)
   frag2 = strmid(in, p+l)
   return, frag1+subin+frag2

end

function choosefiles, selector, prefix0, prefix1, prefix2, prefix3
   if n_elements(prefix0) eq 0 then begin
      return, findfile(selector)
   endif
   if n_elements(prefix0) eq 1 then filenamesf0 = findfile(prefix0+selector)
   if n_elements(prefix1) eq 1 then filenamesf1 = findfile(prefix1+selector) 
   if n_elements(prefix2) eq 1 then filenamesf2 = findfile(prefix2+selector) 
   if n_elements(prefix3) eq 1 then filenamesf3 = findfile(prefix3+selector) 

   if n_elements(prefix3) eq 1 then begin
     for i = 0, n_elements(filenamesf3)-1 do $
       if total(repstr1(filenamesf3[i], prefix3, prefix2) eq filenamesf2) then filenamesf3[i] = repstr1(filenamesf3[i], prefix3, prefix2) ;else filenamesf3[i] = ''
     filenamesf2 = unique([filenamesf2, filenamesf3])
   endif
   if n_elements(prefix2) eq 1 then begin
     for i = 0, n_elements(filenamesf2)-1 do $
       if total(repstr1(filenamesf2[i], prefix2, prefix1) eq filenamesf1) then filenamesf2[i] = repstr1(filenamesf2[i], prefix2, prefix1); else filenamesf2[i] = ''
     filenamesf1 = unique([filenamesf1, filenamesf2])
   endif
   if n_elements(prefix1) eq 1 then begin
     for i = 0, n_elements(filenamesf1)-1 do $
       if total(repstr1(filenamesf1[i], prefix1, prefix0) eq filenamesf0) then filenamesf1[i] = repstr1(filenamesf1[i], prefix1, prefix0) ;else filenamesf1[i] = ''
     filenamesf0 = unique([filenamesf0, filenamesf1])
   endif

   dum = where(filenamesf0 ne '', ct)
   if ct eq 0 then return, ''
   filenamesf0 = filenamesf0[where(filenamesf0 ne '')]
   filenum = strarr(n_elements(filenamesf0))
   for i = 0, n_elements(filenamesf0)-1 do begin
     filenum[i] = strmid(filenamesf0[i],strpos(filenamesf0[i],'_')+1)
   endfor
   order = sort(filenum)
   return, filenamesf0[order]
end
