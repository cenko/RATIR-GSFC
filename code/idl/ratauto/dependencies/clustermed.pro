function clustermed, l
   if n_elements(l) eq 1 then return, l
   if n_elements(l) le 3 then return, median(l, /even)
   s = l[sort(l)]
   d = s[1:n_elements(s)-1] - s[0:n_elements(s)-2]
   nd = n_elements(d)
   g = nd / 16 > 2  ; sensitive to clusters up to a little less than 1/16 of the data set
   if nd lt 12 then g = 1
   minmean = total(d)
   imean = nd / 2
   for i = 0, nd-1 do begin
       r = [i-g>0,i+g<nd-1]
       m = mean(d[r[0]:r[1]])
       if m lt minmean then begin 
          minmean = m
          imean = i
       endif
   endfor
   md = s[imean] ;+ s[imean+1])/2
   return, md
end


