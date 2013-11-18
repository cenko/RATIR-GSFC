pro findcroplines, arrayin, xcrop, ycrop, threshold=threshold

   if n_elements(threshold) eq 0 then threshold = 0.1   

   xcrop = [0, -1]
   ycrop = [0, -1]

   array = arrayin
   nans = where(finite(array) eq 0, ctnans)
   if ctnans gt 0 then array[nans] = 0.

   for dim = 1, 2 do begin
      ztot = total(array, dim) ; sum over x-axis
      naxis = n_elements(ztot)
      z = findgen(n_elements(ztot))
      ztotmean = mean(ztot)
      ztothizone = where(ztot ge ztotmean*threshold) ; works OK with NaNs
      for i = 1, n_elements(ztothizone)-2 do begin
         ii = ztothizone[i]
         ; basically, demand the pixels on left and right also be "hi" to permit drawing the line here.
         ; this avoids sporadic bad columns, etc.
         if z[ii] eq z[ii-1]+1 and z[ii] eq z[ii+1]-1 then begin
            minz = ii
            break
         endif
      endfor
      for i = n_elements(ztothizone)-2, 1, -1 do begin
         ii = ztothizone[i]
         if z[ii] eq z[ii-1]+1 and z[ii] eq z[ii+1]-1 then begin
            maxz = ii
            break
         endif
      endfor
      if minz eq 1       then minz = 0       else minz += 1  ; crop one more to compensate for our safety
      if maxz eq naxis-2 then maxz = naxis-1 else maxz -= 1
      if dim eq 1 then ycrop = [minz, maxz]
      if dim eq 2 then xcrop = [minz, maxz]

   endfor

end

