function percentile75, array
   m = median(array,/even)
   w = where(array ge m, ct)
   if ct eq 0 then return, m ; shouldn't ever happen
   return, median(array[w], /even)
end

