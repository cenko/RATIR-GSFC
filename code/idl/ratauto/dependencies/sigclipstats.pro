pro sigclipstats, datain, sigmahi=sigmahi, sigmalo=sigmalo, median=median, stdevi=stdevi, meani=meani
   data = datain ; make a copy!
   maxiter = 6
   endclipfrac = 0.01
 
   ndata = n_elements(data)
   if ndata eq 1 then median = data
   if ndata le 3 then median = median(data, /even)
   if ndata le 3 then return
   if n_elements(sigmahi) eq 0 then sigmahi = 3.5
   if n_elements(sigmalo) eq 0 then sigmalo = 3.5

   data = data[where(finite(data))]
   clipdata = data[sort(data)] ; the going theory is that this speeds things up

   nremoved = ndata
   niter = 0
   while nremoved gt (ndata*0.01 < 50) and niter lt maxiter do begin
      runningmedian = median(clipdata)
      runningstdev = stdev(clipdata)
      lolimit = runningmedian - sigmalo * runningstdev
      hilimit = runningmedian + sigmahi * runningstdev
      runningn = n_elements(clipdata)
      pass = where(clipdata ge lolimit and clipdata le hilimit, ct)
      if ct gt 0 then clipdata = clipdata[pass]
      nremoved = runningn - ct
      niter = niter + 1
   endwhile

   median = median(clipdata, /even)
   stdevi = stdev(clipdata)
   meani = mean(clipdata)
end

