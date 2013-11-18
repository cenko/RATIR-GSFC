function polyiterfit, x, y, unc=unc, order=order, niter=niter, sigthresh=sigthresh, yfit=yfit

   if n_elements(niter) eq 0 then niter = order > 2
   
   if n_elements(sigthresh) eq 0 then sigthresh = [6.,3.]
   if n_elements(sigthresh) eq 1 then sigthresh = [6.,sigthresh+fltarr(niter-1)]
   if n_elements(sigthresh) eq 2 and niter gt 2 then sigthresh = max(sigthresh) - findgen(niter)*(findgen(niter)/(niter-1)) * (max(sigthresh)-min(sigthresh))/(niter-1)

   iorder = 0 + (order*indgen(niter)/(niter-1))


   good = where(finite(x)) ; all
   for iter = 0, niter-1 do begin
        par = poly_fit(x[good], y[good], iorder[iter],  measure_errors=unc, yfit=yfit)
        residual = y[good] - yfit
 
        sig = stdev(residual)

        model = poly(x, par) ;  includes flagged indices

        possgood = where(abs(y-model) lt sigthresh[iter]*sig, ctgood)
        if ctgood gt 3 then good = possgood ; else print, 'Warning: flagged too many points.'

   endfor  

   par = poly_fit(x[good], y[good], order, measure_errors=unc)

   yfit = poly(x, par)

   return, par

end

