;+
; NAME:
;	sigclipstats_vt
;
; PURPOSE:
;	Find median, stdev, and mean of sigma clipped data (until max iteration or remaining
;	data is less than 1% of data points or 50)
;
; INPUT:
;	datain - array (must be at least 4 elements for sigma clipping to work)
;
; OPTIONAL KEYWORDS:
;	maxiter - sets maximum iteration, default is 6
;	sigmahi - sets upper sigma limit, default is 3.5
;	sigmalo - sets lower sigma limit, default is -3.5
;
; OPTIONAL OUTPUT KEYWORDS:
;	median  - median of sigma clipped data
;	stdevi  - standard deviation of sigma clipped data
;	meani	- mean of sigma clipped data
;
; EXAMPLE:
;	sigclipmedian, datain, sigmahi=3.5, sigmalo=-3.5, median=medname, stdevi=stdname, meani=meanname
;
; Written by Dan Perley 
; Modified by Vicki Toy 11/18/2013
;-

pro sigclipstats_vt, datain, sigmahi=sigmahi, sigmalo=sigmalo, median=median, stdevi=stdevi, meani=meani, maxiter=maxiter 

	;Sets maximum iteration default and copies data
	data = datain 
	if n_elements(maxiter) eq 0 then maxiter = 6
 
 	;Can't run sigclipmedian on less than 3 items, so return those with medians
 	;Set default sigma values
   	ndata = n_elements(data)
   	if ndata eq 1 then median = data
   	if ndata eq 2 then median = median(data, /even)
   	if ndata le 3 then return
   	if n_elements(sigmahi) eq 0 then sigmahi = 3.5
   	if n_elements(sigmalo) eq 0 then sigmalo = 3.5

	;Find finite values of the data, if there are 0 or 1 end program
	;Rest of program only run on finite data and sort data to make program faster
   	good = where(finite(data), ctgood)
   	if ctgood eq 0 then median = !values.f_nan
   	if ctgood eq 1 then median = data[good[0]]
   	if ctgood lt 2 then return
   	data = data[good]
   	clipdata = data[sort(data)]

   	nremoved = ndata
   	niter = 0
   	
   ;While current iteration is less than the set maximum and number of data points to process 
   ;is greater than 50 or 1% of the data, find values within median sigma limits and iterate
   	while niter lt maxiter and nremoved gt ndata*0.01 < 50 do begin
      	runningmedian = median(clipdata, /even)
      	runningstdev = stdev(clipdata)
      	lolimit = runningmedian - sigmalo * runningstdev
      	hilimit = runningmedian + sigmahi * runningstdev
    	runningn = n_elements(clipdata)
      	pass = where(clipdata ge lolimit and clipdata le hilimit, ct)
      	if ct gt 0 then clipdata = clipdata[pass] else break
      	nremoved = runningn - ct
      	niter = niter + 1
   	endwhile

  	median = median(clipdata, /even)
   	stdevi = stdev(clipdata)
   	meani  = mean(clipdata)
   	
end


