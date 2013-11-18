function airmassline, airmass, starid, p
  ; a - airmass
  ; s - star ID number (integer)
  ; p - parameters    :
  ;     p[0]         - airmass term
  ;     p[1:ns]   - zeropoint for each star

return, p[starid+1] + p[0]*(airmass-1)

end




pro jointairmassfit, airmass, star, obszpt, photerr, caterr, $
                     zpt1, zpt1unc, airmassterm, airmasstermunc, airmasstermrange=airmasstermrange, defaultterm=defaultterm, outliers=outliers, $
                     starplot=starplot


if n_elements(starplot) eq 0 then starplot = 'iterfit.ps'
if n_elements(maxterm) eq 0 then maxterm = 100
nobs = n_elements(airmass)

stars = unique(star)
nstars = n_elements(stars)

starid = intarr(nobs)-1
starcaterr = fltarr(nstars)
for i = 0, nstars-1 do begin
   thisstar = where(star eq stars[i], ct)
   starid[thisstar] = i
   starcaterr[i] = (caterr[thisstar]) [0]
endfor

plotmed = median(obszpt)

npar = 1 + nstars
parinfo = replicate({value:0.0, fixed:0, limited:[0,0], tied:'', limits:[0.0,0.0]}, npar)
parinfo[0].value = 0.1
parinfo[1:*].value = median(obszpt)

good = 1 + intarr(nobs)
fit = where(good) ; all points to start with
oldfit = fit

if keyword_set(starplot) then begin
  !p.multi = [0,2,2]
  if nstars gt 4 then !p.multi = [0,3,3]
  if nstars gt 9 then !p.multi = [0,4,3]
  if nstars gt 12 then !p.multi = [0,5,3]
  if nstars gt 15 then !p.multi = [0,5,4]
  if nstars gt 20 then !p.multi = [0,6,4]
  maxstars = !p.multi[1]*!p.multi[2]
  !p.position = 0
  psopen, starplot, xsize=8, ysize=8, /inches
endif


sigmaclip = [12.0,6.0,4.5,3.0,3.0]
niter = n_elements(sigmaclip) ; number of clipping iterations
for tryiter = 0, niter do begin

  ;print, 'ITER ', tryiter
  ;print, 'n = ', n_elements(airmass[fit])

  parinfo.fixed = 0
  for s = 0, nstars-1 do begin
    dum = where(starid[fit] eq s, ctgoodstarpt) ; star with no good data points
    if ctgoodstarpt eq 0 then parinfo[1+s].fixed = 1
  endfor

  if max(airmass[fit]) gt min(airmass[fit]) + 0.1 and n_elements(fit) ge total(parinfo.fixed eq 0) then begin
     ; fit for airmass term, only if there is an actual range of airmasses
     ; and the same stars have been observed on multiple epochs

     par = mpfit2dfun('airmassline',$
                          airmass[fit], starid[fit],$
                          obszpt[fit], photerr[fit], parinfo=parinfo, $
                          perror=parunc, bestnorm=chisq,dof=dof, yfit=yfit $
                          , /quiet)
    fixedairterm = 0
    tryfixtermfit = 0
    if par[0] lt min(airmasstermrange) or par[0] gt max(airmasstermrange) or dof lt 0 then tryfixtermfit = 1
    if n_elements(parunc) gt 0 then begin
       if min(parunc[1:*]) gt 1 then tryfixtermfit = 1
    endif else begin
       tryfixtermfit = 1
    endelse
  endif else begin
    tryfixtermfit = 1
  endelse

  if tryfixtermfit then begin
   ; if the airmass term is not reasonable, fix it to the default value and redo the fit
    ;;;parinfo = replicate({value:0.0, fixed:0, limited:[0,0], tied:'', limits:[0.0,0.0]}, npar)
    parinfo[0].fixed = 1
    parinfo[0].value = defaultterm

    par = mpfit2dfun('airmassline',$
                     airmass[fit], starid[fit],$
                     obszpt[fit], photerr[fit], parinfo=parinfo, $
                     perror=parunc, bestnorm=chisq,dof=dof, yfit=yfit $
                     , /quiet)
    fixedairterm = 1
  endif

  ;print, starid[fit]
  ;print, obszpt[fit]
  ;print, par
  ;print, parunc

  if keyword_set(starplot) then begin
    for s = 0, maxstars-1 do begin
      if s gt nstars-1 then begin
         megaplot, xrange=[1,2.4], yrange=[-3,-2]+0.6
         continue
      endif

      thisstar = where(starid eq s, scts)
      sairmass = airmass[thisstar]
      sstarid  = starid[thisstar]
      sobszpt =  obszpt[thisstar]
      sphoterr = photerr[thisstar]
      sgood = good[thisstar]
      ;sgood[fit] = 1
      megaplot, sairmass, sobszpt, yerr=sphoterr, title=clip(s)+' :  '+stars[s], xrange=[1,2.4], yrange=plotmed+1*[-0.75,+0.75], fill=sgood, /xstyle, /ystyle, size=0.6
      air = [1,3]
      oplot, air, par[1+s]+(air-1)*par[0]
      oplot, air, (par[1+s]-parunc[1+s])+(air-1)*par[0], linestyle=1
      oplot, air, (par[1+s]+parunc[1+s])+(air-1)*par[0], linestyle=1
      ;print, s
      ;print, sairmass
      ;print, sobszpt
      ;print, par[1+s]+(sairmass-1)*par[0]
      ;print, sphoterr
    endfor
  endif


  if tryiter ne niter then begin
     residualsigma = abs(yfit-obszpt[fit])/photerr
     outliers_f = where(residualsigma gt sigmaclip[tryiter], complement=fit_f, ct)
     if ct gt 0 then good[fit[outliers_f]] = 0
     fit = where(good)
  endif
  if total(good) eq 0 then begin
     print, 'WARNING - flagged all stars as outliers!'
     break ; just give up if we flag everything...
  endif

endfor

if keyword_set(starplot) then psclose


;parunc[1:*] = parunc[1:*]  * 0.07/0.047  ; I don't know why this is necessary:  mpfit seems to be underestimating the error.

;print, clip(chisq), '/', clip(fix(dof))
;print, yfit-obszpt[fit]
;print, (yfit-obszpt[fit])/photerr[fit]
;print, median(abs(yfit-obszpt))

if dof gt 9 then begin
   if chisq / dof gt 1.2 then begin; (dof + 3*sqrt(1.0*dof))/dof then begin
     parunc = parunc * sqrt(chisq/dof)
   endif
endif

airmassterm = par[0]
airmasstermunc = parunc[0]

szpt1 = par[1:*]
szpt1unc = sqrt(parunc[1:*]^2 + starcaterr^2)


; Finally, do a weighted mean to get the average zeropoint
goodstar = szpt1 gt -50 and szpt1unc lt 50 and parinfo[1:*].fixed eq 0

if total(goodstar) gt 1 then begin
  for tryiter = 0, niter do begin

    zpt1 =     total(goodstar * szpt1 / szpt1unc^2) / total(goodstar * 1./szpt1unc^2)
    szpt1residual = abs(szpt1-zpt1)/(szpt1unc)

    if tryiter ne niter then begin
       residualsigma = abs(szpt1-zpt1)/szpt1unc
       badi = where(residualsigma ge sigmaclip[tryiter], complement=goodi, ct)
       if ct gt 0 then goodstar[badi] = 0
    endif
    ;print, fpr(residualsigma,3.1)
    ;print, goodstar

  endfor
endif else begin
  zpt1 = szpt1[(where(goodstar)) [0]]
endelse
for s = 0, nstars-1 do begin
   if goodstar[s] eq 0 then good[where(starid eq s)] = 0
endfor

outliers = where(good eq 0)


chisqstar = total((szpt1 - zpt1)^2 * goodstar / starcaterr^2)
dofstar = total(goodstar)-1

zpt1unc =  1./sqrt(total(goodstar * 1./szpt1unc^2)) ; treating this as sqrt statistics, should look into mars fitting here...

if dofstar gt 6 then begin
  if chisqstar/dofstar gt 1.2 then begin ;(dofstar + 2*sqrt(1.0*dof))/dof then begin
    zpt1unc = zpt1unc * sqrt(chisqstar/dofstar)
  endif
endif


  ;print, clip(chisqstar), '/', clip(dofstar)
  ;print, par[1:*] - zpt1
  ;print, (par[1:*] - zpt1)/starcaterr
  ;print, median(abs(par[1:*] - zpt1))

end



pro testjointairmassfit

  niter = 1000L
  allzpt1 = fltarr(niter)
  allairmassterm = fltarr(niter)
  allzpt1unc = fltarr(niter)
  allairmasstermunc = fltarr(niter)

  airmass = [1.1, 1.1, 1.3, 1.3, 1.3, 1.8, 1.8]  ;  1.35, 1.35, 1.5, 1.6, 1.6,
  starname= ['A', 'B', 'A', 'B', 'C', 'D', 'A']  ; 'B', 'C',  'A', 'A', 'B', 
  photerr = 0.10 + fltarr(n_elements(airmass))
  caterr  = 0.10 + fltarr(n_elements(starname))

  for i = 0, niter-1 do begin

   obszpt  = -1 + 0.1*(airmass-1)                                  ; the actual values (perfect measurement)
   obszpt  = obszpt + photerr*randomn(seed, n_elements(airmass))   ; add in photometric error
   stars = unique(starname)
   for s = 0, n_elements(unique(stars))-1 do begin
      thisstar = where(starname eq stars[s])
      obszpt[thisstar] = obszpt[thisstar] + caterr[thisstar]*randomn(seed)     ; add in catalog error
   endfor

   jointairmassfit, airmass, starname, obszpt, photerr, caterr, zpt1, zpt1unc, airmassterm, airmasstermunc

   if i lt 10 then print, fpr(zpt1,2.3), ' (+/- ', fpr(zpt1unc,1.3), ')  + ', fpr(airmassterm,2.3),  ' (+/- ', fpr(airmasstermunc,1.3), ')  *(airmass-1)'

   allzpt1[i] = zpt1
   allzpt1unc[i] = zpt1unc
   allairmassterm[i] = airmassterm
   allairmasstermunc[i] = airmasstermunc

  endfor

  print, 'Mean over all iterations:'
  print, fpr(mean(allzpt1),2.3), ' (+/- ', clip(mean(allzpt1unc),5), ',  stdev', fpr(stdev(allzpt1),2.3),')  + ', fpr(mean(allairmassterm),2.3),  ' (+/- ', clip(mean(allairmasstermunc),5), ',  stdev ', fpr(stdev(allairmassterm),2.3),')  *(airmass-1)'
  print, 'Median over all iterations:'
  print, fpr(median(allzpt1),2.3), ' (+/- ', fpr(median(allzpt1unc),1.3), ', cstdev', fpr(sigclipstdev(allzpt1),2.3),')  + ', fpr(median(allairmassterm),2.3),  ' (+/- ', fpr(median(allairmasstermunc),1.3), ', cstdev ', fpr(sigclipstdev(allairmassterm),2.3),')  *(airmass-1)'


end

pro testweightmean

  niter = 1000L
  nmeas = 20  
  weightmeanarr = fltarr(niter)
  weightmeanuncarr = fltarr(niter)

  sigma = 1. + fltarr(nmeas)
  for i = 0, niter-1 do begin
     y = sigma * randomn(seed, nmeas)
     
     weightmeanarr[i] = total(y / sigma^2) / total(1./sigma^2)
     weightmeanuncarr[i] = sqrt(1. / total(1./sigma^2))

     if i lt 10 then print, weightmeanarr[i], weightmeanuncarr[i]
  endfor

  print
  print, median(weightmeanarr), stdev(weightmeanarr), median(weightmeanuncarr)

end



