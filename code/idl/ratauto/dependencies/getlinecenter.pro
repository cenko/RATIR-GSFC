function getlinecenter, pixin, fluxin, xin, fitrad=fitrad, unc=unc
      if n_elements(xin) eq 0 then begin
          flux = pixin
          x = fluxin
          pix = findgen(n_elements(flux))
      endif else begin
          pix = pixin
          flux = fluxin
          x = xin
      endelse
   
      zone = where(pix ge x-fitrad and pix le x+fitrad)
      pix = pix[zone]
      flux = flux[zone]

      ; old mpfit method
      ; parinfo = replicate({value:0.0, fixed:0, limited:[0,0], limits:[0.,0.]}, 4)
      ; parinfo[0].value = max(flux)
      ; parinfo[1].value = x
      ; parinfo[2].value = fitrad/3.
      ; parinfo[3].value = 0            ; I turned background fitting off; why?
      ;  parinfo[3].fixed = 1
      ;
      ; sol = mpfitfun('gaussian',  pix, flux, sqrt((flux + 10.) > 10), parinfo=parinfo,$
      ;               perr=perr, bestnorm=chisq, dof=dof, yfit=yfit, /quiet)
      

       estimates = [max(flux), x, fitrad/3.]

       sol = gaussfit(pix, flux, measure_errors=sqrt((flux + 10.) > 10), nterms=3, sigma=perr)

       unc = perr[1]     

     return, sol[1]

end

