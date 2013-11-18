function subtractsky2, spec, ysamplestep=ysamplestep, discontinuity=discontinuity, sky=sky, order=order, debug=debug, sourceprofile=sourceprofile, plotx=plotx

   ; greatly simplified sky subtractor.
  
   inmode =  size(spec, /tn)
   if inmode eq 'STRING' then begin
      filename = spec
      img=mrdfits(filename, 0, header, /silent)
      sxaddpar, header, 'BZERO', 0.
   endif else begin
      img = spec
   endelse

   debug = keyword_set(debug)

   dims = size(img, /dimens)
   nx = dims[0]
   ny = dims[1]

   yref1 = [ny/6, ny/2]
   yref2 = [ny/2, 5*ny/6]
   refspec = fltarr(nx)
   for x = 0, nx-1 do begin
        refspec[x] = 0.5 * sigclipmedian(img[x,yref1[0]:yref1[1]]) + $
                     0.5 * sigclipmedian(img[x,yref2[0]:yref2[1]]) 
   endfor

   xblocksizeapprox = 128
   nxspecblock = round(nx*1.0/xblocksizeapprox)
   xspecblock = ceil(nx*1.0/nxspecblock)
   xspecstart =  findgen(nxspecblock)*xspecblock
   xspecend   = (findgen(nxspecblock)*xspecblock+xspecblock-1) < (nx-1)
   continuum = fltarr(nx) + !values.f_nan
   for xi = 0, nxspecblock-1 do begin
      xrefspec = refspec[xspecstart[xi]:xspecend[xi]]
      cut = 4*abs(median(xrefspec)) + 10
      notlines = where(xrefspec lt cut)
      cut = 3*abs(median(xrefspec[notlines]))  + 10
      notlines = where(xrefspec lt cut)
      cut = 1.5*abs(median(xrefspec[notlines])) + 10
      notlines = where(xrefspec lt cut)
      x = xspecstart[xi]
        if refspec[x] lt cut                         and refspec[x+1] lt cut then continuum[x] = median(img[x,*])
      for x = xspecstart[xi]+1, xspecend[xi]-2 do begin
        if refspec[x] lt cut and refspec[x-1] lt cut and refspec[x+1] lt cut then continuum[x] = median(img[x,*])
      endfor
      x =  xspecend[xi]-1
        if refspec[x] lt cut and refspec[x-1] lt cut                         then continuum[x] = median(img[x,*])

      ;xx = indgen(nx)
      ;fit =  where(finite(continuum) and xx ge xspecstart[xi] and xx le xspecend[xi])
      ;contpar = poly_fit(xx[fit], continuum[fit], 2)
      ;nans = where(finite(continuum) eq 0 and xx ge xspecstart[xi] and xx le xspecend[xi])
      ;continuum[nans] = (poly(xx[nans], contpar))
   endfor
   newcontinuum = continuum ; newcontinuum includes estimates for continuum over the sky lines.  Probably not robust.
   for x = 0, nx-1 do begin
       newcontinuum[x] = median(continuum[x - 11 > 0 : x+11 < nx-1])
       if finite(newcontinuum[x]) eq 0 then newcontinuum[x] = median(continuum[x - 25 > 0 : x+25 < nx-1])
       if finite(newcontinuum[x]) eq 0 then newcontinuum[x] = median(continuum[x - 55 > 0 : x+55 < nx-1])
   endfor 

   ;;print, total(finite(continuum))
   contimg = continuum # (1+fltarr(ny))
   contsubimg = img - contimg
   if debug then writefits, 'continuum.fits', contimg, header
   if debug then writefits, 'contsubimg.fits', contsubimg, header

   sourceimg = img * 0.

   slitfluxx = fltarr(nxspecblock,ny)
   for xi = 0, nxspecblock-1 do begin
      slitfluxx[xi,*] = median( ( (contsubimg[xspecstart[xi]:xspecend[xi],*])), dimension=1) ;* (xspecend[xi]-xspecstart[xi]+1)
      for x = xspecstart[xi], xspecend[xi] do $
         sourceimg[x,*] = slitfluxx[xi,*]
   endfor
   kernel = (1+fltarr(xspecblock))
   sourceimg = convol(sourceimg,kernel/total(kernel), /EDGE_TRUNCATE)
   if debug then writefits, 'sourceimg.fits', sourceimg, header


   ;plot, slitfluxx[xi,*]
   ;print, mean(slitfluxx[xi,*])
   slitflux = total(slitfluxx, 1)*1.0/nxspecblock    ; slitflux is a profile of the avg (over wavelengths) flux under the slit.
                                                     ; only exists as a convenience to the user for identifying the trace.
   ;plot, slitflux, psym=10  
   ;for xi = 0, nxspecblock-1 do begin
   ;    oplot, slitfluxx[xi,*], linestyle=1
   ;endfor
   ;oplot, slitflux, psym=10  


   stdslitflux = sigclipstdev(slitflux)
   medslitflux = sigclipmedian(slitflux)
   
   sx = findgen(n_elements(slitflux))
   ok = where(slitflux gt medslitflux-4*stdslitflux and slitflux lt medslitflux+4*stdslitflux)
   trend = poly_fit(sx[ok], slitflux[ok], 2)
   slitflux = slitflux - poly(sx, trend)      ; remove the linear trend

   stdslitflux = sigclipstdev(slitflux)
   medslitflux = sigclipmedian(slitflux)

   sourcethreshold = medslitflux+3*stdslitflux
   ;print, sourcethreshold
   sourcezone = where(convol(slitflux,[1,1,1]/3.) ge sourcethreshold)
   nonsourcezone = where(convol(slitflux,[1,1,1]/3.) lt sourcethreshold)
   ;plot, slitflux, psym=10, yrange=[-50,150]
   ;oplot, convol(slitflux,[1,1,1]/3.), linestyle=1
   ;oplot, [0.,10000], [1,1]*medslitflux+3.5*stdslitflux

   ; preliminaries done, now proceed with sky subtraction.


   img = img - sourceimg ; note that this is temporary; we'll add it back at the end.
                         ; not subtracting the continuum right now - it shouldn't matter.


   if debug then writefits, 'sourcesubimg.fits', img, header

   sky = fltarr(nx, ny)

   blocksizeapprox = 36 ;50 ;36
   buffer =  30 ;16  ;12

   if n_elements(order) eq 0 then order = 2;3

   if n_elements(discontinuity) eq 0 then discontinuity = [0, ny] $
                                     else discontinuity = [0, discontinuity, ny]

   if n_elements(plotx) gt 0 then begin
       psopen, 'plotx.ps', xsize=11.5, ysize=8, /inches
   endif

   for d = 0, n_elements(discontinuity)-2 do begin ; change this to big b
     nyd = discontinuity[d+1] - discontinuity[d]
     nyblockapprox = round(nyd * 1.0 / blocksizeapprox)
     ysamplestep = ceil(ny*1.0 / nyblockapprox)
     nyblock = ceil(nyd*1.0/ysamplestep)
     offy   =  discontinuity[d] + ysamplestep*indgen(nyblock)

     for x = 0, nx-1 do begin
       for i = 0, nyblock-1 do begin
         ystart = (offy[i] - buffer ) > 0
         ystop = (offy[i]+ysamplestep-1 + buffer)  < (discontinuity[d+1]-1)

         ystartedit = offy[i]
         ystopedit  = (offy[i]+ysamplestep-1)  < (discontinuity[d+1]-1)

         ;if x eq 0 then begin
         ;   print, clip(ystartedit,4), '-', clip(ystopedit,4), '  (',clip(ystart,4),'-', clip(ystop,4),')'
         ;endif

         minsize = ystop-ystart-1
         nssec = where(nonsourcezone ge ystart and nonsourcezone le ystop, ct)
         while ct lt minsize do begin
           if ct lt minsize then begin ; expand to make sure we bridge the gap
              if total(nonsourcezone eq ystart) gt 0 and total(nonsourcezone eq ystop) gt 0 then begin
                  ystart = (ystart - 6) > 0
                  ystop  = (ystop + 6) < ny-1
              endif
              if total(nonsourcezone eq ystart) eq 0 then ystart = (ystart - 3) > 0
              if total(nonsourcezone eq ystop)  eq 0 then ystop  = (ystop + 3) < ny-1
           endif
           nssec = where(nonsourcezone ge ystart and nonsourcezone le ystop, ct)
         endwhile
         ysec = nonsourcezone[nssec]
         ;if x eq 1000 then print, clip(i,4), ': ', ysec, '(,',clip(ct),')'

         section = (img[x,ysec]) [*] ;ystart:ystop  ; the (photon) values
         indices = (indgen(ystop+1)) [ysec] ;  actual y pixel indices
 
         sigthresh = [5.0, 3.5, 3.2, 3.2, 3.0]
         orders = [0, 0, (1 < order), (2 < order), order]
         good = where(indices ge 0 and finite(section)) ; all
         for iter = 0, n_elements(sigthresh)-1 do begin
            par = poly_fit(indices[good], section[good], orders[iter], yfit=yfit)
            residual = section[good] - yfit
 
            sig = stdev(residual)

            model = poly(indices, par) ; model includes even flagged indices

            possgood = where(finite(section) and abs(section-model) lt sigthresh[iter]*sig, ctgood)
            if ctgood gt 3 then good = possgood ; else print, 'Warning: flagged too many points.'

            ;if x eq 3819 then print, clip(ystart,3), '-', clip(ystop,3), '  ', clip(ystartedit,3), '-', clip(ystopedit,3), '  (', clip(iter),')', median(section-model), mean(section-model), clip(n_elements(section)-ctgood,4)
         endfor  
         par = poly_fit(indices[good], section[good], order, yfit=yfit)
         if total(finite(par) eq 0) gt 0 then begin 
              print, x, i, indices[good], section[good]
              stop
         endif
         editindices = (indgen(ystopedit+1)) [ystartedit:ystopedit]
         editmodel = poly(editindices, par)
         sky[x,ystartedit:ystopedit] = editmodel
       endfor
       xsky = sky[x,*]
       sky[x,*] = convol(sky[x,*],transpose(1+fltarr(ysamplestep/3))/(ysamplestep/3), /edge_truncate)

       if n_elements(plotx) gt 0 then if total(x eq plotx) gt 0 then begin
           plot, img[x,0:582], psym=10
           oplot, sky[x,*]
       endif
     endfor
   endfor

  ;plot, img[3819,*], psym=10
  ;oplot, sky[3819,*]

  if debug then writefits, 'sky.fits', sky, header

  sourceprofile=slitflux

  if n_elements(plotx) gt 0 then begin
     psclose
  endif

  if inmode eq 'STRING' then begin
     slashpos = strpos(spec,'/',/reverse_search)
     outname = strmid(spec,0,slashpos+1) + 's' + strmid(spec,slashpos+1)
     print, 'Writing subtracted/sky/source extensions to ', removepath(outname)
     mwrfits, img-sky+sourceimg, outname, header, /create, /silent
     mwrfits, sky, outname, /silent               ; extension 1: sky
     mwrfits, sourceimg, outname, /silent         ; extension 2: source (useful for tracing later; NOT for any sort of data extraction!)
     return, 0
  endif else begin
     return, img-sky+sourceimg
  endelse


end

