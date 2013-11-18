;+
; NAME:
;   specextract.pro
;
; PURPOSE:
;   Fully-automated quick reduction of LRIS red-side spectra
;
; CALLING SEQUENCE:
;   specextract, filename, ycent, $
;   ysearch=ysearch, xinc=xinc, tracefile=tracefile, traceycent=traceycent, fixtrace=fixtrace, $
;   arcfile=arcfile, arclist=arclist, brightarclist=brightarclist, skyfile=skyfile, skylist=skylist, $
;   aprad=aprad, apsharp=apsharp, bakdist=bakdist, bakwidth=bakwidth, yoff=yoff, $
;   cwav=cwav, dispersion=dispersion, uncdisp=uncdisp, uncrefwav=uncrefwav, initsol=initsol, maxnudge=maxnudge, $
;   outfile=outfile, keywords=keywords
;
;
; INPUTS:
;     I/O parameters
;   filename -   The 2D spectrum.
;   outfile  -   Output file for the extracted 1D spectrum.
;
;     Tracing parameters
;   ycent      - The y coordinate (+-few pix, in middle of image) of the object to extract
;   ysearch    - Region (possibly large) to search for the object to trace
;   xinc       - Step to use in tracing.
;   tracefile  - External file to use for trace instead of current file (use yoff for offset)
;   traceycent - The y coordinate (+-few pix) of object to trace (if different from object to extract)
;   fixtrace   - Flag, set to use exact ycent (instead of searching nearest few pixels)
;
;     Extraction parameters
;   aprad    -   Radius from the central trace to extract (half width)
;   apsharp  -   Sharpness parameter for extraction kernel:  2=gaussian, 6=default, infinity=square
;   bakdist  -   Optional background extraction separation.  default=7/4xaprad if not set.
;   bakwidth -   Background separation width (full width).  Set to 0 to suppress background subtraction.
;   yoff     -   Additional extraction center offset in pixels from best-fit trace center
;
;     Wavelength solution parameters
;   arcfile  -   The 2D spectrum of an arc.  
;   arclist  -   List of arc lines to use in the wavelength solution
;   skyfile  -   2D spectrum of the sky for wavelength refinement.  
;                If an integer, use the extension of the 2D spectrum
;   skylist  -   List of skylines to use in the refined wavelength solution
;   dispersion - Guess at the dispersion of the spectrum (Ang/pix)
;   initsol   -  Guess at the full wav(x) polynomial, optional
;   maxnudge  -  Maximum change in central wavelength possible from initsol
;  
;     Other parameters
;   keywords -   List of FITS header keywords to preserve in output file
;   success  -   Returns a negative flag if procedure failed, or 1 on success
;
;
;-
;------------------------------------------------------------------------------



function findobject, img, ysearch=ysearch, apmin=apmin, apmax=apmax

      s = size(img)
      nx = s[1]
      ny = s[2]

      ; blind search algorithm not rigorously tested for faint objects

      if n_elements(ysearch) ge 2 then begin
         ysearchmin = ysearch[0]
         ysearchmax = ysearch[1]
      endif
      if n_elements(ysearch) eq 0 then begin
         ysearchmin = 0  + 2
         ysearchmax = ny - 2
      endif
      ysearchmin = round(ysearchmin)
      ysearchmax = round(ysearchmax)

      xmin = (nx*4)/10
      xmax = (nx*6)/10

      print, '  Searching for trace between y=',clip(ysearchmin),'-',clip(ysearchmax), '  (x=', clip(xmin),'-',clip(xmax),')'

      wavsum = fltarr(ysearchmax-ysearchmin)
      for y = ysearchmin, ysearchmax-1 do begin
         wavsum[y-ysearchmin] = median(img[xmin:xmax,y]) 
      endfor
      kernel = [0.2,1.0,1.5,1.0,0.2]
      if n_elements(kernel) gt n_elements(wavsum) then kernel = [0.3,1.0,0.3]
      if n_elements(kernel) gt n_elements(wavsum) then kernel = [1.0]
      kernel = kernel / total(kernel)
      wavsumconv = convol(wavsum, kernel, /edge_truncate)

      centersigma = max(wavsum)/stdev(wavsum)   ; stdev is a huge overestimate for a bright star being present.
      if centersigma lt 2 then begin
         print, '  WARNING:  Cannot find the source!   Best center is only ', centersigma, '-sigma above background.'
         print, '  Continuing anyway, but you will likely want to rerun after specifying ysearch.'
      endif

      avcenty = ysearchmin + (where(wavsumconv eq max(wavsumconv))) [0] ; max val after convolution


      profmax = max(wavsum)
      baseline = median(wavsum, /even)
      for y = avcenty, ysearchmax-1 do begin
         if wavsum[y-ysearchmin] lt baseline + (profmax-baseline)*0.2 then begin
            apmax = y
            break
         endif
      endfor
      for y = avcenty, ysearchmin-1, -1 do begin
         if wavsum[y-ysearchmin] lt baseline + (profmax-baseline)*0.2 then begin
            apmin = y
            break
         endif
      endfor
     
      if n_elements(apmax) eq 0 then apmax = ysearchmax
      if n_elements(apmin) eq 0 then apmin = ysearchmin

      if (apmax-avcenty) gt 6 then apmax = avcenty + 5
      if (avcenty-apmin) gt 6 then apmin = avcenty - 5

      print, '  Found object at y=', strtrim(avcenty,2), ' (+', clip(apmax-avcenty), ' -', clip(avcenty-apmin), ')'

   return, avcenty

end

function traceobject, sci, avcenty=avcenty, xinc=xinc
      ; add avcentx to allow starting further from middle?

      ; trace the object across the chip
   
      ; march across the image, left to right, in blocks of 100, summing up horizontally for every y to get 
      ; the variable spatial profile

      debug = 0
      order = 3
      tracecheckplot = 0

      imsize = size(sci)
      nx = imsize[1]
      ny = imsize[2]

      xmax = nx-1

      ni = fix(xmax/xinc)
      centy = intarr(ni)
      centx = intarr(ni)
      centsig = fltarr(ni)
      isedge = fltarr(ni)
      isgood = intarr(ni)
      initscansize = 3
      scansize = 2  ; how far away from the most recently-found center to look for trace variation
      wavsum = fltarr(ny)  ; fltarr(scansize*2+1)
      lastgoodtracei = 0
      ntraceloss = 0
      ntracefound = 0
      lastgoodtracex = -1
      lastgoodtracey = avcenty
      kernel = [0.1,0.2,0.4,0.2,0.1]
      kernel = kernel / total(kernel)

      starti = fix(ni/2)
      direction = 1
      i = starti
      for j = 0, ni-1 do begin        ; now : starts in the middle
         if j eq 0 then begin
            lasty = avcenty
            rady = initscansize
         endif else begin
            i += direction
            rady = scansize
         endelse
         if i eq ni then begin
            direction = -direction
            i = starti-1
            lasty = firsty
         endif
         
         ymin = (lasty-rady) > 0
         ymax = (lasty+rady) < (ny-1)

         x1 = fix(i*xinc)
         x2 = fix(x1+xinc-1) < xmax

         ; get the maximum (block-summed) pixel and significance
         for y = 0, ny-1 do begin
            wavsum[y] = total(sci[x1:x2,y])  ; sum across xinc pix within the block
         endfor
         bad = where(finite(wavsum) eq 0, ct)
         if ct gt 0 then wavsum[bad] = 0.
         djs_iterstat, wavsum, mean=fullymean, median=fullymedian, sigma=noiserms
         wavsum -= fullymedian
         wavsumsearch = wavsum[ymin:ymax]
         if n_elements(kernel) lt n_elements(wavsumsearch) then $
            wavsumsearchconv = convol(wavsumsearch, kernel, /edge_truncate) $
         else $
            wavsumsearchconv = wavsumsearch
         maxval = max(wavsumsearchconv)
         significance = maxval/noiserms

         ; get the y-coordinate of maximum pixel
         ymaxpix = ymin + (where(wavsumsearchconv eq maxval)) [0]
         if ymaxpix eq ymin or ymaxpix eq ymax then isedge[i] = 1

         ; get the refined centroid
         rymin = (ymaxpix-2) > 0
         rymax = (ymaxpix+2) < (ny-1)
         rwavsum = wavsum[rymin:rymax]
         ycentroid = total((findgen(rymax-rymin+1)+rymin) * rwavsum) / total(rwavsum)
         centxi = i*xinc + xinc/2
 
         ;oplot, yarr, wavsumconv

         if debug gt 0 then print, '  ',clip(i,3),' ', rclip(x1,4)+'-'+clip(x2,4),rclip(ymin,4),'-',clip(ymax,4),  ymaxpix, ycentroid, significance

         if isedge[i] eq 0 then lasty = ymaxpix
         if j eq 0 then firsty = lasty
         
         centx[i] = centxi
         if ycentroid gt rymin and ycentroid lt rymax then $
            centy[i] = ycentroid $
         else $
            centy[i] = ymaxpix
         centsig[i] = significance
      endfor

      sigthresh = 2.5
      if max(centsig) lt 5  then sigthresh = 2.0
      if max(centsig) gt 30 then sigthresh = 3.5
      if max(centsig) gt 50 then sigthresh = 5.0

      goodtrace = where(centsig gt sigthresh, ct)
      if ct lt 3 then begin
          print, '  Failed to trace - insufficient signal.'
          if tracecheckplot then psclose
          return, -1
      endif
      if ct lt n_elements(centsig)/2 then order = 1

      ; need to make this plot more useful - specific file output, write pixels, etc.

      if tracecheckplot then begin
         psopen, 'traceplot.ps', xsize=7, ysize=5, /inch
         device, /color
         colors = transpose([[0,0,0],$
                          [180,180,210],$ 
                          [64, 64, 192],$ 
                          [64, 136,128],$ 
                          [64, 158, 64],$  
                          [112,112, 64],$ 
                          [192, 64, 64],$
                          [128, 64, 128]$
                       ])
          tvlct, colors

          plot, centx, centy, /ynozero
          ;oplot, centx, centy2, linestyle=1   centy2 is not reliable 
       endif

      ; fit a polynomial to the profile, excluding possible bad traces



      tracepar = polyiterfit(centx[goodtrace], centy[goodtrace],order=order,niter=2)   ; this should be updated
      print, '  Trace function: '      , format='($,A)'
      for o = 0, n_elements(tracepar)-1 do begin
         print, clip(tracepar[o]), format='($,A)'
         if o eq 1 then print, 'x', format='($,A)'
         if o ge 2 then print, 'x^'+clip(o), format='($,A)'
         if o ne n_elements(tracepar)-1 then print, ' + ' , format='($,A)'
      endfor
      print
      if debug gt 0 then begin
         for x = 0, xmax, 250 do begin
            print, x, poly(x, tracepar)
         endfor
      endif

      tracey = poly(findgen(nx), tracepar)

      ; should consider doing the linear extrapolation method for regions beyond the end of the "goodtrace"

      if tracecheckplot then begin
         if tracecheckplot then oplot, centx[goodtrace], centy[goodtrace], psym=7

         oplot, findgen(nx), tracey, linestyle=2

         psclose
      endif
      return, tracey
 
end



function extract, img, tracey, aprad, apsharp, bakdist, bakwidth, fbak=fbak, regfile=regfile, saturate=saturate, nsaturate=nsaturate

      imsize = size(img)
      nx = imsize[1]
      ny = imsize[2]

      baksharp = 12.d

      regstep = 40

      nsaturate = 0

      if n_elements(regfile) gt 0 then begin
         openw, 2, regfile
         printf, 2, '# Region file format: DS9 version 4.1'
         printf, 2, 'global color=green'
         printf, 2, 'image'
      endif

      nans = where(finite(img) eq 0, ct)
      if ct gt 0 then img[nans] = 0.

      ; extract the object (and the background, if desired)
      fobj = fltarr(nx)
      fbak = fltarr(nx)
      for x = 0, nx-1 do begin
         cy = tracey[x] 

         if n_elements(regfile) gt 0 then begin
         if x mod regstep eq 0 then begin
             if x gt 0 then begin
                ; the +1 is to offset IDL/data frame (first cell 0) to DS9 frame (first cell 1)
                printf, 2, 'line('+clip(1+prevx)+','+clip(1+prevy)+','+clip(1+x)+','+clip(1+cy)+') # line=0 0 dash=1'
                printf, 2, 'line('+clip(1+prevx)+','+clip(1+prevy+aprad)+','+clip(1+x)+','+clip(1+cy+aprad)+') # line=0 0'
                printf, 2, 'line('+clip(1+prevx)+','+clip(1+prevy-aprad)+','+clip(1+x)+','+clip(1+cy-aprad)+') # line=0 0'
                if n_elements(bakwidth) gt 0 and n_elements(bakdist) gt 0 then if bakwidth gt 0. then begin
                 printf, 2, 'line('+clip(1+prevx)+','+clip(1+prevy+bakdist)+','+clip(1+x)+','+clip(1+cy+bakdist)+') # line=0 0 color=yellow'
                 printf, 2, 'line('+clip(1+prevx)+','+clip(1+prevy+bakdist+bakwidth)+','+clip(1+x)+','+clip(1+cy+bakdist+bakwidth)+') # line=0 0 color=yellow'
                 printf, 2, 'line('+clip(1+prevx)+','+clip(1+prevy-bakdist)+','+clip(1+x)+','+clip(1+cy-bakdist)+') # line=0 0 color=yellow'
                 printf, 2, 'line('+clip(1+prevx)+','+clip(1+prevy-bakdist-bakwidth)+','+clip(1+x)+','+clip(1+cy-bakdist-bakwidth)+') # line=0 0 color=yellow'
               endif
               prevx = x
               prevy = cy
            endif else begin
               prevx = x
               prevy = cy
            endelse
         endif
         endif


         eky = indgen(ny)
         extractionkernel = exp(-abs(((eky-cy)*1.0d/aprad)^apsharp)/2.)
         if n_elements(bakwidth) gt 0 and n_elements(bakdist) gt 0 then if bakwidth gt 0. then $
            bakkernel = exp(-abs(((eky-(cy-(bakdist+bakwidth/2.)))*1.0d/(bakwidth/2.))^baksharp)/2.)  $
                      + exp(-abs(((eky-(cy+(bakdist+bakwidth/2.)))*1.0d/(bakwidth/2.))^baksharp)/2.)

        ; if x eq 1232 then begin
        ;    window, 2
        ;    plot, extractionkernel, xrange=[cy-bakdist-3*bakwidth/2,cy+bakdist+3*bakwidth/2], thick=3, /xstyle
        ;    if bakwidth gt 0. then oplot, bakkernel, linestyle=1, thick=2
        ;    oplot, sci[x,*]/(median(sci[x,cy-20:cy+20])+4*stdev(sci[x,cy-20:cy+20])), psym=10
        ; endif

         if n_elements(saturate) gt 0 then begin
            saturated = where(img[x,*] gt saturate and extractionkernel gt 0.01*max(extractionkernel), ctsaturated)
            nsaturated += ctsaturated
         endif

         fobj[x] = total(img[x,*]*extractionkernel)
         if n_elements(bakwidth) gt 0 and n_elements(bakdist) gt 0 then if bakwidth gt 0. then $
            fbak[x] = total(img[x,*]*bakkernel) * total(extractionkernel)/total(bakkernel)
      endfor

  if n_elements(regfile) then close, 2

  return, fobj
end


pro specextract, filename, ycent, $
   ysearch=ysearch, xinc=xinc, tracefile=tracefile, traceycent=traceycent, fixtrace=fixtrace, $
   arcfile=arcfile, arclist=arclist, brightarclist=brightarclist, skyfile=skyfile, skylist=skylist, $
   aprad=aprad, apsharp=apsharp, bakdist=bakdist, bakwidth=bakwidth, yoff=yoff, $
   cwav=cwav, dispersion=dispersion, unccwav=unccwav, uncdisp=uncdisp, uncrefwav=uncrefwav, initsol=initsol, $
   outfile=outfile, keywords=keywords, maxorder=maxorder, maxnudge=maxnudge, success=success

   ;if n_elements(arclist) eq 0 then arclist = '~/research/redux/arc/HgArNegood.list'
   ;if n_elements(brightarclist) eq 0 then brightarclist = '~/research/redux/arc/HgNeArbright.list'
   ;if n_elements(skylist) eq 0 then skylist = '~/research/Useful/skyspec/gident_most.dat'

   if n_elements(outfile) eq 0 then outfile = repstr(filename,'.fits','.spec')

   filefail = 0
   checkfiles = [filename]
   if n_elements(arcfile) then checkfiles = [checkfiles, arcfile]
   if n_elements(arclist) eq 1 then checkfiles = [checkfiles, arclist]
   if n_elements(brightarclist) eq 1 then checkfiles = [checkfiles, brightarclist]
   if n_elements(skyfile) gt 0 then begin
      skyfile = clip(skyfile)
      if skylist ne '' then checkfiles = [checkfiles, skylist]
      if strlen(skyfile) gt 2 then begin
         checkfiles = [checkfiles, skyfile]
      endif else begin
         ; need to check that the extension exists, not clear how to do this efficiently
      endelse
   endif
   for c = 0, n_elements(checkfiles)-1 do begin
      if file_test(checkfiles[c]) eq 0 then begin
         print, '  Cannot find ', checkfiles[c]
         filefail += 1
      endif
   endfor
   if filefail gt 0 then begin 
      print
      success = -1 ; couldn't find input file
      return
   endif


   ; --- Read in file ---

   sci = mrdfits(filename,0,h, /silent)

   imsize = size(sci)
   nx = imsize[1]
   ny = imsize[2]
   if nx eq 0 or ny eq 0 then begin
     print, '  Not a valid 2D FITS spectrum.'
     print
     success = -2 ; bad file
     return
   endif

   

   ; --- Read in trace file or trace function ---

   if n_elements(xinc) eq 0 then xinc = nx/16.

   if n_elements(tracefile) gt 0 then begin
      if strpos(tracefile,'.fit') gt 0 then begin
         traceimg = mrdfits(tracefile,0,traceh, /silent)
      endif else begin
         ntl = countlines(tracefile)
         if ntl le 2 then begin
            openr, 3, tracefile
            inline = ''
            readline, 3, inline
            close, 3
            tracepar = float(strsplit(inline,/extract))
            tracey = poly(fltarr(nx), tracepar)
         endif else begin
            rdfloat, tracefile, tracey
         endelse
      endelse
   endif else begin
      traceimg = sci  ; (science image *is* the trace image unless specified)
   endelse


   ; --- Find and trace tracing object, if necessary ---

   if n_elements(traceimg) gt 0 then begin

      if n_elements(traceycent) eq 0 and n_elements(ycent) gt 0 and n_elements(tracefile) eq 0 then traceycent = ycent
      if n_elements(traceycent) eq 0 and n_elements(ysearch) eq 1 then traceycent = ysearch

      if n_elements(traceycent) gt 0 then stop
      if n_elements(traceycent) eq 0 then traceycent = findobject(traceimg,ysearch=ysearch, apmin=yfindmin, apmax=yfindmax)
      if n_elements(ycent) eq 0 then ycent = traceycent

      if n_elements(aprad) eq 0 then begin
          aprad = min([yfindmax-traceycent,traceycent-yfindmin])
          if aprad lt 2. then aprad = 2. ; need to work on improving the aperture determiner in findobject.
          print, '  Using aperture radius ', clip(aprad)
      endif else begin
          print, '  Radius was set to ', (aprad)
      endelse

      tracey = traceobject(traceimg, avcenty=traceycent, xinc=xinc)
      if n_elements(tracey) eq 1 then begin
        if tracey eq -1 then begin
           print
           success = -3 ; trace failure
           return
        endif
      endif
   endif 

   if n_elements(aprad) eq 0 then aprad = 12
   if n_elements(apsharp) eq 0 then apsharp = 6
   if n_elements(bakdist) eq 0 then bakdist = aprad*7./4
   if n_elements(bakwidth) eq 0 then bakwidth = aprad*4.
     ; formerly didn't set these if a skyfile was given, but it's good to subtract background even with sky by default
     ; since there might be galaxy, twilight gradients, etc.


   ; --- Get optimum trace y-offset, if necessary ----

   if n_elements(yoff) eq 0 then yoff = 0

   if keyword_set(fixtrace) eq 0 and (n_elements(tracefile) gt 0 or traceycent ne ycent)  then begin
      ; get the optimum yoff by repeatedly extracting and measuring total flux
      yoffsearch = ycent-tracey[nx/2] + [-3,3]

      yoffarr = min(yoffsearch) + findgen(max(yoffsearch)-min(yoffsearch))
      totfluxy = fltarr(n_elements(yoffarr))
      for ei = 0, n_elements(yoffarr)-1 do begin
         ftry = extract(sci, tracey+yoffarr[ei], aprad/2. > 1.5, apsharp) ; note no background
         totfluxy[ei] = total(ftry)
         ;print, yoffarr[ei], tracey[nx/2]+yoffarr[ei], totfluxy[ei]
      endfor
      yoffmax = yoffarr[(where(totfluxy eq max(totfluxy))) [0]] ; maybe fit a gaussian in the future?
      yofffit = getlinecenter(yoffarr, totfluxy, yoffmax, fitrad=3)
      ;print, 'Optimum yoff = ', yoff
      ;print, yoffmax, yofffit
      print, '  Optimum trace center at x=',clip(nx/2), ': y=', clip(yofffit+tracey[nx/2])
      yoff = yofffit + yoff
   endif
   if keyword_set(fixtrace) and n_elements(ycent) gt 0 then begin
      yoff = ycent-tracey[nx/2]
   endif

   ; --- Extract 1D spectrum ---



   traceregfile = repstr(outfile,'.spec','.trace.reg')

   fobj = extract(sci, tracey+yoff, aprad, apsharp, bakdist, bakwidth, fbak=fbak, regfile=traceregfile, saturate=saturate, nsaturate=nsaturate)

   if nsaturate gt 1 then begin
     print, '  WARNING: '+clip(nsaturate)+' saturated pixels detected on the trace!!!'
     print, '     Spectrum may not be usable.'
   endif




   ; ------------------ Calibrate Wavelength -------------------

    checkplotfile = repstr(outfile,'.spec','.wavsol.ps')

    pix = findgen(nx)

    maxsep = nx*dispersion/2.

    if n_elements(arcfile) gt 0 then begin

      arc = mrdfits(arcfile, 0, arch, /silent) ; use the arcfile

      print, '  Extracting arc spectrum.'

      farc = extract(arc, tracey+yoff, aprad, apsharp)

      if n_elements(brightarclist) eq 1 then $
        brightreflinewav = grabcolumn(brightarclist, 0) $
      else $
        brightreflinewav = brightarclist    

      if n_elements(arclist) eq 1 then $
        morereflinewav = grabcolumn(arclist, 0) $
      else $
        morereflinewav = arclist    

      arcsol = solvearcspec(farc, wav, refwav=morereflinewav, brightrefwav=brightreflinewav, $
       cwav=cwav, unccwav=unccwav, dispersion=dispersion, uncdisp=uncdisp, $
       maxsep=maxsep, maxorder=maxorder, uncrefwav=uncrefwav, initsol=initsol, nbrightline=20, $
       nmatch=nmatch, ngoodmatch=ngoodmatch, checkplotfile=checkplotfile, maxnudge=maxnudge, /leaveplotopen)
      checkplotfile = '' ; don't reopen next time

      sol = arcsol
      initsol = arcsol
    endif
    if n_elements(skyfile) gt 0 then begin
       if strlen(clip(skyfile)) gt 1 then begin
          sky = mrdfits(skyfile, 0, skyh, /silent)
       endif else begin
          sky = mrdfits(filename, fix(skyfile), skyh, /silent) ; use extension
       endelse

       print, '  Extracting sky spectrum.'

       fsky = extract(sky, tracey+yoff, aprad, apsharp)  

       if skylist ne '' then begin
          readcol, skylist, idn, allreflinewav, allreflinestr, format='i,f,f'
          thresh = 100.
          minwav = cwav-unccwav[0] - (dispersion+max(uncdisp))*nx/2. ; estimation only
          maxwav = cwav+unccwav[0] + (dispersion+max(uncdisp))*nx/2. ; estimation only
          w = where(allreflinestr gt thresh and allreflinewav gt minwav and allreflinewav lt maxwav)
          skylinewav = allreflinewav[w] ; ultimately should edit the file?
          skylinestr = allreflinestr[w]

          nline = 100
          if maxwav lt 7000. then nline = 20
          if maxwav lt 6000. then nline = 5

          brightskywav = findbrightest(skylinewav, skylinestr, (nline/5)>3)   
          moreskywav = findbrightest(skylinewav, skylinestr, nline)

          print, '  Solving for wavelength using sky lines...'
          ;skysol = solvearc(linex, brightlinex, moreskywav, brightskywav, $
          ;                  cwavrange=cwavrange, disprange=disrange, maxsep=maxsep, initsol=initsol, nx=nx, uncrefwav=1.0)

          nofitsky = 1
          skysol = solvearcspec(fsky, wav, refwav=moreskywav, brightrefwav=brightskywav, $
           cwav=cwav, unccwav=unccwav, dispersion=dispersion, uncdisp=uncdisp, $
           maxsep=maxsep, uncrefwav=uncrefwav, initsol=initsol, nofit=nofitsky, nbrightline=20, $
           nmatch=nmatch, ngoodmatch=ngoodmatch, checkplotfile=checkplotfile, maxnudge=maxnudge)
          sol = skysol

       endif

    endif else begin  ; alternatively, could use fbak
       fsky = 0
    endelse 

    psclose
    !p.multi = 0

    if n_elements(sol) eq 0 and n_elements(initsol) gt 1 then begin
       print, '  Performing no refinements to the input wavelength calibration.'
       sol = initsol
       wav = poly(pix-nx/2., initsol)
    endif 

    success = 1
    if n_elements(wav) eq 0 then begin
       if n_elements(arcfile) gt 0 or n_elements(skyfile) gt 0 then begin
         print, '  Unable to calculate wavelength solution.'
         success = -4  ; failed wavelength solution
         return
       endif else begin
         print, '  No wavelength calibration information provided; expressing in pixels.'
       endelse
    endif
    if n_elements(nmatch) gt 0 then begin
;      if ngoodmatch lt fix(0.3*nmatch)-2 then begin
    if ngoodmatch lt fix(0.3*nmatch)-2 then begin
         success = 2   ; possible bad wavelength solution
      endif
    endif




   ; Save to disk
   objspec = fobj-fbak
   skyspec = fsky+fbak
   readnoise = 4.6*sqrt(2.0*aprad) ; formerly, total of extractionkernel
   gain = 1.2
   uncspec = sqrt(readnoise^2 + ((objspec + skyspec) > 0)*gain)/gain

   x = indgen(nx)

   ;inststr = strlowcase(inst)
   ;gratstr = (strsplit(grat,'/',/extract)) [0]
   ;outfile =  strlowcase(object)+'_'+gratstr+'.dat'  ;+'_'+clip(nsci+1) ;+clip(fix(cwave))

   openw, 1, outfile
   for k = 0, n_elements(keywords)-1 do begin
      kval = sxpar(h, keywords[k], count=ct)
      if ct ge 1 then printf, 1, '# ',clip(keywords[k],8), ' = ', kval
   endfor
   if n_elements(arcfile) gt 0 then begin
      printf, 1, '# ARCFILE  = ', arcfile
   endif
   printf, 1, '# CWAV     = ', cwav
   printf, 1, '# DISPAV   = ', dispersion
   printf, 1, '# YEXTRACT = ', tracey[nx/2]+yoff
   printf, 1, '# APRAD    = ', aprad
   printf, 1, '# APSHARP  = ', apsharp
   printf, 1, '# NSATURAT = ', nsaturate
   printf, 1, '# NWAV     = ', n_elements(wav)
   if n_elements(bakdist) gt 0 then printf, 1, '# BAKDIST  = ', bakdist
   if n_elements(bakdist) gt 0 then printf, 1, '# BAKWIDTH = ', bakwidth

   if n_elements(wav) gt 0 then begin
      for i = 0, n_elements(wav)-1 do begin
         printf, 1, wav[i], objspec[i], skyspec[i], uncspec[i], pix[i]
      endfor
   endif else begin
      for i = 0, n_elements(pix)-1 do begin
         printf, 1, pix[i], objspec[i], skyspec[i], uncspec[i]
      endfor
   endelse
   close, 1
   print, '  Spectrum saved to ', outfile


end

