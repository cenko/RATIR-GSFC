;+
; NAME:
;   pzapspec.pro
;
; PURPOSE:
;   Remove cosmic rays from a 2-D spectrum using model subtraction / percolation.
;
; CALLING SEQUENCE:
;   pzapspec, name, [outname, outmaskname], boxsize=boxsize, nsigma=nsigma, 
;          subsigma=subsigma, nzap=nzap, mask=mask, writemodel=writemodel, 
;          debug=debug, verbose=verbose, skysubtract=skysubtract, zero=zero, usamp=usamp
;
; INPUTS:
;   name        - 2-D image array, name of input FITS file, or name of a 
;                 text file containing a list of FITS filenames.
;                 (input as array has not been bug-tested)
;
; OPTIONAL INPUTS:
;   outname     - Output image array, or name of output FITS file.
;   outmaskname - Output mask array, or name of output FITS file.
;   boxsize     - Boxsize for computing local median; default to 9.
;                 Should be similar to the seeing in pixels.
;   nsigma      - Rejection threshhold in sigma.  10+ advised.
;   subsigma    - Sigma for neighbors to also be considered cosmic rays 
;   mask        - Write out a mask identifying zapped pixels.
;   writemodel  - Write out the X+Y (source+sky) model-subtracted image.
;   zero        - Replace bad pixels with 0 instead of model value (for easy identification)
;   skysubtract - Subtract sky lines.
;   debug       - Write out other intermediate image products.
;   verbose     - Print out information during processing.
;
; OPTIONAL OUTPUTS:
;   nzap        - Number of pixels zapped.
;
; COMMENTS:
;   A very robust cosmic ray rejection and replacement script optimized
;   for the new LRIS-red, a deep-depletion CCD susceptible to muon-trailed
;   cosmic ray impacts which leave long, snaking patterns across the chip.
;   Replacing nnzapspec, this program creates a detailed model of the source
;   and sky based on multiple median filters (sky [X], source [Y], and
;   line [XY] filters) and identifies cosmic rays as large excesses above
;   this background.  A percolation algorithm is used to flag the entire
;   ray trail.
;
;   Has been tested only two images to date (both faint objects) but is
;   designed to work robustly on brighter sources as well.  Significant 
;   caution should be exercised and the /mask option should always be
;   set to enable the user to verify that unaffected pixels are not being
;   flagged.
;   
;   Typically, the only cosmic rays that are "missed" are long, straight
;   and horizontal/vertical (and therefore difficult to distinguish from a
;   source or sky line).  Typically these are actually identified by the
;   algorithm as cosmic rays and will be masked, but the splotching 
;   algorithm does not replace them well.  Cosmic rays within a few pixels
;   of the edge are also typically missed.  Situations where a cosmic ray
;   passes through a strong source line should be treated with extreme 
;   caution as the line is likely to be flagged.
;
;   It is very important to set the boxsize parameter accurately.  The
;   default value of 9 will likely be too large for good seeing and a
;   narrow slit; in this case reducing it to 7 is advised.
;   (This is still killing some lines in old LRIS spectra...)
;   There is typically no need to alter the sigma parameters unless
;   cosmic rays are being systematically missed or real features are
;   being flagged even after optimizing.
;
;   Last significant update May 19, 2010.
;   Last code modification 2012-01-24
;
;   Issue discovered 2012-04-30: effective sigma is strongly affected by
;   bad cropping, and many CR's are missed.
;
;   Currently cannot remove all CR's on LRIS-red without removing strong
;   narrow lines in any setting.  Should do better modeling (borrow from
;   pzap) or add a gradient check.
;   Replacement algorithm seems to not work well sometimes (see 120522_0199)
;
;-
;------------------------------------------------------------------------------

@modelimage

pro pzapspec, name, outname, outmaskname, boxsize=boxsize, nsigma=nsigma, subsigma=subsigma, $
 nzap=nzap, mask=mask, writemodel=writemodel, debug=debug, verbose=verbose, skysubtract=skysubtract, zero=zero, $
 method=method, usamp=usamp

   ;if n_elements(boxsize) eq 0 then boxsize=7 ;9
   if n_elements(nsigma) eq 0 then nsigma=15 
   if n_elements(subsigma) eq 0 then subsigma=2.8

   if n_elements(method) eq 0 then method = 0

   if (size(name, /tname) EQ 'STRING') then begin
      if n_elements(name) eq 0 then begin
         print, '   A filename must be specified.'
         return
      endif
      if n_elements(name) gt 1 then begin
         for f = 0, n_elements(name)-1 do begin
            print, '   ', name[f]
            pzapspec, name[f], boxsize=boxsize, nsigma=nsigma, subsigma=subsigma, mask=mask, $
              writemodel=writemodel, debug=debug, verbose=verbose, skysubtract=skysubtract, zero=zero, usamp=usamp
         endfor
         return
      endif
      if n_elements(name) eq 1 then begin
         namearr = strsplit(name,'.', /extract)
         nameext = namearr[1]
         if nameext eq 'cat' or nameext eq 'lis' or nameext eq 'list' or nameext eq 'txt' then begin
            ;files = grabcolumn(name,0,/str)
            readcol, name, files, format='a', /silent
            pzapspec, files, boxsize=boxsize, nsigma=nsigma, subsigma=subsigma, mask=mask, $ 
              writemodel=writemodel, debug=debug, verbose=verbose, skysubtract=skysubtract, zero=zero, usamp=usamp
            return
         endif


         outimg=mrdfits(name, 0, header, /silent)

         if n_elements(outimg) eq 1 then begin
            print, '   Error reading ', name
            return
         endif

         if n_elements(outname) eq 0 then begin
            slashpos = strpos(name,'/')
            dir = strmid(name,0,slashpos+1)
            outname = dir + 'z' + strmid(name,slashpos+1)
         endif
      endif
   endif else begin
      outimg=name
   endelse


   if keyword_set(mask) then begin
      namearr = strsplit(outname, '.', /extract)
      outmaskname = namearr[0] + '.mask.' + namearr[1]
   endif





   dims = size(outimg, /dimens)
   nx = dims[0]
   ny = dims[1]
   outmask = bytarr(nx, ny)  ;create an array to hold CR position mask
   ymedimage = fltarr(nx, ny)
   xmedimage = fltarr(nx, ny)
   zapimage = fltarr(nx, ny)

   nzap = 0
   iter = 0
   nbad = 1


   holdtime = systime(1)

   ;while (iter LT maxiter and nbad GT 0) do begin

      ; First, model the photon flux in the image using a sum of 1D, slowly-varying arrays
      ; by calculating medians along each axis to remove sky lines and most source flux.

      if keyword_set(verbose) then print, '   Calculating sky/source model...'

      ; First do a crude median subtraction of the sky lines
      for x = 0, nx-1 do begin
         ymedimage[x,*] = median(outimg[x,*], /even)
      endfor
      ysubimage = outimg - ymedimage

      ; Now subtract traces so we can do a better sky-line subtraction
      ; (otherwise traces may be killed or wings oversubtracted)
      sourceblocksize = 100
      nxblock = ceil(1.0*nx/sourceblocksize)  ; this block needs to be significantly bigger than any sky feature
      realsourcefiltsize = sourceblocksize
      x0 = indgen(nxblock)*realsourcefiltsize
      x1 = [x0[1:*]-1, nx-1]
      xs0 = [x0[0], x0[0:nxblock-2]]               ; use 3x blocksize
      xs1 = [x1[1:nxblock-1], x1[nxblock-1]]
      for b = 0, nxblock-1 do begin ;block iteration
         for y = 0, dims[1]-1 do begin
            xmedimage[x0[b]:x1[b],y] = median(ysubimage[xs0[b]:xs1[b],y], /even)
         endfor
      endfor
      kernel = [findgen(realsourcefiltsize/2), reverse(findgen(realsourcefiltsize/2))]  
      kernel = kernel / total(kernel) ; this is a strange kernel, should probably use boxcar
      xmedimage = convol(xmedimage, kernel, /edge_truncate)
      xsubimage = outimg - xmedimage

      ;print, systime(1)-holdtime
      ;holdtime = systime(1)


      ; Find and subtract sky again more carefully now that sources have been removed

      if method eq 0 then begin
         skyblocksize = 40
         nyblock = ceil(1.0*ny/skyblocksize)
         realskyblocksize = skyblocksize
         y0 = indgen(nyblock)*realskyblocksize
         y1 = [y0[1:*]-1, ny-1]
         ys0 = [y0[0], y0[0:nyblock-2]]          ; use 3x blocksize
         ys1 = [y1[1:nyblock-1], y1[nyblock-1]]
         for b = 0, nyblock-1 do begin ;block iteration
            for x = 0, nx-1 do begin
               ymedimage[x,y0[b]:y1[b]] = sigclipmedian(xsubimage[x,ys0[b]:ys1[b]])
            endfor
         endfor
         kernel = transpose([findgen(realskyblocksize/2), reverse(findgen(realskyblocksize/2))]) ; should probably use boxcar
         kernel = kernel / total(kernel)
         ymedimage = convol(ymedimage, kernel, /edge_truncate)
         ysubimage = outimg - ymedimage
      endif else begin
          ysubimage = fltarr(nx, ny)
          skyblocksize = 15
          nyblock = ceil(1.0*ny/skyblocksize)
          realskyblocksize = ceil(1.0*ny/nyblock)
          y0 = indgen(nyblock)*realskyblocksize
          y1 = [y0[1:*]-1, ny-1]
          ;ys0 = [y0[0], y0[0:nyblock-2]]          ; use 3x blocksize
          ;ys1 = [y1[1:nyblock-1], y1[nyblock-1]]
          ys0 = (y0 - 10) > 0
          ys1 = (y1 + 10) < ny-1
          for b = 0, nyblock-1 do begin ;block iteration
             print, '   ', b, ':    ', clip(ys0[b],4), ' ', clip(y0[b],4), '-  ', clip(y1[b],4), ' ', clip(ys1[b], 4)
             sub = subtractsky(xsubimage[*,ys0[b]:ys1[b]], sky=sky, order=1, ysample=4, debug=(b eq 15))
             ;ysubimage[*,y0[b]:y1[b]] = sub[*,(y0[b]-ys0[b]):(y0[b]-ys0[b])+(y1[b]-y0[b])]
             ymedimage[*,y0[b]:y1[b]] = sky[*,(y0[b]-ys0[b]):(y0[b]-ys0[b])+(y1[b]-y0[b])]
          endfor      
          kernel = transpose([findgen(realskyblocksize/2), reverse(findgen(realskyblocksize/2))]) ; should probably use boxcar
          kernel = kernel / total(kernel)
          ymedimage = convol(ymedimage, kernel, /edge_truncate)
          ysubimage = outimg - ymedimage
      endelse

      ;This simple median subtraction is much slower.
      ;for y = 0, ny-1 do begin ;block iteration
      ;   for x = 0, nx-1 do begin
      ;      ymedimage[x,y] = median(xsubimage[x,(y-skyfiltsize/2)>0:y+skyfiltsize/2<(ny-1)], /even)
      ;   endfor
      ;endfor

      ; One could put a more careful source subtraction here for a better sigma map and
      ; sky subtraction, but for the moment it doesn't seem to be needed.
 

      ;print, systime(1)-holdtime
      ;holdtime = systime(1)


      skysubimage = outimg - ymedimage - xmedimage  ; actually subtracts sky AND sources.

      if keyword_set(debug) then begin
        skyiter = ''
        writefits, 'ymedimage.fits', ymedimage, header
        writefits, 'xmedimage.fits', xmedimage, header
        writefits, 'ysubimage.fits', ysubimage, header
        writefits, 'xsubimage.fits', xsubimage, header
      endif


      if keyword_set(writemodel) then begin
         namearr = strsplit(outname, '.', /extract)
         outskyname = namearr[0] + '.subtract.' + namearr[1]
         writefits, outskyname, skysubimage, header
      endif


      if keyword_set(verbose) then print, '   Subtracting features...'

      ; Now try to remove any 2D photon features (lines, etc.) with a 2D median filter.

      if n_elements(boxsize) gt 0 then begin
         filterimage = median(skysubimage, boxsize) ; unfortunately there is no way to mask out flagged pixels.
      endif else begin
         if keyword_set(usamp) eq 0 then filterimage = model2image(skysubimage) $
                                    else filterimage = usampmodelimage(skysubimage)
      endelse
      residualimage = skysubimage - filterimage

      ; Make the sigma image (noise map) using the geometric sum of column sigma (read noise, sky noise)
      ; and source+filter residuals (source noise)
      sigmaimage = fltarr(nx,ny)
      for x = 0, dims[0]-1 do begin
        s = residualimage[x,*]
        ss = s[(sort(s)) [dims[1]*0.03:dims[1]*0.95]]  ;clip out bad values which otherwise bloat sigma a lot
        sigmaimage[x,*] = stdev(ss)                    ;(CRs contribute a significant fraction of the flux in spectroscopy)
      endfor
      sigmaimage = sqrt(sigmaimage^2 + (abs(filterimage) + abs(xmedimage)))  ; this scaling assumes gain ~ 1
      
      residualnsigmaimage = residualimage / sigmaimage

      if keyword_set(debug) then begin
        writefits, 'xmedimage.fits', xmedimage, header
        writefits, 'ymedimage.fits', ymedimage, header
        writefits, 'skysubimage.fits', skysubimage, header
        writefits, 'filterimage.fits', filterimage
        writefits, 'residualimage.fits', residualimage, header
        writefits, 'sigmaimage.fits', sigmaimage, header
        writefits, 'residualnsigmaimage.fits', residualnsigmaimage, header
      endif


      d0 = dims[0]



      if keyword_set(verbose) then print, '   Identifying cosmic rays and percolating...'

      ; Find severely affected cosmic ray pixels (those grossly in excess of the model)

      crcores = where(residualnsigmaimage gt nsigma, newzaps)
      zapimage[crcores] = 1

      if keyword_set(verbose) then print, '   Flagged ', strtrim(newzaps,2), ' initial affected pixels before percolation.'


      ;  The old 1st order NN method is below.  
      ;  Right now the model seems to work well enough that this is no longer necessary.
      ;coarsesigma = 15     ;zap points >15 sigma
      ;coarsesubsigma = 6   ;comparing to >6 sigma immediate neighbors
      ;coarsenthresh = 3    ;don't zap if it has three or more such neighbors (could be a badly subtracted source)
      ;for x = 2, dims[0]-1-2 do begin
      ;   for y = 2, dims[1]-1-2 do begin
      ;      if residualnsigmaimage[x,y] lt coarsesigma then continue
      ;       coarseblockpos = [(y-1)*d0+(x-1),(y-1)*d0+x,(y-1)*d0+(x+1), $
      ;                         (y  )*d0+(x-1),           (y  )*d0+(x+1), $
      ;                         (y+1)*d0+(x-1),(y+1)*d0+x,(y+1)*d0+(x+1)]
      ;       highneighbor = where(residualnsigmaimage[coarseblockpos] gt coarsesubsigma, ct)
      ;       if ct lt coarsenthresh then begin
      ;         zapimage[x,y] = 1
      ;         nzap = nzap + 1
      ;         ;if ct gt 0 then zapimage[coarseblockpos[highneighbor]] = 1
      ;       endif
      ;   endfor
      ;endfor

     if keyword_set(debug) then writefits, 'zapcoreimage.fits', zapimage, header

      ; Now percolate outward to get all pixels covered by each cosmic ray.
      nperczap = 0
      iter = 1
      while iter lt 32 do begin
         nextperc = where(zapimage eq iter, ct)
         iter = iter + 1
         newzaps = 0        
         if ct le 3 then break
         nrays = n_elements(nextperc)
         for c = 0L, nrays-1 do begin
             ci = nextperc[c]
             if ci lt d0-1 or ci gt nx*ny-d0-2 then continue
             coarseblockpos = [         ci-d0,       $
                               ci-1,           ci+1, $
                                        ci+d0         ]
             ;coarseblockpos = [ci-d0-1, ci-d0, ci-d0+1, $
             ;                  ci-1,           ci+1,    $
             ;                  ci+d0-1, ci+d0, ci+d0+1]
             newzap = where(abs(residualnsigmaimage[coarseblockpos]) gt subsigma and zapimage[coarseblockpos] eq 0, addzap)
             if addzap gt 0 then begin
                zapimage[coarseblockpos[newzap]] = iter
                newzaps = newzaps + addzap
             endif
             nperczap = nperczap + addzap
         endfor
         ;if keyword_set(verbose) then print, '   Iteration ', strtrim(iter,2), ':', newzaps, ' pixel zaps.'
      endwhile
      if keyword_set(verbose) then print, '   Flagged ', strtrim(nperczap,2), ' pixels during ', strtrim(iter,2), ' percolation steps.'

      ; Finally, zap anything hiding in a cosmic ray "corner" (three neighbors are cosmic ray pixels)

      countneighbor = intarr(nx, ny)
      nextperc = where(zapimage gt 0, ct)
      if ct gt 3 then begin
         nrays = n_elements(nextperc)
         for c = 0L, nrays-1 do begin
             ci = nextperc[c]
             coarseblockpos = [ci-d0-1, ci-d0, ci-d0+1, $
                               ci-1,           ci+1,    $
                               ci+d0-1, ci+d0, ci+d0+1]
             countneighbor[coarseblockpos] = countneighbor[coarseblockpos] + 1
         endfor
         newzap = where(countneighbor gt 3 and zapimage eq 0, newzaps)
         if newzaps gt 0 then begin
            if keyword_set(verbose) then print, '   Flagged an additional ', strtrim(newzaps,2), ' neighboring pixels.'
            zapimage[newzap] = iter + 1
         endif
      endif



     if keyword_set(debug) then writefits, 'zapimage.fits', zapimage, header





      ;Actually do the zapping by replacing with the sky plus local median value.
      ibad = where(zapimage NE 0, nbad)
      if (nbad GT 0) then begin
         skysubimage[ibad] = !Values.F_NAN           ; NaNs ignored by median() and derivative image modelsers

         ; filterimage is really the replacement image (without sky and sources)
         if n_elements(boxsize) gt 0 then begin
            filterimage = median(skysubimage, boxsize) ; unfortunately there is no way to mask out flagged pixels.
         endif else begin
            filterimage = median(skysubimage, 7)    ; stick with median replacement
            ;if keyword_set(usamp) eq 0 then filterimage = model2image(skysubimage) $
            ;                           else filterimage = usampmodelimage(skysubimage)
         endelse

         ; see if we actually got all the NaNs or if any remain.
         witer = 0
         ctnan = -1
         while ctnan ne 0 and witer le 2 do begin

           if ctnan gt 0 then begin   ; note not active 0th iteration, since we already did the median above
             if n_elements(boxsize) gt 0 then begin
                filterimage[filterbad] = (median(filterimage, boxsize)) [filterbad] ; unfortunately there is no way to mask out flagged pixels.
             endif else begin
                filterimage[filterbad] = (median(filterimage, 7)) [filterbad]    ; stick with median replacement
               ;  if keyword_set(usamp) eq 0 then filterimage[filterbad] = (model2image(filterimage)) [filterbad]  $
               ;                             else filterimage[filterbad] = (usampmodelimage(filterimage)) [filterbad]
             endelse


           endif
           ;mwrfits, filterimage, 'ofiltim'+clip(witer)+'.fits', h, /create
           filterbad = where(finite(filterimage) eq 0 and finite(xmedimage) eq 1 and finite(ymedimage) eq 1, ctnan)  
                      ; any NaNs left (that we can deal with?)
           
           ;There is some strange behavior here I don't fully understand where it takes 2 iterations to flag pixels near the edges(?)
           ; although nothing has really changed from the 1st iteration.  But, doesn't seem to affect results, so don't worry for now.
           ;print, witer, ctnan
           witer += 1
           if witer eq 2 and ctnan gt 0 then filterimage[filterbad] = 0   ; give up
         endwhile
         ;filterbad = where(finite(filterimage) eq 0, ctfbad)

         outimg[ibad] = filterimage[ibad] + ymedimage[ibad] + xmedimage[ibad]
         outmask[ibad] = 1
      endif


      print, '   Zapped ', strtrim(string(nbad),2), ' total affected pixels (', strtrim(string((nbad)*100./(1.0*nx*ny),format='(F6.3)'),2), '% of total).'

      nzap = nzap + nbad

   ;endwhile

   if (nbad GT 0 and keyword_set(zero)) then begin  ; really should skip the whole above iterators if this is active
      outimg[ibad] = 0
   endif

    ;d = where(finite(outimg) eq 0, ctobad)

   if (size(outname, /tname) EQ 'STRING') then begin
      if n_elements(outmaskname) gt 0 then sxaddpar, header, 'BPM', outmaskname, 'Bad pixel mask'
      writeimg = outimg
      if keyword_set(skysubtract) then writeimg = outimg-ymedimage
      writefits, outname, writeimg, header
   endif else begin
      outname = outimg
   endelse

   if (size(outmaskname, /tname) EQ 'STRING') then begin
      writefits, outmaskname, outmask, header
   endif else begin
      outmaskname = outmask
   endelse

end

