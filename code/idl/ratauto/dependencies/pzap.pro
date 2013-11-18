;+
; NAME:
;   pzap.pro
;
; PURPOSE:
;   Remove cosmic rays from a 2-D image using model subtraction / percolation.
;
; CALLING SEQUENdE:
;  pzap, name, outname, outmaskname, nsigma=nsigma, subsigma=subsigma, 
;  maxiter=maxiter, zealfactor=zealfactor, nrings=nrings, nzap=nzap, 
;  mask=mask, weight=weight, zero=zero, debug=debug, verbose=verbose
;
; INPUTS:
;   name        - 2-D image array, or name of input FITS file.
;   outname     - Output image array, or name of output FITS file.
;
; OPTIONAL INPUTS:
;   outmaskname - Output mask array, or name of output FITS file.
;   nsigma      - Rejection threshhold in sigma
;   subsigma    - Sigma for neighbors to also be considered cosmic rays 
;   maxiter     - Number of zapping iterations
;   zealfactor- Raises the threshold for zapping cosmic rays inside of sources (default 1.2)
;   nrings      - Also zap pixels within X pixels of identified bad pixels
;   mask        - Write out a mask identifying zapped pixels.
;   weight      - Write out weight map (inverse of mask).
;   zero        - Replace bad pixels with 0 instead of model value (for easy identification)
;   debug       - Write some intermediate image products to disk
;   verbose     - Print details of the percolation process;
;
; OPTIONAL OUTPUTS:
;   nzap        - Number of pixels zapped.
;
; COMMENTS:
;   A very robust cosmic ray rejection and replacement script optimized
;   for the new LRIS-red, a deep-depletion CCD susceptible to muon-trailed
;   cosmic ray impacts which leave long, snaking patterns across the chip.
;   Replacing nnzaps, this program creates a model of the image using
;   nearby pixels and and identifies cosmic rays as large excesses above
;   this background.  A percolation algorithm is used to flag the entire
;   ray trail.
;  
;   The routine is effective but not perfect, typically removing about 99% of
;   cosmic ray flux.  Missed pixels tend to be the weakest cosmic rays with
;   wide, diagonal trails which the algorithm has difficulty distinguishing
;   from real sources in images with significant sky noise (usually, some but
;   not all of the ray trail is typically flagged in such cases.)
;
;   In testing, real sources are almost never affected, with the only 
;   exceptions observed ocurring along bad columns.
;
;
; PROCEDURES CALLED:
;   djs_iterstat
;
; REVISION HISTORY:
;          2009 - Replaced qzap's flux comparison with nearest neighbor comparison.
;      Jan 2010 - Numerous changes to selection and replacement algorithm, taking
;                 into account image noise statistics and making it harder to flag
;                 regions of detected source to dramaticall reduce number of false
;                 positives in imaging.
;-
;------------------------------------------------------------------------------

; This has performed more poorly than hoped.  Some definite critical improvements:
;   - handle edges properly (ignore them)
;   - do some tests with adaptive PSF/parameter modeling to optimize sigma parameters.

@modelimage

pro pzap, name, outname, outmaskname, nsigma=nsigma, subsigma=subsigma, preiter=preiter, postiter=postiter, $
 zealfactor=zealfactor, nrings=nrings, nzap=nzap, mask=mask, weight=weight, zero=zero, statsec=statsec, debug=debug, usamp=usamp, $
 verbose=verbose, quiet=quiet, silent=silent, $
saturation=saturation

   ;if (NOT keyword_set(skyfiltsize)) then skyfiltsize=15
   ;if (NOT keyword_set(boxsize)) then boxsize=5
   if (NOT keyword_set(nsigma)) then nsigma=5.0
   if n_elements(clnsigma) eq 0 then clnsigma = 3.7               ;3.5 were the LRIS values
   if n_elements(clsubnsigma) eq 0 then clsubnsigma = 2.9         ;2.8
   if (NOT keyword_set(subsigma)) then subsigma = 2.0
   if (NOT keyword_set(zealfactor)) then zealfactor = 1.0    
   if (NOT keyword_set(preiter)) then preiter=3
   if (NOT keyword_set(postiter)) then postiter=1
   if (NOT keyword_set(nrings)) then nrings=1

   nsigma = nsigma / zealfactor
   clnsigma = clnsigma / zealfactor
   clsubnsigma = clsubnsigma / zealfactor

stojkhlj
   if (size(name, /tname) EQ 'STRING') then begin
      if n_elements(name) eq 0 then begin
         print, 'A filename must be specified.'
         return
      endif
      if n_elements(name) gt 1 then begin
         for f = 0, n_elements(name)-1 do begin
            print, name[f]
            pzap, name[f],  $
               nsigma=nsigma, subsigma=subsigma, preiter=preiter, postiter=postiter, zealfactor=zealfactor, $
               nrings=nrings, nzap=nzap, mask=mask, weight=weight, zero=zero, debug=debug, usamp=usamp, $
               verbose=verbose, quiet=quiet, silent=silent, saturation=saturation
         endfor
         return
      endif
      if n_elements(name) eq 1 then begin
         namearr = strsplit(name,'.', /extract)
         nameext = namearr[1]
         if nameext eq 'cat' or nameext eq 'lis' or nameext eq 'list' or nameext eq 'txt' then begin
            files = grabcolumn(name,0,/str)
            pzap, files,  $
              nsigma=nsigma, subsigma=subsigma, preiter=preiter, postiter=postiter, zealfactor=zealfactor,  $
              nrings=nrings, nzap=nzap, mask=mask, weight=weight, zero=zero, debug=debug, usamp=usamp, $
              verbose=verbose, quiet=quiet, silent=silent, saturation=saturation
         return
         endif


         img=readfits(name, header, /silent)

         if n_elements(img) eq 1 then begin
            print, 'Error reading ', name
            return
         endif
         dum = where(finite(img), ctnotnan)
         if ctnotnan eq 0 then begin
            print, 'No valid pixels in image ', name
            return
         endif

         if n_elements(outname) eq 0 then begin
            slashpos = strpos(name,'/')
            dir = strmid(name,0,slashpos+1)
            outname = dir + 'z' + strmid(name,slashpos+1)
stop
         endif
      endif
   endif else begin
      img=name
   endelse
stop
   if keyword_set(mask) then begin
      if strtrim(string(mask),2) eq '1' then begin
         namearr = strsplit(outname, '.', /extract)
         outmaskname = namearr[0] + '.mask.' + namearr[1]
      endif else begin
         outmaskname = mask
      endelse
   endif
   if keyword_set(weight) then begin
      if strtrim(string(weight),2) eq '1' then begin
         namearr = strsplit(outname, '.', /extract)
         outweightname = namearr[0] + '.weight.' + namearr[1]
      endif else begin
         outweightname = weight
      endelse
   endif

   tstart = systime(/seconds)

   dims = size(img, /dimens)
   nx = dims[0]
   ny = dims[1]
   outmask = bytarr(dims[0], dims[1])  ;create an array to hold CR position mask

   initbad = where(finite(img) eq 0, ninitbad)

   ; Calculate the overall st.dev. of the image sky
   if keyword_set(verbose) then print, '  Computing image sigma...'

   if n_elements(statsec) eq 0 then begin
     if n_elements(img) gt 5e5 then begin
       statsec = [0, nx-1, ny/2-100 > 0, ny/2+100 < ny-1]
       ; speed things up by only using a horizontal stripe of the image to calculate statistics.
     endif
   endif
   if n_elements(statsec) eq 4 then begin
     imgstat = img[statsec[0]:statsec[1],statsec[2]:statsec[3]]
   endif else begin
     if n_elements(statsec) gt 1 then begin
        print, '  statsec format:  [x0, x1, y0, y1]'
        return
     endif
     imgstat = img
   endelse

   ; Apply a median filter to remove any gradients, halos, etc. that would screw up our statistics
   skyfiltsize = 11
   skyimgstat = median(imgstat, skyfiltsize)
   ;Subtract the median to remove the sky
   ;skysubimage = img - skyimage

   skystatinitbad = where(finite(skyimgstat) eq 0, ct)
   if ct gt 0 then skyimgstat[skystatinitbad] = median(skyimgstat)
   statinitbad = where(finite(imgstat) eq 0, ct)
   if ct gt 0 then imgstat[statinitbad] = skyimgstat[statinitbad] ; djs barfs if nans are around.
   djs_iterstat, imgstat-skyimgstat, median=skylevel, sigma=sigval, sigrej=4
   skylevel = skylevel + median(skyimgstat)  ; skylevel is a scalar indicating the average sky value.


   skymodel = modelsky(img) ; skymodel is an image indicating the broad variations in the sky.
                              ; we don't use it for statistics, but will be useful later for
                              ; properly dealing with the safety factor.
    if keyword_set(debug) then mwrfits, skymodel, 'skymodel.fits', header, /silent, /create

   ; Get a model


   nzap = 0
   iter = 0
   prevbad = 0
   ;nbad = 1 wtf?
   maxiter = preiter + postiter

   if n_elements(saturation) gt 0 then begin
      satlevel = saturation
   endif else begin
      satlevel = sxpar(header,'SATURATE')
   endelse
   if satlevel eq 0 then begin
     hipix = where(img gt skylevel + 10*sigval,cthipix)
     if cthipix gt 0 then begin
        imghipix = img[hipix]
        sr = sort(imghipix)
        satlevel = imghipix[sr[n_elements(sr)*0.98]] * 0.88 > 48000.
     endif else begin
        satlevel = 48000.
     endelse 
   endif
   if keyword_set(verbose) then stop
;      print, '  sky = ', strtrim(string(skylevel),2), ' (sigma = ', strtrim(string(sigval),2), ', saturation = ', strtrim(long(satlevel),2),')'
   satpix = where(img gt satlevel, ctsat)
   if ctsat gt n_elements(img)/2 then begin
       print, '  The entire image is saturated.  Cannot cosmic-ray clean.'
       return
   endif
   if ctsat gt 0 then begin
     satpix = [satpix, satpix+1, satpix-1, satpix+nx, satpix-nx, satpix+1+nx, satpix-1+nx, satpix+1-nx, satpix-1-nx] ; grow a little 
     satpix = satpix[where(satpix ge 0 and satpix lt n_elements(img))]
   endif
stop
   ;;model2 = model2image(img) 
   ;;mwrfits, model2, 'model2image.fits', header, /silent, /create
   ;;mwrfits, img-model2, 'model2imageresidual.fits', header, /silent, /create
   ;;mwrfits, (img-model2)/img, 'model2imageresidualnorm.fits', header, /silent, /create

   titer0 = systime(/seconds)

   while (iter LT maxiter) do begin
      iter = iter + 1
      nbad = 0

      titerstart = systime(/seconds)

      if (keyword_set(quiet) eq 0 and keyword_set(silent) eq 0) then begin
        if iter le preiter then print, '   Pre-iteration', iter, format='($,A,I3)'
        if iter gt preiter then print, '  Post-iteration', iter-preiter, format='($,A,I3)'
      endif
      
      ;;Apply a smaller median filter to the subtracted image to get a measure of the local average
      ;;This is a crappy way to model the flux - obsolete
      ;filterimage = median(skysubimage, boxsize)  
   
      imgnosat = img
      if ctsat gt 0 then imgnosat[satpix] = !values.f_nan

      if keyword_set(usamp) eq 0 then begin
        if iter eq 1 then filterimage = model2image(imgnosat)
        if iter eq 2 then filterimage = modelimage(imgnosat)
        if iter eq 3 then filterimage = modelimage(imgnosat)
        if iter eq 4 then filterimage = model2image(imgnosat)
        if iter ge 5 then filterimage = modelimage(imgnosat)
        ; the right number for modeluncertaintyfactor is highly uncertain
        if iter eq 1 then modeluncertaintyfactor = 0.05/zealfactor 
        if iter eq 2 then modeluncertaintyfactor = 0.12/zealfactor
        if iter eq 3 then modeluncertaintyfactor = 0.12/zealfactor
        if iter eq 4 then modeluncertaintyfactor = 0.05/zealfactor
        if iter ge 5 then modeluncertaintyfactor = 0.05/zealfactor
      endif else begin
        if iter eq 1 then filterimage = usampmodelimage(imgnosat)
        if iter eq 2 then filterimage = usampmodelimage(imgnosat)
        if iter eq 3 then filterimage = usampmodelimage(imgnosat)
        if iter eq 4 then filterimage = usampmodelimage(imgnosat)
        if iter ge 5 then filterimage = usampmodelimage(imgnosat)
        if iter eq 1 then modeluncertaintyfactor = 0.1/zealfactor 
        if iter eq 2 then modeluncertaintyfactor = 0.1/zealfactor
        if iter eq 3 then modeluncertaintyfactor = 0.1/zealfactor
        if iter eq 4 then modeluncertaintyfactor = 0.1/zealfactor
        if iter ge 5 then modeluncertaintyfactor = 0.05/zealfactor
      endelse


      ;Subtract the model to remove genuine sources
      residualimage = img - filterimage
      if ctsat gt 0 then residualimage[satpix] = 0
      
      ;if iter eq 1 then djs_iterstat, residualimage, sigma=sigval, sigrej=5
      ;print, 'Sigma = ', sigval   ; RMS of the background sky



      ; sigmaimage = sigval * sqrt(filterimage/median(filterimage) > 0.9)
      sigmaimage = sigval * sqrt(filterimage/median(filterimage) > 0.9) + modeluncertaintyfactor*((filterimage-skymodel)>0)
                   ; counting noise term                                   ; modeling error term

      residualnsigmaimage = residualimage / sigmaimage
      ;boxsize = sqrt(16)
      ;sourcensigmaimage = (filterimage-skymodel) *  boxsize / sigmaimage  ; formerly boxsize is because this is filtered, sqrt(N) = sqrt(box^2)
      ;extrathreshimage = safetyfactor * sqrt((sourcensigmaimage-1.5) > 0)  

      ;residualnsigmaimage = residualnsigmaimage - extrathreshimage

      if keyword_set(debug) then begin ;and iter eq 1 then begin
      mwrfits, filterimage, 'filterimage'+clip(iter)+'.fits', header, /silent, /create
      mwrfits, residualimage, 'residualimage'+clip(iter)+'.fits', header, /silent, /create
      mwrfits, sigmaimage, 'sigmaimage'+clip(iter)+'.fits', header, /silent, /create
      ;mwrfits, sourcensigmaimage, 'sourcensigmaimage.fits', header, /silent, /create
      mwrfits, residualnsigmaimage, 'residualnsigmaimage'+clip(iter)+'.fits', header, /silent, /create
      ;mwrfits, nsigma+extrathreshimage,  'threshimage.fits', header, /silent, /create
      endif

      ;Identify regions where the remaining flux is larger than expected given sky noise
      ;zapimage = residualnsigmaimage GT nsigma

      ;Do a nearest-neighbors comparison to see if there are other high-sigma points nearby.
      ;Examine separately 1st-nearest-neighbors and 2nd-nearest-neighbors.
      ;If there are many 2nd nearest neighbors, this is probably a real source and is not flagged.
      ;If there are few 1st nearest neighbors, this is probably a statistical fluctuation in the noise
      ;  and is not flagged (or a bad pixel).
      ;This is poorly optimized.
      zapimage = intarr(dims[0], dims[1])

      if keyword_set(verbose) then print, '  Identifying cosmic rays and percolating...'

      ; Find severely affected cosmic ray pixels
      crcores = where(residualnsigmaimage gt nsigma, initzaps)

      ; Find slightly less severely cosmic ray pixels if they have neighbors
      potentialsubcrcores = where(residualnsigmaimage gt clnsigma and residualnsigmaimage lt nsigma, ncheck4zap)
      if ncheck4zap gt 0 then begin
        subcrcores = lonarr(ncheck4zap)
        b = 0
        for i = 0L, ncheck4zap-1 do begin
           s = potentialsubcrcores[i]
           coarseblockpos = [         s-nx,       $
                             s-1,           s+1, $
                                      s+nx         ]
           if max(residualnsigmaimage[coarseblockpos]) gt clsubnsigma then begin
             subcrcores[b] = s
             b = b + 1
           endif
        endfor
        if b gt 0 then subcrcores = subcrcores[0:b-1]
        ;print, ''
        ;print, 'Additional', b
        if b gt 0 then zapimage[subcrcores] = 1
      endif else begin
        b = 0 ; hazardous to get here, from recent experience (all NaN image)
      endelse

      if initzaps gt 0 then zapimage[crcores] = 1

      if initzaps eq 0 and b eq 0 then begin                 ; make this the last iteration if no CR's found
         if iter lt preiter then iter = preiter   
         if iter gt preiter then iter = maxiter
      endif

      if n_elements(ibad) gt 0 then begin
         if ibad[0] ne -1 then zapimage[ibad] = 1   ; retain zaps from last cycle for percolation.
      endif

      if keyword_set(debug) then begin
        mwrfits, zapimage, 'zapcoreimage'+clip(iter)+'.fits', header, /silent, /create
      endif

      titera = systime(/seconds)

      if keyword_set(verbose) then print, '    Flagged ', strtrim(initzaps,2), ' initial affected pixels before percolation.'

      nperczap = 0
      subiter = 1
      nx = dims[0]
      ny = dims[1]
      while subiter lt 32 and initzaps gt 0 do begin
         nextperc = where(zapimage eq subiter, ct)
         subiter = subiter + 1
         newzaps = 0        
         if ct le 3 then break
         nrays = n_elements(nextperc)
         for c = 0L, nrays-1 do begin
             ci = nextperc[c]
             if ci lt 2*nx or ci gt nx*ny-2*nx then continue
             if ci mod nx lt 3 or ci mod nx ge nx-3 then continue
             coarseblockpos = [         ci-nx,       $
                               ci-1,           ci+1, $
                                        ci+nx         ]
             ;coarseblockpos = [ci-nx-1, ci-nx, ci-nx+1, $
             ;                  ci-1,           ci+1,    $
             ;                  ci+nx-1, ci+nx, ci+nx+1]
             newzap = where(residualnsigmaimage[coarseblockpos] gt subsigma and zapimage[coarseblockpos] eq 0, addzap)
             if addzap gt 0 then begin
                zapimage[coarseblockpos[newzap]] = subiter
                newzaps = newzaps + addzap
             endif
             nperczap = nperczap + addzap
         endfor
         ;if keyword_set(verbose) then print, '    subiteration ', strtrim(subiter,2), ':', newzaps, ' pixel zaps.'
      endwhile
      if keyword_set(verbose) then print, '    Flagged ', strtrim(nperczap,2), ' pixels during ', strtrim(subiter,2), ' percolation steps.'

      ; Finally, zap anything hiding in a cosmic ray "corner" (three neighbors are cosmic ray pixels)
      ;countneighbor = intarr(nx, ny)
      ;nextperc = where(zapimage gt 0, ct)
      ;if ct gt 3 then begin
      ;   nrays = n_elements(nextperc)
      ;   for c = 0L, nrays-1 do begin
      ;       ci = nextperc[c]
      ;       coarseblockpos = [ci-nx-1, ci-nx, ci-nx+1, $
      ;                         ci-1,           ci+1,    $
      ;                         ci+nx-1, ci+nx, ci+nx+1]
      ;       countneighbor[coarseblockpos] = countneighbor[coarseblockpos] + 1
      ;   endfor
      ;   newzap = where(countneighbor ge 4 and zapimage eq 0, newzaps)
      ;   if newzaps gt 0 then begin
      ;      if keyword_set(verbose) then print, '    Flagged an additional ', strtrim(newzaps,2), ' neighboring pixels.'
      ;      zapimage[newzap] = subiter + 1
      ;   endif
      ;endif

      titerb = systime(/seconds)

      if keyword_set(debug) then begin
      mwrfits, zapimage, 'zapimage'+clip(iter)+'.fits', header, /silent, /create
      endif


      ;if (NOT keyword_set(nofluxcompare)) then begin
      ;   ;For every pixel in high-flux regions (bright stars and crays), calculate the ratio
      ;   ; of the local (small-box) median after sky subtraction and the actual value.  If
      ;   ; it is more than 1/fluxcompare, target that pixel for zapping.
      ;   i = where(zapimage NE 0,icount)
      ;   if icount gt 0 then zapimage[i] = (filterimage[i] / residualimage[i]) LT fluxratio
      ;endif

      ;Apply the growth radius by also zapping surrounding pixels, if requested.
      ; For this to function, put it at the end!!
      ;if nrings gt 0 then $
      ; zapimage = smooth(float(zapimage), 1+2*nrings, /edge) GT 1.0E-6

      ;Actually do the zapping by replacing with NaN (to be replaced later)
      ibad = where(zapimage NE 0, nbad)
      if (nbad GT 0) then begin
         img[ibad] = !values.f_nan ;;skyimage[ibad] + filterimage[ibad]
         outmask[ibad] = 1
        ; Recalculate the median without the cosmic rays and replace again, this time for real.
 
      ;Apply a smaller median filter to the subtracted image to get a measure of the local average
      ;;refilterimage = median(img, boxsize)  
      ;;img[ibad] = refilterimage[ibad]
      endif
     
      nbadnew = nbad - prevbad
      if keyword_set(silent) eq 0 and keyword_set(quiet) eq 0 then $
        print, ': Number zapped = ', strtrim(string(nbadnew),2)
      prevbad = nbad
      nzap = nbad ;nzap + nbad

      if (postiter gt 0 and iter eq preiter) or iter eq maxiter then begin    ; replace NaNs
        ibad = where(outmask gt 0, nibad)
        if nibad gt 0 then begin
          if keyword_set(usamp) eq 0 then replacementmodel = modelimage(img) $ ; note slow step here.
                                     else replacementmodel = usampmodelimage(img)
          img[ibad] = replacementmodel[ibad]
        endif
        nans = where(finite(img) eq 0, ct) ; Anything *still* NaN after replacement?
        naniter = 0
        nel = n_elements(img)
        while ct gt 0 do begin
          holdimg = img
          for c = 0L, n_elements(nans)-1 do begin
             ci = nans[c]
             if naniter lt 5 then begin
             coarseblockpos = [ci-nx-1, ci-nx, ci-nx+1, $
                               ci-1,           ci+1,    $
                               ci+nx-1, ci+nx, ci+nx+1] 
             endif else if naniter lt 7 then begin
             coarseblockpos = [ci-2*nx-2, ci-nx*2, ci-nx*2+2,$
                               ci-2,               ci+2,     $
                               ci+2*nx-2, ci+nx*2, ci+nx*2+2] 
             endif else if naniter lt 9 then begin
             coarseblockpos = [ci-20*nx-20,ci-20*nx-10,ci-20*nx,ci-20*nx+10,ci-20*nx+20,$
                               ci-10*nx-20,ci-10*nx-10,ci-10*nx,ci-10*nx+10,ci-10*nx+20,$
                               ci     - 20,ci -     10,         ci     + 10,ci     + 20,$
                               ci+10*nx-20,ci+10*nx-10,ci+10*nx,ci+10*nx+10,ci+10*nx+20,$
                               ci+20*nx-20,ci+20*nx-10,ci+20*nx,ci+20*nx+10,ci+20*nx+20]
             endif else begin
                 img[nans[c]] = 0
                 continue
             endelse
             coarseblockpos = coarseblockpos[where(coarseblockpos gt 0 and coarseblockpos lt nel)]
             img[nans[c]] = median(holdimg[coarseblockpos], /even)

          endfor
          naniter = naniter + 1
          nans = where(finite(img) eq 0, ct)
        endwhile
      endif

     ;print, 'proc time: ', titerb-titera, titera-titerstart

   endwhile

   ; End of zapping iterations.


   ; Final nearest-neighbor growth:

   countneighbor = intarr(nx, ny)
   nextperc = where(outmask gt 0, ct)
   if ct gt 3 then begin
      nrays = n_elements(nextperc)
      for c = 0L, nrays-1 do begin
          ci = nextperc[c]
          coarseblockpos = [ci-nx-1, ci-nx, ci-nx+1, $
                            ci-1,           ci+1,    $
                            ci+nx-1, ci+nx, ci+nx+1]
          countneighbor[coarseblockpos] = countneighbor[coarseblockpos] + 1
      endfor
      newzap = where(countneighbor ge 4 and outmask eq 0, newzaps)
      if newzaps gt 0 then begin
         if keyword_set(verbose) then print, '    Flagged an additional ', strtrim(newzaps,2), ' neighboring pixels.'
         outmask[newzap] = 1
      endif
   endif
   if keyword_set(quiet) eq 0 and keyword_set(silent) eq 0 then $
     print, '  Zapped ', strtrim(string(newzaps),2), ' additional neighboring pixels.'

   if ninitbad gt 0 then outmask[initbad] = 1

   ibad = where(outmask gt 0, ct)
   ; do the splotching
   if keyword_set(zero) then begin
       img[ibad] = 0
   endif else begin
       if keyword_set(usamp) eq 0 then replacementmodel = modelimage(img) $ ; note slow step here.
                                   else replacementmodel = usampmodelimage(img)
       if ct gt 0 then img[ibad] = replacementmodel[ibad]
   endelse

   if keyword_set(silent) eq 0 then $
   print, '  Zapped ', strtrim(string(n_elements(ibad)-ninitbad),2), ' total affected pixels (', strtrim(string((n_elements(ibad)-ninitbad)*100./(1.0*nx*ny),format='(F6.3)'),2), '% of total).'
   ; Numbers don't add up, but whatever

   sxaddpar, header, 'PZAPZEAL', zealfactor
   if keyword_set(usamp) then modeltype = 'usamp' else modeltype = 'reg'
   sxaddpar, header, 'PZAPMDL', modeltype
   sxaddpar, header, 'NPZAP', n_elements(ibad)-ninitbad, 'Num. of pixels zapped by pzap'

   get_date, now
   hist = 'Processed by pzap '+now
   sxaddhist, hist, header

   if (size(outname, /tname) EQ 'STRING') then begin
      if n_elements(outmaskname) gt 0 then sxaddpar, header, 'BPM', outmaskname, 'Bad pixel mask'
      mwrfits, img, outname, header, /create
   endif else begin
      outname = img
   endelse

   if (size(outmaskname, /tname) EQ 'STRING') then begin
      mwrfits, outmask, outmaskname, header, /create
   endif ;else begin
      ;mask = outmask
   ;endelse
   if (size(outweightname, /tname) EQ 'STRING') then begin
      mwrfits, 1-outmask, outweightname, header, /create
   endif ;else begin
      ;weight = 1-outmask
   ;endelse

end


; Improvements: reduce the number of model calculations to optimize speed.
; To get really fancy, break the image into sections to do this.

