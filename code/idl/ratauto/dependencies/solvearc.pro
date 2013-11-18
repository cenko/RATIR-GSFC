function solvearc, linex, brightlinesx, reflines, brightreflines, cwavrange=cwavrange, disprange=disprange, maxsep=maxsep, uncrefwav=uncrefwav, initsol=initsol, nx=nx, nmatch=nmatch, ngoodmatch=ngoodmatch, nofit=nofit, maxorder=maxorder, maxnudge=maxnudge

         ; linex = detected line list of pixel centers
         ; brightlinesx = brightest ~10 lines
         ; reflines = good reference lines, including faint ones
         ; brightreflines = bright reference lines for coarse matching

         ; This routine works by looking for trios of strong lines in the arc spectrum and in the reference list.
         ; The ratio of the spacing between (center line to furthest line)/(center line to nearer line) is calculated,
         ; which is a scale-independent quanity.  Every possible trio is matched to every other possible trio and this
         ; ratio is tabulated, then a simple routine looks for a cluster of points around a similar value.  
         ; Using lines in matched trios the dispersion and wavelength at pix=0 are then estimated.

      debug = 0

      brightcol = 5

      !p.position = 0
      !p.multi = [0, 1, 4]

      brightlinesx = brightlinesx[sort(brightlinesx)]
      brightreflines = brightreflines[sort(brightreflines)]

      isbright = intarr(n_elements(linex))
      for l = 0, n_elements(linex)-1 do begin
         bmatch = where(linex[l] eq brightlinesx, ct)
         if ct gt 0 then isbright[l] = 1
      endfor

      if debug ge 2 then begin
        for l = 0, n_elements(brightlinesx)-1 do  print, l, brightlinesx[l]
        for l = 0, n_elements(brightreflines)-1 do  print, l, brightreflines[l]
      endif

     if n_elements(nx) eq 0 then $
        plotxr = [0, max(linex)+10] $ ;min(linex)-10
     else $
        plotxr = [0,nx]
     x = plotxr[0] + (plotxr[1]-plotxr[0])*findgen(120)/120.
     dx = x - nx/2.
     xl = 0.03
     xr = 0.99
     ybot = reverse([0.05, 0.29, 0.53, 0.77])
     ytop = ybot + 0.21

     if n_elements(initsol) eq 0 then begin    
               ; if no highly-reliable (within 10 Angstroms everywhere) solution is known, use a
               ; pattern-matching technique to determine it

         automatchlines, brightlinesx, brightreflines, bestxmatches, bestwavmatches, $
                         cwavrange=cwavrange, disprange=disprange, maxsep=maxsep, nx=nx, $
                         debug=debug, nblock=1, outcwav=cwav, outdispersion=dispersion
                          

         brightlinesdx = brightlinesx - nx/2.
         linearfit = poly_fit(brightlinesdx[bestxmatches], brightreflines[bestwavmatches], 1, yfit=yfit)
         !p.position = [xl, ybot[0], xr, ytop[0]]
         ; ----- Plot 1 -------
         plot, brightlinesx[bestxmatches], brightreflines[bestwavmatches]-poly(brightlinesdx[bestxmatches],linearfit), xrange=plotxr, /xstyle, /nodata
         oplot, brightlinesx[bestxmatches], brightreflines[bestwavmatches]-poly(brightlinesdx[bestxmatches],linearfit), psym=7, color=brightcol
         for ll = 0, n_elements(brightlinesx)-1 do begin
           oplot, [brightlinesx[ll], brightlinesx[ll]], [-1, 1], color=brightcol
         endfor
         
         outlieriter = 0
         if outlieriter gt 0 then begin
           for outlieriter = 1, 3 do begin 
           nl = n_elements(bestxmatches)
            localresidual = fltarr(nl)
           localdev = fltarr(nl)
           for l = 0, nl-1 do begin
              checklines = (indgen(nl)) [(l-3-(l eq nl-1)) > 0: (l+3+(l eq 0)) < (nl-1)]   ; e.g.  2,3,4,5,6
              checklines = checklines[where(checklines ne l)]   ; e.g.  2,3, 5,6
              checkfit = poly_fit(brightlinesdx[bestxmatches[checklines]], brightreflines[bestwavmatches[checklines]], 1, yfit=yfit)
              localresidual[l] = poly(brightlinesdx[bestxmatches[l]], checkfit) - brightreflines[bestwavmatches[l]]
              localdev[l] = stdev(yfit - brightreflines[bestwavmatches[checklines]])
              if debug gt 0 then print, l, localresidual[l], localdev[l], abs(localresidual[l]/localdev[l])
           endfor
           if debug gt 0 then print, abs(localresidual/(localdev > 1.))
           good = where(abs(localresidual/(localdev > 1.)) lt 30, ct)
           ;bestmatches = bestmatches[good]
           bestxmatches = bestxmatches[good]
           bestwavmatches = bestwavmatches[good]
           endfor
         endif


         if n_elements(bestxmatches) gt 3 then begin
   
            ;print, 'Matched ', clip(n_elements(bestxmatches)), ' lines'

            ; remove outliers one by one
            ; these settings are for the INITIAL fit to the bright lines!
                              ;dispersion*0.75
            targetdev = 1.0*uncrefwav+dispersion*0.25      ; threshold to increase the order of the fit
            initreject = 10.                        ; threshold to continue rejecting individual outliers
                                 ;;dispersion*2
            finalreject = 2.0*uncrefwav+dispersion*0.5 ; largest discrepancy we will tolerate in the final solution 
                                                ; (rejectthreshold slowly drops to this)
            alwaysreject = 50.

            prelimfit = polyfititerauto(brightlinesdx[bestxmatches], brightreflines[bestwavmatches], dxout, refwavout, $
                            targetdev=targetdev, maxorder=maxorder, $
                            initreject=initreject, finalreject=finalreject, alwaysreject=alwaysreject)
            prelimorder = n_elements(prelimfit)-1

            ;print, '   Initial wavelength solution complete.'
            ;if noutlier gt 0 then print, 'Removed ', clip(noutlier), ' lingering outlier line matches.'


         endif else begin
            print, '   Bright line match failed.'
            print, 'WARNING - The line matcher was unable to match bright lines in the arc/sky '
            print, '          spectrum with those in the reference list.  It is extremely unlikely'
            print, '          that you will be able to get an accurate wavelength calibration for'
            print, '          this setup.   You can type ".continue" to proceed anyway but the '
            print, '          wavelength solution may diverge greatly from reality.  Consider '
            print, '          inspecting the arc file to be sure it was properly warmed up and'
            print, '          more than a few lines were not saturated.'
            stop
            prelimorder = 1
            if n_elements(wav0) eq 1 then prelimfit = [wav0, dispersion] else prelimfit = initsol
            bestxmatches = xmatches
            bestwavmatches = wavmatches
         endelse
    
         if n_elements(bestxmatches) eq 0 or n_elements(bestwavmatches) eq 0 or n_elements(linearfit) eq 0 then $
           stop
         oplot, dxout+nx/2., refwavout-poly(dxout,linearfit), psym=1, color=brightcol
         oplot, x, poly(dx,prelimfit)-poly(dx,linearfit), linestyle=1
    
      endif else begin

         ; Even if solution is known from a previous frame, flexure can cause huge horizontal shifts.
         ; Try to correct the 0th order term using bright lines before proceeding to the fit.

        if n_elements(maxnudge) eq 0 then maxnudge = 50.
        wavnudge = (findgen(41)-20)*maxnudge/20.
        ;xnudge = (findgen(20)-10))*5

        brightlinesdx = brightlinesx-nx/2.
 
        bestmatchdiff = !values.f_nan
        nmostmatches = -1 ;  not used???
        mostmatches = 0
        for n = 0, n_elements(wavnudge)-1 do begin
          nudgeinitsol = initsol
          nudgeinitsol[0] += wavnudge[n]
          predwavbright = poly(brightlinesdx, nudgeinitsol)
          brightmatchoffset = fltarr(n_elements(brightlinesdx)) 
          matchx = intarr(n_elements(brightlinesdx)) 
          matchwav = intarr(n_elements(brightlinesdx)) 
          m = 0
          for i = 0, n_elements(brightlinesdx)-1 do begin
                dev =    (brightreflines - predwavbright[i])
             absdev = abs(brightreflines - predwavbright[i])
             closest = (where(absdev eq min(absdev))) [0]
             if absdev[closest] lt 10. then begin
                ;if debug gt 0 then print, clip(brightlinesx[i]), ' @', clip(predwavbright[i]), ' matches ', clip(brightreflines[closest]), '  ', dev[closest]
                matchx[m] = i
                matchwav[m] = closest
                brightmatchoffset[m] = dev[closest]
                m = m + 1
             endif
          endfor
          avmatchdiff = median(abs(brightmatchoffset-median(brightmatchoffset)))
          ;if m gt 0 then print, wavnudge[n], m, '|', brightmatchoffset[0:m-1]
          if (m eq mostmatches and avmatchdiff lt bestmatchdiff) or (m gt mostmatches) then begin
            mostmatches = m

            bestmatchx = matchx[0:m-1]     ; this is only actually needed for plotting.  Full match comes later.
            bestmatchwav = matchwav[0:m-1]

            brightmatchoffset = brightmatchoffset[0:m-1]
            brightmodeoffset = clustermed(brightmatchoffset)
            prelimfit = nudgeinitsol
            prelimfit[0] += brightmodeoffset

            prelimorder = n_elements(prelimfit)

            dispersion = prelimfit[1]
          endif 

        endfor

        if mostmatches gt 0 then begin
           print, '  Nudging solution by ', prelimfit[0]-initsol[0]
        endif else begin
           print, '  Found no line matches to nudge wavelength solution.'
           return, initsol
        endelse

        !p.position = [xl, ybot[0], xr, ytop[0]]
        ; ----- Plot 1 (initsol was given, no automatching) -------

        linearfit = initsol[0:1]
        plot, brightlinesx[bestmatchx], brightreflines[bestmatchwav]-poly(brightlinesdx[bestmatchx],initsol[0:1]), /nodata, /xstyle, xrange=[0,nx-1]
        oplot, brightlinesx[bestmatchx], brightreflines[bestmatchwav]-poly(brightlinesdx[bestmatchx],initsol[0:1]), psym=7, color=brightcol
        oplot, x, poly(dx,initsol)-poly(dx,linearfit), linestyle=1, color=1
        oplot, x, poly(dx,prelimfit)-poly(dx,linearfit), linestyle=1, color=0
        
      endelse




      ; might want to iterate here in the future


      ; Now use ALL the detected lines, including faint ones.
      ; Match by simply finding the line nearest to the predicted curve.

      linedx = linex-nx/2.

      predwav = poly(linedx, prelimfit)
      matchx = intarr(n_elements(linedx)) 
      matchwav = intarr(n_elements(linedx)) 
      matchoffset = fltarr(n_elements(linedx)) 

      matchsearch = 10.
      m = 0
      for i = 0, n_elements(linedx)-1 do begin
            dev =    (reflines - predwav[i])
         absdev = abs(reflines - predwav[i])
         closest = (where(absdev eq min(absdev))) [0]
         if absdev[closest] lt matchsearch then begin  ; formerly 5, but you should be able to use the mode feature to get first-order corr right away in refinement mode
            matchx[m] = i
            matchwav[m] = closest
            if debug gt 0 then print, clip(linex[i]), ' @', clip(predwav[i]), ' matches ', clip(reflines[closest]), '  ', dev[closest]
            matchoffset[m] = dev[closest]
            m = m + 1
         endif else begin
            if debug gt 0 then print, clip(linex[i]), ' @', clip(predwav[i]), '    -   ', dev[closest], ' from ',clip(reflines[closest])
         endelse
      endfor

      if m lt 3 + n_elements(initsol) and keyword_set(nofit) eq 0 then begin
         if n_elements(initsol) gt 0 then begin
            print, 'WARNING: Insufficient line matches in full line list for precise solution.'
            print, 'Not attempting to provide a refined solution.'
            return, initsol
         endif else begin
            print, 'Cannot verify the solution.'
            print, 'The line-matching routine probably failed.  Check that the arc / sky frame is valid,'
            print, ' and that the correct line list is used.
         endelse
      endif

      matchx = matchx[0:m-1]
      matchwav = matchwav[0:m-1]
      matchoffset = matchoffset[0:m-1]
      modeoffset = clustermed(matchoffset) ; not used, but might be useful for later development


     ; Compare to the linear solution (for plotting)
     ;;;if n_elements(linearfit) eq 0 then linearfit = poly_fit(linedx[matchx], reflines[matchwav], 1, yfit=yfit)
     !p.position = [xl, ybot[1], xr, ytop[1]]
     ; ----- Plot 2 -------
     ran = [linex[matchx],reflines[matchwav]-poly(linedx[matchx],linearfit)]
     yrange = [min(ran),max(ran)]
     plotc, linex[matchx], reflines[matchwav]-poly(linedx[matchx],linearfit), psym=7, xrange=plotxr, /xstyle, color=isbright[matchx]*brightcol

     if keyword_set(nofit) eq 0 then begin
     ; Perform the final fit for the wavelength solution using all matched lines.  
   
        finalfit = polyfititerauto(linedx[matchx], reflines[matchwav], dxfinal, wavfinal, startorder=prelimorder, maxorder=maxorder,$
                         alwaysreject=5., targetdev=0.5*uncrefwav+dispersion*0.3,  $
                         initrejectdev=4., finalrejectdev=uncrefwav)
        finalorder = n_elements(finalfit)-1
        ; this isn't quite the final fit - it still (might) have duplicates

        residual = reflines[matchwav]-poly(linedx[matchx],finalfit) ; get the residual
        removeduplicatematches, matchx, matchwav, residual          ; use it to remove duplicates optimally

        ; now redo without the duplicates
        finalfit = polyfititerauto(linedx[matchx], reflines[matchwav], dxfinal, wavfinal, startorder=finalorder, maxorder=maxorder, $
                         alwaysreject=5., targetdev=0.5*uncrefwav+dispersion*0.3, $
                         initrejectdev=4., finalrejectdev=uncrefwav+dispersion*0.3)
        finalorder = n_elements(finalfit)-1
     endif else begin
        finalfit = prelimfit 
        print, '  Nudging by an additional ', modeoffset  ; not sure I trust modeoffset enough

        ; It would be nice to do a rejection iteration and repeat, but for now it doesn't matter.
        finalfit[0] += modeoffset
        dxfinal = linedx[matchx]
        wavfinal = reflines[matchwav]
     endelse


    residual = reflines[matchwav]-poly(linedx[matchx],finalfit) ; again

    xfinal = dxfinal + nx/2.

    ;oplotc, linex[matchx], reflines[matchwav]-poly(linedx[matchx],linearfit), psym=7, color=isbright*brightcol
    isbrightfinal = intarr(n_elements(dxfinal))
    for l = 0, n_elements(dxfinal)-1 do begin
       bmatch = where(dxfinal[l] eq brightlinesdx, ct)
       if ct gt 0 then isbrightfinal[l] = 1
    endfor
    oplotc, xfinal, wavfinal-poly(dxfinal,linearfit), psym=1, color=isbrightfinal*brightcol
    oplot, x, poly(dx,finalfit)-poly(dx,linearfit), linestyle=1

    residualfinal = wavfinal - poly(dxfinal,finalfit)
    !p.position = [xl, ybot[2], xr, ytop[2]]
    ; ------- Plot 3 ----------
    plotc, linex[matchx], residual, psym=7, xrange=plotxr, /xstyle, yrange=[-1,1]*(1.02*max(abs(residual)) > 1.2*uncrefwav), color=isbright[matchx]*brightcol, /ysty
    oplotc, xfinal, residualfinal, psym=1, color=isbrightfinal*brightcol
    oplot, x, x*0, linestyle=1

    if debug gt 0 then begin
       for i = 0, n_elements(matchx)-1 do begin
          bad = ''
          if abs(residual[i]) gt 2.*uncrefwav+0.2*dispersion then bad = 'x'
          print, linex[matchx[i]], reflines[matchwav[i]], residual[i], '   ', bad
       endfor
    endif

    print, '  Wavelength solution: '      , format='($,A)'
    for o = 0, n_elements(finalfit)-1 do begin
       print, clip(finalfit[o]), format='($,A)'
       if o eq 1 then print, 'dx', format='($,A)'
       if o ge 2 then print, 'dx^'+clip(o), format='($,A)'
       if o ne n_elements(finalfit)-1 then print, ' + ' , format='($,A)'
    endfor
    nmatch = n_elements(matchx)
    goodthresh = 2.0*uncrefwav+0.2*median(disprange)
    w = where(abs(residualfinal) lt goodthresh, ngoodmatch)
    print
    print, '     '+clip(nmatch)+'/'+clip(n_elements(linex))+' matched lines, '$
           +clip(ngoodmatch)+' within '+clip(fpr(goodthresh,2.2))+' Ang tolerance'
    print, '     Median absolute residual ', fpr(median(abs(residualfinal)),3.3), ' Ang '
 
    if ngoodmatch lt fix(0.3*nmatch)-2 then begin
       print, '     WARNING - This looks like a poor (wrong) wavelength solution.'
    endif

    !p.position = [xl, ybot[3], xr, ytop[3]] ; final plot will be external

    return, finalfit
end


