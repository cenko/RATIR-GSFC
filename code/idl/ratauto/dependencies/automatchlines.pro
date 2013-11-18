pro automatchlines, obsxin, refwavin, outxmatches, outwavmatches, cwavrange=cwavrange, disprange=disprange, maxsep=maxsep, nblocks=nblocks, debug=debug, nx=nx, outcwav=outcwav, outdispersion=outdispersion

    if n_elements(nblocks) eq 0 then nblocks = 1 ; use of nblocks > 1 is not recommended.  Use maxsep instead?
    
    sortx = sort(obsxin)
    obsx = obsxin[sortx]
    sortwav = sort(refwavin)
    refwav = refwavin[sortwav]

    ; now assumes that they're sorted, but should not assume this, should sort internally

    obsdx = obsx - nx/2.  ; relative to center (now doing cwav style)

    if n_elements(debug) eq 0 then debug = 0

  ; "Fully" automated line matching.

    if n_elements(cwavrange) eq 0 then cwavrange = [4500,8500]
    if n_elements(disprange) eq 0 then disprange = [0.2, 3.0]
    if n_elements(maxsep) eq 0 then maxsep = 2000.
  
    blockxs =  (findgen(nblocks+1)/nblocks)*nx
    blockwavs = median(cwavrange,/even)-(nx/2.)*mean(disprange) + blockxs*median(disprange,/even)
    wavunc = max(cwavrange)-min(cwavrange)
   
    outxmatches = [-1]
    outwavmatches = [-1]
    dispersions = fltarr(nblocks)
    cwavs = fltarr(nblocks)

    for b = 0, nblocks-1 do begin

         wx = where(obsx ge blockxs[b] and obsx lt blockxs[b+1], ct)
         if ct lt 3 then continue
         obsdxb = obsdx[wx]
         obsxb = obsdxb + nx/2.
         ww = where(refwav ge blockwavs[b]-wavunc and refwav lt blockwavs[b+1]+wavunc, ct)
         if ct lt 3 then continue
         refwavb = refwav[ww]

         if disprange[1] lt disprange[0] then disprange = disprange[[1,0]]

         order = 3

         ndetlines = n_elements(obsdxb)
         xdiffrat = fltarr(ndetlines, ndetlines, ndetlines)
         for i = 1, ndetlines-2 do begin
            for j = 0, i-1 do begin
               for k = i+1, ndetlines-1 do begin
                  rat = (1.0*obsdxb[j]-obsdxb[i])/(1.0*obsdxb[k]-obsdxb[i])    ; this is biased against small k-i
                  xdiffrat[i,j,k] = rat
               endfor
            endfor
         endfor

         nreflines = n_elements(refwavb)
         wavdiffrat = fltarr(nreflines, nreflines, nreflines)
         for l = 1, nreflines-2 do begin
            for m = 0, l-1 do begin
               for n = m+1, nreflines-1 do begin
                  rat = (refwavb[m]-refwavb[l])/(refwavb[n]-refwavb[l])
                  if abs(refwavb[m]-refwavb[n]) lt maxsep then $
                     wavdiffrat[l,m,n] = rat  ; ensures that lines are within maxsep Angstroms
               endfor
            endfor
         endfor

         maxmatch = 1024
         candidatematchi = intarr(maxmatch)
         candidatematchj = intarr(maxmatch)
         candidatematchk = intarr(maxmatch)
         candidatematchl = intarr(maxmatch)
         candidatematchm = intarr(maxmatch)
         candidatematchn = intarr(maxmatch)
         candidatematchdispersion = fltarr(maxmatch)
         candidatematchcwav = fltarr(maxmatch)      
         candidatematchsoft = fltarr(maxmatch) ; soft is the discrepancy between ratios in x and wav

         softthresh = 0.06  ;0.01, then 0.025
         bestsoft = !values.f_infinity

         if debug gt 0 then openw, 4, 'matchlog.txt'

         alldisp = fltarr(ndetlines^3/4.*nreflines^3/4.)
         allcwav = fltarr(ndetlines^3/4.*nreflines^3/4.)
         allsoft = fltarr(ndetlines^3/4.*nreflines^3/4.)
         a = 0L

         nmatch = 0
         for i = 1, ndetlines-2 do begin
            for j = 0, i-1 do begin
               for k = i+1, ndetlines-1 do begin         ; that's right, 6 nested loops in IDL
                  if xdiffrat[i,j,k] eq 0 then continue
                  for l = 1, nreflines-2 do begin
                     for m = 0, l-1 do begin
                        for n = l+1, nreflines-1 do begin
                           if wavdiffrat[l,m,n] eq 0 then continue
                           ; would make much more sense to measure (dispersion_ji - dispersion_ki) / dispersion_jk
                           soft = abs(((xdiffrat[i,j,k]-wavdiffrat[l,m,n])/ $
                                   ((xdiffrat[i,j,k]+wavdiffrat[l,m,n])/2.)))
                           dispersion = (refwavb[n]-refwavb[m]) / $
                                         (obsdxb[k]-obsdxb[j])
                           cwav = refwavb[l] - obsdxb[i]*dispersion
                           alldisp[a] = dispersion
                           allcwav[a] = cwav
                           allsoft[a] = soft 
                           a += 1
                           ijkstr = strmid(rclip(j,2)+'-'+clip(i)+'-'+clip(k,2)+'     ',0,9)
                           lmnstr = strmid(rclip(m,2)+'-'+clip(l)+'-'+clip(n,2)+'     ',0,9)
                           if soft lt bestsoft then bestsoft = soft
                           if soft lt softthresh then begin
                                ; possible match

                              if dispersion gt disprange[0] and dispersion lt disprange[1] and $
                                    cwav gt cwavrange[0] and cwav lt cwavrange[1] then begin
                                 candidatematchdispersion[nmatch] = dispersion
                                 candidatematchcwav[nmatch] = cwav
                                 candidatematchi[nmatch] = i
                                 candidatematchj[nmatch] = j
                                 candidatematchk[nmatch] = k
                                 candidatematchl[nmatch] = l
                                 candidatematchm[nmatch] = m
                                 candidatematchn[nmatch] = n
                                 candidatematchsoft[nmatch] = soft

                                 nmatch = nmatch + 1

                                 if keyword_set(debug) then printf,4,  ijkstr+' '+lmnstr+' '+'|'+$
                                                                       rclip(fix(obsxb[j]),4)+' '+rclip(fix(obsxb[i]),4)+' '+rclip(fix(obsxb[k]),4)+'|'+$
                                                                       clip(fix(refwavb[m]),5)+clip(fix(refwavb[l]),5)+clip(fix(refwavb[n]),5)+$
                                                                   '|'+ fpr(soft,2.4)+ ' '+ fpr(dispersion,3.4)+ fpr(cwav,6.2)+ ' * '
                                                                    ;clip(i,3)+clip(j,3)+clip(k,3)+clip(l,3)+clip(m,3)+clip(n,3)+$
                                 if nmatch ge maxmatch then begin 
                                    print, '  Maximum possible matches reached!! Cycling back to zero (may result in bad solution)'
                                    nmatch = 0
                                 endif
                              endif else begin
                                 if debug ge 2 then printf,4,  ijkstr+' '+lmnstr+' '+'|'+$
                                                               rclip(fix(obsxb[j]),4)+' '+rclip(fix(obsxb[i]),4)+' '+rclip(fix(obsxb[k]),4)+'|'+$
                                                               clip(fix(refwavb[m]),5)+clip(fix(refwavb[l]),5)+clip(fix(refwavb[n]),5)+$
                                                             '|'+ fpr(soft,2.4)+ ' '+ fpr(dispersion,2.4)+ fpr(cwav,6.2)+ ' ~'
                              endelse 
                           endif else begin
                              if debug ge 3 then printf,4,  ijkstr+' '+lmnstr+' '+'|'+$
                                                            rclip(fix(obsxb[j]),4)+' '+rclip(fix(obsxb[i]),4)+' '+rclip(fix(obsxb[k]),4)+'|'+$
                                                             clip(fix(refwavb[m]),5)+clip(fix(refwavb[l]),5)+clip(fix(refwavb[n]),5)+$
                                                             '|'+fpr(soft,2.4)
                           endelse
                        endfor
                     endfor
                  endfor
               endfor
            endfor
         endfor

         allsoft = allsoft[0:a-1]
         allcwav = allcwav[0:a-1]
         alldisp = alldisp[0:a-1]

         if debug gt 0 then close, 4
         if debug gt 9 then begin
            s2 = where(allsoft ge softthresh and allsoft lt softthresh*3)
            s1 = where(allsoft lt softthresh and allsoft gt 0.)
             psclose
             window, 0, xsize=800, ysize=800
             !p.position = 0
             !p.multi = 0
             plot, [0], [0], xrange=cwavrange+[-100,100], yrange=disprange+[-0.05,0.05], /xsty, /ysty
             oplot, allcwav[s2], alldisp[s2], psym=7
             oplot, allcwav[s1], alldisp[s1], psym=1
             stop
         endif

         if nmatch eq 0 then begin
           continue
           ;print, 'NO candidate line matches found!!'
           ;;print, 'Best match has softness ', bestsoft
           ;print, 'Check arc file or input dispersion range.'
         endif

         candidatematchdispersion = candidatematchdispersion[0:nmatch-1]
         candidatematchcwav = candidatematchcwav[0:nmatch-1]

         dispersion = clustermed(candidatematchdispersion)
         cwav = clustermed(candidatematchcwav)               ; even better would be to do this in 2D

         if debug gt 0 then begin
            print, 'Guessing the dispersion is ', clip(dispersion)
            print, 'Guessing cwav is ', clip(cwav)
         endif

         goodmatch = where(abs(candidatematchdispersion - dispersion)/abs(dispersion) lt 0.2 and $
                           abs(candidatematchcwav - cwav)/abs(cwav) lt 0.2)
         ; this step is unnecessary if disprange and cwavrange were actually set.


         ; collect all the indices together.  There are probably lots of repeats.
         xmatchindex = [candidatematchi[goodmatch], candidatematchj[goodmatch], candidatematchk[goodmatch]]
         wavmatchindex = [candidatematchl[goodmatch], candidatematchm[goodmatch], candidatematchn[goodmatch]]

         ; Find the "real" matches by majority vote (which wav-index corresponds most often to each x-index)
         xmatches = unique(xmatchindex) ; list of x indices that actually match lines, no repeats
         nm = n_elements(xmatches)
         wavmatches = -1 + intarr(nm)  ; list of wav indices corresponding to above
         nmatches = -1 + intarr(nm)
         nothermatches = -1 + intarr(nm)
         victorymargin = -1 + intarr(nm)
         for iii = 0, n_elements(xmatches)-1 do begin
            wm = where(xmatchindex eq xmatches[iii], ct)  ; wm: xmatchindex subscripts for this unique line in x
            if ct eq 0 then continue
            wmw = wavmatchindex[wm]                       ; wmw: wavmatchindex subscripts " "   "     "     " wav
            wmwhist = histogram(wmw,min=0)                ; count of nmatches for each index [e.g. 000010007002]
            wmwhist = [wmwhist, 0, 0]                     ; enables 2nd place finder below
            wmwmax = max(wmwhist)                         ; most matches
            wmw2nd = wmwhist[(reverse(sort(wmwhist))) [1]] ; 2nd most matches
            if max(wmwhist) ge 2 then begin
               wavmatches[iii] = (where(wmwhist eq wmwmax)) [0]
               nmatches[iii] = wmwmax
               nothermatches[iii] = total(wmwhist)-wmwmax
               victorymargin[iii] = wmwmax-wmw2nd

               if debug gt 0 then $
                   print, iii, obsdxb[xmatches[iii]],  refwavb[wavmatches[iii]], $
                               nmatches[iii], nothermatches[iii], victorymargin[iii], '|', wmwhist[where(wmwhist gt 0)]
            endif
         endfor


         bettermatches = where(wavmatches ge 0, nbtm)
         betterxmatches = xmatches[bettermatches]
         betterwavmatches = wavmatches[bettermatches]
         bettervictorymargin = victorymargin[bettermatches] ; use this one for now
         

         ; search for conflicts the other way (multiple x's match to the same wav)
         uwavmatches = unique(betterwavmatches)       
         bestmatchi = intarr(n_elements(uwavmatches))
         for iu = 0, n_elements(uwavmatches)-1 do begin
           wm = where(betterwavmatches eq uwavmatches[iu], ct)
           if ct ge 2 then begin
                bestwm = (where(bettervictorymargin[wm] eq max(bettervictorymargin[wm]))) [0]
           endif else begin
                bestwm = wm
           endelse
           bestmatchi[iu] = wm[bestwm]
         endfor
         bestorder = sort(obsdxb[betterxmatches[bestmatchi]])
         bestmatches = (bettermatches[bestmatchi]) [bestorder]
         bestxmatches = (betterxmatches[bestmatchi]) [bestorder]
         bestwavmatches = (betterwavmatches[bestmatchi]) [bestorder]

         outxmatches = [outxmatches, wx[bestxmatches]]
         outwavmatches = [outwavmatches, ww[bestwavmatches]]
         dispersions[b] = dispersion
         cwavs[b] = cwav

  endfor

  outxmatches = outxmatches[1:*]
  outwavmatches = outwavmatches[1:*]

  outdispersion = median(dispersions,/even)
  outcwav = median(cwavs, /even)

  if debug gt 0 then begin
       print, 'Final bright match list:'
       for iv = 0, n_elements(bestmatches)-1 do begin
          print, obsdxb[bestxmatches[iv]]+nx/2., '=', refwavb[bestwavmatches[iv]]
       endfor
    endif

  outxmatches = sortx[outxmatches]
  outwavmatches = sortwav[outwavmatches]   ; un-sort

  debug = 0
  
end


