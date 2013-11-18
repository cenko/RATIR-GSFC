; +
; NAME:
;   solvearcspec.pro
;
; PURPOSE:
;   Fully-automated solution of an arc spectrum.
;
; CALLING SEQUENCE:
;
; INPUTS:
;   farc               The arc/sky spectrum flux(x)
;   reflinewav         List of reference line wavelengths (array)
;   brightrefwav   List of bright reference line wavelengths (array)
;   x0                 Central x value for polynomial fitting, usually nx/2
;   cwav               Estimated wavelength at x0
;                       (Return value will be improved estimate)
;   unccwav            Uncertainty in estimate; scalar or 2-element vector
;   disp               Estimated dispersion
;                       (Return value will be improved estimate)
;   uncdisp            Uncertainty in estimate; scalar or 2-element vector
;                       WARNING: unccwav and uncdisp include physical variation
;                       across chip, so don't make too small!
;   maxsep             Maximum separation between wavelength triplets during
;                       automated solving
;   uncrefwav          Uncertainty in listed line positions.  Does not take
;                       centroiding error in account, that is dealt with 
;                       internally right now (potentially not well?)
;   initsol            Initial solution.  Will skip the automated line search
;                       if this is set and identify lines from this solution
;                       plus a wavelength shift. Set only if you are confident
;                       the >=1st order terms are close to accurate.
;   nofit              Only solve the offset term (e.g. in sky line refinement)
;   nbrightline        Number of lines to use in automated matching.
;   checkplotfile      File to print out verification plots.  Plots to active
;                       device if this is not set.
;   leaveplotopen      Don't close the checkplot, allowing future runs to
;                       add new pages.
;
; OUTPUTS:
;   wav                Solved wavelength of each pixel.
;   nmatch             Number of "likely" line matches in final solution.
;   ngoodmatch         Number of these matches that meet the uncrefwav tolerance


function solvearcspec, farc, wav, refwav=refwav, brightrefwav=brightrefwav, $
       x0=x0, cwav=cwav, unccwav=unccwav, dispersion=dispersion, uncdisp=uncdisp, $
       maxsep=maxsep, uncrefwav=uncrefwav, initsol=initsol, nofit=nofit, nbrightline=nbrightline, $
       nmatch=nmatch, ngoodmatch=ngoodmatch, checkplotfile=checkplotfile, $
       leaveplotopen=leaveplotopen, maxorder=maxorder, maxnudge=maxnudge
 
      nx = n_elements(farc)
      pix = findgen(nx)

      if n_elements(x0) eq 0 then x0 = nx/2. ;

      disprange = dispersion+[-1.,1.]*abs(uncdisp)
      cwavrange  = cwav + [-1.,1.]*abs(unccwav)
      minwav = cwavrange[0]-(x0- 0)*max(disprange)
      maxwav = cwavrange[1]+(nx-x0)*max(disprange)

      wavinrange = where(brightrefwav gt minwav and brightrefwav lt maxwav, ct)
      brightrefwav = brightrefwav[wavinrange]
      wavinrange = where(refwav gt minwav and refwav lt maxwav, ct)
      refwav = refwav[wavinrange]


      if n_elements(checkplotfile) gt 0 then begin
      if checkplotfile ne '' then begin
         !p.position = 0
         !p.multi = [0, 2, 1]
         psopen, checkplotfile, xsize=7, ysize=5, /inch
         device, /color
         colors = transpose([[0,0,0],$
                       [170,170,200],$  ; 1 = grey          (reference line)
                       [210,210,220],$  ; 2 = light grey    (faint reference line)
                       [64, 158, 64],$  ; 3 = green         (faint detected line)
                       [112,112, 64],$ 
                       [0,    0,192],$  ; 5 = blue          (strong detected line)
                       [192, 64, 64],$  ; 6 = reddish       (strong reference line)
                       [128, 64, 128]$
                    ])
         tvlct, colors
      endif
      endif


      linex = findlines(farc, farcsub, linestrength=linestrength)
      if n_elements(linex) lt 3 then begin
          ; not many lines -> sky in twilight?  try again.
          linex = findlines(farc, farcsub, linestrength=linestrength, thresh=5.5)
      endif
      if n_elements(linex) eq 1 and linex[0] eq -1 then begin
         ; try again if no lines found
         linex = findlines(farc, farcsub, linestrength=linestrength, thresh=4)
         if n_elements(linex) lt 3 and n_elements(initsol) le 1 then begin
            print, '  ERROR: Insufficient lines to solve the arc.'  
            return, [0]
         endif
         if n_elements(linex) eq 1 and linex[0] eq -1 then begin
            print, '  WARNING: Spectrum has no lines.  Cannot refine initial wavelength solution.'
            arcsol = initsol
         endif
      endif

      if n_elements(linex) gt 0 then begin
         brightlinex = findbrightest(linex, linestrength, nbrightline)
         arcsol = solvearc(linex, brightlinex, refwav, brightrefwav, $
                            cwavrange=cwavrange, disprange=disprange, maxsep=maxsep, initsol=initsol, nx=nx, uncrefwav=uncrefwav, $
                            nmatch=nmatch, ngoodmatch=ngoodmatch, nofit=nofit, maxorder=maxorder, maxnudge=maxnudge)
         initsol = arcsol
      endif else begin
         arcsol = initsol
         xl = 0.03
         xr = 0.99
         !p.multi = 0
         !p.position = [xl, 0.1, xr, 0.5]
      endelse

      plot, pix, farcsub, /xstyle, /ylog, yrange=[10,max(farc)*1.02], /ystyle
      for l = 0, n_elements(brightlinex)-1 do begin
        oplot, [brightlinex[l],brightlinex[l]], [10,1000], color=5, thick=3
      endfor
      for l = 0, n_elements(linex)-1 do begin
        oplot, [linex[l],linex[l]], [10,50], color=3
      endfor

      wav = poly(pix-x0,arcsol)
      cwav = arcsol[0]
      disp = arcsol[1]

      ; temporary kludge to prevent awkward double-valued or negative solution in far-UV
      if min(wav) lt 2950 then begin
        dwav3000 = 0
        ; Freeze dispersion at 3000 and extrapolate linearly
        for i = n_elements(wav)-1, 0, -1 do begin
           if wav[i] lt 3000 and dwav3000 eq 0 then begin
             dwav3000 = wav[i]-wav[i-1]
           endif
           if dwav3000 gt 0 then begin
             wav[i] = wav[i+1] - dwav3000
           endif
        endfor
      endif



      plotarccheck, pix, wav, arcsol, farcsub, brightrefwav, refwav, nx=nx

      if n_elements(checkplotfile) gt 0 and keyword_set(leaveplotopen) eq 0 then begin
        psclose
      endif

      return, arcsol

end

