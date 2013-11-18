;function doublegaussian, x, p
;    return, gaussian(x,p[0:3]) + gaussian(x,p[4:7])
;end

function findlines, f, fbaksub, linestrength=linestrength, centroidrad=centroidrad, threshold=threshold

   ; Input a 1D spectrum (usually an arc or sky); returns a list of line locations (in pixels)
   ; also background-subtracts the spectrum which may or may not be desirable...
   ; strengths are returned in linestrength

   ; could be improved - absorption lines (solar calcium, usually) cause fake detections on the edges

      if n_elements(threshold) eq 0 then threshold = 7.0

      farc = f
      nx = n_elements(farc)

      forig = farc ; farc IS altered permanently

      if n_elements(centroidrad) eq 0 then cr = 5 else cr = centroidrad
      ; do a really rough background continuum subtraction
      for iter = 0, 3 do begin
         bakkernel =  [0.5,1,1,1,1,1,0.5]
         bakkernel = bakkernel / total(bakkernel)
         fbakconv = convol(farc, bakkernel, /edge_truncate)
         for x = 0, nx-1 do begin
            farc[x] = farc[x] - min(fbakconv[x-15>1:x+15<nx-1])
         endfor
      endfor
      
      farcx = findgen(n_elements(farc))
      arckernel = [0.2,0.4,0.9,1,0.9,0.4,0.2]
      arckernel = arckernel / total(arckernel)
      arcconv = convol(farc, arckernel, /edge_truncate)
      ; identify lines

      diff = abs(farc - shift(farc,-1))
      readnoise = sigclipmedian(diff)
      noise = sqrt((forig>0) + readnoise^2)

      arcthreshold = noise*threshold
      ;print, 'arc threshold ', arcthreshold

      linex = fltarr(1000) ; up to 1000 lines...
      linestrength = fltarr(1000)
      l = 0
      for x = cr, nx-cr-1 do begin
         if arcconv[x] lt arcthreshold[x] then continue
         if arcconv[x] gt arcconv[x-1] and arcconv[x] gt arcconv[x+1] then begin

            linex[l] = total(farc[x-cr:x+cr] * farcx[x-cr:x+cr])/total(farc[x-cr:x+cr])
            linestrength[l] = arcconv[x]

            ;print, x, linex[l], arcconv[x]

            improvedlinex = getlinecenter(forig, linex[l], fitrad=7)  ; seems to produce better results,
            ;print, linex[l], improvedlinex                          ; but monitor this in future

            if abs(linex[l]-improvedlinex) lt 0.7 then linex[l] = improvedlinex
                ; if larger, is probably a blend or something and gaussian fit not reliable
                ; we should actually do a new iteration to look for blends and refit those as doubles

            l = l + 1

         endif
      endfor
      if l ge 1 then linex = linex[0:l-1]   else linex = -1
      if l ge 1 then linestrength = linestrength[0:l-1]   else linestrength = -1

      farc = farc > 1e-5 ; for plotting purposes...

  fbaksub = farc
  return, linex

end




