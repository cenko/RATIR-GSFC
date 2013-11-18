function findbrightest, linex, linestrength, nfind, xrange=xrange

   ; find the "brightest" nfind lines (as ranked by isolation from next-brightest line)

   if nfind ge n_elements(linex) then return, linex

   if n_elements(xrange) gt 0 then xmin = xrange[0] else xmin = min(linex)-0.1
   if n_elements(xrange) gt 1 then xmax = xrange[1] else xmax = max(linex)+0.1
   minstrength = 0.0

   w = where(linex ge xmin and linex le xmax and linestrength ge minstrength)
   x = linex[w]   
   y = linestrength[w]
 
   n = n_elements(x)
   prom = fltarr(n)
   method = 'sum'
   if method eq 'abs' then begin
      for i = 0, n-1 do begin
         dist = abs(x-x[i])
         brighter = where(y gt y[i], ct)
         if ct eq 0 then begin
            prom[i] = !values.f_infinity
         endif else begin
            prom[i] = min(dist[brighter])
         endelse
      endfor
   endif 
   if method eq 'sum' then begin
      for i = 0, n-1 do begin
         dist = abs(x-x[i])
         brighterleft = where(y gt y[i] and x lt x[i], ctleft)
         brighterright = where(y gt y[i] and x gt x[i], ctright)
         if ctleft  gt 0 then lprom = min(dist[brighterleft])  else lprom = x[i]-min(x)
         if ctright gt 0 then rprom = min(dist[brighterright]) else rprom = max(x)-x[i]
         prom[i] = lprom + rprom
         ;if n lt 200 then print, x[i], y[i], ' | ', ctleft, ctright, ' | ',  lprom, rprom, '|', prom[i]
      endfor
   endif

   sprom = reverse(sort(prom))

   xout = x[sprom[0:nfind-1]]
   yout = y[sprom[0:nfind-1]]

   s = sort(xout)
   xout = xout[s]
   yout = yout[s]

;psclose
;window, 0
;!p.multi = 0
;!p.position = 0
;plot, x, y, psym=1, yrange=[1.0,2.1]*1.02e6, /ystyle
;oplot, xout, yout, psym=7
;stop

   return, xout

end

