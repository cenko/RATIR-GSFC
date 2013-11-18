function poly2d, x, y, p

  return, p[0] + $    
          p[1]*y + $  ; no x term
          p[2]*y^2 + p[3]*x*y + $
          p[4]*y^3 + p[5]*x*y^2 + p[6]*x^2*y + $
          p[7]*y^4 + p[8]*x*y^3 + p[9]*x^2*y^2 + p[10]*x^3*y + $
          p[11]*y^5 + $
          p[12]*y^6 + $
          p[13]*y^7
  
end

function subtractsky, spec, yref=yref, xsamplestep=xsamplestep, ysamplestep=ysamplestep, sky=sky, order=order, debug=debug

  
  ; solves for the shift in the wavelength solution across the chip in the y direction.

  ; right now this is only a stand-in, ultimately I want to work into pzapspec

   inmode =  size(spec, /tn)
   if inmode eq 'STRING' then begin
      filename = spec
      img=mrdfits(filename, 0, header, /silent)
   endif else begin
      img = spec
   endelse

   dims = size(img, /dimens)
   nx = dims[0]
   ny = dims[1]

   if n_elements(yref) eq 0 then begin
      if ny le 100 then yref = [0, ny-1] $
                   else yref = [ny/2-50, ny/2+50]
   endif
   yrefmed = (yref[0] + yref[1]) / 2.
   if debug then print, 'med', yrefmed

   ;refspec = median(img[*,yref[0]:yref[1]], dimension=2)
   refspec = fltarr(nx)
   for x = 0, nx-1 do begin
        refspec[x] = sigclipmedian(img[x,yref[0]:yref[1]])
   endfor


   if n_elements(ysamplestep) eq 0 then ysamplestep = 8
   if n_elements(xsamplestep) eq 0 then xsamplestep = 512
   nycheck = ny/ysamplestep
   nxblock = nx/xsamplestep
   offy   = ysamplestep*indgen(nycheck)
   gridx = fltarr(nxblock*nycheck)
   gridy = fltarr(nxblock*nycheck)
   gridf = fltarr(nxblock*nycheck)
   d = 0
   for xb = 0, nxblock-1 do begin
      x0 = xb * xsamplestep
      xm =   x0 + xsamplestep/2.
      x1 =   x0 + xsamplestep - 1
      for i = 0, nycheck-1 do begin
         shift = maxcrosscorrelate(refspec[x0:x1], img[x0:x1,offy[i]], -5, 5, 0.125)
         gridx[d] = xm
         gridy[d] = offy[i]-yrefmed
         gridf[d] = shift
         ;print, gridx[d], gridy[d], gridf[d]
         d += 1
      endfor
   endfor

  ;plot, gridy[where(gridx eq 2304.)], gridf[where(gridx eq 2304)]

  npar = 14
  parinfo = replicate({value:0.0, fixed:0, limited:[0,0], tied:'', limits:[0.0,0.0]}, npar)
  parinfo[[0,6,8,9,19]].fixed = 1

  if n_elements(order) gt 0 then begin
    if order eq 1 then parinfo[2:13].fixed = 1
    if order eq 2 then parinfo[4:13].fixed = 1
    if order eq 3 then parinfo[7:13].fixed = 1
  endif

  for i = 0,0 do begin
    ;if i eq 0 then parinfo[6].fixed = 0
    ;if i eq 1 then parinfo[8].fixed = 0
    ;if i eq 2 then parinfo[9].fixed = 0
    ;if i eq 3 then parinfo[10].fixed = 0

    par = mpfit2dfun('poly2d',gridx,gridy,gridf,1+fltarr(d), parinfo=parinfo, bestnorm=chisq, dof=dof, yfit=yfit, /quiet)
    ;writefits, 'shiftresidual'+clip(i)+'.fits', reform(yfit-gridf,nycheck,nxblock)
  

    yy = findgen(ny)
    ;print, min(yy-yrefmed), max(yy-yrefmed)
    ;oplot, yy-yrefmed, poly2d(fltarr(ny)+2304., yy-yrefmed, par)


    offsetmap = findgen(nx,ny)
    sky = fltarr(nx, ny)
    for y = 0, ny-1 do begin
      shift = poly2d(findgen(nx), y-yrefmed, par)
      offsetmap[*,y] = shift
      sky[*,y] =  interpol(refspec, findgen(nx), findgen(nx)-shift)
      ;if y mod 100 eq 0 then print, y
      if debug then print, y, shift[2460]
    endfor
 
    if debug then writefits, 'in.fits', spec, header
    if debug then writefits, 'sky'+clip(i)+'.fits', sky, header
    if debug then writefits, 'residual'+clip(i)+'.fits', img-sky, header
  endfor

  if inmode eq 'STRING' then begin
     writefits, 's'+spec, img-sky, header
     return, 0
  endif else begin
     return, img-sky
  endelse

end

