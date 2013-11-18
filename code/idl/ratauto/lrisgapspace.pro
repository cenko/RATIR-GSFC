pro lrisgapspace, filename, mode=mode

  ; artificially insert the gap into an un-split LRIS spectrum.

  ; LRIS R-3: gap appears to be 29.9" (222 pix) - wrong
  ; in imaging it is clearly 210 pix (28.35")


  im = readfits(filename,h)
  
  ampl2 = fix(strsplit(sxpar(h,'AMPL2'),',',/extract))
  ampr1 = fix(strsplit(sxpar(h,'AMPR1'),',',/extract))
  
  spacepix = 210 / sxpar(h,'YBIN')
  s = size(im)
  nx = s[1]
  ny = s[2]

  if mode eq 'sp' then begin
     nny = ny + spacepix
     newim = fltarr(nx, nny)

     newim[*,0:ampl2[1]] = im[*,0:ampl2[1]]
     newim[*,ampr1[0]+spacepix:nny-1]  = im[*,ampr1[0]:*]
  endif
  if mode eq 'im' then begin
     nnx = nx + spacepix
     newim = fltarr(nnx, ny)

     newim[0:ampl2[1],*] = im[0:ampl2[1],*]
     newim[ampr1[0]+spacepix:nnx-1,*]  = im[ampr1[0]:*,*]
  endif


  newname = repstr(filename,'.fits','space.fits')

  writefits, newname, newim, h
 

  ; LRIS-B is tilted, don't worry about it now


end

