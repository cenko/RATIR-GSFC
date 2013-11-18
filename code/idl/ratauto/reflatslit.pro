pro reflatslit, filename, slitflatname, dy, outname=outname

  ; did the slit move?  Use this to divide out the incorrect slit response in the flatfield and then
  ; re-add it.

  im = readfits(filename, h)
  slitflat = readfits(slitflatname, hflat)

  s = size(im)  
  nx = s[1]
  ny = s[2]
  x = findgen(nx)
  y = findgen(ny)

  dy0 = -sxpar(h,'LTV2')+sxpar(hflat,'LTV2')
  slitflat = slitflat[*,dy0:*] ; assumes flat is larger

  newx = x
  newy = y - dy

  shiftflat = interpolate(slitflat, newx, newy, /grid)

  im *= slitflat
  im /= shiftflat

  if n_elements(outname) eq 0 then outname = repstr(filename,'fp','cfp')
  writefits, outname, im, h

end

