function addspec, spec1, spec2, w1=w1, w2=w2

  if n_elements(w1) eq 0 then w1 = 1.
  if n_elements(w2) eq 0 then w2 = 1.
  wtot = w1+w2

  n1 = n_elements(spec1)
  n2 = n_elements(spec2)
  nmin = n1 < n2

  disp = abs(median(spec1.wav - shift(spec1.wav,1)))
  wavshift = median(spec2[0:nmin-1].wav-spec1[0:nmin-1].wav,/even)
  dpix = round(wavshift / disp)

  ord1 = indgen(n1)
  ord2 = indgen(n2)
  nnew = nmin - dpix
  if dpix gt 0 then begin
     ord1 = ord1[dpix:*]
  endif
  if dpix lt 0 then begin
     ord2 = ord2[-dpix:*]
  endif
  n = n_elements(ord1) < n_elements(ord2)
  ord1 = ord1[0:n-1]
  ord2 = ord2[0:n-1]
  sumspec = spec1[ord1]

  sumspec.wav = (w1*spec1[ord1].wav + w2*spec2[ord2].wav) / wtot
  sumspec.flux= (w1*spec1[ord1].flux+ w2*spec2[ord2].flux) / wtot
  if tag_exist(spec1,'sky') then $
  sumspec.sky = (w1*spec1[ord1].sky + w2*spec2[ord2].sky) / wtot
  if tag_exist(spec1,'unc') and tag_exist(spec2,'unc') then $
  sumspec.unc = sqrt( (w1*spec1[ord1].unc^2 + w2*spec2[ord2].unc^2)/wtot )
  if tag_exist(spec1,'pix') and tag_exist(spec2,'pix') then $
  sumspec.pix = (w1*spec1[ord1].pix + w2*spec2[ord2].pix) / wtot
  if tag_exist(spec1,'norm') and tag_exist(spec2,'norm') then $
  sumspec.norm = (w1*spec1[ord1].norm + w2*spec2[ord2].norm) / wtot

  return, sumspec

end


pro testaddspec
  s1 = replicate({wav:0.,flux:0.,sky:0.,unc:0.,pix:0L,norm:0.},10)
  s2 = replicate({wav:0.,flux:0.,sky:0.,unc:0.,pix:0L,norm:0.},12)

  s1.wav = 100+findgen(10)
  s2.wav = 100+findgen(12)-4.1

  s1.flux = 10.
  s2.flux = 11.

  print, s1.wav
  print, s2.wav

  add = addspec(s1, s2)

  print, add.wav

end

