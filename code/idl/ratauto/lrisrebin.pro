; Rebin an LRIS image.  (could generalize...?)

; Use this if you took calibrations in 1x1 binning mode but switched to 2x2 (or 2x1, 1x2) mode later.

pro lrisrebin, imagename, newbinning, outimagename=outimagename, fluxsum=fluxsum

  
  if strpos(newbinning,',') gt 0 then splitchar = ','
  if strpos(newbinning,'x') gt 0 then splitchar = 'x'
  newbin = fix(strsplit(newbinning,splitchar,/extract))

  if n_elements(outimagename) eq 0 then outimagename = repstr(imagename,'.fits',repstr(newbinning,',','x')+'.fits')

  data = mrdfits(imagename,0,h, /silent)
  ; shift(data, 0, 1) ; Is this needed?  seems to be in the blue (2 pix, actually)

  oldbinning = clip(sxpar(h,'BINNING'))
  oldbin = fix(strsplit(oldbinning,',',/extract))

  dim = (size(data))[1:*]

  binfactor = newbin*1.0/oldbin

  extrapix0 = dim[0] mod fix(dim[0]/binfactor[0])
  if extrapix0 gt 0 then data = data[0:dim[0]-1-extrapix0,*] ; perhaps better to expand?
  extrapix1 = dim[1] mod fix(dim[1]/binfactor[1] )
  if extrapix1 gt 0 then data = data[*,0:dim[1]-1-extrapix1]

  
  newdata = rebin(data, fix(dim[0]/binfactor[0]), fix(dim[1]/binfactor[1]))

  if keyword_set(fluxsum) then newdata *= binfactor[0]*binfactor[1]
  if keyword_set(fluxsum) then sxaddpar, h,'SATURATE',sxpar('SATURATE')*binfactor[1]*binfactor[2]

  newdim = (size(newdata)) [1:*]

  print,    imagename + ' ('+clip(   dim[0])+'x'+clip(   dim[1])+ ' @ '+repstr(oldbinning,',','x')+') --> ', $
         outimagename + ' ('+clip(newdim[0])+'x'+clip(newdim[1])+ ' @ '+repstr(newbinning,',','x')+')'

  sxaddpar, h, 'BINNING', repstr(newbinning,'x',',')
  sxaddpar, h, 'LTV1', sxpar(h,'LTV1')/binfactor[0]
  sxaddpar, h, 'LTV2', sxpar(h,'LTV2')/binfactor[1]
  sxaddpar, h, 'LTM1_1', sxpar(h,'LTM1_1')/binfactor[0]
  sxaddpar, h, 'LTM2_1', sxpar(h,'LTM2_1')/binfactor[1]
  sxaddpar, h, 'LTM1_2', sxpar(h,'LTM1_2')/binfactor[0]
  sxaddpar, h, 'LTM2_2', sxpar(h,'LTM2_2')/binfactor[1]

  sxaddpar, h, 'REBIN', imagename, 'Original file before rebinning'

  ; leaving window alone for now
  ; leaving AMP?? alone for now

  mwrfits, newdata, outimagename, h, /create


end

