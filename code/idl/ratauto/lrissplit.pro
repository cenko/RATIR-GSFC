pro lrissplit, file, rightonly=rightonly, leftonly=leftonly

   arr = mrdfits(file,0,h, /silent)

   get_date, now
   hist = 'Processed by lrissplit '+now
   sxaddhist, hist, h

   rfile = repstr(file,'.fits','r.fits')
   lfile = repstr(file,'.fits','l.fits')

   outstr =  removepath(file) + ' -> '
   if  keyword_set(leftonly) eq 0 then outstr += removepath(rfile)
   if  keyword_set(rightonly) eq 0 and keyword_set(leftonly) eq 0 then outstr += ','
   if  keyword_set(rightonly) eq 0  then outstr += removepath(lfile)
   print, outstr

   slit = clip(sxpar(h,'SLITNAME'))
   if strpos(slit,'_') ge 0 then mode = 'sp'
   if strlowcase(slit) eq 'direct' then mode = 'im'
   if n_elements(mode) eq 0 then begin
      print, '  Cannot identify mode.  Assuming imaging.'
      mode = 'im'
   endif

   inst = clip(sxpar(h,'INSTRUME'))
   if inst eq 'LRIS' then camchar = 'r'
   if inst eq 'LRISBLUE' then camchar = 'b'

   lrisversion = 0
   mjd = float(sxpar(h,'MJD-OBS'))
   if mjd gt 0 and mjd lt 54952 then lrisversion = 1
   if mjd ge 54952 and mjd lt 55562 then lrisversion = 2
   if mjd ge 55562 then lrisversion = 3
   if lrisversion eq 0 then begin
     print, '  Cannot identify LRIS version - assuming 3'
     lrisversion = 3
   endif

   lh = h
   rh = h
   s = size(arr)
   if mode eq 'im' then nx = s[1]
   if mode eq 'sp' then nx = s[2]


   ltv1 = sxpar(h,'LTV1') ; always negative
   ltv2 = sxpar(h,'LTV2')
   if mode eq 'im' then ltvx = ltv1
   if mode eq 'sp' then ltvx = ltv2

   ;ampl1 = fix(strsplit(sxpar(h,'AMPL1'),',',/extract))
   ampl2 = fix(strsplit(sxpar(h,'AMPL2'),',',/extract))
   ampr1 = fix(strsplit(sxpar(h,'AMPR1'),',',/extract))
   ;ampr2 = fix(strsplit(sxpar(h,'AMPR2'),',',/extract))

   cbuf = [0,0]
   if camchar eq 'b' then imx = [360, 3790]
   if camchar eq 'r' and lrisversion eq 2 then imx = [430, 3778]   ; this is NOT CONFIRMED
   if camchar eq 'r' and lrisversion eq 2 then cbuf = [44, 12]     ; this is NOT CONFIRMED
   if camchar eq 'r' and lrisversion eq 3 then imx = [430, 3778]
   if camchar eq 'r' and lrisversion eq 3 then cbuf = [44, 12]    ; buffer around gap

   lcrop = [imx[0]+ltvx > 0 , ampl2[1]-cbuf[0]]
   rcrop = [ampr1[0]+cbuf[1]  , imx[1] < nx-1 ]

   if mode eq 'im' then begin
     if keyword_set(rightonly) eq 0 then arrl = arr[lcrop[0]:lcrop[1],*]
     if keyword_set(leftonly) eq 0 then  arrr = arr[rcrop[0]:rcrop[1],*]

     crpix1 = sxpar(h,'CRPIX1')
     sxaddpar, lh, 'CRPIX1', crpix1-lcrop[0]
     sxaddpar, rh, 'CRPIX1', crpix1-rcrop[0]
     sxaddpar, lh, 'LTV1', ltv1-lcrop[0]
     sxaddpar, rh, 'LTV1', ltv1-rcrop[0]
   endif
   if mode eq 'sp' then begin

     if keyword_set(rightonly) eq 0 then arrl = arr[*,lcrop[0]:lcrop[1]]
     if keyword_set(leftonly) eq 0 then  arrr = arr[*,rcrop[0]:rcrop[1]]
     sxaddpar, lh, 'LTV2', ltv2-lcrop[0]
     sxaddpar, rh, 'LTV2', ltv2-rcrop[0]
   endif

   sxdelpar, h, ['AMP1','AMP2','AMP3','AMP4','AMPL1','AMPL2','AMPR1','AMPR2']
   sxaddpar, lh, 'CHIP', 'left'
   sxaddpar, rh, 'CHIP', 'right'

   dotpos = strpos(file,'.',/reverse_search)
   fileroot = strmid(file,0,dotpos)
   if keyword_set(rightonly) eq 0 then mwrfits, arrl, fileroot+'l.fits', lh, /silent, /create
   if keyword_set(leftonly) eq 0 then  mwrfits, arrr, fileroot+'r.fits', rh, /silent, /create

end

