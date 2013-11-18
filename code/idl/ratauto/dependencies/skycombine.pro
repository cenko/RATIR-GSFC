pro skycombine, filelist, outfile, removeobjects=removeobjects, response=response, indresp=indresp, yresp=yresp, lines=lines, objthresh=objthresh, objbuffer=objbuffer, verbose=verbose, mediani=mediani, meani=meani, trimlo=trimlo, trimhi=trimhi, mincounts=mincounts, maxcounts=maxcounts, type=type, normzone=normzone, linewidth=linewidth, lineblocks=lineblocks, satlevel=satlevel, satradius=satradius, xnormkey=xnormkey, slitnorm=slitnorm

; combine a bunch of flats with optional response division

algorithm = 'median' ; default
if keyword_set(mediani) then algorithm = 'median'
if keyword_set(meani)   then algorithm = 'mean'  ; the 'i' is because otherwise it overwrites mean()?
if n_elements(trimlo) eq 0 then trimlo = (n_elements(filelist)+1)/4
if n_elements(trimhi) eq 0 then trimhi = trimlo
if n_elements(mincounts) eq 0 then mincounts = 1
if n_elements(maxcounts) eq 0 then maxcounts = 55000
if n_elements(satradius) eq 0 then satradius = 50
if n_elements(satlevel) eq 0 then satlevel=30000.

if n_elements(outfile) eq 0 then begin
  filearr = strsplit(filelist,'.', /extract)
  fileroot = filearr[0]
  outfile = fileroot + '.fits'
endif

if n_elements(filelist) eq 1 then begin
  files = grabcolumn(filelist,0,/str)
  files = files[where(files ne '')]
endif else begin
  files = filelist
endelse

nfiles = n_elements(files)
nmid = nfiles/2

refdata = mrdfits(files[nmid], 0, h, /silent, /fscale)
s = size(refdata)
nx = s[1]
ny = s[2]

data = fltarr(nx, ny, nfiles) + !values.f_nan
inmeds = fltarr(nfiles)
usefiles = strarr(nfiles)
skymed=fltarr(nfiles)

z = 0
for f = 0, nfiles-1 do begin
   indata = mrdfits(files[f], 0, h, /silent, /fscale)
   ins = size(indata)
   inx = ins[1]
   iny = ins[2]
   if inx ne nx or ny ne ny then begin
      print, 'File ', files[f], ' has wrong dimensions ('+clip(inx)+'x'+clip(iny)+'; should have '+clip(nx)+'x'+clip(ny)+')'
      continue
   endif

   data[*,*,z] = float(indata)
   usefiles[z] = files[f]

   ;indata = indata + sxpar(h,'BZERO')
   if n_elements(xnormkey) gt 0 then begin
     xnormstr = sxpar(h,xnormkey)
     xnormzone = fix(strsplit(xnormstr,',',/extract))
     normzone = [xnormzone[0], 0, xnormzone[1], -1]
   endif
   if n_elements(normzone) eq 4 then begin
      if normzone[2] lt 0 then normzone[2] = nx+normzone[2]
      if normzone[3] lt 0 then normzone[3] = ny+normzone[3]
      inmed = sigclipmedian(indata[normzone[0]:normzone[2],normzone[1]:normzone[3]])
   endif else begin
      inmed = sigclipmedian(indata)
   endelse
   inmeds[f] = inmed
   if inmed ge mincounts and inmed le maxcounts then begin
      if keyword_set(indresp) then begin
         iresp = median(indata) ; if each flat has a changing sky background
         for x = 0, nx-1 do indata[x,*] = inmed * indata[x,*] / iresp[x] ; inmed will immediately get normalized out, which is what we want
      endif
      print, '  ' + removepath(files[f]) + ' (' + strtrim(long(inmed),2) + ' counts/pix)'
      skymed[z] = inmed
      data[*,*,z] = float(indata); / inmed
      usefiles[z] = files[f]
      z = z + 1
   endif else begin
      if inmed lt mincounts then $
      print, '  ' + removepath(files[f]) + ' (' + strtrim(long(inmed),2) + ' counts/pix) - too few counts; excluding'
      if inmed gt maxcounts then $
      print, '  ' + removepath(files[f]) + ' (' + strtrim(long(inmed),2) + ' counts/pix) - too many counts; excluding'
   endelse
;   z = z+1
endfor

medsky=median(skymed)
for f = 0, z-1 do begin
   inmed = median(data[*,*,f]) ; if each flat has a changing sky background
   factor=medsky/inmed
;print,factor
   data[*,*,f]=data[*,*,f]*factor(0)
endfor

if z lt 2 then begin
    print, 'ERROR - Not enough counts to make a flat with these data!'
    return
endif

if z ne nfiles then data = data[*,*,0:z-1]

;flat = fltarr(nx, ny)
;;medmax=0
;;sigmax=0
;if algorithm eq 'median' then begin
;   print, '  Median-combining...'
;   i = 0l
;   for y = 0, ny-1 do begin    
;     for x = 0, nx-1 do begin  ; x is the fast array and should be iterated inside
;       i = i + 1
;;       flat[x, y] = median(data[x, y, *], /even) ; mediansigclip?
;       djs_iterstat, data[x,y,*], sigrej=3, mean=m, median=me, sigma=s
;       flat[x,y] = me
;;       IF me GT medmax THEN medmax = me
;;       IF s GT sigmax THEN sigmax = s
;;print,me,' ',s
;;print,medmax,' ',sigmax
;       if i mod 157557 eq 0 then begin
;;          st = stddev(data[x, y, *])
;          val = flat[x, y]
;          if keyword_set(verbose) then print, clip(x,5), clip(y,5), 'stdev=',fpr(st,1.4), '  val=',fpr(val,1.4);, '  stde;;/sqrt(val)=',fpr(st/sqrt(val),3.2) 
;       endif
;     endfor
;   endfor
;endif
;;print,medmax,sigmax
;if algorithm eq 'mean' then begin
;   print,  '  Combining via trimmed mean...'
;   i = 0l
;   nslice = (size(data)) [3]
;   for y = 0, ny-1 do begin    
;     for x = 0, nx-1 do begin  ; x is the fast array and should be iterated inside
;       slice = data[x,y,*]
;       slice = slice[sort(slice)]
;       slice = slice[trimlo:nslice-1-trimhi]
;       flat[x,y] = mean(slice)
;     endfor
;   endfor
;endif
;device,decomposed=0
;imSize=SIZE(flat)
;flat=congrid(flat,imSize[1]/2.,imSize[2]/2.)
;imSize=SIZE(flat)
;WINDOW,0,XSIZE=imSize[1],YSIZE=imSize[2]
;scale=bytscl(flat,min=330,max=340,top=255)
;tv,scale

;stop
;if keyword_set(indresp) then begin
;   ; Remove this again so the more refined response method can be applied later
;   ; for x = 0, nx-1 do data[x,*,*] = data[x,*,*] * iresp[x]   ;  data not used again for spec
;   for x = 0, nx-1 do flat[x,*]   = flat[x,*]   * iresp[x]   ; 
;endif
;
if keyword_set(removeobjects) then begin
   print, '  Identifying objects...'

   if n_elements(objthresh) eq 0 then objthresh = 6
   if n_elements(objbuffer) eq 0 then objbuffer = 12
   for f = 0, z-1 do begin
      if keyword_set(verbose) then print, '  ', removepath(files[f])
      indata = data[*,*,f] ;/ flat
      mask = fltarr(nx, ny)
      satmask = fltarr(nx, ny)
      sigclipstats, indata, sigmahi=5, sigmalo=5, median=datamed, stdev=datastdev
      ;datamed = sigclipmedian(indata, sigmahi = 5, sigmalo = 5)
      ;datastdev = sigclipstdev(indata, sigmahi = 5, sigmalo = 5)

      sourcepixels = where(indata ge datamed + objthresh * datastdev, ctsourcepix)
      if ctsourcepix gt 0 then begin
         mask[sourcepixels] = 1
         smoothmask = smooth(mask, 1+2*objbuffer, /edge_truncate)
      endif
      ctsatpix = 0
      if n_elements(satlevel) gt 0 then begin
         satpixels = where(indata ge satlevel, ctsatpix)
         if ctsatpix gt 0 then begin
            satmask[satpixels] = 1
            smoothsatmask = smooth(mask, 1+2*satradius, /edge_truncate)
         endif
      endif
      if ctsourcepix gt 0 then begin
         sourcepixels = where(smoothmask gt 1.0E-6, ct)
         sourcepixels = where(mask EQ 1., ct)
         indata[sourcepixels] = !Values.F_NAN
      endif
      if ctsatpix gt 0 then begin
         nearsatpixels = where(smoothsatmask gt 1.0E-6, ct)
         nearsatpixels = where(satmask EQ 1., ct)
         indata[nearsatpixels] = !Values.F_NAN
      endif


      data[*,*,f] = indata  ; inefficient to copy over unaffected elements, but oh well
   endfor

   reflat = fltarr(nx, ny)

   if algorithm eq 'median' then begin
     print, '  Median-combining...'
     i = 0l
     for y = 0, ny-1 do begin
       for x = 0, nx-1 do begin
         i = i + 1
         vector=data[x,y,*]
         tmp=where(finite(vector))
         djs_iterstat, vector(tmp), sigrej=3, mean=m, median=me, sigma=s
;         djs_iterstat, data[x,y,*], sigrej=3, mean=m, median=me, sigma=s
;         IF (y GT 10 AND finite(me) NE 1) THEN stop
         reflat[x,y] = me
;         reflat[x, y] = median(data[x, y, *], /even)   ;sigclipping this is REALLY slow...
         if i mod 157557 eq 0 and keyword_set(verbose) then begin
            slice = data[x, y, *]
            good = where(finite(slice), ct)
            if ct gt 3 then st = stddev(slice[good]) else st = 0
            val = flat[x,y] * reflat[x, y]
            print, '  ', clip(x,5), clip(y,5), 'n=', clip(ct,2), ' stdev=',fpr(st,1.4), '  val=',fpr(val,1.4);, '  stdev/sqrt(val)=',fpr(st/sqrt(val),3.2) 
         endif
       endfor
     endfor
     bad = where(finite(reflat) eq 0, ct)
     if ct gt 0 then reflat[bad] = 1
   endif

   if algorithm eq 'mean' then begin
     print,  '  Combining via trimmed mean...'
     i = 0l
     for y = 0, ny-1 do begin
       for x = 0, nx-1 do begin
         slice = data[x,y,*]
         good = where(finite(slice) eq 1, ctgood)
         if ctgood eq 0 then begin
             reflat[x,y] = 1
         endif else begin
             slice = slice[good]
             itrimlo = trimlo
             itrimhi = trimhi
             while ctgood-itrimlo-itrimhi lt 1 do begin
               itrimlo = (itrimlo - 1) > 0
               itrimhi = (itrimhi - 1) > 0
             endwhile
             slice = slice[sort(slice)]
             slice = slice[itrimlo:ctgood-1-itrimhi]
             reflat[x,y] = mean(slice)
         endelse
       endfor
     endfor
   endif

   flat = reflat
;   flat = flat * reflat    ; want to divide by flat, then divide by reflat, so divide by (flat * reflat)
;device,decomposed=0
;imSize=SIZE(flat)
;flat=congrid(flat,imSize[1]/2.,imSize[2]/2.)
;imSize=SIZE(flat)
;WINDOW,0,XSIZE=imSize[1],YSIZE=imSize[2]
;scale=bytscl(flat,min=330,max=340,top=255)
;tv,scale
;stop
endif
;
;
;
;if keyword_set(response) then begin
;  if n_elements(linewidth) eq 0 then linewidth = 12
;  if n_elements(lineblocks) eq 0 then lineblocks = 5
;
;  if n_elements(yresp) eq 0 then begin
;    yresp = indgen(ny)
;  endif
;  if n_elements(lines) eq 55 then begin ; this is not used
;    if lines eq 'find' then begin
;        lines = [-1]
;        nyr = n_elements(yresp)
;        lineyresp = yresp[(nyr/2-20) > 0 : (nyr/2+20) < (nyr-1)]
;        medcolx = fltarr(nx)
;        for x = 0, nx-1 do begin
;           medcolx[x] = median(flat[x,lineyresp], /even)
;        endfor
;        for x = 50, nx-50-1, 100 do begin
;           xb     = (indgen(nx)) [(x-75)>0:(x+75)<(nx-1)]
;           par = poly_fit(xb, medcolx[xb], 3, yfit=yfit)
;           ydiff = medcolx[xb]-yfit
;           posslines = where(ydiff gt 3.5*stddev(ydiff), ct) ; positive deviation only
;           if ct gt 0 then print, xb[posslines]
;           if ct gt 0 then lines = [lines,xb[posslines]]
;        endfor
;    endif
;  endif
;
;  medcolx = fltarr(nx)
;  for x = 0, nx-1 do begin
;     medcolx[x] = median(flat[x,yresp], /even)
;  endfor
;  ;plot, medcolx
;  
;  responsex = medcolx
;  for x = 0, nx-1 do begin
;     xi = (findgen(nx)) [(x-10)>0:(x+10)<(nx-1)]
;     par = polyiterfit(xi, medcolx[xi], order=1, yfit=yfit)
;     responsex[x] = yfit[where(xi eq x)]
;  endfor
;  for x = 0, nx-1 do begin
;     flat[x,*] = flat[x,*] / responsex[x]
;  endfor
;
;  ; mwrfits, flat, 'flatresp1.fits', /create
;
;  ystart = fix(ny*(findgen(lineblocks)/lineblocks))
;  ystop =  fix(ny*((1+findgen(lineblocks))/lineblocks))-1
;  for l = 0, n_elements(lines)-1 do begin
;     linestart = lines[l] - linewidth/2
;     linestop = lines[l] + linewidth/2
;     for b = 0, lineblocks-1 do begin
;        medleft = median(flat[linestart-9:linestart-3,ystart[b]:ystop[b]], /even)
;        medright = median(flat[linestop+3:linestop+9,ystart[b]:ystop[b]], /even)
;        medtarget = (medleft + medright)/2.
;        for x = linestart, linestop do begin
;           medzone = median(flat[x,ystart[b]:ystop[b]], /even)
;           flat[x,ystart[b]:ystop[b]] = medtarget * flat[x,ystart[b]:ystop[b]] / medzone
;        endfor
;     endfor
;     medleft = median(flat[linestart-9:linestart-3,*], /even)
;     medright = median(flat[linestop+3:linestop+9,*], /even)
;     medtarget = (medleft + medright)/2.
;     for x = linestart, linestop do begin
;        flat[x,*] = medtarget * flat[x,*] / median(flat[x,*])
;     endfor
;  endfor
;
;  ; don't even bother to flatfield for counts < ~1000*gain: even with 20 flatfields
;  ; error introduced is 1% which exceeds typical flat variation
;
;  reqcounts = 1000 < (max(medcolx)/10.)
;
;  lowcounts = where(medcolx lt reqcounts, ct, complement=goodcounts) ; formerly medcolx*median(inmeds)
;  maxx = (where(medcolx eq max(medcolx))) [0]
;  if ct gt 0 then begin
;    if min(lowcounts) lt maxx then $
;       print, '  Low counts for x=', $
;           clip(min(lowcounts[where(lowcounts lt maxx)])),'-',clip(max(lowcounts[where(lowcounts lt maxx)]))
;    if max(lowcounts) gt maxx then $
;       print, '  Low counts for x=', $
;           clip(min(lowcounts[where(lowcounts gt maxx)])),'-',clip(max(lowcounts[where(lowcounts gt maxx)]))
;    print, '  Using only a slit-averaged flat in this region.'
;    for y = 0, ny-1 do begin
;      flat[lowcounts,y] = median(flat[goodcounts[0:100<(n_elements(goodcounts)-1)],y])   
;        ; this does assume directionality of the bad region...
;    endfor
;  endif
;
;endif
;
;if keyword_set(slitnorm) then begin  ; halogen flats have an artificial slit response superimposed; remove it here.
;  if n_elements(goodcounts) eq 0 then goodcounts = indgen(nx)
;  for y= 0, ny-1 do begin
;     medrowy = median(flat[goodcounts,y], /even)
;     if medrowy gt 0.25 then flat[*,y] /= medrowy ; less than this is assumed to not be a false slit variation but genuinely 
;                                                  ; dimmed part of the chip (off slit), which should be kept for later cropping
;  endfor
;  
;endif


;badval = where(flat le 0.01, ctbad)  ; zero values will create problems (nans) later and negative values are nonsense.
                                     ; also flagging very low values in advance.
;if ctbad gt 0 then begin
;  flat[badval] = 0.01  ; Prevents noise from being amplified too ridiculously if this flat is actually used.  Could set to Nan...
;endif

;sxaddpar, h, 'BZERO', 0
for f = 0, z-1 do begin
  sxaddpar, h, 'SKY'+strtrim(string(f),2), removepath(usefiles[f])
endfor

if keyword_set(type) then sxaddpar,  h, 'SKYTYPE', type

get_date, now
hist = 'Processed by flatcombine '+now
sxaddhist, hist, h
;sxaddpar, h, 'BZERO', 0
;sxaddpar, h, 'BSCALE', 1
mwrfits, flat, outfile, h, /create

print, '  Written to ', outfile

end
