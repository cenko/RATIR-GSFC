function nndist, x, y
   allmindist = fltarr(n_elements(x))
   for i = 0L, n_elements(x)-1 do begin
      mindist = 1e9
      for j = 0, n_elements(x)-1 do begin
         if i eq j then continue
         dist = sqrt((x[i]-x[j])^2+(y[i]-y[j])^2)
         if dist lt mindist then mindist = dist
      endfor
      allmindist[i] = mindist
   endfor
   return, allmindist
end

function findneighbors, x, y, i, mindist=mindist, count=count

      neighbors = [-1]
      for j = 0L, n_elements(x)-1 do begin
         if i eq j then continue
         dist = sqrt((x[i]-x[j])^2+(y[i]-y[j])^2)
         if dist lt mindist then neighbors = [neighbors,j]
      endfor
      count = n_elements(neighbors)-1
      if count gt 0 then neighbors = neighbors[1:*]
      return, neighbors
end

function nearest, x, y, xarr, yarr, mindist=mindist, count=count

      count = 0
      index = -1
      imindist = mindist

      for j = 0L, n_elements(xarr)-1 do begin
         dist = sqrt((x-xarr[j])^2+(y-yarr[j])^2)
         if dist lt imindist then begin
             imindist = dist
             index = j
             count = count + 1
         endif
     endfor

     return, index
end

function maxpixel, imagename, x, y, r=r
   data = mrdfits(imagename, /silent)
   if n_elements(r) eq 0 then r = 2
   nx = (size(data)) [1]
   ny = (size(data)) [2]
   m = fltarr(n_elements(x))
   for i = 0, n_elements(x)-1 do begin
     xi = fix(x[i])
     yi = fix(y[i])
     if xi-r ge nx or yi-r ge ny then m[i] = !values.f_nan else $
     m[i] = max(data[(xi-r)>0:(xi+r)<(nx-1),  (yi-r)>0:(yi+r)<(ny-1)])
;stop
   endfor

   return, m
end





pro fieldphot, imagename, ra, dec, mag, magerr, zpt1=zpt1, airtermi=airtermi, rinner=rinner, router=router, seeing=seeing, relative=relative, getseeing=getseeing

; r1 is the optimum aperture.
; r2 is a large aperture that can be treated as infinite.

common lrisconfig, lrisautopath, refdatapath, defaultspath

if file_test('temp.param') eq 0 then spawn, 'cp '+defaultspath+'/temp.param ./'
if file_test('sex.config') eq 0 then spawn, 'cp '+defaultspath+'/sex.config ./'
if file_test('sex.conv') eq 0 then spawn, 'cp '+defaultspath+'/sex.conv ./'

if airtermi lt 0 then airterm = abs(airtermi[0]) else airterm = airtermi[0]

data = mrdfits(imagename,0,h,/silent)
;filter = sxpar(h,'FILTNAM')
filter = sxpar(h,'FILTER')
;airmass = sxpar(h,'AIRMASS')
airmass = 1
exptime = sxpar(h,'EXPTIME')
;saturate = sxpar(h,'SATURATE')
saturate = 100000.

sexpath = getsexpath()

if file_test('temp.cat') ge 1 then spawn, 'rm -f temp.cat'
cmd = sexpath+'sex '+imagename+' -c sex.config'
spawn, cmd

if file_test('temp.cat') eq 0 then begin
   print, 'Sextractor failed on '+ imagename
   print, 'Command was ', cmd
   stop
   return
endif

readcol, 'temp.cat', x, y, ra, dec, sexmag, sexmagerr, ellip, fwhm, /silent
; proper way to do this would be to use a struct.

;megaplot, fwhm, ellip, xrange=[0,20], size=10.0^(-(25+sexmag-16.)/10)

mindist = nndist(x,y)

;for i = 0, n_elements(x)-1 do print, x[i], y[i], ra[i], dec[i], sexmag[i], sexmagerr[i], ellip[i], fwhm[i], mindist[i]

skypix = sxpar(h,'COUNTS')
maxpix = maxpixel(imagename, x,y, r=2)

satval = 48000 < saturate
ct = 0
IF strcmp(filter,'H') THEN abovesky = 6000 ELSE IF strcmp(filter,'J') THEN abovesky = 4500 ELSE IF strcmp(filter,'Y') THEN abovesky = 800 ELSE abovesky = 100
while ct lt 7 do begin
   if abovesky lt 0 then break
   pstars = where(fwhm gt seeing[0] and fwhm lt seeing[1]*airmass and maxpix lt satval and maxpix gt abovesky+skypix, ct)
   if ct lt 7 and abovesky gt 1600 then abovesky = abovesky - 1000
   if ct lt 7 and abovesky le 1600 then abovesky = abovesky - 100
endwhile
if ct eq 0 then begin
   print, "Couldn't identify any calibration stars in the image."
   print, "Check input seeing and try again."
   stop ; crash deliberately here
endif
medstarfwhm = median(fwhm[pstars], /even)
fwhmthresh = medstarfwhm*1.2 < airmass*seeing[1]

getseeing = medstarfwhm
; could do the others too, but no need - watch out for this however

if n_elements(rinner) eq 0 then r1 = medstarfwhm*1.2 else r1 = rinner
if keyword_set(relative) then begin
  r2 = r1
endif else begin
  ;if n_elements(router) eq 0 then r2 = 35.         else    this is now too dangerous.
  r2 = router
endelse

print, 'Using aperture radius ', clip(r1)

rcheck1 = medstarfwhm     
rcheck2 = medstarfwhm*1.5 

; trust sex's centers
aper, data, x, y, apmagarr, apmagerrarr, sky, skyerr, 1.0, [r1,r2,rcheck1,rcheck2], [40,50], [-100,satval], /silent
apmagarr = apmagarr + 2.5*alog10(exptime) ; adjust for exposure time!!!
dmfloor = 0.03 ; some danger here if the fwhm-estimated aperture is way too large or small

dm = (apmagarr[2,*] - apmagarr[3,*]) [*]; flux difference between two apertures
;dmunc = apmagerrarr[3,*]
;print, dm

stars = where(fwhm lt fwhmthresh and dm gt dmfloor)
;IF clip(r1) EQ 6.98400 THEN stop
bstars = where(fwhm lt fwhmthresh and dm gt dmfloor and maxpix ge median(maxpix[stars], /even) )
dmthresh = (median(dm[bstars],/even) * [0.8, 1.2]) + [-0.1, 0.1]
dmthresh[0] = dmthresh[0] > 0.04 ; some danger here if the fwhm-estimated aperture is way too large or small

apmag1 = (apmagarr[0,*]) [*]
apmag2 = (apmagarr[1,*]) [*]
apmagerr1 = (apmagerrarr[0,*]) [*]
apmagerr2 = (apmagerrarr[1,*]) [*]

;for i = 0, n_elements(x)-1 do print, x[i], y[i], ra[i], dec[i], apmag1[i], apmagerr1[i], apmag2[i], apmagerr2[i], mindist[i]

isgood1mag = apmag1 gt -50 and apmag1 lt 50 and dm gt 0
isgood2mag = apmag2 gt -50 and apmag2 lt 50 and dm gt 0

if keyword_set(relative) eq 0 then begin
  brightthresh = 0.02-0.0025
  minsep = r2+2*r1
  ct = 0
  while ct lt 5 do begin ; used to be 3...
     brightthresh = brightthresh + 0.0025
     brightstars = where(apmagerr2 lt brightthresh and mindist gt minsep and fwhm lt fwhmthresh and $
                         dm gt dmthresh[0] and dm lt dmthresh[1] and isgood2mag, ct)
     if brightthresh gt 0.09 then begin ; not willing to accept any more than a 10 percent error, at this point start degrading the search radius
        if minsep gt r2 then begin 
           minsep = minsep - r1
           brightthresh = 0.02
           print, 'WARNING - not enough isolated bright stars ('+clip(ct)+') found; decreasing minsp to ', fpr(minsep,2.2), ' pixels.'
        endif else begin
           if ct ge 1 then begin
              print, 'WARNING - very few isolated bright stars ('+clip(ct)+') found.'
              break
           endif else begin
              print, 'ERROR - found no isolated stars.'
              ;print, '        Could only find '+clip(ct)+' stars separated by > '+fpr(minsep,2.1)+'pixels.'
              ;stop
              return
           endelse
        endelse
     endif
  endwhile

  for i = 0, n_elements(brightstars)-1 do begin
     ;print, x[brightstars[i]], y[brightstars[i]], apmag1[brightstars[i]], apmag2[brightstars[i]], apmag1[brightstars[i]]-apmag2[brightstars[i]]
  endfor
  
  apcorr = median(apmag1[brightstars]-apmag2[brightstars], /even)
  if ct gt 1 then begin
    apunc  = stdev(apmag1[brightstars]-apmag2[brightstars])
    apdev  = median(abs(apmag1[brightstars]-apmag2[brightstars] - apcorr),/even)
  endif else begin
    apunc = 1.0
    apdev = 1.0
  endelse
  napstars = n_elements(brightstars)

  print, 'Aperture correction:',fpr(apcorr,2.3), ' +/-', fpr(apunc,2.3), '(', fpr(apdev,1.3), '):  ', clip(napstars), ' stars'
endif else begin
  apcorr = 0.
  apunc = 0.
  apdev = 0.
  napstars = 0
endelse

apmag1corr = apmag1 - apcorr

zpteff = zpt1[0] - airterm*(airmass-1)

printf, 3, fpr(medstarfwhm,2.3), ' ', fpr(r1,2.3), ' ', fpr(dmthresh[0],2.3), ' ', fpr(dmthresh[1],1.3), ' | ', fpr(zpteff,2.3), ' ', fpr(apcorr,2.3), ' ', fpr(apunc,1.3), ' ', fpr(apdev,1.3), ' ', clip(napstars,3)
;printf, 3, clip('file',25), 'r_ap', 'fwhm', 'dmmin', 'dmmax', 'apcor', 'apunc', 'apdev', 'nap', 'zpt'


calthresh = 0.1
calstars = where(apmagerr1 lt calthresh and mindist gt 3*r1 and fwhm lt fwhmthresh and dm gt dmthresh[0] and dm lt dmthresh[1] and isgood2mag, ct)
ra = ra[calstars]
dec = dec[calstars]
mag = apmag1corr[calstars] + zpteff
magerr = apmagerr1[calstars]

end





pro multifieldphot, imdata, stardata, zptdata, seeing=seeing, router=router

  allra = [0.D]
  alldec = [0.D]
  allmag = [0.]
  allmagerr = [0.]
  allfilt = ['']
  allchip = ['']
  alli = [-1]

  for i = 0, n_elements(imdata)-1 do begin
     print, imdata[i].filename
     filt = imdata[i].filt
;     chip = imdata[i].chip
;     dich = imdata[i].dich
     chip = 0
     dich = 0
     printf, 3, clip(imdata[i].filename,23)+' '+clip(fix(imdata[i].exp),3)+' '+clip(filt,2)+' '+clip(dich,3)+' '+fpr(imdata[i].air,1.2)+' | ', format='($,A)'
     
;     fi = where(zptdata.filt eq filt and zptdata.dich eq dich, ct)
     fi = where(zptdata.filt eq filt, ct)
     if ct gt 0 then begin
        zpt1 = zptdata[fi].zpt1
        airterm = zptdata[fi].airterm
        fieldphot, imdata[i].filename, ra, dec, mag, magerr, zpt1=zpt1, airterm=airterm, seeing=seeing, router=router, getseeing=getseeing
     endif else begin
        fieldphot, imdata[i].filename, ra, dec, mag, magerr, zpt1=0, airterm=0, seeing=seeing, router=router, /relative, getseeing=getseeing
     endelse
;stop
     if n_elements(mag) eq 0 then begin
        print, 'Unable to determine field photometry for ', imdata[i].filename
        continue
     endif

     imdata[i].seeing = getseeing

     good = where(mag lt 50 and mag gt -50)
     good = good[sort(mag[good])]
     allra  = [allra, ra[good]]   ; append good stars to megalist
     alldec = [alldec, dec[good]]
     allmag = [allmag, mag[good]]
     allmagerr = [allmagerr, magerr[good]]
     allfilt = [allfilt, replicate(imdata[i].filt,n_elements(good))]
     allchip = [allchip, replicate(imdata[i].chip,n_elements(good))]
     alli = [alli, replicate(i,n_elements(good))]
  endfor

  nstardata = n_elements(allra)-1
  stardata = replicate({ra:0.D, dec:0.D, mag:0., magerr:0., filt:'', chip:'', image:0, starid:0}, nstardata)
 
  stardata.ra = allra[1:*]
  stardata.dec = alldec[1:*]
  stardata.mag = allmag[1:*]
  stardata.magerr = allmagerr[1:*]
  stardata.filt = allfilt[1:*]
  stardata.chip = allchip[1:*]
  stardata.image = alli[1:*]
  stardata.starid = intarr(nstardata)-1

  cosdec = cos(median(alldec, /even)*!pi/180.)
  s = 0
  for a = 0, nstardata-1 do begin
     if stardata[a].starid eq -1 then begin
        stardata[a].starid = s
        s = s + 1
     endif
;     neighbors = findneighbors(stardata.ra*cosdec, stardata.dec, a, mindist=1./3668200., count=ct)
     neighbors = findneighbors(stardata.ra*cosdec, stardata.dec, a, mindist=1./3600., count=ct) ;1./12000. = 0.3 arcsec
     if ct gt 0 then stardata[neighbors].starid = stardata[a].starid
  endfor
  ns = s
return

end

pro makecal, cat, stardata, imdata, calname, object, uncthresh
;imdata should be imdata[thisobj]

    common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite
    common lrisfail, flatfail, catastrofail, relastrofail, fullastrofail, extractfail, wavsolfail, wavsolwarn

    filters = unique(stardata.filt)
    nfilt = n_elements(filters)
    maxs = max(stardata.starid)
    nexp = n_elements(unique(stardata.image))
    cosdec = cos(median(stardata.dec, /even)*!pi/180.)


    ; DIRECT cal - just get stars straight from the catalog that are in all images / not saturated

    caloutfile = imworkingdir + strlowcase(object) + '.'+strlowcase(calname)+'.direct.cal'

    openw, 1, caloutfile
    tstr = '#'
    tstr = tstr +      clip('RA',10-1)
    tstr = tstr + ' '+ ' ' + clip('dec',10-1)
    for f = 0, n_elements(filters)-1 do begin
       tstr = tstr + ' '+ clip(filters[f],6)
       tstr = tstr + ' '+ 'unc  '
    endfor
    printf, 1, tstr

    nminf = intarr(n_elements(filters)) ; star must be catalogued in this many images to be included
    for f = 0, n_elements(filters)-1 do begin
       maxoverlap = 0
       for s = 0, maxs do begin
          w = where(stardata.starid eq s and stardata.filt eq filters[f], n)
          if n gt maxoverlap then maxoverlap = n  ; new code block
       endfor
       nminf[f] = maxoverlap*0.5
    endfor

    nprinted = 0
    for s = 0, maxs do begin
       w = where(stardata.starid eq s, n)
       if n lt maxoverlap or n lt 1 then continue   ; n lt nexp-nfilt-1 
       if n eq 1 then begin
         avgra = stardata[w[0]].ra
         avgdec = stardata[w[0]].dec
       endif else begin
         avgra = median(stardata[w].ra, /even)
         avgdec = median(stardata[w].dec, /even)
      endelse

       smatch = nearest(avgra*cosdec, avgdec, cat.ra*cosdec, cat.dec, mindist=1./3600., count=ct)
       print,'SMATCH = ',smatch,s
;348.072833, 10.770431
       if smatch eq -1 then continue

       sstr = ''
       sstr = sstr +       fpr(cat[smatch].ra,3.6)
       sstr = sstr + ' ' + fpr(cat[smatch].dec,3.6)
       usestar = 0
       for f = 0, n_elements(filters)-1 do begin
          unc = 99.
          filt = filters[f]
          mf = (where(filt eq cat[smatch].filts, ct)) [0]
          if ct gt 0 then begin
            mag = cat[smatch].mags[mf]
            unc = cat[smatch].maguncs[mf]
          endif else begin
            mag = 999
            unc = 999
          endelse

          nexpfilt = n_elements(unique(stardata[where(stardata.filt eq filt)].image)) ; could use imdata again
          dum = where(stardata.starid eq s and stardata.filt eq filters[f], ngoodfilt)
          
print,s,f,mag,unc,filters[f],ngoodfilt
;stop
                                               ; formerly nexpfilt-1.  need to come up with a better solution for RL chips.
          if unc ge -1. and unc lt uncthresh and ngoodfilt ge nminf[f] and mag ge 0. and mag lt 99. then begin  ; print the mag
              usestar = 1
              sstr = sstr + ' ' + fpr(mag,2.3) + ' ' + fpr(unc,1.3)
              nprinted = nprinted + 1
          endif else begin
              sstr = sstr + ' ' + '-     '      + ' ' + '-    '
           endelse
;print,usestar
;stop
          ;print, s, ct, mag, unc, ngoodfilt, nexpfilt, usestar
          ;print, s, ct, mag, unc, ngoodfilt, nminf[f], usestar

       endfor
       if usestar then printf, 1, sstr
    endfor
    close, 1
    if nprinted gt 0 then begin
       print, 'Catalog printed to ', caloutfile, ' (',clip(nprinted),' stars)'
    endif else begin
       print, 'WARNING - no stars overlap all fields!  Cannot make direct catalog.'
       ;return  - no, keep going, although this might be perilous
    endelse

  
    ; EXTENDED cal using the shortest exposure

    ; need to make this look at both left and right, as usual.  
    ; loop (right, left, concatenate)?

    if strlowcase(object) eq 'grb051008' or object eq '051008' then return ; skip this due to bright star
    if (strlowcase(object) eq 'grb070810b' or strlowcase(object) eq '070810b') and strpos(calname,'sdss') ge 0 then return ; sloan is apparently near the field but not near enough
    if (strlowcase(object) eq 'grb080210' or object eq '080210') and strpos(calname,'sdss') ge 0 then return ; sloan is apparently near the field but
    if (strlowcase(object) eq 'grb080603a' or strlowcase(object) eq '080603a') and strpos(calname,'sdss') ge 0 then return ; sloan is apparently near the field but

    noshort = 0
    fdzpt = fltarr(n_elements(filters))
    fdzptstdev = fltarr(n_elements(filters))
    for f = 0, n_elements(filters)-1 do begin
       filt = filters[f]

       checkcatfilt = where(filt eq cat[0].filts, ctf)
       if ctf eq 0 then begin
          noshort = noshort+1
          print, 'Filter ', filt, ' not in catalog.'
          continue
       endif

       ff = where(imdata.filt eq filt)

       shortestexp = ff[(where(imdata[ff].exp eq min(imdata[ff].exp))) [0]]  ; the i index of the shortest file, hopefully
       if min(imdata[ff].exp) gt 1000 then begin
           noshort = noshort + 1
           continue
       endif

       print, 'Calibrating from ', imdata[shortestexp].filename

       dzpt = [-1]
       for s = 0, maxs do begin
           w = where(stardata.starid eq s and stardata.image eq shortestexp, n)
           if n eq 0 then continue
           if n gt 1 then continue ; n > 1 means two matching sources - sketchy (shouldn't happen?)
           avgra = stardata[w].ra
           avgdec = stardata[w].dec
           smatch = nearest(avgra*cosdec, avgdec, cat.ra*cosdec, cat.dec, mindist=1./3600, count=ct)

           if ct eq 0 then continue

           mag = -1
           mf = where(filt eq cat[smatch].filts, ctf)
           if ctf gt 0 then begin
             mag = cat[smatch].mags[mf]
             unc = cat[smatch].maguncs[mf]
           endif

           if ctf gt 0 and mag gt 0. and unc ge 0. and unc lt 99. then dzpt = [dzpt, mag-stardata[w].mag]
       endfor
       if n_elements(dzpt) eq 1 then continue ; some sort of failure has occurred
       dzpt = dzpt[1:*]
       ntot = n_elements(dzpt)
       meddzpt = median(dzpt, /even)
       if ntot gt 1 then stddzpt = stdev(dzpt) else stddzpt = 2
       clipdzpt = dzpt[where(abs(dzpt-meddzpt) lt 3*stddzpt, ngood)]
       nout = ntot - ngood
       fdzpt[f] = median(clipdzpt, /even)
       if ngood gt 1 then fdzptstdev[f] = stdev(clipdzpt) else  fdzptstdev[f] = 2
       absdev = median(abs(clipdzpt - fdzpt[f]), /even)

       print, filt, ' catalog zeropoint adjustment:', fpr(fdzpt[f],3.3), ' +/- ', fpr(fdzptstdev[f],1.3), $
                    ' (', clip(ntot), ' stars, ', clip(nout), ' outliers)'
       printf, 3, '# dzpt_'+strmid(calname,5)+'_'+filt+' =', fpr(fdzpt[f],2.3), ' ',fpr(fdzptstdev[f],1.3), ' ', fpr(absdev,1.3), ' ', clip(ngood,3), clip(nout,3)

    endfor

    if noshort lt n_elements(filters) then begin
      caloutfile = imworkingdir + strlowcase(object) + '.'+strlowcase(calname)+'.extend.cal'
      openw, 1, caloutfile
      tstr = '#'
      tstr = tstr +      clip('RA',10-1)
      tstr = tstr + ' '+ ' ' + clip('dec',10-1)
      for f = 0, n_elements(filters)-1 do begin
         tstr = tstr + ' '+ clip(filters[f],6)
         tstr = tstr + ' '+ 'unc  '
      endfor
      printf, 1, tstr

      for s = 0, maxs do begin
         usestar = 0
         w = where(stardata.starid eq s, n)
         if n eq 0 then continue ; this shouldn't happen but it does?
         avgra = median([stardata[w].ra], /even)
         avgdec = median([stardata[w].dec], /even)
         sstr = ''
         sstr = sstr +       fpr(avgra,3.6)
         sstr = sstr + ' ' + fpr(avgdec,3.6)
         for f = 0, n_elements(filters)-1 do begin
            

            ;shortexp = ff[(where(imdata[ff].exp lt 100)) [0]]  ; the i index of the shortest file, hopefully
            shortestexp = ff[(where(imdata[ff].exp eq min(imdata[ff].exp))) [0]]
            w = where(stardata.starid eq s and stardata.image eq shortestexp, n)

            nexpfilt = n_elements(unique(stardata[where(stardata.filt eq filt)].image)) ; could use imdata again
            dum = where(stardata.starid eq s and stardata.filt eq filters[f], ngoodfilt)

            ; n = 1:  good in this reference frame
            ; ngood < nexp-1 : good in all (but 1) other frames too

            if n eq 1 then begin
              mag = stardata[w].mag + fdzpt[f]
              unc = stardata[w].magerr
              if ngoodfilt ge nminf[f] and unc lt 0.05 then begin   ; again, formerly nexpfilt-1
                 sstr = sstr + ' ' + fpr(mag,2.3) + ' ' + fpr(unc,1.3)
                 usestar = usestar + 1
              endif else begin
                 sstr = sstr + ' ' + '-    ' + ' ' + '-    '
              endelse
            endif

          endfor
          if usestar then printf, 1, sstr
      endfor
      close, 1
    endif
    print, 'Catalog printed to ', caloutfile

return

end


pro lrisautofieldphot, red=red, blue=blue, seeingarcsec=seeingarcsec, forcephotometric=forcephotometric, objcat=objcat, chip=chip

common lrisconfig, lrisautopath, refdatapath, defaultspath
common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite
common lrisfail, flatfail, catastrofail, relastrofail, fullastrofail, extractfail, wavsolfail, wavsolwarn


if n_elements(imworkingdir) eq 0 then imworkingdir = '' ; can run outside of lrisautoproc

;if n_elements(chip) gt 0 then chipchar = strmid(chip,0,1) else chipchar = '[rl]' ; this actually works!
;if keyword_set(blue) then imagefiles =  findfile(imworkingdir+'a*f*b??????_????'+chipchar+'.fits')
;if keyword_set(red)  then imagefiles =  findfile(imworkingdir+'a*f*r??????_????'+chipchar+'.fits')
wildcharimg = '?????????????????_img_?'
imagefiles = choosefiles(imworkingdir+'a*f*p2'+wildcharimg+'.fits')

;if keyword_set(blue) eq 0 and keyword_set(red) eq 0 then begin
;   lrisautofieldphot, /red,  seeing=seeingarcsec, forcephotometric=forcephotometric, objcat=objcat
;   lrisautofieldphot, /blue,  seeing=seeingarcsec, forcephotometric=forcephotometric, objcat=objcat
;   return
;endif

if n_elements(unique(imagefiles)) eq 1 and imagefiles[0] eq '' then return
imagefiles = imagefiles[where(imagefiles ne '')]
close, /all

h = headfits(imagefiles[0])
mjd = float(sxpar(h,'MJD-OBS'))

;if keyword_set(blue) then camera = 'blue'
;if keyword_set(red) and mjd ge 55562 then camera = 'red3'
;if keyword_set(red) and mjd ge 54983 and  mjd lt 55562 then camera = 'red2' ; after 2009-06-01
;if keyword_set(red) and mjd lt 54983 then camera = 'red1'


if keyword_set(seeingarcsec) eq 0 then begin
  if file_test('seeing.txt') then begin
    rdfloat, 'seeing.txt', see1, see2
    seeingarcsec = [see1[0], see2[0]]
  endif else begin
    seeingarcsec = [0.4, 3.0]
  endelse
endif
;if camera eq 'red1' then seeing = seeingarcsec/0.210
;if camera eq 'red2' then seeing = seeingarcsec/0.135
;if camera eq 'red3' then seeing = seeingarcsec/0.135
;if camera eq 'blue' then seeing = seeingarcsec/0.135

seeing = seeingarcsec/0.135
zptfile = imworkingdir+'zeropoints.dat'
if file_test(zptfile) eq 0 then begin
   zptfile = defaultspath+'/'+'zeropoints.dat'
   photometric = 0
endif else begin
   photometric = 1
endelse
openr, 1, zptfile
utdatestr = ''
aperstr = ''
readf, 1, utdatestr
readf, 1, aperstr
close, 1
fullaperrad = float((strsplit(aperstr,/extract)) [0])
;zptdata = replicate({filt:'', dich:'', zpt1:0., airterm:0.},countlines(zptfile)-1)
;readcol, zptfile, zfilt, zdich, zzpt1, zzpt1unc, zairterm, zairtermunc, format='a,a,f,f,f,f', /silent
zptdata = replicate({filt:'', zpt1:0., airterm:0.},countlines(zptfile)-1)
readcol, zptfile, zfilt, zzpt1, zzpt1unc, zairterm, zairtermunc, format='a,a,f,f,f', /silent
zptdata = zptdata[0:n_elements(zfilt)-1]
zptdata.filt = zfilt
;zptdata.dich = zdich
zptdata.zpt1 = zzpt1
zptdata.airterm = zairterm

;grbsposfile = '~/research/redux/allgrbpos.txt'
if n_elements(objcat) gt 0 then begin
   grbsposfile = objcat
   readcol, grbsposfile, grbname, grbra, grbdec, format='a,f,f', /silent
   ;grbname = 'GRB'+grbname
   ngrb = n_elements(grbra)
endif

nim = n_elements(imagefiles)
imdata = replicate({filename:'', isgrb:0, object:'', filt:'', dich:'', chip:'', exp:0., air:0., seeing:0., fluxratio:0.}, nim)
for i = 0, nim-1 do begin
  h = headfits(imagefiles[i])
  ;dum = where(sxpar(h,'OBJECT') eq stdname, ctmatch)
  ;if ctmatch eq 0 then imagefiles[i] = ''
  ra = sxpar(h,'CRVAL1')  ; use the WCS header keywords to get the real center
  dec = sxpar(h,'CRVAL2')

  imdata[i].filename = imagefiles[i]
  imdata[i].filt = clip(sxpar(h,'FILTER'))
;  if imdata[i].filt eq 'RG850' then imdata[i].filt = 'z'
;  if imdata[i].filt eq 'G' then imdata[i].filt = 'g'
;  if imdata[i].filt eq 'U' then imdata[i].filt = 'u'
;  imdata[i].chip = clip(sxpar(h,'CHIP'))
;  if imdata[i].chip eq '' or imdata[i].chip eq '0' then begin
;     dotpos = strpos(imagefiles[i],'.',/reverse_search)
;     camchar = strmid(imagefiles[i],dotpos-1,1)
;     if camchar eq 'l' then imdata[i].chip = 'left'
;     if camchar eq 'r' then imdata[i].chip = 'right'
;     if camchar eq 'o' then imdata[i].chip = 'both'
;  endif
  imdata[i].exp = clip(sxpar(h,'ELAPTIME'))
;  imdata[i].exp = 
;  imdata[i].dich = clip(sxpar(h,'DICHNAME'))
  imdata[i].air = sxpar(h,'AIRMASS')

  if n_elements(objcat) gt 0 then begin
    distance = fltarr(ngrb)
    for r = 0, ngrb-1 do begin
       dist = 3600*(abs(dec-grbdec[r]) + abs(ra-grbra[r]))
       if dist lt 300 then gcirc, 2, ra, dec, grbra[r], grbdec[r], dist
       distance[r] = dist
    endfor
    mindist = min(distance)
    minr = (where(distance eq mindist)) [0]
    if mindist lt 300. then begin
      imdata[i].isgrb = 1
      print, clip(sxpar(h,'OBJECT')), ' -> ',  grbname[minr]
      imdata[i].object = grbname[minr]     ;clip(sxpar(h,'OBJECT'))
    endif 
  endif else begin
    imdata[i].isgrb = 1
    imdata[i].object = clip(sxpar(h,'TARGNAME'))
  endelse

endfor
nim = i ; ignore non-GRBs?
imdata = imdata[where(imdata.isgrb)]
imdata = imdata[where(imdata.filt ne 'OG570')]
objects = unique(imdata.object)

openw, 3, imworkingdir+'autophotsummary.txt'
openw, 4, imworkingdir+'autophotsummaryflux.txt'

for b = 0, n_elements(objects)-1 do begin  ; in this context objects means object keyword of filename
  object = objects[b]
  if object eq '' then continue

  imfiles = where(imdata.object eq object, nobjim)  
  objimdata = imdata[imfiles]                           ; only look at images of the relevant field

  print
  print, object
  printf, 3, '# ', object
  printf, 4, '# ', object

  multifieldphot, objimdata, stardata, zptdata, seeing=seeing, router=fullaperrad
     ; stardata is the only *output*: list of ra, dec, mag, magerr, filt:'', image (#), starid (#)

  ; Compare exposures to tell if this was a photometric observation

  filters = unique(stardata.filt)
;  chips = unique(stardata.chip)
  for f = 0, n_elements(filters)-1 do begin
     print, filters[f], ' transmission'
;  for c = 0, n_elements(chips)-1 do begin
;     thisfilt = where(stardata.filt eq filters[f] and stardata.chip eq chips[c], ct)
     thisfilt = where(stardata.filt eq filters[f], ct)

     if ct eq 0 then continue
     fstardata = stardata[thisfilt]
     
     fimagelist = unique(fstardata.image)   ; fimagelist is f subscripts, fstardata.image is o subscripts
     fimagedata = objimdata[fimagelist]     ; fimagedata is f subscripts
     fstarlist = unique(fstardata.starid)

     nfi = n_elements(fimagelist) ; number of filter/chip images
     dmagarr = fltarr(nfi) - 999
     dmagij = fltarr(nfi,nfi)
     ncompij = intarr(nfi,nfi)

     ;j = ni/2 ; just use the middle image to compare to

     ; probably should do something here when nstars gets extremely large (>1000) to boost speed
     if nfi gt 1000 then begin
      print, 'Crowded field, may run slowly...'
     endif

     for fi = 0, nfi-1 do begin ; primary image loop          ; fi and fj indices are for fimagelist
       for fj = 0, nfi-1 do begin ; compare to this image     
         if fi eq fj then continue
         oi = fimagelist[fi]                           ; oi and oi indices are for the object/starlist (all filters, chips)
         oj = fimagelist[fj]
         dmagarr = [-1]
         for fs = 0, n_elements(fstarlist)-1 do begin ; collect the stars
            os = fstarlist[fs]
            wi = where(fstardata.starid eq os and fstardata.image eq oi, cti)
            wj = where(fstardata.starid eq os and fstardata.image eq oj, ctj)
            if cti eq 1 and ctj eq 1 then begin
               sdm = (fstardata[wi].mag - fstardata[wj].mag) [0]
               dmagarr = [dmagarr, sdm]
            endif
         endfor
         ncompstars = n_elements(dmagarr)-1
         ncompij[fi,fj] = ncompstars
         if ncompstars gt 0 then dmagij[fi,fj] = median(dmagarr[1:*], /even) else dmagij[fi,fj] = !values.f_nan
         ;print, fimagedata[fi].filename, ' vs. ', fimagedata[fj].filename, ' :  dmag = ', fpr(dmagij[fi,fj],3.3), '  (', clip(ncompstars), ' stars)'

       endfor
     endfor

     dmag = fltarr(nfi)
     for fi = 0, nfi-1 do begin 
       for fj = fi+1, nfi-1 do begin 
          if dmagij[fi,fj] lt -10 or finite(dmagij[fi,fj]) eq 0 then continue
          newdmagij = dmagij[fi,fj] - dmag[fi] + dmag[fj]
          dmag[fj] = dmag[fj] - newdmagij
       endfor
     endfor
     fluxr = 10.0^(-dmag/2.5)

     for fi = 0, nfi-1 do begin
       objimdata[fimagelist[fi]].fluxratio = fluxr[fi] ; hopefully assigns correctly
       print, fimagedata[fi].filename, ' :  dmag = ', fpr(dmag[fi],2.3), '  relflux = ', fpr(fluxr[fi],2.3), '  seeing = ', fpr(fimagedata[fi].seeing,2.3)
       printf, 4, clip(fimagedata[fi].filename,23)+' '+clip(fix(fimagedata[fi].exp),3)+' '+clip(filters[f],2)+' '+clip(fimagedata[fi].dich,3)+' '+fpr(fimagedata[fi].air,1.2)+' | ', fpr(dmag[fi],2.3), ' ', fpr(fluxr[fi],2.3), ' ', fpr(fimagedata[fi].seeing,2.3)
     endfor

     maxabsdmag = max(abs(dmag))
print,stdev(dmag)
;stop
     if n_elements(dmag) gt 1 then stdevdmag = stdev(dmag) else stdevdmag = 0
     if stdevdmag gt 0.05 then begin
       print, 'Not photometric!'  
       photometricobs = 0
     endif else begin
       photometricobs = 1
     endelse
 ; endfor
  endfor

  ; Look up 2MASS for this field and find matches

  cra = median(stardata.ra, /even)
  cdec =  median(stardata.dec, /even)
;  sdss = getsdss(cra, cdec, 360.)
  craref=cra
  cdecref=cdec
  tmpsc = gettmpsc(cra, cdec, 360.)

  mindist = 3600. * min(sqrt(((stardata.ra-craref)*cos(cdecref*!pi/180.))^2 + (stardata.dec-cdecref)^2))
  uncthresh = 0.12
  if n_elements(tmpsc) gt 10 and mindist lt 120 then begin
    ;tmpsc = tmpsc[where(tmpsc.type eq 6)] ; stars only
    tmpsccat = replicate({ra:0.D, dec:0.D, filts:strarr(3), mags:fltarr(3), maguncs:fltarr(3)}, n_elements(tmpsc))
    tmpsccat.ra = tmpsc.ra
    tmpsccat.dec = tmpsc.dec
    tmpsccat.filts = ['J','H','K'] ; kind of stupidly inefficient, but whatever
    tmpsccat.mags[0] = tmpsc.j
    tmpsccat.mags[1] = tmpsc.h
    tmpsccat.mags[2] = tmpsc.k
    tmpsccat.maguncs[0] = 0
    tmpsccat.maguncs[1] = 0
    tmpsccat.maguncs[2] = 0
;stop
    makecal, tmpsccat, stardata, objimdata, 'tmpsc', strlowcase(object), uncthresh

  endif else begin
    print, 'Not a 2MASS field.'
  endelse



  ; Look up the Nickel-based catalog and use that as a calibrator

  ; in the future, point to a general catalog and look at everything

  gipos = strpos(object,'0') < strpos(object,'1')
  grb = strlowcase(strmid(object,gipos))
  grb = (strsplit(grb,'_',/extract)) [0]

  grbcalfile = '~/research/redux/nickelcal/grb' + grb + '.cat' 
  if strlen(grb) eq 6 and file_test(grbcalfile) eq 0 then $
    grbcalfile = '~/research/redux/nickelcal/grb' + grb + 'a.cat' 

  if file_test(grbcalfile) then begin
     ns = countlines(grbcalfile)-1
     nickelcat = replicate({ra:0.D, dec:0.D, filts:strarr(5), mags:(fltarr(5)-1), maguncs:(fltarr(5)-1)}, ns)

     inline = ''
     openr, 1, grbcalfile
     readf, 1, inline
     headings = strsplit(strmid(inline,1),/extract)
     nh = n_elements(headings)
     f = 0
     for i = 5, nh-1, 4 do begin
        nickelcat[*].filts[f] = headings[i]
        f = f + 1
     endfor
     nfilt = f
     for s = 0, ns-1 do begin
        readf, 1, inline
        inarr = strsplit(inline,/extract)
        nickelcat[s].ra = double(inarr[0])
        nickelcat[s].dec = double(inarr[1])
        for f = 0, nfilt-1 do begin
           if clip(inarr[5+f*4]) ne '-' then nickelcat[s].mags[f] = float(inarr[5+f*4])
           if clip(inarr[5+f*4+2]) ne '-' and inarr[5+f*4+1] ne '-' then nickelcat[s].maguncs[f] = float(inarr[5+f*4+1])/3. > float(inarr[5+f*4+2]) 
        endfor
     endfor
     close, 1

     makecal, nickelcat, stardata, objimdata, camera + '.nickel', strlowcase(object), uncthresh

  endif else begin
     print, 'Not a Nickel field.'
  endelse


  ; Do an absolute calibration
  if (keyword_set(photometric) and photometricobs) or keyword_set(forcephotometric) then begin

    landoltoutfile = imworkingdir + strlowcase(object)+'.'+camera + '.landolt.cal'
    openw, 1, landoltoutfile


    ns = max(stardata.starid)+1
    filters = unique(stardata.filt)
    cosdec = cos(median(stardata.dec, /even)*!pi/180.)
    tstr = '#'
    tstr = tstr +      clip('RA',10-1)
    tstr = tstr + ' '+ ' ' + clip('dec',10-1)
    tstr = tstr + ' '+ clip('unc_RA',8)
    tstr = tstr + ' '+ clip('unc_dec',8)
    tstr = tstr + ' '+ clip('n',2)
    for f = 0, n_elements(filters)-1 do begin
       tstr = tstr + ' '+ clip(filters[f],6)
       tstr = tstr + ' '+ clip('var',5)
       tstr = tstr + ' '+ clip('unc',5)
       tstr = tstr + ' '+ clip('n',2)
    endfor
    printf, 1, tstr
    for s = 0, ns-1 do begin
       w = where(stardata.starid eq s, n)
       if n eq 0 then continue ; should almost never happen, but occasionally a flpt error...?
       if n gt 1 then begin
          avgra = median(stardata[w].ra, /even)
          avgdec = median(stardata[w].dec, /even)
          uncra = stdev(1.0D*stardata[w].ra)*3600.*cosdec
          uncdec = stdev(1.0D*stardata[w].dec)*3600.
       endif else begin
          avgra = stardata[w[0]].ra
          avgdec = stardata[w[0]].dec
          uncra = 1
          uncdec = 1
       endelse
   
       sstr = ''
       sstr = sstr +       fpr(avgra,3.6)
       sstr = sstr + ' ' + fpr(avgdec,3.6)
       sstr = sstr + ' ' + fpr(uncra,1.6)
       sstr = sstr + ' ' + fpr(uncdec,1.6)
       sstr = sstr + ' ' + clip(n,2)
       for f = 0, n_elements(filters)-1 do begin
          thisfilt = where(stardata[w].filt eq filters[f], ct)
          if ct eq 0 then begin
             sstr = sstr + ' ' + clip('-',6)
             sstr = sstr + ' ' + clip('-',5)
             sstr = sstr + ' ' + clip('-',5)
             sstr = sstr + ' ' + clip(0,2)
          endif else begin
             if ct gt 1 then begin
                avgmag = median(stardata[w[thisfilt]].mag, /even)
                avgmagerr = median(stardata[w[thisfilt]].magerr, /even)
                varmag = stdev(stardata[w[thisfilt]].mag)
             endif else begin
                avgmag = stardata[w[thisfilt[0]]].mag
                avgmagerr = stardata[w[thisfilt[0]]].magerr
                varmag = 1
             endelse
             sstr = sstr + ' ' + fpr(avgmag,2.3)
             sstr = sstr + ' ' + fpr(varmag,1.3)
             sstr = sstr + ' ' + fpr(avgmagerr,1.3)
             sstr = sstr + ' ' + clip(ct,2)
          endelse
       endfor
       ;if n ge 3 then 
       printf, 1, sstr
    endfor
    close, 1
    print, 'Catalog printed to ', landoltoutfile
  endif else begin
    print, 'Not photometric.'
  endelse
  
  print

endfor
close, 3
close, 4


end



