@readlandoltlris

; Calculate airmass solutions and zeropoints as a function of airmass.

pro lrisautostdphot, blue=blue, red=red, photaperrad=photaperrad, fullaperrad=fullaperrad

common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite
common lrisconfig, lrisautopath, refdatapath, defaultspath

if n_elements(refdatapath) eq 0 then refdatapath = '~/progs/idl/lris/refdata/'
if n_elements(imworkingdir) eq 0 then imworkingdir = ''
if n_elements(overwrite) eq 0 then overwrite = 0

close, /all

if keyword_set(blue) then imagefiles = findfile(imworkingdir+'a*f*b??????_????'+'*'+'.fits')
if keyword_set(red)  then imagefiles = findfile(imworkingdir+'a*f*r??????_????'+'*'+'.fits')
if keyword_set(blue) eq 0 and keyword_set(red) eq 0 then begin
   lrisautostdphot, /red, photaperrad=photaperrad, fullaperrad=fullaperrad
   lrisautostdphot, /blue, photaperrad=photaperrad, fullaperrad=fullaperrad
   return
endif

validfiles = where(imagefiles ne '', ct)
if ct eq 0 then begin
   cd, current=pwd
   print, 'No LRIS files found in directory '+pwd
   return
endif


h = headfits(imagefiles[0])
mjd = float(sxpar(h,'MJD-OBS'))

camera = ''
if keyword_set(blue) then camera = 'blue'
if keyword_set(red) and mjd ge 55562 then camera = 'red3'
if keyword_set(red) and mjd ge 54983 and  mjd lt 55562 then camera = 'red2' ; after 2009-06-01
if keyword_set(red) and mjd lt 54983 then camera = 'red1'

if n_elements(photaperrad) eq 0 then begin
  if camera eq 'red1' then photaperrad = 6.5  else photaperrad = 10.
endif
if n_elements(fullaperrad) eq 0 then begin
  if camera eq 'red1' then fullaperrad = 22.5 else fullaperrad = 35.
endif
minedgedist = fullaperrad

if file_test(imworkingdir+'zeropoints.'+camera+'.dat') and overwrite eq 0 then begin
   print, 'File exists: ', imworkingdir+'zeropoints.'+camera+'.dat'
   return
endif


landolt = readlandoltlris(file=refdatapath+'/landolt2009.txt')
nref = n_elements(landolt)



isstd = intarr(n_elements(imagefiles))

for i = 0, n_elements(imagefiles)-1 do begin
  h = headfits(imagefiles[i])
  ;dum = where(sxpar(h,'OBJECT') eq stdname, ctmatch)
  ;if ctmatch eq 0 then imagefiles[i] = ''
  ra = sxpar(h,'CRVAL1')  ; use the WCS header keywords to get the real center
  dec = sxpar(h,'CRVAL2')
  trapdoor = clip(sxpar(h,'TRAPDOOR'))
  exptime = sxpar(h,'ELAPTIME')

  ;filter = clip(sxpar(h,'FILTER'))
  ;if filter ne 'R' then continue

  if trapdoor eq 'open' and exptime lt 30 then begin
     distance = fltarr(nref)
     for r = 0, nref-1 do begin
        dist = 3600*(abs(dec-landolt[r].dec) + abs(ra-landolt[r].ra));
        if dist lt 300 then gcirc, 2, ra, dec, landolt[r].ra, landolt[r].dec, dist
        distance[r] = dist
     endfor
     mindist = min(distance)
     minr = (where(distance eq mindist)) [0]
     if mindist lt 300. then isstd[i] = 1
  endif
endfor

if total(isstd) eq 0 then begin
   print, 'No Landolt standards found.'
   return
endif

imagefiles = imagefiles[where(isstd)]



; now we just have images that match known standards.


sobs = replicate({imagename:'', filter:'', starname:'', date:'', instmag:0., catmag:0., obsmag:0., caterr:0., photoerr:0., airmass:0., good:1},1000)

landolt = readlandoltlris()
nref = n_elements(landolt)

; iterate through the images 
obsn = 0
for f = 0, n_elements(imagefiles)-1 do begin

     data = float(mrdfits(imagefiles[f], 0, h, /silent, /fscale))
     ra = sxpar(h,'CRVAL1')  ; use the WCS header keywords to get the real center
     dec = sxpar(h,'CRVAL2')
     filter = clip(sxpar(h,'FILTER'))
     if filter eq 'G' then filter = 'g'
     if filter eq 'RG850' then filter = 'z'
     if filter eq 'U' then filter = 'u'
     dichroic = clip(sxpar(h,'DICHNAME'))
     exptime = sxpar(h, 'ELAPTIME')
     airmass = sxpar(h,'AIRMASS')
     gain = sxpar(h,'GAIN')
     date = strmid(sxpar(h,'DATE',/silent),0,10)

     nx = (size(data)) [1]
     ny = (size(data)) [2]

     ; find nearby stars in the full Landolt catalog
     distance = fltarr(nref)
     for r = 0, nref-1 do begin
        dist = 3600*abs(dec-landolt[r].dec) + 3600*abs(ra-landolt[r].ra)
        if dist lt 600 then gcirc, 2, ra, dec, landolt[r].ra, landolt[r].dec, dist
        distance[r] = dist
     endfor
     nearby = where(distance lt 500, nnearby)
     if nnearby eq 0 then continue
     nearbylandolt = landolt[nearby]


     print, imagefiles[f], ' (', filter, '/', dichroic, ')', '  ', string(exptime,format='(I4)'), 's, air=',fpr(airmass,1.3)

     adxy, h, nearbylandolt.ra, nearbylandolt.dec, nearbylandoltx, nearbylandolty

     apdiff = fltarr(nnearby) + !values.f_NaN
     apdiffunc = fltarr(nnearby) + !values.f_NaN
     iobs2obs = intarr(nnearby) - 1

     ; iterate through stars, add them (w/photometry) if they're in good places in the field
     iobsn = 0 
     for istar = 0, nnearby-1 do begin
        if nearbylandoltx[istar] lt 0 or nearbylandoltx[istar] gt nx-1 or $
           nearbylandolty[istar] lt 0 or nearbylandolty[istar] gt ny-1 then begin
           print, clip(nearbylandolt[istar].name,12), ' out of field.'
           continue
        endif
        if nearbylandoltx[istar] lt minedgedist    or nearbylandolty[istar] lt minedgedist   or $
           nearbylandoltx[istar] ge nx-minedgedist or nearbylandolty[istar] ge ny-minedgedist then begin
           print, clip(nearbylandolt[istar].name,12), ' near edge of chip.'
           continue
        endif
 ;      if  abs(nearbylandoltx[istar]-255) gt 2        and abs(nearbylandoltx[istar]-783) gt 2  $             ;  this was for old LRIS-2 I guess....???
 ;      endif

           sobs[obsn].imagename = imagefiles[f]
           sobs[obsn].filter    = filter+'_'+dichroic
           sobs[obsn].airmass   = airmass
           sobs[obsn].starname  = nearbylandolt[istar].name
           sobs[obsn].catmag    = getstarmag(nearbylandolt[istar], filter, unc=unc)
           sobs[obsn].caterr    = sqrt(unc^2 + 0.01^2)   ; 0.01 is an estimate of additional uncertainty due to color terms, etc.
           sobs[obsn].date      = date

           ; do the aperture photometry
           cntrd, data, nearbylandoltx[istar], nearbylandolty[istar], xcen, ycen, 15  ; get the exact position
           aper,  data, xcen, ycen, mag, errap, sky, skyerr, gain, [photaperrad,fullaperrad], [40,50], [-100,50000], /silent
           mag = mag + 2.5*alog10(exptime) ; adjust for exposure time!!!
           ; don't bother to scale sky

           ;if mag gt 50 then continue
           if mag[0] gt 50 or mag[1] gt 50 then begin
              print, clip(nearbylandolt[istar].name,12), ' produced a bad photometry value'
              continue
           endif

           print, clip(nearbylandolt[istar].name,12), '[', fpr(nearbylandoltx[istar],4.2),',',fpr(nearbylandolty[istar],4.2),']->[',fpr(xcen,4.2),',',fpr(ycen,4.2),']', $
                 fpr(sobs[obsn].catmag,4.3),  ' -', fpr(mag[1],3.3),  ' =', fpr(sobs[obsn].catmag-mag[1],4.3);, mag[1]


           sobs[obsn].instmag = mag[0]
           sobs[obsn].photoerr = sqrt((errap[1])^2 + 0.02^2) ; use the outer-aperture radius (being conservative)

           apdiff[iobsn] = mag[0]-mag[1]
           apdiffunc[iobsn] = errap[1]
           iobs2obs[iobsn] = obsn

           obsn = obsn + 1
           iobsn = iobsn + 1
     endfor

     if iobsn gt 0 then begin
        apdiff = apdiff[0:iobsn-1]
        apdiffunc = apdiffunc[0:iobsn-1]
        iobs2obs = iobs2obs[0:iobsn-1]

        ;print, apdiff
        ;print, apdiffunc

        apcorr = median(apdiff[where(apdiffunc le median(apdiffunc))], /even)

        ;print, apcorr

        sobs[iobs2obs].instmag = sobs[iobs2obs].instmag - apcorr

        ;for io = 0, iobsn-1 do begin

           ;sobs[iobs2obs[io]].instmag = sobs[iobs2obs[io]].instmag - apcorr
           ;sobs[iobs2obs[io]].photoerr = sqrt((errap [1])^2 + 0.02^2)  ; 0.02 is an estimate of flatfielding, pixel, etc. errors
                                                                       ; be conservative by using the outer-aperture uncertainty and the inner-aperture magnitude

           ;oobs[obsn].obszpt = 
           ;oobs[obsn].obszptunc = 
        ;endfor

     endif 
endfor
sobs = sobs[0:obsn-1]



filts = unique(sobs.filter)
filts = filts[sort(filtwv(filts))]

fzpt1 = fltarr(n_elements(filts))
fzpt1unc = fltarr(n_elements(filts))
fairmassterm = fltarr(n_elements(filts)) 
fairmasstermunc = fltarr(n_elements(filts)) 

openw, 1, imworkingdir+'zeropoints.'+camera+'.dat'
datestr = ''
dates = unique(sobs.date)
for i = 0, n_elements(datestr)-1 do datestr = datestr + dates[i] + ' '
printf, 1, datestr
printf, 1, clip(fullaperrad), ' ', clip(photaperrad)

for f = 0, n_elements(filts)-1 do begin
   filter = filts[f]

   thisfilt = where(sobs.filter eq filts[f])
   filtsobs = sobs[thisfilt]
   zpt = filtsobs.catmag - filtsobs.instmag
   zptunc = sqrt(filtsobs.caterr^2 + filtsobs.photoerr^2)
   ntot = n_elements(zpt)

   bad = where(filtsobs.instmag gt 50 or filtsobs.instmag lt -20 or $
               filtsobs.catmag gt 50  or filtsobs.catmag lt -20, ctbad)
   if ctbad gt 0 then begin
      filtsobs[bad].good = -1
      ;print, ctbad, filtsobs[bad].instmag
   endif

   faint = where(filtsobs.good gt 0 and zptunc gt median(zptunc)*3, ctfaint, complement=bright)
   if ctfaint gt 0 then begin
      filtsobs[faint].good = -1  
   endif
   
   cloudy = where(filtsobs.good gt 0 and zpt-(filtsobs.airmass-1.)*(-0.3) lt min(zpt[bright]-3*zptunc[bright])-0.25, ctcloudy, complement=ok)
   if ctcloudy gt 0 then begin
     print,  clip(ctcloudy), '/', clip(ntot), ' cloudy observations (', fpr(100.*ctcloudy/ntot,2.1), '% of total)'
     filtsobs[cloudy].good = -1
   endif

   wacky =  where(filtsobs.good gt 0 and abs(zpt - median(zpt[ok])) gt 1.5, ctwacky)
   if ctwacky gt 0 then begin
     print,  clip(ctwacky), '/', clip(ntot), ' wacky observations (', fpr(100.*ctwacky/ntot,2.1), '% of total)'
     filtsobs[wacky].good = -1
   endif



   ; solve the airmass term by going star by star
  
   gfiltsobs = filtsobs[where(filtsobs.good gt 0)]
   flt = (strsplit(filter,'_',/extract)) [0]
   if flt eq 'z' then airmasstermrange=-[0,0.15]
   if flt eq 'I' then airmasstermrange=-[0,0.15]
   if flt eq 'R' then airmasstermrange=-[0,0.2]
   if flt eq 'V' then airmasstermrange=-[0,0.25]
   if flt eq 'g' then airmasstermrange=-[0,0.3]
   if flt eq 'B' then airmasstermrange=-[0,0.4]
   if flt eq 'u' then airmasstermrange=-[0,0.7]
   if flt eq 'z' then defaultterm = -0.06
   if flt eq 'I' then defaultterm = -0.07
   if flt eq 'R' then defaultterm = -0.11
   if flt eq 'V' then defaultterm = -0.12
   if flt eq 'g' then defaultterm = -0.14
   if flt eq 'B' then defaultterm = -0.17
   if flt eq 'u' then defaultterm = -0.37
   jointairmassfit, gfiltsobs.airmass, gfiltsobs.starname, gfiltsobs.catmag-gfiltsobs.instmag, gfiltsobs.photoerr, gfiltsobs.caterr, $
                     zpt1, zpt1unc, airmassterm, airmasstermunc, airmasstermrange=airmasstermrange, defaultterm=defaultterm, outliers=outliers, starplot='plotairstar'+camera+filts[f]+'.ps'
   if min(outliers) ge 0 then begin; make sure there were some outliers
      gfiltsobs[outliers].good = 0
      filtsobs[where(filtsobs.good gt 0)].good = gfiltsobs.good ; push this back to the main filtsobs struct array
      sobs[thisfilt].good = filtsobs.good                   ;                to the main sobs struct array
   endif

   fzpt1[f] = zpt1
   fzpt1unc[f] = zpt1unc
   fairmassterm[f] = airmassterm
   fairmasstermunc[f] = airmasstermunc

   print, filter+'-band zeropoint = ', fpr(zpt1,2.3), ' (+/- ', fpr(zpt1unc,1.3), ')  - ', fpr(-airmassterm,2.3),  ' (+/- ', fpr(airmasstermunc+0.0001,2.3), ')  *(airmass-1)'
   ;print, filter+'-band zeropoint = ', fpr(zpt1,2.3), ' (+/- ', fpr(zpt1stdev,1.3), ')  - ', fpr(-airmassterm,1.3), '*(airmass-1)'

   fstr = (strsplit(filter,'_',/extract)) [0]
   dstr = (strsplit(filter,'_',/extract)) [1]

   printf, 1, fstr, ' ', dstr, ' ', fpr(zpt1,2.4), ' ', fpr(zpt1unc,1.4), ' ',  fpr(airmassterm,2.4), ' ', fpr(airmasstermunc,2.4)

endfor

close, 1

!p.multi = [0,1,2]
xp = 0
yp = 0
x1bars = [0.06,  0.56]
x2bars = [0.48, 0.98]
y1bars = [0.55, 0.05]
y2bars = [0.95, 0.45]
psopen, imworkingdir+'airmasscurve'+camera+'.ps', xsize=7, ysize=7, /inches
!p.font=0
device, /helvetica, font_index = 17
device, /color
colors = transpose([[0,0,0],$
                    [172,172,172], $
                    [255,255,255]])

tvlct, colors



for f = 0, n_elements(filts)-1  do begin
   !p.position = [x1bars[xp], y1bars[yp], x2bars[xp], y2bars[yp]]
   if f mod 4 ne 0 then !p.multi[0] = 1
   xp = xp + 1
   if xp eq 2 then begin
      xp = 0
      yp = yp + 1
   endif
   if yp ge 2 then yp = 0

   filtsobs = sobs[where(sobs.filter eq filts[f])]

   color = intarr(n_elements(filtsobs)) + 1
   color[where(filtsobs.good ge 0)] = 0 

  fyr = 0
  if strlowcase(strmid(filts[f],0,1)) eq 'u' then fyr = -1
  if strlowcase(strmid(filts[f],0,1)) eq 'B' then fyr = 0.5
  if strlowcase(strmid(filts[f],0,1)) eq 'z' then fyr = 1.0
  if camera eq 'red2' then fyr = 0.5

   megaplot, title='!17'+filts[f], xrange=[1,2.4], yrange=[1.5,3.5]+fyr, /xstyle, /ystyle
   megaplot, filtsobs.airmass, filtsobs.catmag-filtsobs.instmag, yerr=sqrt(filtsobs.caterr^2 + filtsobs.photoerr^2), fill=(filtsobs.good eq 1), color=color, size=0.5,/noerase
   air = [1,3]
   oplot, air, fzpt1[f]+(air-1)*fairmassterm[f]

endfor
psclose


;gobs = sobs[where(sobs.filter eq 'g')]
;
;!p.multi = [0,2,2]
;psopen, 'gairmasscurveind.ps'
;for i = 0, n_elements(indstars)-1 do begin
;   gobss = gobs[where(gobs.starname eq indstars[i])]
;   megaplot, gobss.airmass, gobss.catmag-gobss.instmag, yerr=sqrt(gobss.caterr^2 + gobss.photoerr^2), title=indstars[i], xrange=[1,2.4], yrange=[-3,0]
;endfor
;psclose


end


