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

;******

function circle, xcenter, ycenter, radius
   points = (2 * !PI / 99.0) * FINDGEN(100)
   x = xcenter + radius * COS(points )
   y = ycenter + radius * SIN(points )
   return, TRANSPOSE([[x],[y]])
end

;******

PRO icoords

;Define a number of global variables
filters=['r','i','z','y','J','H']

;Identify files
prefchar = 'coadd'
wildcharimg = '?????-????-?????_?'
zffiles = choosefiles(prefchar+wildcharimg+'.fits')
weightfiles = choosefiles(prefchar+wildcharimg+'.weight.fits')

;perform initial crop to remove noisy edges
FOR i =0,n_elements(zffiles)-1 DO BEGIN
   img=readfits(zffiles(i),h)
   tmp=strsplit(zffiles(i),'_',/extract)
   filter=tmp(1)
   tmp=strsplit(filter,'.',/extract)
   filter=tmp(0)  

;   IF strcmp(filter,'i') NE 1 THEN BEGIN
      s=size(img)
      x1=0.1*s(1) 
      x2=0.9*s(1)
      y1=0.1*s(2)
      y2=0.9*s(2)
;   ENDIF
;   IF strcmp(filter,'H') THEN BEGIN
;      x1=550.
;      x2=1250.
;      y1=1500.
;      y2=3100.
;   ENDIF
;   IF strcmp(filter,'i') THEN BEGIN
;      x1=200.
;      x2=1700.
;      y1=380.
;      y2=1250.
;   ENDIF
   IF strcmp(filter,'r') THEN BEGIN
      x1=105.
      x2=1580.
      y1=500.
      y2=1500.
   ENDIF
;   IF strcmp(filter,'J') THEN BEGIN
;      x1=700.
;      x2=1500.
;      y1=200.
;      y2=1800.
;   ENDIF
;   IF strcmp(filter,'Y') THEN BEGIN
;      x1=170.
;      x2=900.
;      y1=200.
;      y2=1800.
;   ENDIF
;   IF strcmp(filter,'Z') THEN BEGIN
;      x1=125.
;      x2=975.
;      y1=150.
;      y2=1950.
;   ENDIF

   hextract,img,h,x1,x2,y1,y2
   img=img
   tmp=strsplit(zffiles(i),'.',/extract)
   ofile=strcompress([tmp(0)+'.crop.fits'],/remove_all)
   writefits,ofile,img,h

   wimg=readfits(weightfiles(i),wh)
   hextract,wimg,wh,x1,x2,y1,y2
   tmp=strsplit(weightfiles(i),'.',/extract)
   ofile=strcompress([tmp(0)+'.crop.weight.fits'],/remove_all)
   writefits,ofile,wimg,wh
ENDFOR

;Resample all images using SWarp to a reference image called multicolor
FOR i =0,n_elements(zffiles)-1 DO BEGIN
   tmp=strsplit(zffiles(i),'.',/extract)
   ifile=strcompress([tmp(0)+'.crop.fits'],/remove_all)
   IF i EQ 0 THEN swarpstring=['"'+ifile+'"'] ELSE swarpstring=[swarpstring+' "'+ifile+'" ']
ENDFOR
stackcmd = ['swarp '+swarpstring+' -DELETE_TMPFILES N']
stackcmd = stackcmd + ' -IMAGEOUT_NAME multicolor.fits -WEIGHTOUT_NAME multicolor.weight.fits'
spawn,stackcmd

;Rename all the resampled files
FOR i =0,n_elements(zffiles)-1 DO BEGIN
   tmp=strsplit(zffiles(i),'.',/extract)
   ifile=strcompress([tmp(0)+'.crop.resamp.fits'],/remove_all)
   ofile=strcompress([tmp(0)+'.crop.fits'],/remove_all)
   mvcmd = ['mv -f '+ifile+' '+ofile]
   spawn, mvcmd
   ifile=strcompress([tmp(0)+'.crop.resamp.weight.fits'],/remove_all)
   ofile=strcompress([tmp(0)+'.crop.weight.fits'],/remove_all)
   mvcmd = ['mv -f '+ifile+' '+ofile]
   spawn, mvcmd
ENDFOR

;run sextractor on pipeline reduced files to identify point sources
FOR i = 0,n_elements(zffiles)-1 DO BEGIN
   tmp=strsplit(zffiles(i),'_',/extract)
   filter=tmp(1)
   tmp=strsplit(filter,'.',/extract)
   filter=tmp(0)                                        ;define filter for working file
   IF strcmp(filter,'Z') THEN filter=strlowcase(filter) ;format filter for input into later code
   IF strcmp(filter,'Y') THEN filter=strlowcase(filter) ;

   ;run SExtractor on each image
   tmp=strsplit(zffiles(i),'.',/extract)
   ifile=strcompress([tmp(0)+'.crop.fits'],/remove_all)
   cmd=['sex '+ifile+' -c "ratir_nir.sex"'] 
   print, cmd
   spawn, cmd

   ;read in results from sextractor and produce IRAF coordinate file
   readcol,'temp.cat',format='f,f,d,d,f,f,f,f',x,y,ra,dec,mag,magerr,e,fwhm
   close,1
   openw,1,strcompress(['coords'+filter],/remove_all)
   FOR j = 0,n_elements(ra)-1 DO BEGIN
      printf,1,format='(f15.6,f15.6)',ra(j),dec(j)
   ENDFOR
   close,1

   cmd=['mv -f temp.cat fluxes1_'+filter+'.txt']
   print, cmd
   spawn, cmd

   fwhm=median(fwhm)
   close,1
   openw,1,strcompress([filter+".aannulus"],/remove_all)
   printf,1,f='(f15.6)',fwhm
   close,1

ENDFOR

END
