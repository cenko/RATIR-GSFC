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

PRO assoc

filters=['r','i','z','y','J','H']

skipastrometry:
prefchar = 'coadd'
wildcharimg = '?????-????-?????_?'
zffiles = choosefiles(prefchar+wildcharimg+'.crop.fits')

;RE-crop images to produce overlapping images
FOR i =0,n_elements(zffiles)-1 DO BEGIN
   tmp=strsplit(zffiles(i),'.',/extract)
   reffile=strcompress([tmp(0)+'.crop.fits'],/remove_all)
   img=readfits(reffile,h)
   tmp=strsplit(zffiles(i),'_',/extract)
   filter=tmp(1)
   tmp=strsplit(filter,'.',/extract)
   filter=tmp(0)  
   imSize=size(img)
   xyad,h,0,0,ra1,dec1
   xyad,h,imSize(1)-1,imSize(2)-1,ra2,dec2
   IF i EQ 0 THEN ra1arr=ra1 ELSE ra1arr=[ra1arr,ra1]
   IF i EQ 0 THEN dec1arr=dec1 ELSE dec1arr=[dec1arr,dec1]
   IF i EQ 0 THEN ra2arr=ra2 ELSE ra2arr=[ra2arr,ra2]
   IF i EQ 0 THEN dec2arr=dec2 ELSE dec2arr=[dec2arr,dec2]
ENDFOR

raleft=min(ra1arr)
raright=max(ra2arr)
decbot=max(dec1arr)
dectop=min(dec2arr)

FOR i =0,n_elements(zffiles)-1 DO BEGIN
   tmp=strsplit(zffiles(i),'.',/extract)
   reffile=strcompress([tmp(0)+'.crop.fits'],/remove_all)
   newfile=strcompress([tmp(0)+'.crop.multi.fits'],/remove_all)
   img=readfits(reffile,h)
   adxy,h,raleft,decbot,x1,y1
   adxy,h,raright,dectop,x2,y2
   print,x1+1,x2,y1+1,y2
   hextract,img,h,x1+1,x2,y1+1,y2
   writefits,newfile,img,h
ENDFOR

;Crop master image
reffile='multicolor.fits'
img=readfits(reffile,h)
adxy,h,raleft,decbot,x1,y1
adxy,h,raright,dectop,x2,y2
print,x1+1,x2,y1+1,y2
hextract,img,h,x1+1,x2,y1+1,y2
writefits,reffile,img,h

;run sextractor on pipeline reduced files 
;create a unique photometry list defined by multicolor.fits

reffile='multicolor.fits'

FOR i = 0,n_elements(zffiles)-1 DO BEGIN
   tmp=strsplit(zffiles(i),'.',/extract)
   compfile=strcompress([tmp(0)+'.crop.multi.fits'],/remove_all)
   amfile=strcompress([tmp(0)+'.am'],/remove_all)

   tmp=strsplit(zffiles(i),'_',/extract)
   compfilter=tmp(1)
   tmp=strsplit(compfilter,'.',/extract)
   compfilter=tmp(0)            ;define filter for working file
    
   IF strcmp(compfilter,'Z') THEN compfilter=strlowcase(compfilter) ;format filter for input into later code
   IF strcmp(compfilter,'Y') THEN compfilter=strlowcase(compfilter) ;

   cmd=['sex ' + reffile + ',' + compfile + ' -c "ratir_nir.sex"'] 
   print, cmd
   spawn, cmd

   readcol,'temp.cat',format='f,f,d,d,f,f,f,f',x,y,ra,dec,mag,magerr,e,fwhm

   cmd=['mv -f temp.cat fluxes2_'+compfilter+'.txt']
   print, cmd
   spawn, cmd
ENDFOR

;;Print ultimate photometry list back to a coords file to be run in IRAF
;close,1
;openw,1,'coords'
;FOR i = 0,n_elements(ra)-1 DO BEGIN
;   printf,1,format='(f15.6,f15.6)',ra(i),dec(i)
;ENDFOR
;close,1

END
