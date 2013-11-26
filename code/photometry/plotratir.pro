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

PRO plotratir

filters=['r','i','z','y','J','H']
plotra=dblarr(10000)
plotdec=dblarr(10000)
plotrmag=dblarr(10000)
plotrmagerr=dblarr(10000)
plotimag=dblarr(10000)
plotimagerr=dblarr(10000)
plotzmag=dblarr(10000)
plotzmagerr=dblarr(10000)
plotymag=dblarr(10000)
plotymagerr=dblarr(10000)
plotJmag=dblarr(10000)
plotJmagerr=dblarr(10000)
plotHmag=dblarr(10000)
plotHmagerr=dblarr(10000)
plotcat=dblarr(10000)
plotfilt=dblarr(10000)

;make a global list of detections
prefchar = 'coadd'
wildcharimg = '?????-????-?????_?'
zffiles = choosefiles(prefchar+wildcharimg+'.crop.multi.fits')
loadct,0
;find overlapping stars
FOR i = 0,n_elements(zffiles)-1 DO BEGIN
   ifile=zffiles(i)
   tmp=strsplit(ifile,'.',/extract)
   tmp=strsplit(zffiles(i),'_',/extract)
   filter=tmp(1)
   tmp=strsplit(filter,'.',/extract)
   filter=tmp(0)                                        ;define filter for working file
   IF strcmp(filter,'Z') THEN filter=strlowcase(filter) ;format filter for input into later code
   IF strcmp(filter,'Y') THEN filter=strlowcase(filter) ;
   
   ifile=strcompress(['finalphot'+filter+'.am'],/remove_all)
   readcol,ifile,f='i,f,f,d,d,f,f,f',id,x,y,ra,dec,mag,magerr

   IF i EQ 0 THEN BEGIN
      plotra(0:n_elements(ra)-1)=ra
      plotdec(0:n_elements(dec)-1)=dec
      cmd=strcompress(["plot"+filter+"mag(0:n_elements(mag)-1)=mag"],/remove_all)
      void = execute(cmd[0])
      cmd=strcompress(["plot"+filter+"magerr(0:n_elements(magerr)-1)=magerr"],/remove_all)
      void = execute(cmd[0])
;      plotcat(0:n_elements(cat)-1)=cat
      nstars=n_elements(ra)
   ENDIF ELSE BEGIN
      compra=ra
      compdec=dec
      compmag=mag
      compmagerr=magerr
;      compcat=cat

      count=0
      FOR k = 0,n_elements(compra)-1 DO BEGIN
         smatch = nearest(compra(k)*cos(compdec(k)*!pi/180.),compdec(k),plotra*cos(plotdec*!pi/180.),plotdec,mindist=1./3600.,count=ct)
         IF smatch(0) GE 0 THEN BEGIN
            cmd=strcompress(["plot"+filter+"mag(smatch(0))=compmag(k)"],/remove_all)
            void = execute(cmd[0])
            cmd=strcompress(["plot"+filter+"magerr(smatch(0))=compmagerr(k)"],/remove_all)
            void = execute(cmd[0])
;            plotcat(smatch(0))=1
            count=count+1
         ENDIF ELSE BEGIN

            plotra(nstars)=compra(k)
            plotdec(nstars)=compdec(k)
;            IF compcat(k) EQ 1 THEN plotcat(nstars) = 1
            cmd=strcompress(["plot"+filter+"mag(nstars)=compmag(k)"],/remove_all)
            void = execute(cmd[0])
            cmd=strcompress(["plot"+filter+"magerr(nstars)=compmagerr(k)"],/remove_all)
            void = execute(cmd[0])
            nstars=nstars+1
            count=count+1
         ENDELSE
      ENDFOR
      skiploop:
   ENDELSE
ENDFOR

close,1
openw,1,'finalmags.txt'
FOR i = 0,n_elements(plotra)-1 DO BEGIN
   printf,1,format='(f15.6,f15.6,f15.6,f15.6,f15.6,f15.6,f15.6,f15.6,f15.6,f15.6,f15.6,f15.6,f15.6,f15.6)',plotra(i),plotdec(i),plotrmag(i),plotrmagerr(i),plotimag(i),plotimagerr(i),plotzmag(i),plotzmagerr(i),plotymag(i),plotymagerr(i),plotJmag(i),plotJmagerr(i),plotHmag(i),plotHmagerr(i);,plotcat(i)
ENDFOR
close,1

set_plot,'x'

;Make plot of each image with circles on star identification
prefchar = 'coadd'
wildcharimg = '?????-????-?????_?'
zffiles = choosefiles(prefchar+wildcharimg+'.crop.multi.fits')
loadct,0
readcol,'finalmags.txt',f='d,d,f,f,f,f,f,f,f,f,f,f,f,f',ra,dec,rmag,rmagerr,imag,imagerr,zmag,zmagerr,ymag,ymagerr,Jmag,Jmagerr,Hmag,Hmagerr;,cat

realdetections=where(rmagerr GT 0 AND imagerr GT 0)

FOR i = 0,n_elements(zffiles)-1 DO BEGIN
;   ifile=zffiles(i)
   tmp=strsplit(zffiles(i),'.',/extract)
   ifile=strcompress([tmp(0)+'.crop.multi.fits'],/remove_all)
   fits_read,ifile,img,h
   IF i EQ 0 THEN refh=h
   hastrom,img,h,refh
   imSize=SIZE(img)
   scalefactor=1
   img=congrid(img,imSize[1]/scalefactor,imSize[2]/scalefactor)
   IF i EQ 0 THEN imgarr=img ELSE imgarr=[[[imgarr]],[[img]]]
   tmp=strsplit(zffiles(i),'_',/extract)
   filter=tmp(1)
   tmp=strsplit(filter,'.',/extract)
   filter=tmp(0)                                        ;define filter for working file
   IF strcmp(filter,'Z') THEN filter=strlowcase(filter) ;format filter for input into later code
   IF strcmp(filter,'Y') THEN filter=strlowcase(filter) ;
   IF i EQ 0 THEN filterarr=filter ELSE filterarr=[filterarr,filter]
ENDFOR

cr=where(strcmp(filterarr,'r')) 
ci=where(strcmp(filterarr,'i'))
cz=where(strcmp(filterarr,'z'))
cy=where(strcmp(filterarr,'y'))
cH=where(strcmp(filterarr,'H'))
cJ=where(strcmp(filterarr,'J'))

IF cJ GE 0 and cH GE 0 then r=imgarr[*,*,cJ]*0.5+imgarr[*,*,ch]*0.5
IF cJ GE 0 and cH LT 0 then r=imgarr[*,*,cJ]
IF cH GE 0 and cJ LT 0 then r=imgarr[*,*,cH]
IF cH LT 0 and cJ LT 0 THEN r=0.

IF cz GE 0 and cy GE 0 then g=imgarr[*,*,cz]*0.5+imgarr[*,*,cy]*0.5
IF cz GE 0 and cy LT 0 then g=imgarr[*,*,cz]
IF cy GE 0 and cz LT 0 then g=imgarr[*,*,cy]
IF cy LT 0 and cz LT 0 THEN g=0.

IF cr GE 0 and ci GE 0 then b=imgarr[*,*,cr]*0.5+imgarr[*,*,ci]*0.5
IF cr GE 0 and ci LT 0 then b=imgarr[*,*,cr]
IF ci GE 0 and cr LT 0 then b=imgarr[*,*,ci]
IF ci LT 0 and cr LT 0 THEN b=0.

red=r
green=g
blue=b

mi=min(blue)
ma=max(blue)*0.01+mi
IF n_elements(blue) GT 1 THEN blue=bytscl(blue,min=0,max=8,top=250)

mi=min(green)
ma=max(green)*0.005+mi
IF n_elements(green) GT 1 THEN green=bytscl(green,min=0,max=8,top=250)

mi=min(red)
ma=max(red)*0.004
IF n_elements(red) GT 1 THEN red=bytscl(red,min=0,max=8,top=250) 

IF n_elements(red) GT 1 THEN imSize=size(red) ELSE IF n_elements(green) GT 1 THEN imSize=size(green) ELSE IF n_elements(blue) GT 1 THEN imSize=size(blue)
WINDOW,0,XSIZE=imSize[1],YSIZE=imSize[2]

color=bytarr(3,imSize[1],imSize[2])
IF n_elements(red) GT 1 THEN color(0,*,*)=red*0.5
IF n_elements(green) GT 1 THEN color(1,*,*)=green*0.5
IF n_elements(blue) GT 1 THEN color(2,*,*)=blue*0.5

;[[[red]],[[green]],[[blue]]]
;ring=reform(ring)
device,decomposed=1
;loadct,2
loadct,1
tv,color,true=1
ofile='color.bmp'
write_bmp,ofile,color,color(0,*,*),color(1,*,*),color(2,*,*),/rgb
wdelete

loadct,0
FOR i = 0,n_elements(zffiles)-1 DO BEGIN
;   ifile=zffiles(i)
   tmp=strsplit(zffiles(i),'.',/extract)
   ifile=strcompress([tmp(0)+'.crop.multi.fits'],/remove_all)
   fits_read,ifile,img,h

   tmp=strsplit(zffiles(i),'_',/extract)
   filter=tmp(1)
   tmp=strsplit(filter,'.',/extract)
   filter=tmp(0)                                        ;define filter for working file
   IF strcmp(filter,'Z') THEN filter=strlowcase(filter) ;format filter for input into later code
   IF strcmp(filter,'Y') THEN filter=strlowcase(filter) ;

   device,decomposed=0
   imSize=SIZE(img)
   img=congrid(img,imSize[1]/scalefactor,imSize[2]/scalefactor)
   imSize=SIZE(img)
   device,retain=2
   WINDOW,0,XSIZE=imSize[1],YSIZE=imSize[2]
   scale=bytscl(img,min=0,max=10,top=255)
   tv,scale

   tmp=strsplit(ifile,'.',/extract)
;   ofile=strcompress([tmp(0)+'.am'],/remove_all)
;   readcol,ofile,f='f,f,d,d,f,f,f',x,y,ra,dec,mag,magerr,cat
   ADXY,h,ra,dec,x,y

   ;overplot the apertures used for extraction
   green=getcolor('green',!D.Table_Size-2)
   FOR j=0,n_elements(x)-1 DO BEGIN
      tmp=where(realdetections EQ j)
      IF tmp(0) GE 0 THEN BEGIN
         plots,circle(x(j)/scalefactor,y(j)/scalefactor,10),color=green,/device
         xyouts,x(j)/scalefactor,y(j)/scalefactor,string(j),color=green,/device,charsize=1.1,charthick=2,alignment=0.5
      ENDIF
   ENDFOR

   red=getcolor('red',!D.Table_Size-3)
   ostring=strcompress([filter+'-Band'],/remove_all)
   xyouts,0.2,0.9,ostring,/normal,charsize=5,charthick=5,color=red

   ;make .jpegs
   tmp=strsplit(ifile,'.',/extract)
   ofile=strcompress([tmp(0)],/remove_all)
   WSHOW, 0, ICONIC=0
   device,retain=2
   tmp=tvread(filename=ofile,/bmp,/nodialog)
   wdelete
ENDFOR

webpage:

;create an easy to view webpage
printratirhtml

END
