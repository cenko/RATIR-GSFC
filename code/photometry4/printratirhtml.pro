PRO printratirhtml

plotps,'photcomp.eps'
device,/portrait;,xsize=8.,ysize=8,/inches
!p.font=0
;multiplot,/reset
!p.multi=0

xtit='AB Magnitude'
ytit='!9D!X Mag'

purple=getcolor('purple',!D.Table_Size-2) ;rband
blue=getcolor('blue',!D.Table_Size-3)     ;iband
aqua=getcolor('aqua',!D.Table_Size-4)     ;Zband
green=getcolor('green',!D.Table_Size-5)   ;Yband
orange=getcolor('orange',!D.Table_Size-6) ;Jband
red=getcolor('red',!D.Table_Size-7)       ;Hband
black=getcolor('black',!D.Table_Size-8)
white=getcolor('white',!D.Table_Size-9)

readcol,'finalmags.txt',f='d,d,f,f,f,f,f,f,f,f,f,f,f,f',ra,dec,rmag,rmagerr,imag,imagerr,zmag,zmagerr,ymag,ymagerr,Jmag,Jmagerr,Hmag,Hmagerr;,cat

syms,1,1,1
tmp=sort(rmag)
plot,rmag(tmp),rmagerr(tmp),psym=8,thick=4,xthick=4,ythick=4,xrange=[15,22],yrange=[-0.01,0.2],ystyle=1,xtitle=xtit,ytitle=ytit,title='SN2000ch',charsize=1.2,/nodata;,color=black,background=white,/nodata
tmp=sort(rmag)
oplot,rmag(tmp),rmagerr(tmp),psym=8,color=purple
tmp=sort(imag)
oplot,imag(tmp),imagerr(tmp),psym=8,color=blue
tmp=sort(zmag)
oplot,zmag(tmp),zmagerr(tmp),psym=8,color=aqua
tmp=sort(ymag)
oplot,ymag(tmp),ymagerr(tmp),psym=8,color=green
tmp=sort(jmag)
oplot,jmag(tmp),Jmagerr(tmp),psym=8,color=orange
tmp=sort(hmag)
oplot,hmag(tmp),Hmagerr(tmp),psym=8,color=red

items=['r','i','z','y','J','H']
colors=[purple,blue,aqua,green,orange,red]
legend,items,colors=colors,psym=8,charsize=1.2

;ofile='photometry'
;WSHOW, 0, ICONIC=0
;device,retain=2
;tmp=tvread(filename=ofile,/bmp,/nodialog)
device,/close
set_plot,'x'

close,1
openw,1,'ratir.html'

printf,1,'<HTML><HEAD><TITLE>RATIR DATA</TITLE></HEAD>'
printf,1,'<BODY BGCOLOR="#FFFFFF" TEXT="#003300">'
;<FONT SIZE="+2" COLOR="#006600">RATIR Sources for RA, DEC: 145.7222083333 , 9.494944444444 </FONT><P>
;<IMG SRC="thumb_full_stack.fits.jpg"> <IMG SRC="thumb_dss.fits.jpg"><BR>
;<A HREF="index_nocircles.html"> Version of this page </A> with no
;sources marked in postage stamp images.<BR> <HR><P>
prefchar = 'coadd'
wildcharimg = '?????-????-?????_?'
zffiles = choosefiles(prefchar+wildcharimg+'.bmp')
ostring=strcompress(['<IMG SRC="./color.bmp">'])
printf,1,ostring
;stop
FOR i = 0,n_elements(zffiles)-1 DO BEGIN
   tmp=strsplit(zffiles(i),'.',/extract)
   bmpfile=strcompress([tmp(0)+'.bmp'],/remove_all)
   ostring=strcompress(['<IMG SRC="./'+bmpfile+'">'])
;   IF i EQ 0 OR i EQ 2 OR i EQ 4 THEN printf,1,'<br />'
   printf,1,ostring
ENDFOR
printf,1,'<BR><HR><FONT SIZE="+2" COLOR="#006600">AB System Photometry (sources within 1 arcmin):</FONT><BR>'
printf,1,'Notes: Non-zero magnitudes with uncertainty of zero are 3-sigma upper limits.  Sources with magnitudes of 0.0000 are unobserved.<BR>'
;printf,1,'<A HREF="phot.html">Full sourcelist.</A><BR> <PRE>'
printf,1,'# RA DEC phot,phot_err [riZYJH]'
readcol,'finalmags.txt',f='d,d,f,f,f,f,f,f,f,f,f,f,f,f',ra,dec,rmag,rmagerr,imag,imagerr,zmag,zmagerr,ymag,ymagerr,Jmag,Jmagerr,Hmag,Hmagerr;,cat
realdetections=where(rmagerr GT 0 AND imagerr GT 0)
printf,1,'<br />'
FOR i = 0,n_elements(ra)-1 DO BEGIN
   tmp=where(realdetections EQ i)
   IF tmp(0) GE 0 THEN BEGIN
      printf,1,format='(f15.6,f15.6,"&nbsp",f15.6,"&nbsp",f15.2,"&nbsp",f15.2,"&nbsp",f15.2,"&nbsp",f15.2,"&nbsp",f15.2,"&nbsp",f15.2,"&nbsp",f15.2,"&nbsp",f15.2,"&nbsp",f15.2,"&nbsp",f15.2,"&nbsp",f15.2,"&nbsp",f15.2)',i,ra(i),dec(i),rmag(i),rmagerr(i),imag(i),imagerr(i),zmag(i),zmagerr(i),ymag(i),ymagerr(i),Jmag(i),Jmagerr(i),Hmag(i),Hmagerr(i);,cat(i)
      printf,1,'<br />'
   ENDIF
ENDFOR
printf,1,'</PRE><BR><HR>'
printf,1,'<IMG SRC="photcomp.jpg">'
printf,1,'</PRE><BR><HR>'
;<FONT SIZE="+2" COLOR="#006600">Exposure Info:</FONT><BR> <PRE>
;Camera C0 (filter=r) exposure information:
;  Timespan:  
;  A-side exposure: 
;  B-side exposure: 
;Camera C1 (filter=i) exposure information:
;  Timespan: 20121211T093420 20121211T104356
;  A-side exposure: 1440
;  B-side exposure: 1440
;Camera C2 (filter=ZY) exposure information:
;  Timespan:  
;  A-side exposure: 
;  B-side exposure: 
;Camera C3 (filter=JH) exposure information:
;  Timespan:  
;  A-side exposure: 
;  B-side exposure: 
;</PRE><HR><P>
;<FONT SIZE="+2" COLOR="#006600">Photometric Calibration Info:</FONT><BR>
;<PRE>

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
   ifile=strcompress(['irafoffset'+filter],/remove_all)
   readcol,ifile,ABoffset,err

;   ostring=['Photometry for '+strcompress(filter)+'-band:']
;   printf,1,ostring
;   ostring=['Magnitude Offset '+strcompress(string(ABoffset))'+/-'strcompress(string(err))]
;   printf,1,ostring
ENDFOR
;Photometry for C0 r band:
;Photometry for C1 i band:
; Matches 287
; Median Sextractor FWHM 10.0500
; Position Offset (arcsec): 0.31, 0.21 (0.38)
; Magnitude Offset 4.8764+/-0.0009
; 10-sigma limiting magnitude 21.66 (20.40--22.14)
;Photometry for C2 Z band:
; Matches 59
; Median Sextractor FWHM 11.2200
; Position Offset (arcsec): -0.30, -0.35 (0.46)
; Magnitude Offset 2.9577+/-0.0131
; 10-sigma limiting magnitude 20.61 (19.32--21.07)
;Photometry for C2 Y band:
; Matches 52
; Median Sextractor FWHM 10.6200
; Position Offset (arcsec): -0.65, -0.67 (0.94)
; Magnitude Offset 3.2189+/-0.0097
; 10-sigma limiting magnitude 20.22 (19.19--20.65)
;Photometry for C3 J band:
; Matches 61
; Median Sextractor FWHM 9.9650
; Position Offset (arcsec): -0.66, 0.16 (0.68)
; Magnitude Offset 3.4502+/-0.0055
; 10-sigma limiting magnitude 19.86 (19.02--20.23)
;Photometry for C3 H band:
; Matches 55
; Median Sextractor FWHM 9.3200
; Position Offset (arcsec): 0.10, -0.04 (0.11)
; Magnitude Offset 3.8150+/-0.0080
; 10-sigma limiting magnitude 19.20 (18.58--19.66)
;</PRE><HR><P>
printf,1,'<FONT SIZE="+2" COLOR="#006600">Fits thumbnails:</FONT><BR>'
;<A HREF="thumb_full_stack.fits.bz2">thumb_full_stack.fits.bz2</A><BR>
;<A HREF="thumb_dss.fits.bz2">thumb_dss.fits.bz2</A><BR>
FOR i = 0,n_elements(zffiles)-1 DO BEGIN
   tmp=strsplit(zffiles(i),'.',/extract)
   fitsfile=strcompress([tmp(0)+'.crop.fits'],/remove_all)
   ostring=strcompress(['<A HREF="'+fitsfile+'">'+fitsfile+'</A><BR>'])
   printf,1,ostring
ENDFOR
;printf,1,'<P><FONT SIZE="+2" COLOR="#006600">Full-frame Images:</FONT><BR>'
;printf,1,'<A HREF=""></A><BR>'
;<A HREF="stack_20121211T093420_C1_i.jpg">stack_20121211T093420_C1_i.jpg</A><BR>
;<A HREF="stack_20121211T093420_C2_ZY.jpg">stack_20121211T093420_C2_ZY.jpg</A>
;<A HREF="stackA_20121211T093420_C2_ZY.jpg">stackA_20121211T093420_C2_ZY.jpg</A>
;<A HREF="stackB_20121211T093420_C2_ZY.jpg">stackB_20121211T093420_C2_ZY.jpg</A><BR>
;<A HREF="stack_20121211T093420_C3_JH.jpg">stack_20121211T093420_C3_JH.jpg</A>
;<A HREF="stackA_20121211T093420_C3_JH.jpg">stackA_20121211T093420_C3_JH.jpg</A>
;<A HREF="stackB_20121211T093420_C3_JH.jpg">stackB_20121211T093420_C3_JH.jpg</A><BR>
;<HR WIDTH="100%"> Last Updated: Tue Jan 15 04:40:10 UTC 2013 <P>
;<ADDRESS> Nat Butler (natbutler@asu.edu)</ADDRESS> </BODY></HTML>
close,1
END
