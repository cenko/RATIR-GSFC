PRO prepfiles

;Make imagelist files
prefchar = 'coadd'
wildcharimg = '?????-????-?????_?'
zffiles = choosefiles(prefchar+wildcharimg+'.crop.fits')

FOR i = 0,n_elements(zffiles)-1 DO BEGIN
   tmp=strsplit(zffiles(i),'.',/extract)
   reffile=strcompress([tmp(0)+'.crop'],/remove_all)
;   img=readfits(reffile,h)
   tmp=strsplit(zffiles(i),'_',/extract)
   filter=tmp(1)
   tmp=strsplit(filter,'.',/extract)
   filter=tmp(0)  
   ofile=strcompress(['imagelist'+filter],/remove_all)
print,ofile
   close,1
   openw,1,ofile
   printf,1,reffile
   close,1
ENDFOR

FOR i = 0,n_elements(zffiles)-1 DO BEGIN
   tmp=strsplit(zffiles(i),'.',/extract)
   reffile=strcompress([tmp(0)+'.crop.multi'],/remove_all)
;   img=readfits(reffile,h)
   tmp=strsplit(zffiles(i),'_',/extract)
   filter=tmp(1)
   tmp=strsplit(filter,'.',/extract)
   filter=tmp(0)  
   ofile=strcompress(['imagelist_multi'+filter],/remove_all)
print,ofile
   close,1
   openw,1,ofile
   printf,1,reffile
   close,1
ENDFOR

END
