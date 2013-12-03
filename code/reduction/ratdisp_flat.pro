pro ratdisp_flat,camera=camera

;dirn='/Volumes/Apps_and_Docs/ofox/data/ratir/20120427/1.7micron/'
;dirout = '/Volumes/Apps_and_Docs/ofox/data/ratir/20120427/1.7micron/out/'
device, decomp=0
loadct, 1

;get list of directories
target_dirs=find_all_dir('.')

FOR i = 1,n_elements(target_dirs)-1 DO BEGIN
   target=strsplit(target_dirs(i),'/',/extract)
   IF i EQ 1 THEN target_arr=target(1) ELSE target_arr=[target_arr,target(1)]
;   stop
ENDFOR

target_arr=target_arr(uniq(target_arr))

;yl=[100,1799] & xl=[100,799] & xs=xl(1)-xl(0)+1 & ys=yl(1)-yl(0)+1 
fc=3

window, xs=512, ys=512
direction='';initialize the variable
close,1
close,2
close,3
close,4

FOR i=1,1 DO BEGIN
gocam='C0'
IF keyword_set(camera) THEN gocam = strcompress(['C'+string(camera)],/remove_all)

IF strcmp(gocam,'C0') THEN BEGIN
;goto,C2
C0:
;CAMERA 0

   filelist0=strcompress('C0_'+target_arr(i)+'.list')
   readcol,filelist0,f='a',files0
   openw,1,strcompress(target_arr(i)+'/im_0.list')
stop
   FOR k = 0,n_elements(files0)-1 DO BEGIN
      erase
      filename=strsplit(files0(k),'.',/extract)
      filename=filename(0)
      ifile=strcompress(target_arr(i)+'/'+filename+'.fits',/remove_all)
      print,ifile
      fxread, ifile, im, h

      ;no rotation needed
      im=rotate(im,3)

      ;Let user determine if the image is good or not
      implot=rebin(im, 512, 512) 
      m=median(im) & s=ROBUST_SIGMA(im) ; mean and stdev of image 1
      tv, bytscl(implot, min=m-2*s,max=m+2*s) , 0, 0
      print,'Median is ',m,' counts'
      read, direction, prompt='Type Y for YES or N for NO:'
      IF strcmp(direction,'N',/fold_case) THEN BEGIN
         goto,SKIP0
      ENDIF ELSE IF strcmp(direction,'Y',/fold_case) THEN BEGIN
         printf,1,filename
      ENDIF
      SKIP0:
   ENDFOR

   close,1
stop
ENDIF 
IF strcmp(gocam,'C0') or strcmp(gocam,'C1') THEN BEGIN
C1:
;CAMERA 1

   filelist1=strcompress('C1_'+target_arr(i)+'.list')
   readcol,filelist1,f='a',files1
   openw,1,strcompress(target_arr(i)+'/im_1.list')
   FOR k = 0,n_elements(files1)-1 DO BEGIN
      erase
      filename=strsplit(files1(k),'.',/extract)
      filename=filename(0)
      ifile=strcompress(target_arr(i)+'/'+filename+'.fits',/remove_all)
      print,ifile
      fxread, ifile, im, h

      ;no rotation needed
      ;Let user determine if the image is good or not
      implot=rebin(im, 512, 512) 
      m=median(im) & s=ROBUST_SIGMA(im) ; mean and stdev of image 1
      print,'Median is ',m,' counts'
      tv, bytscl(implot, min=m-fc*s,max=m+fc*s) , 0, 0
      read, direction, prompt='Type Y for YES or N for NO:'
      IF strcmp(direction,'N',/fold_case) THEN BEGIN
         goto,SKIP1
      ENDIF ELSE IF strcmp(direction,'Y',/fold_case) THEN BEGIN
         printf,1,filename
      ENDIF
      SKIP1:
   ENDFOR

   close,1
stop
ENDIF
IF strcmp(gocam,'C0') or strcmp(gocam,'C2') THEN BEGIN
;CAMERA 2
C2:
   filelist2=strcompress('C2_'+target_arr(i)+'.list')
   readcol,filelist2,f='a',files2
   openw,1,strcompress(target_arr(i)+'/im_z.list')
   openw,2,strcompress(target_arr(i)+'/im_y.list')

   FOR k = 0,n_elements(files2)-1 DO BEGIN
      erase
      filename=strsplit(files2(k),'.',/extract)
      filename=filename(0)
      ifile=strcompress(target_arr(i)+'/'+filename+'.fits',/remove_all)
      print,ifile
      fxread, ifile, im, h
      asic=sxpar(h,'ASIC_NUM')

      ;no rotation needed
      im=rotate(im,3)

      ;Let user determine if the target is oriented on the E or W side of image
      implot=rebin(im, 512, 512) 
      m=median(im) & s=ROBUST_SIGMA(im) ; mean and stdev of image 1
      mright=median(im[1144:2043,1:1700])
      mleft=median(im[1:900,1:1700])
      print,'Median of Right Side is ',mright,' counts'
      print,'Median of Left Side is ',mleft,' counts'
      tv, bytscl(implot, min=m-fc*s,max=m+fc*s) , 0, 0
      ;Z & J face the EASTERN part of the field
      read, direction, prompt='Type Y for YES or N for NO or Q for Quit:'
      IF strcmp(direction,'Q',/fold_case) THEN BEGIN
         goto,GETOUT2
      ENDIF ELSE IF strcmp(direction,'N',/fold_case) THEN BEGIN
         goto,SKIP2
      ENDIF ELSE IF strcmp(direction,'Y',/fold_case) THEN BEGIN
         imfits=strcompress(target_arr(i)+'/'+filename+'_im_Y.fits',/remove_all)
;         sxaddpar,h,'NAXIS1','900'
         im_y = im[1144:2043,1:1700]
         sxaddpar, h, 'FILTER', 'Y'
         writefits, imfits, im_y, h
         printf,2,strcompress(filename+'_im_Y')

         imfits=strcompress(target_arr(i)+'/'+filename+'_im_Z.fits',/remove_all)
;         sxaddpar,h,'NAXIS1',900
         im_z = im[1:900,1:1700]
         sxaddpar, h, 'FILTER', 'Z'
         writefits, imfits, im_z, h
         printf,1,strcompress(filename+'_im_Z')
      ENDIF
      SKIP2:
   ENDFOR
   GETOUT2:

   close,1
   close,2
stop
ENDIF
IF strcmp(gocam,'C0') or strcmp(gocam,'C3') THEN BEGIN
;CAMERA 3
C3:

   filelist3=strcompress('C3_'+target_arr(i)+'.list')
   readcol,filelist3,f='a',files3
   openw,1,strcompress(target_arr(i)+'/im_j.list')
   openw,2,strcompress(target_arr(i)+'/im_h.list')
   FOR k = 0,n_elements(files3)-1 DO BEGIN
      erase
      filename=strsplit(files3(k),'.',/extract)
      filename=filename(0)
      ifile=strcompress(target_arr(i)+'/'+filename+'.fits',/remove_all)
      print,ifile
      fxread, ifile, im, h
      asic=sxpar(h,'ASIC_NUM')

      ;rotation needed (flip in X)
      im=rotate(im,4)

      ;Let user determine if the target is oriented on the E or W side of image
      implot=rebin(im, 512, 512) 
      m=median(im) & s=ROBUST_SIGMA(im) ; mean and stdev of image 1
      tv, bytscl(implot, min=m-10*s,max=m+15*s) , 0, 0

      mright=median(im[1144:2043,1:1700])
      mleft=median(im[1:900,1:1700])
      print,'Median of Right Side is ',mright,' counts'
      print,'Median of Left Side is ',mleft,' counts'

      ;Z & J face the EASTERN part of the field
      read, direction, prompt='Type Y for YES or N for NO or Q for Quit:'
      IF strcmp(direction,'Q',/fold_case) THEN BEGIN
         goto,GETOUT3
      ENDIF ELSE IF strcmp(direction,'N',/fold_case) THEN BEGIN
         goto,SKIP3
      ENDIF ELSE IF strcmp(direction,'Y',/fold_case) THEN BEGIN
         imfits=strcompress(target_arr(i)+'/'+filename+'_im_H.fits',/remove_all)
;         sxaddpar,h,'NAXIS1','900'
         im_h = im[1144:2043,1:1700]
         sxaddpar, h, 'FILTER', 'H'
         writefits, imfits, im_h, h
         printf,2,strcompress(filename+'_im_H')

         imfits=strcompress(target_arr(i)+'/'+filename+'_im_J.fits',/remove_all)
;         sxaddpar,h,'NAXIS1','900'
         sxaddpar, h, 'FILTER', 'J'
         im_j = im[1:900,1:1700]
         writefits, imfits, im_j, h
         printf,1,strcompress(filename+'_im_J')
      ENDIF
      SKIP3:
   ENDFOR
   GETOUT3:
   close,1
   close,2
stop
ENDIF
ENDFOR

;FOR k=0, n_elements(filenames)-1 DO BEGIN

;   ir12=stregex(asic, '2-32') GE 0 
;   if stregex(asic, '2-32') GE 0 then im=reverse((transpose(im)),2) ; orient them in N-E direction 
;   if stregex(asic, '1-29') GE 0 then im=reverse((transpose(im)),1)
;   im1=im(xl(0):xl(1),yl(0):yl(1))
;   im1=rebin(im1, xs/2, ys/2)
;   im2=im(xl(0)+1024:xl(1)+1024,yl(0):yl(1))
;   im2=rebin(im2, xs/2, ys/2)   
;   m1=median(im1) & s1=ROBUST_SIGMA(im1) ; mean and stdev of image 1
;   m2=median(im2) & s2=ROBUST_SIGMA(im2)
;   tv, bytscl(im1, min=m1-fc*s1,max=m1+fc*s1) , 50, 50
;   tv, bytscl(im2, min=m2-fc*s2,max=m2+fc*s2) , 100+xs/2, 50
;   xyouts, 20, 20, asic , col=0, /dev
;    img=tvrd(true=3)
;    impng=transpose([[[reform(img(*,*,0))]],[[reform(img(*,*,1))]],[[reform(img(*,*,2))]]],[2,0,1])
;    SKIP: 
;ENDFOR

;return
end
