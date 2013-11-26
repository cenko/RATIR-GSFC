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

PRO calcoff
filters=['r','i','z','y','J','H']
;ABoffsetarr=[0.146,0.366,0.528,0.634,0.938,1.379]
ABoffsetarr=[0.,0.,0.,0.,0.938,1.379] ;Isaac's code yields rizy in AB and JH in Vega (2mass)

wildcharimg = 'fluxes1_?.txt'
zffiles = choosefiles(wildcharimg)
print,zffiles

FOR i = 0,n_elements(zffiles)-1 DO BEGIN

   tmp=strsplit(zffiles(i),'_',/extract)
   filter=tmp(1)
   tmp=strsplit(filter,'.',/extract)
   filter=tmp(0)  
   IF strcmp(filter,'Z') THEN filter=strlowcase(filter) ;format filter for input into later code
   IF strcmp(filter,'Y') THEN filter=strlowcase(filter) ;

   readcol,zffiles(i),format='f,f,d,d,f,f,f,d',x,y,ra,dec,mag_aper,magerr_aper,e,fwhm;,flux_aper,fluxerr_aper,mag_auto,magerr_auto,flux_auto,fluxerr_auto
   IF strcmp(filter,'J') OR strcmp(filter,'H') THEN $
      tmp=where(mag_aper GT 10 and mag_aper LT 20 and magerr_aper LT 0.4 AND FWHM LT 9) $   
   else tmp=where(mag_aper GT 10 and mag_aper LT 20 and magerr_aper LT 0.4 AND FWHM LT 9)

   x=x(tmp)
   y=y(tmp)
   ra=ra(tmp)
   dec=dec(tmp)
   e=e(tmp)
   fwhm=fwhm(tmp)
   aper_mag=mag_aper(tmp)
   aper_err=magerr_aper(tmp)
;   flux=flux_auto(tmp)
;   fluxerr=fluxerr_auto(tmp)

   n=n_elements(aper_mag)
;   aper_mag=aper_mag(0:n/2.)
;   aper_err=aper_err(0:n/2)
;   x=x(0:n/2)
;   y=y(0:n/2)
;   ra=ra(0:n/2)
;   dec=dec(0:n/2)

   filelist=strcompress('imagelist'+filter,/remove_all)
   readcol,filelist,f='a',filename

   ;convert all x,y positions to ra and dec to calculate zero-point
   fits_read,strcompress([filename[0]+'.fits'],/remove_all),img,h
   x=float(x)
   y=float(y)
   aper_mag=float(aper_mag)
   xyad,h,x,y,ra_out,dec_out

   ;CALCULATE ZERO POINTS
   imfile=strcompress([filename(0)+'.im'],/remove_all)
   close,1
   openw,1,imfile
   printf,1,'#RA ',' DEC ',' INST_MAG'
   FOR j = 0,n_elements(ra_out)-1 DO BEGIN
      printf,1,format='(f15.6,f15.6,f15.6)',ra_out(j),dec_out(j),aper_mag(j)
   ENDFOR
   close,1

   ra=ra_out
   dec=dec_out
   mag=aper_mag
   magerr=aper_err

   ;run Isaac Shivers routines to obtain catalogued magnitudes
   catfile=strcompress([filename(0)+'.cat'],/remove_all)
   amfile=strcompress([filename(0)+'.am'],/remove_all)

   cmd=['python zeropoint.py '+imfile+' '+filter+' '+catfile]
   print, cmd
   spawn, cmd

   ;read catalog magnitude results
;   readcol,catfile,f='d,d,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,i',refra,refdec,u_mag,g_mag,r_mag,i_mag,z_mag,y_mag,bigB_mag,bigV_mag,bigR_mag,bigI_mag,J_mag,H_mag,K_mag,u_err,g_err,r_err,i_err,z_err,y_err,bigB_err,bigV_err,bigR_err,bigI_err,J_err,H_err,K_err,mode
;   readcol,'../../ptf13abc_calib_sdss.txt',f='f,f,f,f,f,f,f,f,f,f,f,f,f',starra,stardec,inst_mag,search_id,matched_id,refra,refdec,g_mag,r_mag,i_mag,g_err,r_err,i_err
;   IF strcmp(filter,'r') OR strcmp(filter,'i') OR strcmp(filter,'Z') OR strcmp(filter,'Y') OR strcmp(filter,'J') OR strcmp(filter,'H') THEN BEGIN
;      print,filter,' for ptf13s_calib_star.lst'
;      readcol,'../../ptf13s_calib_star.lst',f='f,f,f,f,f,f,f,f,f,f,f,f',refra,refdec,r_mag,r_err,i_mag,i_err,z_mag,z_err,J_mag,J_err,H_mag,H_err
;   ENDIF ELSE BEGIN
      readcol,catfile,f='d,d,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,i',refra,refdec,u_mag,g_mag,r_mag,i_mag,z_mag,y_mag,bigB_mag,bigV_mag,bigR_mag,bigI_mag,J_mag,H_mag,K_mag,u_err,g_err,r_err,i_err,z_err,y_err,bigB_err,bigV_err,bigR_err,bigI_err,J_err,H_err,K_err,mode
;   ENDELSE


   ;find stars in catalogues
;   IF strcmp(filter,'r') OR strcmp(filter,'i') OR strcmp(filter,'Z') OR strcmp(filter,'Y') OR strcmp(filter,'J') OR strcmp(filter,'H') THEN BEGIN
;      tmp=where(r_mag GT 0)
;   ENDIF ELSE BEGIN
      tmp=where(mode EQ 0)
;   ENDELSE
   refra=refra(tmp)
   refdec=refdec(tmp)
   cmd=strcompress(['refmag = '+filter+'_mag(tmp)'],/remove_all)
   void = execute(cmd[0])
   cmd=strcompress(['referr = '+filter+'_err(tmp)'],/remove_all)
   void = execute(cmd[0])

   ;Calculate zero-point offset
   count=0
   FOR j = 0,n_elements(ra)-1 DO BEGIN ;for all stars in sextractor catalog...match to corresponding catalog stars
      smatch = nearest(ra(j)*cos(dec(j)*!pi/180.),dec(j),refra*cos(refdec*!pi/180.),refdec,mindist=1./3600.,count=ct)
print,smatch
      IF smatch(0) GE 0 THEN BEGIN
         IF refmag(smatch) LE 0 THEN continue
         dif=mag(j)-refmag(smatch)
         err=sqrt(magerr(j)^2.+referr(smatch)^2.)
         IF count EQ 0 THEN difarr=dif ELSE difarr=[difarr,dif]
         IF count EQ 0 THEN errarr=err ELSE errarr=[errarr,err]
         IF count EQ 0 THEN catarr=1 ELSE catarr=[catarr,1]
         count=count+1
      ENDIF
      IF smatch(0) EQ -1 THEN BEGIN
         IF count EQ 0 THEN catarr=0 ELSE catarr=[catarr,0]
;         count=count+1
      ENDIF
   ENDFOR

   ;take the 3-sigma-clip median
   djs_iterstat,difarr,median=difmag,sigma=sig
   errarr=difarr-difmag
   djs_iterstat,errarr,median=err,sigma=sigerr
   err=sig/sqrt(N)
print,filter,difarr
print,filter,err
stop
;   IF n_elements(difarr) GT 1 THEN difmag=median(difarr) ELSE difmag = difarr[0]
   print,'Filter ',filter,' offset:', difmag, ' mag ',' error:', sig/sqrt(N)
   ;apply the offset and the VEGA to AB conversion
   ofile=strcompress(['irafoffset'+filter],/remove_all)
   close,1
   openw,1,ofile
   printf,1,format='(f15.6,f15.6)',difmag,err
   close,1
ENDFOR

END
