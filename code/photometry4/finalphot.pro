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

PRO finalphot

filters=['r','i','z','y','J','H']

wildcharimg = 'fluxes2_?.txt'
zffiles = choosefiles(wildcharimg)
print,zffiles
FOR i = 0,n_elements(zffiles)-1 DO BEGIN

   filter=filters(i)
   tmp=strsplit(zffiles(i),'_',/extract)
   filter=tmp(1)
   tmp=strsplit(filter,'.',/extract)
   filter=tmp(0)  
   IF strcmp(filter,'Z') THEN filter=strlowcase(filter) ;format filter for input into later code
   IF strcmp(filter,'Y') THEN filter=strlowcase(filter) ;

   amfile=strcompress(['finalphot'+filter+'.am'],/remove_all)
   ifile=strcompress(['irafoffset'+filter],/remove_all)
   readcol,ifile,ABoffset,err

;   readcol,'coords',ra,dec

   filelist=strcompress('imagelist_multi'+filter,/remove_all)
   readcol,filelist,f='a',filename

   ;APERTURE INSTRUMENT MAGS
;   ifile=strcompress('sn_'+filename[0]+'.sky.aper.als',/remove_all)
;   readcol,ifile,f='a,f,i,a,a',fname,y,aperid,coords,lid,comment='#'
;   index=where(strcmp(coords,'coords')) 
;   x=fname(index+1)
;   y=y(index+1)
;   fname=fname(index)
;   aperid=aperid(index)
;   readcol,ifile,f='f,f,a,a,a,a,a,a',aper,sum,area,flux,aper_mag,aper_err,pier,perror,comment='#'
;   index=(indgen(n_elements(aper)/3.)+1.)*3.-1
;   aper_mag=aper_mag(index)
;   aper_err=aper_err(index)
;   tmp=where(strcmp(aper_mag,'INDEF') EQ 1)
;   aper_mag(tmp)=0
;   aper_err(tmp)=0
;   tmp = sort(aperid)
;   aperid = aperid(tmp)
;   aper_mag = aper_mag(tmp)
;   aper_err = aper_err(tmp)
;   x=x(tmp)
;   y=y(tmp)
;   ra=ra(tmp)
;   dec=dec(tmp)
;   aperid=aperid-1              ;subtract 1 off the index to comply
;   with IDL array indexing

   readcol,zffiles(i),format='f,f,d,d,f,f,f,d',x,y,ra,dec,mag,magerr,e,fwhm
   aperid=findgen(n_elements(mag))+1
   aper_mag=mag
   aper_err=sqrt(magerr^2.+err[0]^2.)
;   tmp=sort(mag)
;   aper_mag=mag(tmp)
;   aper_err=magerr(tmp)
;   id=aperid(tmp)
;   x=x(tmp)
;   y=y(tmp)
;   ra=ra(tmp)
;   dec=dec(tmp)
;   e=e(tmp)
;   fwhm=fwhm(tmp)

   ;apply previously calculated ABoffset
   aper_mag = aper_mag-ABoffset(0)
   close,1
   openw,1,amfile               ;amfile defined above (am = absolute magnitudes)
   printf,1,'#ID ','X ','Y ','RA ',' DEC ',' CAL_MAG ',' CAL_MAG_ERR'
   FOR k=0,n_elements(ra)-1 DO BEGIN
      printf,1,format='(f15.6,f15.6,f15.6,f15.6,f15.6,f15.6,f15.6)',aperid(k),x(k),y(k),ra(k),dec(k),aper_mag(k),aper_err(k)
   ENDFOR
   close,1 
;stop
ENDFOR

END
