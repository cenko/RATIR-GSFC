function gettmpsc, ra, dec, radius ;, catalog=catalog

   ; boxsize in arcsec, apparently

   catalog = 'tmpsc'

   rainp = clip(string(ra))
   decinp = clip(string(dec))
   radinp = clip(radius/60.) ; sdss direct query only

   if strpos(rainp,':') gt 0 then rainp = clip(15*ten(ra))
   if strpos(decinp,':') gt 0 then decinp = clip(ten(dec))

;   queryurl = '"http://cas.sdss.org/astro/en/tools/search/x_radial.asp?ra='+rainp+'&dec='+decinp+'&radius='+radinp+'&min_u=0&max_u=30&min_g=0&max_g=30&min_r=0&max_r=30&min_i=0&max_i=30&min_z=0&max_z=30&entries=top&topnum=1000&format=csv"'
   boxsize = -500                  ;arcsec
   queryurl = strcompress("http://tdc-www.harvard.edu/cgi-bin/scat?catalog=" + catalog +  "&ra=" + string(ra) + "&dec=" + string(dec) + "&system=J2000&rad=" + string(boxsize) + "&sort=mag&epoch=2000.00000&nstar=6400",/remove_all)
   usercat = 0
   racolumn = 1
   deccolumn = 2
   jmagcolumn = 3
   hmagcolumn = 4
   kmagcolumn = 5
   print,queryurl
;   cat = urllib.urlopen(queryurl)
;   catlines = cat.readlines()
   command = 'wget "' + queryurl + '" -O tempcat.dat -q'
   spawn, command
   n = countlines('tempcat.dat')-1 ;  -12
   if n le 0 then return, -1

;   readcol, 'tempcat.dat', id, run, rerun, camcol, field, obj, type,
;   ra, dec, u, g, r, i, z, u_unc, g_unc, r_unc, i_unc, z_unc,
;   format='ll,i,i,i,i,i,i,d,d,f,f,f,f,f,f,f,f,f,f' , skipline=1,
;   /silent   ; formerly skipline 12
   readcol,'tempcat.dat',format='ll,a,a,f,f,f,f',id,ra,dec,j,h,k,arcsec   ; formerly skipline 12
   n = n_elements(id)
   FOR i = 0,n_elements(id)-1 DO BEGIN
      rasplit=strsplit(ra(i),':',/extract)
      tmpra = 15.*(rasplit(0)+(rasplit(1)+(rasplit(2)/60.))/60.)
      decsplit = strsplit(dec(i),':',/extract)
      tmpdec = (decsplit(0)+(decsplit(1)+(decsplit(2)/60.))/60.)
      IF i EQ 0 THEN raarr = tmpra ELSE raarr=[raarr,tmpra]
      IF i EQ 0 THEN decarr = tmpdec ELSE decarr=[decarr,tmpdec]
   ENDFOR
   ra=raarr
   dec=decarr
   ;if n le 1 then return, -1

   tmpsc = replicate({id:0LL, sra:'',sdec:'',ra:0.D, dec:0.D, $
                     j:0.,h:0.,k:0.,Jj:0.,Hj:0.,Kj:0., $
                     j_unc:0.,h_unc:0.,k_unc:0.,Jj_unc:0.,Hj_unc:0.,Kj_unc:0., $
                     type:0},n)
   tmpsc.id = id
   tmpsc.j = j
   tmpsc.h = h
   tmpsc.k = k
;   sdss.j_unc = j_unc
;   sdss.h_unc = h_unc
;   sdss.k_unc = k_unc
;   sdss.type = type
   tmpsc.ra = ra
   tmpsc.dec = dec

;   sdss.Uj =  u - 0.0216*(u - g) - 0.799 ; combination of 2 equations, fortunately not needed
;
;   sdss.Bj = (u - 0.8116*(u - g) + 0.1313)*0+$;  sigma = 0.0095
;             (g + 0.3130*(g - r) + 0.2271)*1;    sigma = 0.0107
;
;   sdss.Vj = (g - 0.2906*(u - g) + 0.0885)*0+$;  sigma = 0.0129
;             (g - 0.5784*(g - r) - 0.0038)*1;    sigma = 0.0054
;
;   sdss.Rc = (r - 0.1837*(g - r) - 0.0971)/2+$;  sigma = 0.0106
;             (r - 0.2936*(r - i) - 0.1439)/2;    sigma = 0.0072
;
;   sdss.Ic = (r - 1.2444*(r - i) - 0.3820)*1+$;  sigma = 0.0078
;             (i - 0.3780*(i - z) - 0.3974)*0;    sigma = 0.0063
;
;   sdss.Bj_unc = sqrt((g_unc * 1.3130)^2 + (r_unc*0.2271)^2 + 0.0107^2)
;
;   sdss.Vj_unc = sqrt((g_unc*(1-0.5784))^2 + (r_unc*0.5784)^2 + 0.0054^2)
;
;   sdss.Rc_unc = sqrt(r_unc^2 + 0.0106^2)
;
;   sdss.Ic_unc = sqrt((r_unc*0.2444)^2 + (i_unc*1.2444)^2 + 0.0078^2)

   return, tmpsc

end
