function getsdss, ra, dec, radius ;, catalog=catalog

   ; boxsize in arcsec, apparently

   catalog = 'sdss'

   rainp = clip(string(ra))
   decinp = clip(string(dec))
   radinp = clip(radius/60.) ; sdss direct query only

   if strpos(rainp,':') gt 0 then rainp = clip(15*ten(ra))
   if strpos(decinp,':') gt 0 then decinp = clip(ten(dec))

   ;;queryurl = "http://tdc-www.harvard.edu/cgi-bin/scat\?catalog=" + catalog +  "\&ra=" + rainp + "\&dec=" + decinp + "\&system=J2000\&rad=" + clip(-boxsize) + "\&sort=mag\&epoch=2000.00000\&nstar=1600"

   queryurl = '"http://cas.sdss.org/astro/en/tools/search/x_radial.asp?ra='+rainp+'&dec='+decinp+'&radius='+radinp+'&min_u=0&max_u=30&min_g=0&max_g=30&min_r=0&max_r=30&min_i=0&max_i=30&min_z=0&max_z=30&entries=top&topnum=1000&format=csv"'

   ; &min_u=0&max_u=30&min_g=0&max_g=30&min_r=0&max_r=30&min_i=0&max_i=30&min_z=0&max_z=30

   command = 'wget ' + queryurl + ' -O tempcat.dat -q'
   spawn, command

   n = countlines('tempcat.dat')-1 ;  -12
   if n le 0 then return, -1

   ;readcol, 'tempcat.dat', id, sra, sdec, u, g, r, i, z, type, dist, format='l,a,a,f,f,f,f,f,i,f', skipline=12, /silent
   readcol, 'tempcat.dat', id, run, rerun, camcol, field, obj, type, ra, dec, u, g, r, i, z, u_unc, g_unc, r_unc, i_unc, z_unc, format='ll,i,i,i,i,i,i,d,d,f,f,f,f,f,f,f,f,f,f' , skipline=1, /silent   ; formerly skipline 12

   n = n_elements(id)
stop
   ;if n le 1 then return, -1

   sdss = replicate({id:0LL, sra:'',sdec:'',ra:0.D, dec:0.D, $
                     u:0.,g:0.,r:0.,i:0.,z:0.,                     Uj:0., Bj:0.,Vj:0.,Rc:0.,Ic:0., $
                     u_unc:0.,g_unc:0.,r_unc:0.,i_unc:0.,z_unc:0., Uj_unc:0.,Bj_unc:0.,Vj_unc:0.,Rc_unc:0.,Ic_unc:0., $
                     type:0},n)
   sdss.id = id
   ;sdss.sra = sra
   ;sdss.sdec = sdec
   sdss.u = u
   sdss.g = g
   sdss.r = r
   sdss.i = i
   sdss.z = z
   sdss.u_unc = u_unc
   sdss.g_unc = g_unc
   sdss.r_unc = r_unc
   sdss.i_unc = i_unc
   sdss.z_unc = z_unc
   sdss.type = type
   ;for i = 0, n-1 do begin
   ;  sdss[i].ra = 15*ten(sra[i])
   ;  sdss[i].dec = ten(sdec[i])
   ;endfor
   sdss.ra = ra
   sdss.dec = dec

   sdss.Uj =  u - 0.0216*(u - g) - 0.799 ; combination of 2 equations, fortunately not needed

   sdss.Bj = (u - 0.8116*(u - g) + 0.1313)*0+$;  sigma = 0.0095
             (g + 0.3130*(g - r) + 0.2271)*1;    sigma = 0.0107

   sdss.Vj = (g - 0.2906*(u - g) + 0.0885)*0+$;  sigma = 0.0129
             (g - 0.5784*(g - r) - 0.0038)*1;    sigma = 0.0054

   sdss.Rc = (r - 0.1837*(g - r) - 0.0971)/2+$;  sigma = 0.0106
             (r - 0.2936*(r - i) - 0.1439)/2;    sigma = 0.0072

   sdss.Ic = (r - 1.2444*(r - i) - 0.3820)*1+$;  sigma = 0.0078
             (i - 0.3780*(i - z) - 0.3974)*0;    sigma = 0.0063

   sdss.Bj_unc = sqrt((g_unc * 1.3130)^2 + (r_unc*0.2271)^2 + 0.0107^2)

   sdss.Vj_unc = sqrt((g_unc*(1-0.5784))^2 + (r_unc*0.5784)^2 + 0.0054^2)

   sdss.Rc_unc = sqrt(r_unc^2 + 0.0106^2)

   sdss.Ic_unc = sqrt((r_unc*0.2444)^2 + (i_unc*1.2444)^2 + 0.0078^2)

   return, sdss

end
