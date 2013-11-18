pro getcatalog, ra, dec, boxsize, catalog=catalog

   

   queryurl = "http://tdc-www.harvard.edu/cgi-bin/scat?catalog=" + catalog +  "&ra=" + str(ra) + "&dec=" + str(dec) + "&system=J2000&rad=" + clip(-boxsize) + "&sort=mag&epoch=2000.00000&nstar=1600"

       
   command = 'wget ' + queryurl + ' -O tempcat.dat'
   spawn, command
   

   rdfloat, 'tempcat.dat', ra, dec, u, g, r, i, z, u_unc, g_unc, r_unc, i_unc, z_unc
        cat = urllib.urlopen(queryurl)
        catlines = cat.readlines()
        cat.close()




end
