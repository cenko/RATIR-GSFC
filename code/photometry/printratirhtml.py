"""
NAME:
	printratirhtml
PURPOSE:
	Create photcomp.png comparing magnitude and errors as well as create HTML page that has
	all of the data easily displayed
OUTPUT:
	photcomp.png - shows magnitude vs. error for each filter
	ratir.html   - html page showing information about sources

Translated from printratirhtml.pro by John Capone (jicapone@astro.umd.edu).
Modified 12/10/2013 by Vicki Toy (vtoy@astro.umd.edu)
"""

import numpy as np
import pylab as pl
import photprocesslibrary as pplib

def printratirhtml():
	
	#Reads in final magnitudes
	plotra, plotdec, plotrmag, plotrmagerr, plotimag, plotimagerr, plotzmag, plotzmagerr, plotymag, plotymagerr, plotJmag, plotJmagerr, plotHmag, plotHmagerr = np.loadtxt('./finalmags.txt', unpack=True)

	pl.figure()
	pl.xlim([15,22])
	pl.ylim([-0.01,0.2])
	pl.plot(plotrmag, plotrmagerr, marker='o',linestyle='None', label='r', color='purple')
	pl.plot(plotimag, plotimagerr, marker='o',linestyle='None', label='i', color='blue')
	pl.plot(plotzmag, plotzmagerr, marker='o',linestyle='None', label='z', color='aqua')
	pl.plot(plotymag, plotymagerr, marker='o',linestyle='None', label='y', color='green')
	pl.plot(plotJmag, plotJmagerr, marker='o',linestyle='None', label='J', color='orange')
	pl.plot(plotHmag, plotHmagerr, marker='o',linestyle='None', label='H', color='red')
	
	pl.xlabel('AB Magnitude')
	pl.ylabel(r"$\Delta$ Mag")
	pl.legend(loc='lower right')
	pl.savefig( 'photcomp.png', bbox_inches='tight')
	pl.clf()

	f = open( './ratir.html', 'w' )
	f.write( '<!DOCTYPE HTML>\n' )
	f.write( '<HTML>\n' )
	f.write( '<HEAD>\n' )
	f.write( '<TITLE>RATIR DATA</TITLE>\n' )
	f.write( '</HEAD>\n' )
	f.write( '<BODY BGCOLOR="#FFFFFF" TEXT="#003300">\n' )

	#Finds filter image files with green circles over sources and writes to HTML
	prefchar = 'coadd'
	zffiles = pplib.choosefiles( prefchar + '*_?.png' )

	im_wid = 400

	f.write( '<IMG SRC="./color.png" width="' + `im_wid` + '"><BR>\n' )
	for i in range(len(zffiles)):
	    f.write( '<IMG SRC="./' + zffiles[i] + '" width="' + `im_wid` + '">\n' )
	    if (i+1)%3 == 0:
	    	f.write( '<BR>\n' )

	f.write( '<BR><HR><FONT SIZE="+2" COLOR="#006600">AB System Photometry (sources within 1 arcmin):</FONT><BR>\n' )
	f.write( 'Notes: Non-zero magnitudes with uncertainty of zero are 3-sigma upper limits.  Sources with magnitudes of 0.0000 are unobserved.<BR>\n' )
	f.write( '		 Circles in images above match aperture size used in sextractor.<BR>\n' )
	f.write( '<br />\n' )

	#Writes table with all magnitudes
	f.write( '<table border="1" width="100%">\n' )
	f.write( '<tr><th>#</th><th>RA</th><th>DEC</th><th>r MAG</th><th>r MAG ERR</th><th>i MAG</th><th>i MAG ERR</th><th>z MAG</th><th>z MAG ERR</th><th>y MAG</th><th>y MAG ERR</th><th>J MAG</th><th>J MAG ERR</th><th>H MAG</th><th>H MAG ERR</th><tr>\n' )
	for i in range(len(plotra)):
	    f.write( '<tr><td>{:.0f}</td><td>{:.6f}</td><td>{:.6f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td></tr>\n'.format(i,plotra[i],plotdec[i],plotrmag[i],plotrmagerr[i],plotimag[i],plotimagerr[i],plotzmag[i],plotzmagerr[i],plotymag[i],plotymagerr[i],plotJmag[i],plotJmagerr[i],plotHmag[i],plotHmagerr[i]) )
	f.write( '</table>\n' )

	f.write( '</PRE><BR><HR>\n' )
	f.write( '<IMG SRC="photcomp.png">\n' )
	f.write( '</PRE><BR><HR>\n' )

	f.write( '</BODY>\n' )
	f.write( '</HTML>\n' )
	f.close()