"""
Translated from printratirhtml.pro by John Capone (jicapone@astro.umd.edu).
"""

import numpy as np
import pylab as pl
import scipy as sp
import fnmatch
import os

def printratirhtml():
	
	plotra, plotdec, plotrmag, plotrmagerr, plotimag, plotimagerr, plotzmag, plotzmagerr, plotymag, plotymagerr, plotJmag, plotJmagerr, plotHmag, plotHmagerr = np.loadtxt('./finalmags.txt', unpack=True)

	f = open( './ratir.html', 'w' )
	f.write( '<!DOCTYPE HTML>\n' )
	f.write( '<HTML>\n' )
	f.write( '<HEAD>\n' )
	f.write( '<TITLE>RATIR DATA</TITLE>\n' )
	f.write( '</HEAD>\n' )
	f.write( '<BODY BGCOLOR="#FFFFFF" TEXT="#003300">\n' )

	# returns files in directory "loc" which start with prefix and end with postfix
	def get_files( selection, loc='.' ):
	    matches = []
	    for files in os.listdir(loc):
	        if fnmatch.fnmatch( files, selection ):
	            matches.append(files)
	    return matches
	prefchar = 'coadd'
	wildcharimg = '?????-????-?????_?'
	zffiles = get_files( prefchar + wildcharimg + '.png' )

	im_wid = 400

	f.write( '<IMG SRC="./color.png" width="' + `im_wid*3` + '"><BR>\n' )
	for i in range(len(zffiles)):
	    f.write( '<IMG SRC="./' + zffiles[i] + '" width="' + `im_wid` + '">\n' )
	    if (i+1)%3 == 0:
	    	f.write( '<BR>\n' )

	f.write( '<BR><HR><FONT SIZE="+2" COLOR="#006600">AB System Photometry (sources within 1 arcmin):</FONT><BR>\n' )
	f.write( 'Notes: Non-zero magnitudes with uncertainty of zero are 3-sigma upper limits.  Sources with magnitudes of 0.0000 are unobserved.<BR>\n' )

	realdetections = np.where( (plotrmagerr > 0) & (plotimagerr > 0) )[0]
	f.write( '<br />\n' )

	f.write( '<table border="1" width="100%">\n' )
	f.write( '<tr><th>#</th><th>RA</th><th>DEC</th><th>r MAG</th><th>r MAG ERR</th><th>i MAG</th><th>i MAG ERR</th><th>z MAG</th><th>z MAG ERR</th><th>y MAG</th><th>y MAG ERR</th><th>J MAG</th><th>J MAG ERR</th><th>H MAG</th><th>H MAG ERR</th><tr>\n' )
	for i in realdetections:
	    f.write( '<tr><td>{:.0f}</td><td>{:.6f}</td><td>{:.6f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td></tr>\n'.format(i,plotra[i],plotdec[i],plotrmag[i],plotrmagerr[i],plotimag[i],plotimagerr[i],plotzmag[i],plotzmagerr[i],plotymag[i],plotymagerr[i],plotJmag[i],plotJmagerr[i],plotHmag[i],plotHmagerr[i]) )
	f.write( '</table>\n' )

	f.write( '</PRE><BR><HR>\n' )
	f.write( '<IMG SRC="photcomp.jpg">\n' )
	f.write( '</PRE><BR><HR>\n' )

	f.write( '</BODY>\n' )
	f.write( '</HTML>\n' )
	f.close()