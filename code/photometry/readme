This is a photometric zeropoint server, using cataloged
 magnitudes from SDSS, USNOB1, APASS and 2-MASS to either 
 lookup or calculate the SEDs for sources detected in an
 image, then using their cataloged/calculated magnitudes
 in a passband to calculate an estimate of the photometric
 zeropoint for a set of instrumental magnitudes.
 (E.G. SExtractor output.)

The main zeropoint.py script assumes that input files are readable
 with numpy.loadtxt --- i.e. #-commented, tab-or-space-delimited, and uniform.
 The first column should be RA, the second DEC, the third the instrumental mag,
 and the rest are arbitrary and unused.

REQUIRES commmand line tools findsdss8, find2mass, findusnob1
 - available from http://cdsarc.u-strasbg.fr/doc/cdsclient.html


++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

EXAMPLES:

To calculate the zeropoint for a field of y-band instrumental
 magnitudes (in the properly-formatted file `00146952_z_DECam.16.p.w.txt'):
$ python zeropoint.py 00146952_z_DECam.16.p.w.txt y
To save the result as a catalog:
$ python zeropoint.py 00146952_z_DECam.16.p.w.txt y catalog.txt

To produce a catalog of all objects possible in a field, run:
> import get_SEDs as gs
> c = gs.catalog( (314.136483, -6.081352), 900. )
*** catalog of objects in a 1/4 degree (900") field centered at 314.136483, -6.081352 ***
*** can run anywhere in the sky, in or out of SDSS footprint ***
Catalog entries are accessible as:
> c.coords
> c.SEDs
And you can save the catalog to an ascii file:
c.save_catalog( 'catalog.txt' )


++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CURRENT NOTES:




+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

