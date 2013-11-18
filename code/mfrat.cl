## #####################################################################
## mfrat.cl ; procedure that generates flat fields and darks ;
## #####################################################################
#  make dark list; dark-ks.list ; dark-files for k-filter reduction
#  make flat list; flat-ks.list ; flat is exptime-independent
#  make object list; g135n-ks-30.list
#  make sky list; sky-g135n-ks-30.list
#  copy bad-pixel map; badmap.pl or badpix.dat
#  epar combine; check options
## #####################################################################
#  noao.imred.ccdred
#  noao.imred.generic
#  proto; imscale, fixpix

   procedure mfrat (filter,exptime)
      string filter
      string exptime
      string *flist
   begin
      string fname
      string frun
      string frout
      string fout
      string lista
      string darkmaster
      string flatmaster
      string flatdark
      string filterarr[6]
      int    i

      filterarr[1]='0'
      filterarr[2]='1'
      filterarr[3]='z'
      filterarr[4]='y'
      filterarr[5]='j'
      filterarr[6]='h'

      FOR(i=1; i<=6; i=i+1){
         flatmaster = "flat_"//filterarr[i]

###CREATE DARK
#No Darks

###CREATE FLAT
	 lista = "im_"//filterarr[i]//".list"
         flist = lista
         print (lista)
         imdel (flatmaster)
         flist = lista
         combine ("@"//lista,output=flatmaster,combine="median",reject="avsigclip",nlow=2,nhigh=2,scale='median',zero='median',lthreshold=-1000,hthreshold=25000)
         normflat(flatmaster,flatmaster)
         while (fscan(flist,fname)==1){print (fname)}
      }
end
