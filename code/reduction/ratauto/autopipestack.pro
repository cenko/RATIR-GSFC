;+
; NAME:
;	autopipestack
;
; PURPOSE:
;	Does zeropoint correction on each individual frame using sextractor and get_SEDs. 
;   Creates flux scale (newflxsc) from how close to median of zeropoint values.  Uses
;   flux scale to stack images in Swarp (has moved bad zeropoint values and bad newflxsc
;   values to marked folders - badzptfit/ and badflxsc/) and calculates absolute zeropoint 
;   correction of coadd.  Saves zeropoint plot as zpt_(FILTER).ps
;
; OPTIONAL KEYWORDS:
;	outpipevar - output pipeline parameters
;	inpipevar  - input pipeline parameters (typically set in autopipedefaults.pro or ratautoproc.pro, but can be set to default) 
;
; EXAMPLE:
;	autopipestack, outpipevar=pipevar, inpipevar=pipevar
;
; DEPENDENCIES:
;	SWarp, get_SEDs, calc_zpt, findsexobj (sextractor)
;
; Written by Dan Perley 
; Modified by Vicki Toy 9/29/2014
;
; FUTURE IMPROVEMENTS:
;	header keywords to keep
;-

pro autopipestack, outpipevar=outpipevar, inpipevar=inpipevar

	print, 'STACK'
	spawn, 'export CDSCLIENT=http' ;Fix for problem with timeout with CDSCLIENT
	
	;Setup pipeline variables that carry throughout the pipeline
	if keyword_set(inpipevar) then begin
		pipevar = inpipevar
		if pipevar.verbose gt 0 then print, 'Using provided pipevar'
	endif else begin
	    pipevar = {autoastrocommand:'autoastrometry', getsedcommand:'get_SEDs', $
					sexcommand:'sex' , swarpcommand:'swarp' , $
					prefix:'', datadir:'' , imworkingdir:'' , overwrite:0 , verbose:0, rmifiles:0,$
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:'' }
	endelse
  	
  	if pipevar.verbose gt 0 then quiet = 0 else quiet = 1
  
  	; if swarp configuration file is not present as 'default.swarp', have swarp output default configuration to this file name
  	if file_test('default.swarp') eq 0 then $
    	spawn, pipevar.swarpcommand+' -d > default.swarp'
 	
 	;Find files that have had astrometry performed on them, stop program if don't exist
  	azffiles = findfile(pipevar.imworkingdir+'a*'+pipevar.prefix+'*_img_?.fits')
  	if azffiles[0] eq '' then return

  	filetargets   = strarr(n_elements(azffiles))
  	fileexposures = fltarr(n_elements(azffiles))
  	filefilters   = strarr(n_elements(azffiles))
  	filesatval    = fltarr(n_elements(azffiles))
  	fileskyval 	  = fltarr(n_elements(azffiles))
  	fileairval    = fltarr(n_elements(azffiles))
  	fileastrrms1  = fltarr(n_elements(azffiles))
  	fileastrrms2  = fltarr(n_elements(azffiles))
  	datestr = ''

	;Grab information in the headers of astrometry corrected file and save to array
  	for f = 0, n_elements(azffiles)-1 do begin
    	h = headfits(azffiles[f], /silent)
    	if f eq 0 then datestr = sxpar(h,'DATE-OBS')
    	
    	filetargets[f]   = repstr(strtrim(sxpar(h,'TARGNAME'),2),' ', '_')  ; spaces cause barfing in filenames
    	fileexposures[f] = sxpar(h, 'EXPOSURE')
    	filefilters[f]   = string(sxpar(h,'FILTER'))
    	filesatval[f]    = sxpar(h,'SATURATE')
    	fileskyval[f]    = sxpar(h,'COUNTS')
    	fileairval[f]    = sxpar(h,'AIRMASS')
    	fileastrrms1[f]  = sxpar(h,'ASTRRMS1')
    	fileastrrms2[f]  = sxpar(h,'ASTRRMS2')
    	    
  	endfor
  	
  	targets = unique(filetargets)
  	for t = 0, n_elements(targets)-1 do begin
  	
  		;Finds files with same target and the filters associated with this target
    	target     = targets[t]
    	thistarget = where(filetargets eq target, cttarg)
    	if cttarg eq 0 then continue
    	thistargetfilts = unique(filefilters[thistarget])
    
    	;For each filter find files that have the same target and same filter 
    	;and store information on the exposure times and airmass
    	;Puts constraints on how good of an astrometric fit we have from Scamp
    	for l = 0, n_elements(thistargetfilts)-1 do begin
      		filter = thistargetfilts[l]
      		stacki = where(filetargets eq target and filefilters eq filter and (fileastrrms1 lt 2.0e-4 and fileastrrms1 gt 5.0e-6) and (fileastrrms2 lt 2.0e-4 and fileastrrms2 gt 5.0e-6), ctstack)

      		if ctstack eq 0 then continue
      		
      		stacklist = azffiles[stacki]
      		stackexps = fileexposures[stacki]
      		
      		medianexp = median(stackexps)
      		medair    = median(fileairval[stacki])
      		minair    = min(fileairval[stacki])
      		maxair    = max(fileairval[stacki])
      		totalexp  = total(stackexps)
      		nstack    = n_elements(stacklist)
      		
      		filter = strcompress(filter, /REMOVE_ALL)
      		textslist = strjoin(stacklist, ' ')
      						
			zpts = []			
			
			;Find stars for each individual frame and try to find matches with coadded frame with each
			;frame optimized with PSF size
			for i = 0, n_elements(stacklist)-1 do begin
				im = stacklist[i]
				h  = headfits(im)
				ipixscl = sxpar(h, 'PIXSCALE')
				findsexobj, im, 3.0, pipevar, pix=ipixscl, aperture=20.0, quiet=quiet
				starfile = im + '.stars'
				
				readcol, starfile, num, xim, yim, magaper, magerraper, flag, aim, bim, elon, fwhmim, class,xwor,ywor, fluxaper, fluxerraper
				xyad, h, xim, yim, imra, imdec

				;Save image file
				imfile  = im + '.im'
				catfile = im + '.cat'
				writecol, imfile, imra, imdec, magaper
				
				;Filter name correction
				if filter eq 'Z' or filter eq 'Y' then filter = strlowcase(filter)
			
				;Create catalog star file (python get_SEDs.py imfile filter catfile USNOB_THRESH alloptstars
				if pipevar.verbose gt 0 then qtcmd = 'True' else qtcmd = 'False'
				sedcmd = pipevar.getsedcommand + ' ' + imfile + ' ' + filter + ' ' + catfile + " 15 True "+ qtcmd
				if pipevar.verbose gt 0 then print, sedcmd
				spawn, sedcmd
				print, sedcmd
				
				if not(file_test(catfile)) then begin
					zpts = [zpts,!Values.F_NAN]
					continue
				endif			
				readcol, catfile, refra,refdec,u_mag,g_mag,r_mag,i_mag,z_mag,y_mag,bigB_mag,bigV_mag,bigR_mag,$
					bigI_mag,J_mag,H_mag,K_mag,u_err,g_err,r_err,i_err,z_err,y_err,bigB_err,bigV_err, $
					bigR_err,bigI_err,J_err,H_err,K_err,mode
			
				;Hash table/dictionary to use various filters
				maghash = hash('g_mag', g_mag, 'r_mag', r_mag, 'i_mag', i_mag, 'z_mag', z_mag, 'y_mag', y_mag, 'J_mag', J_mag, 'H_mag', H_mag, 'K_mag', K_mag)
				errhash = hash('g_err', g_err, 'r_err', r_err, 'i_err', i_err, 'z_err', z_err, 'y_err', y_err, 'J_err', J_err, 'H_err', H_err, 'K_err', K_err)
				
				;Find relevant catalog filter values and only use values or actual detections
				refmag  = maghash[filter+'_mag']
				goodind = where(mode ne -1 and refmag lt 90.0 and flag lt 8 and elon le 1.3)
			 
				refmag = refmag[goodind]			
				obsmag = magaper[goodind]
				obserr = magerraper[goodind]
				obswts = fltarr(n_elements(obserr))
				obskeepmag = fltarr(n_elements(obsmag))
				
				;Store magnitudes and weights (with minimum magnitude error of 0.01)
				for j = 0, n_elements(obserr)-1 do begin
					if obserr[j] lt 0.1 then begin
						obskeepmag[j] = obsmag[j]
						obswts[j] = 1.0/(max([obserr[j], 0.01])^2)
					endif
				endfor
				mix   = calc_zpt(refmag, obskeepmag, obswts, sigma=3.0)
				zpt   = mix[0]
				scats = mix[1]
				rmss  = mix[2]
				
				sxaddpar, h, 'ABSZPT', zpt + 25.0, 'Relative zeropoint from calc_zpt'
				sxaddpar, h, 'ABSZPTSC', scats, 'Robust scatter of relative zeropoint'
				sxaddpar, h, 'ABSZPRMS', rmss, 'RMS of relative zeropoint'
				
				modfits, im, 0, h
				zpts = [zpts,zpt]
		
			endfor
				
			;Move files with bad zeropoint calculations to folder 'badzptfit' and do not use those frames
			badframes  = where(finite(zpts) eq 0, nfinct, complement=goodframes)
			removedframes = []
			newstacklist=stacklist
				
			if nfinct ne 0 then begin
				if dir_exist(pipevar.imworkingdir+'/badzptfit') eq 0 then spawn, 'mkdir '+pipevar.imworkingdir+'/badzptfit'

				for b = 0, n_elements(badframes)-1 do begin
					removedframes = [removedframes, stacklist[badframes[b]]]
					spawn, 'mv ' + stacklist[badframes[b]] +' '+ pipevar.imworkingdir+'badzptfit/'
				endfor
					
				zpts  = zpts[goodframes]
				newstacklist = stacklist[goodframes]
			endif
		
			badnewflxsc = []
			;Add relative zeropoint values to headers.  Calculate flux scale.  Remove unphysical fluxscale files
			medzp = median(zpts)
			for i = 0, n_elements(newstacklist)-1 do begin
				im = newstacklist[i]
				h  = headfits(im)
				sxaddpar, h, 'NEWFLXSC', 1.0/(10.0^( (zpts[i]-medzp)/2.5 ) ), 'Flux scaling based on median zp'
					
				if 1.0/(10.0^( (zpts[i]-medzp)/2.5 ) ) lt 0.1 then badnewflxsc = [badnewflxsc, im]
				modfits, im, 0, h
			endfor
			    
			;Removes files that have bad newflxsc values and removes these from stack list as well
			if n_elements(badnewflxsc) gt 0 then begin
				if dir_exist(pipevar.imworkingdir+'/badflxsc') eq 0 then spawn, 'mkdir '+pipevar.imworkingdir+'/badflxsc'

				spawn, 'mv ' + badnewflxsc +' '+ pipevar.imworkingdir+'badflxsc/'
                
                removedframes = [removedframes, badnewflxsc]
                ;Remove files that have bad newflxsc values from list of stack
                match, badnewflxsc, newstacklist, subbad, substack
                remove, substack, newstacklist
			endif			    
			    
			newtextslist = strjoin(newstacklist, ' ')

      		stackcmd      = pipevar.swarpcommand+' '     		
      		
      		;Keywords to carry through, change for RIMAS
      		stackcmd = stackcmd + ' -COPY_KEYWORDS OBJECT,TARGNAME,TELESCOP,FILTER,'+$
                             	'INSTRUME,OBSERVAT,ORIGIN,CCD_TYPE,JD,SOFTGAIN,'+$
                             	'PIXSCALE,WAVELENG,DATE-OBS,AIRMASS,FLATFLD,FLATTYPE '
			
			;Create output variables that will be used by SWarp
      		outfile       = pipevar.imworkingdir + 'coadd' + strtrim(target,2) +'_'+ strtrim(filter,2) + '.fits'
      		outweightfile = pipevar.imworkingdir + 'coadd' + strtrim(target,2) +'_'+ strtrim(filter,2) + '.weight.fits'
      
            if pipevar.verbose gt 0 then stackcmd = stackcmd + ' -VERBOSE_TYPE NORMAL ' else stackcmd = stackcmd + ' -VERBOSE_TYPE QUIET ' 
                
			;Coadd with flux scale
			istackcmd = stackcmd + ' -SUBTRACT_BACK N -WRITE_XML N -IMAGEOUT_NAME ' + outfile + ' -WEIGHTOUT_NAME ' + outweightfile + ' -FSCALE_KEYWORD NEWFLXSC ' + newtextslist
			
			if pipevar.verbose gt 0 then print, istackcmd
      		spawn, istackcmd
      			
      		h = headfits(outfile)
      		pixscl = sxpar(h, 'PIXSCALE')
      		findsexobj, outfile, 10.0, pipevar, pix=pixscl, aperture=20.0, wtimage=outweightfile, quiet=quiet

      		h = headfits(outfile)
      		cpsfdiam = 1.34 * float(sxpar(h, 'SEEPIX'))
		
			;Run sextractor again on new coadd file
      		findsexobj, outfile, 3.0, pipevar, pix=pixscl, aperture=cpsfdiam, wtimage=outweightfile, quiet=quiet      		
      		h = headfits(outfile)

			readcol, outfile + '.stars', num, xim, yim, magaper, magerraper, flag, aim, bim, elon, fwhmim, class,xwor,ywor, fluxaper, fluxerraper
			xyad, h, xim, yim, imra, imdec
			
			;Save image file
			imfile  = outfile + '.im'
			catfile = outfile + '.cat'
			writecol, imfile, imra, imdec, magaper
				
			;Filter name correction
			if filter eq 'Z' or filter eq 'Y' then filter = strlowcase(filter)
			
			;Create catalog star file (python get_SEDs.py imfile filter catfile USNOB_THRESH alloptstars
			if pipevar.verbose gt 0 then qtcmd = 'True' else qtcmd = 'False'
			sedcmd = pipevar.getsedcommand + ' ' + imfile + ' ' + filter + ' ' + catfile + " 15 False "+ qtcmd
			if pipevar.verbose gt 0 then print, sedcmd
			spawn, sedcmd
			
			readcol, catfile, refra,refdec,u_mag,g_mag,r_mag,i_mag,z_mag,y_mag,bigB_mag,bigV_mag,bigR_mag,$
				bigI_mag,J_mag,H_mag,K_mag,u_err,g_err,r_err,i_err,z_err,y_err,bigB_err,bigV_err, $
				bigR_err,bigI_err,J_err,H_err,K_err,mode
			
			;Hash table/dictionary to use various filters
			maghash = hash('g_mag', g_mag, 'r_mag', r_mag, 'i_mag', i_mag, 'z_mag', z_mag, 'y_mag', y_mag, 'J_mag', J_mag, 'H_mag', H_mag, 'K_mag', K_mag)
			errhash = hash('g_err', g_err, 'r_err', r_err, 'i_err', i_err, 'z_err', z_err, 'y_err', y_err, 'J_err', J_err, 'H_err', H_err, 'K_err', K_err)
				
			;Find relevant catalog filter values and only use values or actual detections
			refmag = maghash[filter+'_mag']
			goodind = where(mode eq 0 and refmag lt 90.0 and flag eq 0 and elon le 1.3)
			 
			refmag = refmag[goodind]			
			obsmag = magaper[goodind]
			obserr = magerraper[goodind]
			obswts = fltarr(n_elements(obserr))
			obskeepmag = fltarr(n_elements(obsmag))
				
			;Store magnitudes and weights (with minimum magnitude error of 0.01)
			for j = 0, n_elements(obserr)-1 do begin
				if obserr[j] lt 0.1 then begin
					obskeepmag[j] = obsmag[j]
					obswts[j] = 1.0/(max([obserr[j], 0.01])^2)
				endif
			endfor
				
			;Calculate zeropoint of coadded frame with catalog and observed
			mix2 = calc_zpt([refmag], [obskeepmag], [obswts], sigma=3.0, plotter=pipevar.imworkingdir+'zpt_'+filter+'.ps')
			czpts  = mix2[0]
			cscats = mix2[1]			
			crmss  = mix2[2]
			
			hc = headfits(outfile)
				
			;Add zeropoint keywords to header
			sxaddpar, hc, 'SPIX'    , cpsfdiam, 'Final aperture size'
			sxaddpar, hc, 'ABSZPT'  , czpts+25.0, 'Absolute zeropoint from calc_zpt'
			sxaddpar, hc, 'ABSZPTSC', cscats, 'Robust scatter of absolute zeropoint'
			sxaddpar, hc, 'ABSZPRMS', crmss, 'RMS of absolute zeropoint'
			
			;Add summary of stack information to header
			sxaddpar, hc, 'DATE'    , datestr
        	sxaddpar, hc, 'NSTACK'  , nstack
        	sxaddpar, hc, 'AIRMASS' , medair, 'Median exposure airmass'
        	sxaddpar, hc, 'AIRMIN'  , minair, 'Minimum exposure airmass'
        	sxaddpar, hc, 'AIRMAX'  , maxair, 'Maximum exposure airmass'
        	sxaddpar, hc, 'EXPOSURE', medianexp, 'Effective rescaled exposure time'
        	sxaddpar, hc, 'TOTALEXP', totalexp, 'Total summed integration time'
        	sxaddpar, hc, 'MAXEXP'  , max(stackexps), 'Length of longest exposure'
        	sxaddpar, hc, 'MINEXP'  , min(stackexps), 'Length of shortest exposure'
        	sxaddpar, hc, 'SATURATE', min(filesatval[stacki]-fileskyval[stacki])
        	sxaddpar, hc, 'MEDSKY'  , median(fileskyval[stacki], /even)
        		
        	for f = 0, n_elements(newstacklist)-1 do begin
  				sxaddpar, hc, 'STACK'+strtrim(string(f),2), removepath(newstacklist[f])
			endfor
        		
			modfits, outfile, 0, hc
			if removedframes ne [] then begin
				print, 'Removed frames with bad zeropoint fits: ' 
				print, removedframes
			endif
    	endfor
  	endfor

	if pipevar.rmifiles then begin
	
	    ;If remove intermediate files keyword set, delete p(PREFIX)*.fits files
	    pfiles   = findfile(pipevar.imworkingdir+'p'+pipevar.prefix+'*.fits')
	    ffiles   = findfile(pipevar.imworkingdir+'fp'+pipevar.prefix+'*.fits')
	    skyfiles = findfile(pipevar.imworkingdir+'sky-*.fits')
	    sfiles   = findfile(pipevar.imworkingdir+'sfp'+pipevar.prefix+'*.fits')
	    zfiles   = findfile(pipevar.imworkingdir+'zsfp*'+pipevar.prefix+'*.fits')
	    aimfiles = findfile(pipevar.imworkingdir+'a*fp'+pipevar.prefix+'*.im')
	    astfiles = findfile(pipevar.imworkingdir+'a*fp'+pipevar.prefix+'*.stars')
	    actfiles = findfile(pipevar.imworkingdir+'a*fp'+pipevar.prefix+'*.cat')
	    rfiles  = [pfiles,ffiles,skyfiles,sfiles,zfiles, aimfiles, astfiles,actfiles]
	
	    good = where(rfiles ne '', ngood)
	    if ngood gt 0 then begin
	        rfiles = rfiles[good]
	        file_delete, rfiles
	    endif
	    
	endif

	outpipevar = pipevar
	
end

