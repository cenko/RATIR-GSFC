;+
; NAME:
;	autopipestack
;
; PURPOSE:
;	Stacks images with same target and filter using SWarp.  Save files as coadd*_(FILTER).fits.
;	Does zeropoint correction on each individual frame then does a stacked coadd with flux scale 
;	for SWarp and calculates absolute zeropoint correction of coadd
;
; OPTIONAL KEYWORDS:
;	outpipevar - output pipeline parameters
;	inpipevar  - input pipeline parameters (typically set in autopipedefaults.pro or ratautoproc.pro, but can be set to default) 
;
; EXAMPLE:
;	autopipestack, outpipevar=pipevar, inpipevar=pipevar
;
; DEPENDENCIES:
;	SWarp, get_SEDs, calc_zpt, findsexobj
;
; Written by Dan Perley 
; Modified by Vicki Toy 12/08/2013
;
; FUTURE IMPROVEMENTS:
;	prefchar in variable structure?, header keywords to keep
;-

pro autopipestack, outpipevar=outpipevar, inpipevar=inpipevar

	;Setup pipeline variables that carry throughout the pipeline
	if keyword_set(inpipevar) then begin
		pipevar = inpipevar
		print, 'Using provided pipevar'
	endif else begin
		pipevar = {autoastrocommand:'autoastrometry' , sexcommand:'sex' , swarpcommand:'swarp' , $
					datadir:'' , imworkingdir:'' , overwrite:0 , $
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:'' }
	endelse
  
  	; if swarp configuration file is not present as 'default.swarp', have swarp output default configuration to this file name
  	if file_test('default.swarp') eq 0 then $
    	spawn, pipevar.swarpcommand+' -d > default.swarp'
 	
 	;Find files that have had astrometry performed on them, stop program if don't exist
 	;VLT CHANGE FOR RIMAS
  	prefchar = '2'
  	azffiles = findfile(pipevar.imworkingdir+'a*'+prefchar+'*_img_?.fits')
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
    	if f eq 0 then datestr = sxpar(h,'DATE')
    	
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
      		stacki = where(filetargets eq target and filefilters eq filter and (fileastrrms1 lt 1.0e-4 and fileastrrms1 gt 5.0e-6) and (fileastrrms2 lt 1.0e-4 and fileastrrms2 gt 5.0e-6), ctstack)

      		if ctstack eq 0 then continue
      		
      		stacklist = azffiles[stacki]
      		stackexps = fileexposures[stacki]
      		
      		medianexp = median(stackexps)
      		medair    = median(fileairval[stacki])
      		minair    = min(fileairval[stacki])
      		maxair    = max(fileairval[stacki])
      		totalexp  = total(stackexps)
      		nstack    = n_elements(stacklist)

			;Create output variables that will be used by SWarp
      		outfile       = pipevar.imworkingdir + 'coadd' + strtrim(target,2) +'_'+ strtrim(filter,2) + '.fits'
      		outweightfile = pipevar.imworkingdir + 'coadd' + strtrim(target,2) +'_'+ strtrim(filter,2) + '.weight.fits'
      		stackcmd      = pipevar.swarpcommand+' '      		
      		
      		;Keywords to carry through, change for RIMAS
      		stackcmd = stackcmd + ' -COPY_KEYWORDS OBJECT,TARGNAME,TELESCOP,FILTER,'+$
                             	'INSTRUME,OBSERVAT,ORIGIN,CCD_TYPE,JD,SOFTGAIN,'+$
                             	'PIXSCALE,WAVELENG,DATE-OBS,AIRMASS,FLATFLD,FLATTYPE '
      		
      		 if (file_test(outfile) eq 0) or pipevar.overwrite then begin
      			filter = strcompress(filter, /REMOVE_ALL)
      			textslist = strjoin(stacklist, ' ')
      			
      			;Run initial weighted stack with all files to make coadded file
      			istackcmd = stackcmd + ' -WRITE_XML N -IMAGEOUT_NAME ' + outfile + ' -WEIGHTOUT_NAME ' + outweightfile + ' ' + textslist
      			print, istackcmd
      			spawn, istackcmd
      			h = headfits(outfile)
      			pixscl = sxpar(h, 'PIXSCALE')
      			
      			;Run sextractor twice on coadded frame (outfile) and find stars with good PSF aperture 
      			;(first run to find aperture to use, 1.34*seeing)
      			findsexobj, outfile, 10.0, pipevar, pix=pixscl, aperture=20.0, wtimage=outweightfile
      			h = headfits(outfile)
      			cpsfdiam = 1.34 * float(sxpar(h, 'SEEPIX'))
      			findsexobj, outfile, 10.0, pipevar, pix=pixscl, aperture=cpsfdiam, wtimage=outweightfile
				refstars1 = outfile+'.stars'
				readcol, refstars1, refnum, refxim, refyim, refmagaper, refmagerraper, refflag, refaim, refbim, refelon, reffwhmim, refclass,refxwor,refywor, reffluxaper, reffluxerraper
				xyad, h, refxim, refyim, refra, refdec
			
				;Only keep reference information where there is an actual sextractor detection
				goodref = where(refmagaper lt 99)
				refdec = refdec[goodref]
				refra  = refra[goodref]
				refmagaper = refmagaper[goodref]
			
				reflist = []
				maglist = []
				errlist = []
			
				;Find stars for each individual frame and try to find matches with coadded frame with each
				;frame optimized with PSF size
				for i = 0, n_elements(stacklist)-1 do begin
					im = stacklist[i]
					h  = headfits(im)
					ipixscl = sxpar(h, 'PIXSCALE')
					findsexobj, im, 3.0, pipevar, pix=ipixscl, aperture=20.0
					h = headfits(im)
					psfdiam = 1.34 * float(sxpar(h, 'SEEPIX'))
					findsexobj, im, 3.0, pipevar, pix=ipixscl, aperture=psfdiam
					starfile = im + '.stars'
				
					newmags = fltarr(n_elements(refmagaper))
					newwts  = fltarr(n_elements(refmagaper))
				
					readcol, starfile, num, xim, yim, magaper, magerraper, flag, aim, bim, elon, fwhmim, class,xwor,ywor, fluxaper, fluxerraper
					xyad, h, xim, yim, imra, imdec
					for u = 0, n_elements(xim)-1 do begin
						matchind = nearest(imra[u] * cos(imdec[u]*!pi/180.), imdec[u], refra * cos(refdec*!pi/180.), refdec, 10.0/3600.0) 
						if matchind eq -1 then continue
						newmags[matchind] = magaper[u]
						newwts[matchind]  = 1.0/(max([magerraper[u], 0.01])^2)
					endfor
				
					maglist = [ [maglist], [newmags] ]
					errlist = [ [errlist], [newwts] ]
					reflist = [ [reflist], [refmagaper] ]
				endfor
			
				;Calculate zeropoint information for each frame and store in header
				mix = calc_zpt(reflist, maglist, errlist, sigma=3.0)

				zpts  = mix[*,0]
				scats = mix[*,1]			
				rmss  = mix[*,2]	
				
				;Move files with bad zeropoint calculations to folder 'badzptfit' and do not use those frames
				badframes = where(finite(zpts) eq 0,nfinct)
				goodframes = where(finite(zpts) eq 1)
				removedframes = []
				newstacklist=stacklist
				
				if nfinct ne 0 then begin
					if dir_exist(pipevar.imworkingdir+'/badzptfit') eq 0 then spawn, 'mkdir '+pipevar.imworkingdir+'/badzptfit'

					for b = 0, n_elements(badframes)-1 do begin
						removedframes = [removedframes, stacklist[badframes[b]]]
						spawn, 'mv ' + stacklist[badframes[b]] +' '+ pipevar.imworkingdir+'badzptfit/'
					endfor
					
					zpts  = zpts[goodframes]
					scats = scats[goodframes]
					rmss  = rmss[goodframes]
					newstacklist = stacklist[goodframes]
				endif

				newtextslist = strjoin(newstacklist, ' ')
				
				;Add relative zeropoint values to headers.  Calculate flux scale
				medzp = median(zpts)
				for i = 0, n_elements(newstacklist)-1 do begin
					im = newstacklist[i]
					h  = headfits(im)
					sxaddpar, h, 'RELZPT', zpts[i], 'Relative zeropoint from calc_zpt'
					sxaddpar, h, 'RELZPTSC', scats[i], 'Robust scatter of relative zeropoint'
					sxaddpar, h, 'RELZPRMS', rmss[i], 'RMS of relative zeropoint'
					sxaddpar, h, 'NEWFLXSCALE', 1.0/(10.0^( (zpts[i]-medzp)/2.5 ) ), 'Flux scaling based on median zp'
					modfits, im, 0, h
				endfor
			
				;Second stacked coadd but now with FSCALE
				istackcmd2 = stackcmd + ' -WRITE_XML N -IMAGEOUT_NAME ' + outfile + ' -WEIGHTOUT_NAME ' + outweightfile + ' -FSCALE_KEYWORD NEWFLXSCALE ' + newtextslist

				print, istackcmd2
      			spawn, istackcmd2
		
				;Run sextractor again on new coadd file
      			findsexobj, outfile, 10.0, pipevar, pix=pixscl, aperture=cpsfdiam, wtimage=outweightfile      		
      			h = headfits(outfile)

				readcol, outfile + '.stars', num, xim, yim, magaper, magerraper, flag, aim, bim, elon, fwhmim, class,xwor,ywor, fluxaper, fluxerraper
				xyad, h, xim, yim, imra, imdec
			
				;Save image file
				imfile  = outfile + '.im'
				catfile = outfile + '.cat'
				writecol, imfile, imra, imdec, magaper
				
				;Filter name correction
				if filter eq 'Z' or filter eq 'Y' then filter = strlowcase(filter)
			
				;Create catalog star file
				sedcmd = "python /Users/vickitoy/Research/RATIR-GSFC/code/photometry/dependencies/get_SEDs.py " + imfile + ' ' + filter + ' ' + catfile + " 15 True"
				print, sedcmd
				spawn, sedcmd
			
				readcol, catfile, refra,refdec,u_mag,g_mag,r_mag,i_mag,z_mag,y_mag,bigB_mag,bigV_mag,bigR_mag,$
					bigI_mag,J_mag,H_mag,K_mag,u_err,g_err,r_err,i_err,z_err,y_err,bigB_err,bigV_err, $
					bigR_err,bigI_err,J_err,H_err,K_err,mode
			
				;Hash table/dictionary to use various filters
				maghash = hash('g_mag', g_mag, 'r_mag', r_mag, 'i_mag', i_mag, 'z_mag', z_mag, 'y_mag', y_mag, 'J_mag', J_mag, 'H_mag', H_mag, 'K_mag', K_mag)
				errhash = hash('g_err', g_err, 'r_err', r_err, 'i_err', i_err, 'z_err', z_err, 'y_err', y_err, 'J_err', J_err, 'H_err', H_err, 'K_err', K_err)
				
				;Find relevant catalog filter values and only use values or actual detections
				refmag = maghash[filter+'_mag']
				goodind = where(mode ne -1 and refmag lt 90.0)
			 
				refmag = refmag[goodind]			
				obsmag = magaper[goodind]
				obserr = magerraper[goodind]
				obswts = fltarr(n_elements(obserr))
				obskeepmag = fltarr(n_elements(obsmag))
				
				;Store magnitudes and weights (with minimum magnitude error of 0.01)
				for j = 0, n_elements(obserr)-1 do begin
					if obserr[j] lt 90.0 then begin
						obskeepmag[j] = obsmag[j]
						obswts[j] = 1.0/(max([obserr[j], 0.01])^2)
					endif
				endfor
			
				;Calculate zeropoint of coadded frame with catalog and observed
				mix2 = calc_zpt([refmag], [obsmag], [obswts], sigma=3.0)

				czpts  = mix2[0]
				cscats = mix2[1]			
				crmss  = mix2[2]
			
				hc = headfits(outfile)
				
				;Add zeropoint keywords to header
				sxaddpar, hc, 'ABSZPT', czpts+25.0, 'Absolute zeropoint from calc_zpt'
				sxaddpar, hc, 'ABSZPTSC', cscats, 'Robust scatter of absolute zeropoint'
				sxaddpar, hc, 'ABSZPRMS', crmss, 'RMS of absolute zeropoint'
			
				;Add summary of stack information to header
				sxaddpar, h, 'DATE'    , datestr
        		sxaddpar, h, 'NSTACK'  , nstack
        		sxaddpar, h, 'AIRMASS' , medair, 'Median exposure airmass'
        		sxaddpar, h, 'AIRMIN'  , minair, 'Minimum exposure airmass'
        		sxaddpar, h, 'AIRMAX'  , maxair, 'Maximum exposure airmass'
        		sxaddpar, h, 'EXPOSURE', medianexp, 'Effective rescaled exposure time'
        		sxaddpar, h, 'TOTALEXP', totalexp, 'Total summed integration time'
        		sxaddpar, h, 'MAXEXP'  , max(stackexps), 'Length of longest exposure'
        		sxaddpar, h, 'MINEXP'  , min(stackexps), 'Length of shortest exposure'
        		sxaddpar, h, 'SATURATE', min(filesatval[stacki]-fileskyval[stacki])
        		sxaddpar, h, 'MEDSKY'  , median(fileskyval[stacki], /even)
        		
				modfits, outfile, 0, hc
				
				print, 'Removed frames with no overlapping stars: ' 
				print, removedframes
			endif
    	endfor
  	endfor

	outpipevar = pipevar
	
end

