;+
; NAME:
;	findsexobj
;
; PURPOSE:
;	Finds sextractor objects with optional input parameters and default values.  Saves
;	object mask, starfile, and estimates the seeing from the stars found. 
;
; INPUTS:
;	inlist  - input list of fits files or single fits file to run sextractor on
;	sigma   - detection threshold and analysis threshold for sextractor
;	pipevar - pipeline parameters (typically set in autopipedefaults.pro or ratautoproc.pro)
;
; OPTIONAL KEYWORDS:
;	masksfx   - text string identifier for sextractor CHECKIMAGE_NAME (no default: ex. masksfx='mask' saves file as filename.mask.fits)
;	zeropt    - input value for sextractor MAG_ZEROPOINT (default is 25.0)
;	wtimage   - file for input for sextractor WEIGHT_IMAGE (uses weight map, no default)
;	fwhm      - input value for sextractor SEEING_FWHM (default is 1.5)
;	pix       - input value for sextractor PIXEL_SCALE (default is 0.3787)
;	aperture  - input value for sextractor PHOT_APERTURES (default is 5.0)
;	elong_cut - cutoff limit for use in FWHM calculation of elongation to eliminate non-stars (default is 1.3)
;	quiet     - no output from sextractor if set
;
; EXAMPLE:
;	findsexobj, outfile, 10.0, pipevar, pix=pixscl, aperture=20.0, wtimage=outweightfile
;
; DEPENDENCIES:
;	Sextractor
;
; Written by Brad Cenko 
; Abridged version written by Vicki Toy 7/7/2014
;
; FUTURE IMPROVEMENTS:
;	More keywords to sextractor?
;-

pro findsexobj, inlist, sigma, pipevar, masksfx=masksfx, zeropt=zeropt, $
	wtimage=wtimage, fwhm=fwhm, pix=pix, aperture=aperture, elong_cut=elong_cut, quiet=quiet

	;Set default values if keywords not set
	if (not keyword_set(zeropt)) 	then zeropt 	= 25.0
	if (not keyword_set(wtcut)) 	then wtcut 		= 0.1	
	if (not keyword_set(fwhm))		then fwhm 		= 1.5	
	if (not keyword_set(pix)) 		then pix 		= 0.3787	
	if (not keyword_set(aperture)) 	then aperture 	= 5.0	
	if (not keyword_set(elong_cut)) then elong_cut 	= 1.30	
	if (not keyword_set(quiet))		then quiet		= 0
	
	;Move necessary sextractor configuration files if they are not in current directory
	if file_test('coadd.param') eq 0 then spawn, 'cp '+ pipevar.defaultspath +'/coadd.param .'
	if file_test('coadd.conv') eq 0 then spawn, 'cp '+ pipevar.defaultspath +'/coadd.conv .'
	if file_test('coadd.config') eq 0 then spawn, 'cp '+ pipevar.defaultspath +'/coadd.config .'
	if file_test('default.nnw') eq 0 then spawn, 'cp '+ pipevar.defaultspath +'/default.nnw .'

	if quiet gt 0 then verbosetype = 'QUIET' else verbosetype = 'NORMAL'
	
	;For each fits file run sextractor with given input parameters. Saves temp.cat as starfile, 
	;saves starmask, and calculates seeing from starlike objects. Saves necessary parameters to header
	for i = 0, n_elements(inlist)-1 do begin
		if inlist[i] eq '' then continue
		image = inlist[i]
		
		if (not file_test(image)) then continue
		starfile = image + '.stars'
		
		extpos = strpos(image, '.')
		trunim = strmid(image, 0, extpos)
		
		sexcommand = pipevar.sexcommand + ' -c coadd.config -DETECT_THRESH ' + strcompress(sigma, /REMOVE_ALL) + ' -ANALYSIS_THRESH ' + strcompress(sigma, /REMOVE_ALL) + $
			' -PHOT_APERTURES ' + strcompress(aperture, /REMOVE_ALL) + ' -MAG_ZEROPOINT ' + strcompress(zeropt, /REMOVE_ALL) + ' -PIXEL_SCALE ' + strcompress(pix, /REMOVE_ALL) + $
			' -SEEING_FWHM ' + strcompress(fwhm, /REMOVE_ALL) + ' -VERBOSE_TYPE ' +verbosetype
		
		if keyword_set(masksfx) then begin
			mskimg = trunim + '_' + masksfx + '.fits'	
			sexcommand = sexcommand + ' -CHECKIMAGE_TYPE OBJECTS' + ' -CHECKIMAGE_NAME ' + mskimg
		endif
			
		if keyword_set(wtimage) then begin
			if n_elements(wtimage) eq 0 then iwtimage = wtimage else iwtimage = wtimage[i]
			sexcommand = sexcommand + ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE ' + iwtimage + ' '
		endif
		
		sexcommand = sexcommand + ' ' + image
		if quiet eq 0 then print, sexcommand
		spawn, sexcommand
		
		if quiet eq 0 then print, 'mv -f test.cat ' + starfile
		spawn, 'mv -f test.cat ' + starfile
		
		;Calculates seeing with starlike objects
		if file_test(starfile) then begin
			readcol, starfile, num, xim, yim, magaper, magerraper, flag, aim, bim, elon, fwhmim, class,xwor,ywor, fluxaper, fluxerraper
			keep = where( (flag eq 0) and (elon lt elong_cut) and (fwhmim gt 0.25) and (fwhmim lt 20.0), keepct )
			if keepct le 1 then seepix=!values.f_nan else seepix = median(fwhmim[keep])
		endif else begin
			print, 'Failed to find Sextractor output file!'
			seepix=!values.f_nan
		endelse
		
		;Writes to header	
		h = headfits(image)
		
		if keyword_set(masksfx) then begin
			sxaddpar, h, 'MASKNAME', mskimg, "Object mask image from Sextractor"
		endif
		
		sxaddpar, h, "STARFILE", starfile, "Objects file from Sextractor" 
        sxaddpar, h, "ZEROPT", zeropt,   "Photometric zero-point used for Sextractor"
        sxaddpar, h, "SEEPIX", seepix,   "Estimated seeing from Sextractor objects (pix)"
        sxaddpar, h, "NSTARS", n_elements(num), "Estimated number of objects from Sextractor"
        modfits, image, 0, h
						
	endfor
	
	;Removes config files after done
	spawn, 'rm -f coadd.param'
	spawn, 'rm -f coadd.conv'
	spawn, 'rm -f coadd.config'	
	spawn, 'rm -f default.nnw'
	
end
		