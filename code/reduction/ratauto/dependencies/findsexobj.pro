pro findsexobj, inlist, sigma, pipevar, skyval=skyval, masksfx=masksfx, zeropt=zeropt, $
	wtimage=wtimage, wtcut=wtcut, fwhm=fwhm, pix=pix, gain=gain, aperture=aperture, $
	minlim=minlim, minval=minval, class_cut=class_cut, elong_cut=elong_cut

	;Set default values if keywords not set
	if (not keyword_set(skyval)) 	then skyval		= 0.0
	if (not keyword_set(masksfx)) 	then masksfx 	= 'mask'
	if (not keyword_set(zeropt)) 	then zeropt 	= 25.0
	if (not keyword_set(wtcut)) 	then wtcut 		= 0.1	
	if (not keyword_set(fwhm))		then fwhm 		= 1.5	
	if (not keyword_set(pix)) 		then pix 		= 0.3787	
	if (not keyword_set(gain)) 		then gain 		= 25.0	
	if (not keyword_set(aperture)) 	then aperture 	= 5.0	
	if (not keyword_set(class_cut)) then class_cut 	= 0.80	
	if (not keyword_set(elong_cut)) then elong_cut 	= 1.30	
	
	redo = pipevar.overwrite
	
	if file_test('coadd.param') eq 0 then spawn, 'cp '+ pipevar.defaultspath +'/coadd.param .'
	if file_test('coadd.conv') eq 0 then spawn, 'cp '+ pipevar.defaultspath +'/coadd.conv .'
	if file_test('coadd.config') eq 0 then spawn, 'cp '+ pipevar.defaultspath +'/coadd.config .'
	if file_test('default.nnw') eq 0 then spawn, 'cp '+ pipevar.defaultspath +'/default.nnw .'

	
	for i = 0, n_elements(inlist)-1 do begin
		if inlist[i] eq '' then continue
		image = inlist[i]
		
		if (not file_test(image)) then continue
		starfile = image + '.stars'
		
		extpos = strpos(image, '.')
		trunim = strmid(image, 0, extpos)
		mskimg = trunim + '_' + masksfx + '.fits'
		
		;if file_test(starfile) and redo eq 0 then continue
		;if file_test(mskimg)   and redo eq 0 then continue	
		
		;if keyword_set(minlim) then begin
			
		;endif
		sexcommand = pipevar.sexcommand + ' -c coadd.config -DETECT_THRESH ' + strcompress(sigma, /REMOVE_ALL) + ' -ANALYSIS_THRESH ' + strcompress(sigma, /REMOVE_ALL) + $
			' -PHOT_APERTURES ' + strcompress(aperture, /REMOVE_ALL) + ' -MAG_ZEROPOINT ' + strcompress(zeropt, /REMOVE_ALL) + ' -PIXEL_SCALE ' + strcompress(pix, /REMOVE_ALL) + $
			' -SEEING_FWHM ' + strcompress(fwhm, /REMOVE_ALL) + ' -CHECKIMAGE_TYPE OBJECTS' + ' -CHECKIMAGE_NAME ' + mskimg
			
		if keyword_set(wtimage) then begin
			sexcommand = sexcommand + ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE ' + wtimage + ' '
		endif
		
		sexcommand = sexcommand + ' ' + image
		print, sexcommand
		spawn, sexcommand
		
		print, 'mv -f test.cat ' + starfile
		spawn, 'mv -f test.cat ' + starfile
		
		if file_test(starfile) then begin
			readcol, starfile, num, xim, yim, magaper, magerraper, flag, aim, bim, elon, fwhmim, class,xwor,ywor, fluxaper, fluxerraper
			keep = where( (flag eq 0) and (elon lt elong_cut) and (fwhmim gt 0.25) and (fwhmim lt 20.0), keepct )
			if keepct le 1 then seepix=!values.f_nan else seepix = median(fwhmim[keep])
			regfile = trunim + '.reg'
		endif else begin
			print, 'Failed to find Sextractor output file!'
			seepix=!values.f_nan
		endelse
			
		h = headfits(image)
		sxaddpar, h, 'MASKNAME', mskimg, "Object mask image from Sextractor"
		sxaddpar, h, "STARFILE", starfile, "Objects file from Sextractor" 
        sxaddpar, h, "ZEROPT", zeropt,   "Photometric zero-point used for Sextractor"
        sxaddpar, h, "SEEPIX", seepix,   "Estimated seeing from Sextractor objects (pix)"
        sxaddpar, h, "NSTARS", n_elements(num), "Estimated number of objects from Sextractor"
        modfits, image, 0, h
						
	endfor
	
	spawn, 'rm -f coadd.param'
	spawn, 'rm -f coadd.conv'
	spawn, 'rm -f coadd.config'	
	spawn, 'rm -f default.nnw'
	
end
		