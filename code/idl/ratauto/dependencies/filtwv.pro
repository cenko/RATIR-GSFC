function filtwv, filter, fwhm=fwhm
  ; return the effective wavelength in units of Angstroms

  if n_elements(filter) gt 1 then begin
    warr = fltarr(n_elements(filter))
    fwhm = fltarr(n_elements(filter))
    for i = 0, n_elements(filter)-1 do begin
      warr[i] = filtwv(filter[i], fwhm=fwhmi)
      fwhm[i] = fwhmi
    endfor
    return, warr
  endif

  w = 1.
  fwhm = 1.
  
  f = (strsplit(filter[0],'_',/extract)) [0]
  f = (strsplit(f,"'",/extract)) [0]

  if f eq 'RC' or f eq 'Rc' then f = 'R'
  if f eq 'IC' or f eq 'Ic' then f = 'I'
  if f eq 'z0' then f = 'z'

  if f eq 'none' then w = 6588.

  if f eq 'White' then w = 5000.
  if f eq 'white' then w = 6588.
  if f eq 'clear' then w = 6588.

  if f eq 'UVW2' then w = 1800.
  if f eq 'UVM2' then w = 2200.
  if f eq 'UVW1' then w = 2600.
  if f eq 'UVU' then w = 3465.
  if f eq 'UVB' then w = 4392.
  if f eq 'UVV' then w = 5468.

  if f eq 'U' then w = 3652.
  if f eq 'B' then w = 4458.
  if f eq 'V' then w = 5505.   ;Fukujita et al. 1995
  if f eq 'R' then w = 6588.
  if f eq 'I' then w = 8060.

  if f eq 'u' then w = 3585.
  if f eq 'g' then w = 4858.
  if f eq 'r' then w = 6290.   ;Fukujita et al. 1995
  if f eq 'i' then w = 7706.
  if f eq 'z' then w = 9222.
  if f eq 'Z' then w = 9222.
  if f eq 'Y' then w = 10200.

  if f eq 'J' then w = 12350.
  if f eq 'H' then w = 16620. ;Cohen et al. 2003
  if f eq 'K' then w = 21590.
  if f eq 'Ks' then w = 21590.
  if f eq 'Kp' or filter[0] eq "K'" then w = 21240.
  
  if f eq 'FUV' then w = 1528.
  if f eq 'NUV' then w = 2271.  ; GALEX filters, from Dale et al.

  if f eq 'keV' then w = 12.398
  if f eq 'X-ray' then w = 12.398

  ;VLA bands use units of 10^-5 m, not Angstroms!
  ; do NOT use simultaneously with optical
  if f eq 'X' then w = 3570 ;
  if f eq 'C' then w = 6310 ;
  
  if strmid(f,0,4) eq 'F222' then w = 2220.
  if strmid(f,0,4) eq 'F435' then w = 4350.
  if strmid(f,0,5) eq 'F450W' then w = 4500.
  if strmid(f,0,5) eq 'F606W' then w = 6000.
  if strmid(f,0,5) eq 'F625W' then w = 6317.
  if strmid(f,0,5) eq 'F702W' then w = 7020.
  if strmid(f,0,5) eq 'F775W' then w = 7764.
  if strmid(f,0,5) eq 'F814W' then w = 8140.
  if strmid(f,0,6) eq 'F850LP' then w = 8690.
  if strmid(f,0,4) eq 'F110' then w = 11000.  
  if strmid(f,0,4) eq 'F150' then w = 15000.    ; These are not precise.  I really need to update them.
  if strmid(f,0,4) eq 'F160' then w = 16000.  
  if strmid(f,0,5) eq 'F222M' then w = 22200.

  if strmid(f,0,3) eq '3.6' then w = 36000.
  if strmid(f,0,3) eq '4.5' then w = 45000. 
  if strmid(f,0,3) eq '5.8' then w = 58000. 
  if strmid(f,0,3) eq '8.0' then w = 80000.
  if strmid(f,0,2) eq '24'  then w =240000.

  if strmid(f,0,5) eq 'Wise1' or strmid(f,0,2) eq 'W1' then w = 34000.
  if strmid(f,0,5) eq 'Wise2' or strmid(f,0,2) eq 'W2' then w = 46000. 
  if strmid(f,0,5) eq 'Wise3' or strmid(f,0,2) eq 'W3' then w = 120000. 
  if strmid(f,0,5) eq 'Wise4' or strmid(f,0,2) eq 'W4' then w = 220000.

  if strlen(f) ge 2 then if strmid(f,0,2) eq 'LO' then w = float(strmid(f,2))



  if f eq 'none' then fwhm = 3000. ; approx

  if f eq 'White' then fwhm = 1500.
  if f eq 'white' then fwhm = 3000. ;approx, obviously
  if f eq 'clear' then fwhm = 5600. ;3000.

  if f eq 'UVW2' then fwhm = 760.
  if f eq 'UVM2' then fwhm = 510. ; Swift UVOT digest
  if f eq 'UVW1' then fwhm = 700. ; heasarc.nasa.gov/docs/swift/analysis/
  if f eq 'UVU' then fwhm = 785. 
  if f eq 'UVB' then fwhm = 975.
  if f eq 'UVV' then fwhm = 769.

  if f eq 'U' then fwhm = 526.
  if f eq 'B' then fwhm = 1008.
  if f eq 'V' then fwhm = 827.   ;Fukujita et al. 1995
  if f eq 'R' then fwhm = 1568.
  if f eq 'I' then fwhm = 1542.
  if f eq 'Y' then fwhm = 1000.

  if f eq 'u' then fwhm = 556.
  if f eq 'g' then fwhm = 1297.
  if f eq 'r' then fwhm = 1358.   ;Fukujita et al. 1995
  if f eq 'i' then fwhm = 1547.
  if f eq 'z' then fwhm = 1530.
  if f eq 'Z' then fwhm = 1530.

  if f eq 'J' then fwhm = 1620.
  if f eq 'H' then fwhm = 2510. ;Cohen et al. 2003
  if f eq 'K' then fwhm = 2620.
  if f eq 'Ks' then fwhm = 2620.
  if f eq 'Kp' or filter[0] eq "K'" then fwhm = 3510.

  
  if f eq 'F450W' then fwhm = 925.
  if f eq 'F606W' then fwhm = 1579.
  if f eq 'F814W' then fwhm = 1756.  
  if f eq 'F702W' then fwhm = 2100.
  if f eq 'F222M' then fwhm = 500.  ; guess 
  if f eq 'F110W' then fwhm = 3000. ; guess
  if strmid(f,0,4) eq 'F150' then fwhm = 4000.  
  if strmid(f,0,4) eq 'F160' then fwhm = 5000.   ; guess

  if strmid(f,0,3) eq '3.6' then fwhm = 4000.
  if strmid(f,0,3) eq '4.5' then fwhm = 6000. 
  if strmid(f,0,3) eq '5.8' then fwhm = 10000.  ; these are guesses
  if strmid(f,0,3) eq '8.0' then fwhm = 12000.
  if strmid(f,0,2) eq '24'  then fwhm = 80000. ; guess ><

  if strmid(f,0,5) eq 'Wise1' or strmid(f,0,2) eq 'W1' then fwhm = 12000.
  if strmid(f,0,5) eq 'Wise2' or strmid(f,0,2) eq 'W2' then fwhm = 14000. 
  if strmid(f,0,5) eq 'Wise3' or strmid(f,0,2) eq 'W3' then fwhm = 100000. 
  if strmid(f,0,5) eq 'Wise4' or strmid(f,0,2) eq 'W4' then fwhm = 90000.

  if strlen(f) ge 2 then if strmid(f,0,2) eq 'LO' then fwhm = float(strmid(f,2))*0.01

  return, w
  
end
