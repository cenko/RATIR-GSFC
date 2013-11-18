function zeropt, filter

  if n_elements(filter) gt 1 then begin
    karr = strarr(n_elements(filter))
    for i = 0, n_elements(filter)-1 do begin
      karr[i] = zeropt(filter[i])
    endfor
    return, karr
  endif

  ; F = 10^(k-mag/2.5)
  ; k = log10(F0) 

  f = (strsplit(filter[0],'_',/extract)) [0]
  f = (strsplit(filter[0],"'",/extract)) [0]
  if f eq 'RC' or f eq 'Rc' then f = 'R'
  if f eq 'IC' or f eq 'Ic' then f = 'I'
  if f eq 'z0' then f = 'z'

  k = 6.0 ;default

  if f eq 'U' then k = 29-19.716
  if f eq 'B' then k = 29-19.384
  if f eq 'V' then k = 29-19.433 ;Fukujita et al. 1995
  if f eq 'R' then k = 29-19.508
  if f eq 'I' then k = 29-19.614

  if f eq 'J' then k = 9.202
  if f eq 'H' then k = 9.010 ;Cohen et al. 2003
  if f eq 'K' or f eq 'Ks' then k = 8.824

  if f eq 'u' then k = 9.54
  if f eq 'g' then k = 9.56
  if f eq 'r' then k = 9.56
  if f eq 'i' then k = 9.56
  if f eq 'z' then k = 9.57    ; or f eq 'Z'

  if f eq 'Z' then k = 9.356 ; Haislip
  if f eq 'Y' then k = 9.314 ; Haislip

  if f eq 'AB' then k = 9.56

  if strmid(f,0,1) eq 'F' then k = 9.56

  if strlen(f) ge 2 then if strmid(f,0,2) eq 'LO' then k = 9.56
  if strlen(f) gt 2 then if strmid(f,0,2) eq 'LR' then k = 9.56

  if f eq 'none' then k = 29-19.508  ; treat as R band
  if f eq 'clear' then k = 29-19.508 ;
  if f eq 'white' then k = 29-19.508 ;

  if f eq 'UVW2' then k = 8.8258
  if f eq 'UVM2' then k = 8.8834 ;calckuvot3.pro
  if f eq 'UVW1' then k = 8.9512
  if f eq 'White' then k = 9.3783
  if f eq 'UVU' then k = 9.1508
  if f eq 'UVB' then k = 9.6198  ; B and V are consistent with Johnson
  if f eq 'UVV' then k = 9.5715

  ;if f eq 'UVW2' then k = 8.8563
  ;if f eq 'UVM2' then k = 8.8232 ;calckuvot2.pro
  ;if f eq 'UVW1' then k = 8.9939
  ;if f eq 'White' then k = 9.151
  ; UVOT U and B are also discrepant from JC system by ~0.05 mag.

  ;if f eq 'UVW2' then k = 8.866
  ;if f eq 'UVM2' then k = 8.905 ;calckuvot.pro
  ;if f eq 'UVW1' then k = 8.861
  ;if f eq 'White' then k = 29-19.433 ; use V-band?

  if f eq '3.6' then k = 8.448
  if f eq '4.5' then k = 8.254
  if f eq '5.8' then k = 8.061
  if f eq '8.0' then k = 7.812
  ;if f eq '3.6' then k = 9.56  ; Spitzer no longer as AB mags
  ;if f eq '4.5' then k = 9.56
  ;if f eq '5.8' then k = 9.56
  ;if f eq '8.0' then k = 9.56
  if f eq '24' then  k = 9.56

  if f eq 'W1' then k = 8.49 
  if f eq 'W2' then k = 8.23
  if f eq 'W3' then k = 7.50
  if f eq 'W4' then k = 6.92

  return, k                          
end
