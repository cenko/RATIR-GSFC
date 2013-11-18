pro measuregain, filenames, oldblue=oldblue

; use a raw file (flat field) to measure the gain and read noise of each chip individually.

for i = 0, n_elements(filenames)-1 do begin

if keyword_set(oldblue) then begin
   amp1x = [205-1, 1228-1]
   amp2x = [1229-1,2252-1]
   amp3x = [2253-1,3276-1]
   amp4x = [3277-1,4300-1]

   bias1x = [4302,4379]
   bias2x = [4382,4457]
   bias3x = [4462,4537]
   bias4x = [4542,4618]

   data=float(readfits(filenames[i],fhead,/silent))
   avbias1 = median(data[bias1x[0]:bias1x[1],*])
   avbias2 = median(data[bias2x[0]:bias2x[1],*])
   avbias3 = median(data[bias3x[0]:bias3x[1],*])
   avbias4 = median(data[bias4x[0]:bias4x[1],*])

   nx = (size(data)) [1]
   ny = (size(data)) [2]

   for y=0,ny-1 do begin
      data[amp1x[0]:amp1x[1],y] = ( data[amp1x[0]:amp1x[1],y] - median(data[bias1x[0]:bias1x[1],y]) ) 
      data[amp2x[0]:amp2x[1],y] = ( data[amp2x[0]:amp2x[1],y] - median(data[bias2x[0]:bias2x[1],y]) ) 
      data[amp3x[0]:amp3x[1],y] = ( data[amp3x[0]:amp3x[1],y] - median(data[bias3x[0]:bias3x[1],y]) ) 
      data[amp4x[0]:amp4x[1],y] = ( data[amp4x[0]:amp4x[1],y] - median(data[bias4x[0]:bias4x[1],y]) ) 
      data[bias1x[0]:bias1x[1],y] = ( data[bias1x[0]:bias1x[1],y] - median(data[bias1x[0]:bias1x[1],y]) ) 
      data[bias2x[0]:bias2x[1],y] = ( data[bias2x[0]:bias2x[1],y] - median(data[bias2x[0]:bias2x[1],y]) ) 
      data[bias3x[0]:bias3x[1],y] = ( data[bias3x[0]:bias3x[1],y] - median(data[bias3x[0]:bias3x[1],y]) ) 
      data[bias4x[0]:bias4x[1],y] = ( data[bias4x[0]:bias4x[1],y] - median(data[bias4x[0]:bias4x[1],y]) ) 
   endfor

    mwrfits, data, 'bluegain.fits', /create
    med1 = median(data[amp1x[1]-50:amp1x[1],1000:3000])
    med2 = median(data[amp2x[0]:amp2x[0]+50,1000:3000])
    medl = median(data[amp2x[0]:amp2x[0]+1000,1000:3000])
    medr = median(data[amp3x[1]-1000:amp3x[1],1000:3000])
    med3 = median(data[amp3x[1]-50:amp3x[1],1500:2500])
    med4 = median(data[amp4x[0]:amp4x[0]+50,1500:2500])

    print, med1, med2, med3, med4
    invgain = [(med1/med2)*(medl/medr), medl/medr, med3/med3, med4/med3]
    print, invgain
    print, 1./invgain

    openw, 1, 'bluegain.dat'
    printf, 1, 1./invgain
    close, 1

   return
endif





fullarray = readmhdufits(filenames[i], header=header, /notrim)

if sxpar(header,'INSTRUME') eq 'LRISBLUE' then begin
    array = fullarray[204:4299,0:4095]

    amplx = 1025
    amprx = 3072

    med1 = median(array[amplx-50:amplx-1,1000:3000])
    med2 = median(array[amplx+1:amplx+50,1000:3000])
    medl = median(array[amplx+1:amplx+1000,1000:3000])
    medr = median(array[amprx-1000:amprx-1,1000:3000])
    med3 = median(array[amprx-50:amprx-1,1500:2500])
    med4 = median(array[amprx+1:amprx+50,1500:2500])

    print, med1, med2, med3, med4
    invgain = [(med1/med2)*(medl/medr), medl/medr, med3/med3, med4/med3]
    print, invgain
    print, 1./invgain

    openw, 1, 'bluegain.dat'
    printf, 1, 1./invgain
    close, 1

    ;array[0:1024,*] = array[0:1024,*] *gainout/gain1
    ;array[1025:2048-1,*] = array[1025:2048-1,*] *gainout/gain1
    ;array[2049-1:3072-1,*] = array[2049-1:3072-1,*]*gainout/gain3
    ;array[3073-1:*] = array[3073-1:*]*gainout/gain4
endif
if strtrim(sxpar(header,'INSTRUME'),2) eq 'LRIS' then begin
    array = fullarray[48:3447,0:2499]

    amplx = 624
    amprx = 2672

    med1 = median(array[amplx-50:amplx-1,600:2300])
    med2 = median(array[amplx+1:amplx+50,600:2300])
    medl = median(array[amplx+1:amplx+800,600:2300])
    medr = median(array[amprx-800:amprx-1,600:2300])
    med3 = median(array[amprx-50:amprx-1,600:2300])
    med4 = median(array[amprx+1:amprx+50,600:2300])

    print, med1, med2, med3, med4
    invgain = [(med1/med2)*(medl/medr), medl/medr, med3/med3, med4/med3]
    print, invgain
    print, 1./invgain

    openw, 1, 'redgain.dat'
    printf, 1, 1./invgain
    close, 1

    ; gain correction
   ; array[0:623,*] = array[0:623,*] / 0.85714
   ; array[624:1647,*] = array[624:1647,*] / 0.796040
   ;array[1648:2671,*] = array[1648:2671,*] 
   ; array[2672:*,*] = array[2672:*,*] / 1.02458

endif



endfor

;chipmed = fltarr(4)
;
;for c=1,4 do begin
;  chip = mrdfits, infile, c, h, /fscale
;  chip = mrdfits, infile, c, h, /fscale
;  chip = mrdfits, infile, c, h, /fscale
;  chip = mrdfits, infile, c, h, /fscale
;
;  illum = where(chip1 gt median(chip1)/2.)
;  illum = where(chip1 gt median(chip1[illum1])/2.)
;
;  chipmed = median(chipmed)



end
