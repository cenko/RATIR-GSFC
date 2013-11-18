function readstandardtable, filename

   ; reads in a structure with information about the spectroscopic standard stars

   common specstdblock, standardtable

   if n_elements(standardtable) eq 0 then begin
     nl = countlines(filename)
     standardtable = replicate({name:'', rastr:'', decstr:'', ra:0.D, dec:0.D, mag:0., type:'', file:'', altnames:''},nl)
     openr, 11, filename
     inline = ''
     i = 0L
     for l = 0L, nl-1 do begin
       readf, 11, inline
       inline = clip(inline)
       if strlen(inline) eq 0 then continue
       if strmid(inline,0,1) eq '#' then continue
       indata = strsplit(inline,' ',/extract)
       if n_elements(indata) le 6 then continue       
     
       standardtable[i].name = indata[0]
       standardtable[i].rastr = indata[1]
       standardtable[i].decstr = indata[2]
       standardtable[i].ra = 15*ten(indata[1])
       standardtable[i].dec = ten(indata[2])
       if valid_num(indata[3]) eq 0 then standardtable[i].mag = !values.f_nan $
                                    else standardtable[i].mag = float(indata[3])
       standardtable[i].type = indata[4]
       standardtable[i].file = indata[5]
       if n_elements(indata) gt 6 then $
          standardtable[i].altnames = indata[6]
       i += 1
     endfor
     standardtable = standardtable[0:i-1]
   endif

   return, standardtable

end

