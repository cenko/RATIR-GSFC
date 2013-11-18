pro namefix, header, tables, radius=radius

  ; Correct incorrect object names in a FITS header using a positional catalog (specified in tables).
  ; If the RA/DEC position is within radius arcsec of anything in the catalog, the catalog name
  ; becomes the object name.

  common storeref, prevtables, ref

  if n_elements(radius) eq 0 then radius = 300.

  reload = 0
  if n_elements(prevtables) ne n_elements(tables) then reload = 1
  if n_elements(prevtables) eq n_elements(tables) then if total(prevtables ne tables) gt 0 then reload = 1
  if reload then begin
     for t = 0, n_elements(tables)-1 do begin
       if file_test(tables[t]) eq 0 then begin
         print, tables[t], ' not found!'
       endif
     endfor
     refname = ['']
     refra = [0.]
     refdec = [0.]
     reftype = ['']
     for t = 0, n_elements(tables)-1 do begin
       if file_test(tables[t]) eq 0 then continue
       openr, 1, tables[t]
       inline = ''
       while not eof(1) do begin
          readf, 1, inline
          inline = clip(inline)
          if strlen(inline) lt 2 then continue
          if strmid(inline,0,1) eq '#' then continue
          inarr = strsplit(inline,/extract)
          if n_elements(inarr) lt 3 then continue
          refname = [refname, inarr[0]]
          if strpos(inarr[1],':') gt 0 then inra = 15*ten(inarr[1]) else inra  = float(inarr[1])
          if strpos(inarr[2],':') gt 0 then indec =   ten(inarr[2]) else indec = float(inarr[2])
          refra   = [refra,  inra]
          refdec  = [refdec, indec]
          if n_elements(inarr) gt 3 then  inreftype = inarr[3] else inreftype = ''
          reftype = [reftype, inreftype]
       endwhile
       close, 1
     endfor     
     nref = n_elements(refra)-1
     ref = replicate({name:'', ra:0., dec:0., type:''},nref)
     ref.name = refname[1:*]
     ref.ra = refra[1:*]
     ref.dec = refdec[1:*]
     ref.type = reftype[1:*]

     prevtables = tables
  endif else begin
     nref = n_elements(ref)
  endelse

  obj = ''
  objtype = ''

  rastr = strtrim(sxpar(header,'RA'),2)   ; should do something to allow for decimal input
  decstr = strtrim(sxpar(header,'DEC'),2)
  if strpos(rastr,':') gt 0 then ra = ten(rastr)*15  else ra = float(rastr) ; hours to degrees
  if strpos(decstr,':') gt 0 then dec = ten(decstr)  else dec = float(decstr)

  distance = fltarr(nref)

  for r = 0, nref-1 do begin
    dist = 3600*(abs(dec-ref[r].dec));
    if dist lt radius then gcirc, 2, ra, dec, ref[r].ra, ref[r].dec, dist
    distance[r] = dist
  endfor
  mindist = min(distance)
  minr = (where(distance eq mindist)) [0]
  if mindist lt radius then begin
      obj = ref[minr].name
      objtype = ref[minr].type
  endif

  trapdoor = clip(sxpar(header,'TRAPDOOR'))
  lampstr = clip(sxpar(header,'LAMPS'))
  lamps = fix(strsplit(lampstr,',',/extract))
  if obj ne '' then begin
      arc = total(lamps[0:4])
      hal = lamps[5]
      if arc gt 0 then obj += ' (arc)'
      if hal gt 0 then obj += ' (hal)'
      if arc eq 0 and hal eq 0 then begin
        if trapdoor eq 'closed' then obj += ' (closed)'
      endif
  endif


  origobj = sxpar(header, 'OBJECT', /silent)
  sxaddpar, header, 'ORIGOBJ', origobj
  if obj ne origobj and obj ne '' then begin
     print, ' ['+clip(origobj)+' -> '+obj+'] ', format='($,A)'
     sxaddpar, header, 'OBJECT', obj 
  endif
  ;if objtype ne '' then sxaddpar, header, 'OBJTYPE', objtype  ; not sure what this ever did

end

