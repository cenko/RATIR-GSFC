pro joinlrisspace, filename, nobias=nobias, linebias=linebias, gaindata=gaindata, skysubtract=skysubtract

;convert a multi-HDU blue-side file into single-HDU, like the old format, without any processing.

if n_elements(filename) gt 1 then begin
  for f = 0, n_elements(filename)-1 do begin
    print, filename[f]
    joinlris, filename[f]
  endfor
  return
endif

filearr = strsplit(filename,'.', /extract)
fileroot = filearr[0]
fileext = filearr[1]

if fileext eq 'cat' or fileext eq 'lis' or fileext eq 'list' or fileext eq 'txt' then begin
  files = grabcolumn(filename,0,/str)
  joinlris, files
  return
endif

header = headfits(filename, /silent)
sxaddpar, header, 'BZERO', 0

array1 = transpose(lris_read_amp(filename, 1, nobias=nobias, linebias=linebias, gaindata=gaindata))
array2 = transpose(lris_read_amp(filename, 2, nobias=nobias, linebias=linebias, gaindata=gaindata))

s1 = size(array1)
nx1 = s1[1]
ny1 = s1[2]
s2 = size(array2)
nx2 = s2[1]
ny2 = s2[2]

if nx1 ne nx2 then print, 'Arrays do not match!' 

if keyword_set(skysubtract) then begin
   nblock = 5
   blocksize = fix((ny2-20)/nblock)
   for x = 0, nx2-1 do begin
    for block = 0, nblock-1 do begin
      ymin = block*blocksize
      ymax = block*blocksize + blocksize-1
      array2[x,ymin:ymax] = array2[x,ymin:ymax] - median(array2[x,ymin:ymax])
      ;if x eq 1445 then print, ymin, ymax, median(array2[x,ymin:ymax])
    endfor
   endfor
   nblock = 5
   blocksize = fix((ny1-37)/nblock)
   for x = 0, nx1-1 do begin
    for block = 0, nblock-1 do begin
      ymin = 37 + block*blocksize
      ymax = 37 + block*blocksize + blocksize-1
      array1[x,ymin:ymax] = array1[x,ymin:ymax] - median(array1[x,ymin:ymax])
    endfor
   endfor
endif


; "The red CCDs are slightly tilted with respect to one another. The gap, specifically the distance between the first column of the two red CCDs, at pixel (0,4096) is 24.96 arcsec and at pixel (0,0) is 27.91 arcsec. "
if ny1 lt 400 then pxscale = 0.135*2 else pxscale = 0.135

gap = fix(round(   ((24.96 + 27.91)/2) / pxscale  ))  ; pixels

array = fltarr(nx1, ny1+gap+ny2)
array[*,0:ny1-1]   = array1
array[*,ny1+gap:*] = array2


mwrfits, array, 'js'+fileroot+'.fits',  header, /create

end
