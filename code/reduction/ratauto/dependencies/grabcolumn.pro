function grabcolumn, infile, colnum, extstr=extstr, str=str, skiplines=skiplines
  
nlines = countlines(infile)
str = keyword_set(str)
if str then arr = strarr(nlines) else arr = dblarr(nlines)
if n_elements(skiplines) eq 0 then skiplines = 0

iline = '' ; I hate IDL
ndat = 0

openr, 1, infile
for i = 0, nlines-1 do begin
  readf, 1, iline
  if i lt skiplines then continue
  if strmid(iline, 0, 1) eq '#' then continue
  if n_elements(extstr) eq 0 then $
     idata = strtrim(strsplit(iline, /extract),2) $
  else $
     idata = strtrim(strsplit(iline, extsr, /extract),2)
  if n_elements(idata) le colnum then continue
  if str then arr[ndat] = idata[colnum] else arr[ndat] = double(idata[colnum])
  ndat = ndat + 1
endfor  
close, 1

return, arr[0:ndat-1]

end
