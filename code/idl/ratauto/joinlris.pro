pro joinlris, filename, notrim=notrim, nobias=nobias, linebias=linebias, gaindata=gaindata, bspec=bspec

;convert a multi-HDU blue-side file into single-HDU, like the old format, without any processing.

if keyword_set(linebias) eq 0 and keyword_set(nobias) eq 0 then print, 'Warning: use /linebias option for correct bias subtraction'

if n_elements(filename) gt 1 then begin
  for f = 0, n_elements(filename)-1 do begin
    print, filename[f]
    joinlris, filename[f], notrim=notrim, nobias=nobias, linebias=linebias, gaindata=gaindata, bspec=bspec
  endfor
  return
endif

filearr = strsplit(filename,'.', /extract)
fileroot = filearr[0]
fileext = filearr[1]

if fileext eq 'cat' or fileext eq 'lis' or fileext eq 'list' or fileext eq 'txt' then begin
  files = grabcolumn(filename,0,/str)
  joinlris, files, notrim=notrim, nobias=nobias, linebias=linebias, gaindata=gaindata, bspec=bspec
  return
endif

array = readmhdufits(filename, header=header, notrim=notrim, nobias=nobias, linebias=linebias, gaindata=gaindata)

if keyword_set(bspec) then begin
  slitlen = 1184
  array = rotate( reverse(array[1460:1460+slitlen,*],2) ,1)
endif

mwrfits, array, 'j'+fileroot+'.fits',  header, /create

end
