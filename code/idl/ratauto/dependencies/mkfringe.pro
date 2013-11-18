pro mkfringe, filenames, outfile=outfile, scale=scale, removeobjects=removeobjects, objthresh=objthresh, objbuffer=objbuffer
; removeobjects no longer works - consider unifying this system with flatcombine...?

  ; scale specifies the expected fringe amplitude explicitly - else use median sky level

  ; Make a fringe frame

  if n_elements(objthresh) eq 0 then objthresh = 6
  if n_elements(objbuffer) eq 0 then objbuffer = 12

  if n_elements(filenames) eq 1 then begin
    filename = filenames

    if strpos(filename,'.fit') eq -1 then begin
      if strpos(filename, '.txt') eq -1 and $
         strpos(filename,'.lis') eq -1 then begin
        print, 'File must be fits (.fits) or a list (.txt, .lis)'
        return
      endif
      filenames = strarr(200)
      i = 0
      openr, 1, filename
      iline = ''
      while not eof(1) do begin
        readf, 1, iline
        if strmid(iline,0,1) ne '#' then  filenames[i] = iline
        i = i + 1
      endwhile
      close, 1
      filenames = filenames[0:i-1]
    endif
  endif

   nim = n_elements(filenames)
   im = mrdfits(filenames[0], 0, h, /silent, /fscale)
   nx = (size(im)) [1]
   ny = (size(im)) [2]
   imstack = fltarr(nx,ny,nim)

   maskobjects, im, /nan

   fringeim = im*0

   ; older sigma clipping system
   
   for f = 0, nim-1 do begin
     if f gt 0 then im = mrdfits(filenames[f], 0, h, /silent, /fscale)
        immed = median(im)
        imclipmed = median(im[where(im lt 2*immed,ctclip)])

        if n_elements(scale) gt 0 then imscale=scale[f] else imscale=imclipmed


        imstack[*,*,f] = (im-imclipmed)/imscale
     ;print, im[f], immed, imclipmed, nx*ny-ctclip
   endfor
 
   for x = 0, nx-1 do begin
     for y = 0, ny-1 do begin
        pixstack = imstack[x,y,*]
        pixmed = median(pixstack)
        pixstdev = stdev(pixstack[where(finite(pixstack))])
        pixclipmed = median(pixstack[where(pixstack lt pixmed + 3*pixstdev)])
        fringeim[x,y] = pixclipmed
     endfor
   endfor

   for f = 0, n_elements(filenames)-1 do begin
      sxaddpar, h, 'FRINGE'+strtrim(string(f),2), filenames[f]
   endfor


   if n_elements(outfile) eq 0 then outfile = 'fringe.fits'
   writefits, outfile, fringeim, h

end

