pro rewindowspec, filename

  ; window a full-frame LRIS-R1 image to be like a spectrum-mode image.
  ; Not relevant to upgraded versions of LRIS.
  ; edited 2011-03-20

nfile = n_elements(filename)
for f = 0, nfile-1 do begin
   data = mrdfits(filename[f],0,h, /silent)

   cropdata = data[*,498:1498-1]

   ;sxaddpar, h, 'WINDOW','0,0,0,550,1000'
   sxaddpar, h, 'WINDOW','0,0,550,2048,1000' ; they changed this?!
   filearr = strsplit(filename[f],'/',/extract)
   filenam = filearr[n_elements(filearr)-1]
   outfilenam = 'w'+filenam
   mwrfits, cropdata, outfilenam, h, /create
endfor



end
