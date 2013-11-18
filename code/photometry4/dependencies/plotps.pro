PRO PLOTPS, fileout

;set up postscript plotting
set_plot,'ps'
ofile= fileout
device,/inches,/color,bits_per_pixel=8,filename = ofile,/encapsulated

;load colors
;loadct,13,ncolors=15,bottom=0

END
