pro plotarccheck, pix, wav, arcsol, farc, brightreflinewav, morereflinewav, faintlinewav, nx=nx

      ; pix not actually used, but could be for right axis etc.

      xmin = 0
      xmax = n_elements(wav)-1

      wavrange = poly([xmin,xmax]-nx/2.,arcsol)
      print,   '  Wavelength range: ', clip(wavrange[0]), '-', clip(wavrange[1]), $
                ' (average dispersion ', clip((wavrange[1]-wavrange[0])/(xmax-xmin)), ' A/pix)' 


      !p.position = 0
      !p.multi = [0,1,3]
      ;
     
      nplot = 3
      wavlines = min(wav) + findgen(nplot+1)*(max(wav)-min(wav))/nplot
      ybot = reverse(0.1 + (0.98-0.1)*findgen(nplot)/nplot)
      ytop = ybot + abs(ybot[1]-ybot[0])*0.8

      ; main plot

      yr = [median(farc)/3. > 10., 1.02*max(farc)]
      for pi = 0, nplot-1 do begin
         xr = wavlines[[pi,pi+1]]
         !p.position = [0.05, ybot[pi], 0.98, ytop[pi]]
         plot, [0], [0], xrange=xr, yrange=yr, /xstyle, /ystyle, /ylog
         for l = 0, n_elements(faintlinewav)-1 do oplot, [faintlinewav[l], faintlinewav[l]], yr*0.5, color=2
         for l = 0, n_elements(morereflinewav)-1 do oplot, [morereflinewav[l], morereflinewav[l]], yr, color=1
         for l = 0, n_elements(brightreflinewav)-1 do oplot, [brightreflinewav[l], brightreflinewav[l]], yr, color=6, thick=2
         oplot, wav, farc
      endfor

      ; line zoom-ins

      !p.position = 0
      !p.multi=[0,ceil(n_elements(brightreflinewav)/4.),4]
      for l = 0, n_elements(brightreflinewav)-1 do begin
         plot, wav, farc, yrange=yr, /xstyle, /ystyle, /ylog, xrange=brightreflinewav[l]+[-10,+10], psym=10, xtickformat='(I)'
         for ll = 0, n_elements(faintlinewav)-1 do oplot, [faintlinewav[ll], faintlinewav[ll]], yr, color=2
         for ll = 0, n_elements(morereflinewav)-1 do oplot, [morereflinewav[ll], morereflinewav[ll]], yr, color=1
         oplot,  [brightreflinewav[l], brightreflinewav[l]], yr, color=6, thick=2
      endfor

      !p.position = 0
end


