pro lrisplotall

  all = findfiles('lris*.spec')
  all = all[sort(all)]

  psopen, 'lrisspectra.ps', xsize=11, ysize=8.5, /inches
  !p.multi = [0,1,4]
  !p.position = 0
  device, /color
  colors = transpose([[0,0,0],$
                      [180,180,210],$ 
                      [64, 64, 192],$ 
                      [64, 136,128],$ 
                      [64, 158, 64],$  
                      [112,112, 64],$ 
                      [192, 64, 64],$
                      [128, 64, 128]$
                       ])
  tvlct, colors

  for f = 0, n_elements(all)-1 do begin
    filename = all[f]
    s = readspec(filename, h)

    targname = sxpar(h,'TARGNAME')
    targname = clip((strsplit(targname,'_',/extract)) [0])
    print, targname, ' (', removepath(filename), ')'
    if strpos(targname,'BD') ge 0 then continue
    if strpos(targname,'HD') ge 0 then continue
    if strpos(targname,'Feige') ge 0 then continue
    if strpos(targname,'Kopff') ge 0 then continue

    nx = n_elements(s)
    wn = where(s.wav gt 3300 and s.wav lt 10000 and s.wav lt max(s.wav)*0.95 and s.wav gt min(s.wav) * 1.05)
    xrange = [3100, 10350]
    yrange = [0,1.05*max(s[wn].flux)]
    plot, s.wav, s.flux, xrange=xrange, yrange=yrange, /xstyle, /ystyle
    oplot, [5000,12000], [0,0]

    xyouts, 8500, yrange[1]*0.85, targname
  endfor

  psclose

end
