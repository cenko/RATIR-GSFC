pro copydependencies

  readcol, 'dependencies.list', deps, format='a', /silent

  for d = 0, n_elements(deps)-1 do begin
     print, 'cp -pf ../'+deps[d]+'.pro dependencies/'
  endfor

end

