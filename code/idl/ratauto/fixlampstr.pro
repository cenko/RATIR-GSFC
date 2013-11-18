pro fixlampstr, lampstr, badarclist

      if badarclist[0] ne '' then begin
         lampstr = sxpar(header,'LAMPS')
         lamparr = fix(strsplit(lampstr,',',/extract))
         arclamps[f] = lampstr
         if total(badarclist eq 'hg' or badarclist eq 'mercury')  then lampstr =                     '0'+strmid(lampstr,1)
         if total(badarclist eq 'ne' or badarclist eq 'neon')     then lampstr = strmid(lampstr,0,2)+'0'+strmid(lampstr,3)
         if total(badarclist eq 'ar' or badarclist eq 'argon')    then lampstr = strmid(lampstr,0,4)+'0'+strmid(lampstr,5)
         if total(badarclist eq 'cd' or badarclist eq 'cadmium')  then lampstr = strmid(lampstr,0,6)+'0'+strmid(lampstr,7)
         if total(badarclist eq 'zn' or badarclist eq 'zinc')     then lampstr = strmid(lampstr,0,8)+'0'+strmid(lampstr,9)
         if total(badarclist eq 'hal' or badarclist eq 'halogen') then lampstr = strmid(lampstr,0,10)+'0'
      endif


end

