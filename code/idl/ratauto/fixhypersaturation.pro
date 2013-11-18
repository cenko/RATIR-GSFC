function fixhypersaturation, array, deadthreshold=deadthreshold

    truesaturation = max(array)

   ; correct for hypersaturation of an amplifier.


     ;mwrfits, array, 'orig.fits', /create
     s = size(array)
     nx = s[1]
     ny = s[2]
     for iter = 0, 5 do begin
       ; use y since not rotated yet
       diffl = array - shift(array,0,+1)  ; difference from pixel to the left (down)
       diffl[*,0] = 0
       diffr = array - shift(array,0,-1)  ; difference from pixel to the right (up)
       diffr[*,ny-1] = 0 

       ; if the pixel value is dropped by more than ~10000 from the left or the right and
       ; has minimal flux it is probably oversaturated.
                              ; ~ zero on the right chip, 1430 on the left chip.
       hypersaturated = where(array lt deadthreshold and $
                             (diffl lt -truesaturation*0.15 $
                           or diffr lt -truesaturation*0.15), ct)
                       ; this might trigger on cosmic rays... should probably check there are lots in a 
                       ; column

       if ct gt 0 then array[hypersaturated] = truesaturation+1
     endfor

     ;mwrfits, array, 'ocor1.fits', /create

     for y = 0, ny-1 do begin
        col = array[*,y]
        ww = where(col eq truesaturation+1, nhypersat)
        if nhypersat gt 5 then begin ; need at least a few, otherwise could just be bad pixels
                                    ; or cosmic rays or something
          sat = where(col gt truesaturation * 0.95 and col le truesaturation, ctsat)
          if ctsat eq 0 then continue
          hypersat = where(col eq truesaturation+1, cthypersat)
          nx = n_elements(col)
          for s = 0, 1 do begin
             if s eq 0 then x = min(hypersat)
             if s eq 1 then x = max(hypersat)
             ; start at the edge of a hypersaturated block (this assumes there are only two - 
             ;   min will always be in the first and max will always be in the second.)
             ; Then run from the hypersaturated block towards the edges until a regular-saturated
             ;   block is encountered, setting all the pixels to hypersaturation.
             while 1 do begin
               col[x] = truesaturation+1
               x += 1
               if x eq nx then break   ;or x eq chipboundaryl
               if col[x] gt truesaturation*0.95 and col[x] le truesaturation then break
               ; stop when it reaches the edge of the chip or the regularly-saturated block.
             endwhile
             if s eq 0 then x = min(hypersat)
             if s eq 1 then x = max(hypersat)
             while 1 do begin
               col[x] = truesaturation+1
               x -= 1
               if x lt 0 then break  ;or x eq chipboundaryr
               if col[x] gt truesaturation*0.95 and col[x] le truesaturation then break
               ; stop when it reaches the edge of the chip or the regularly-saturated block.
             endwhile
          endfor
          array[*,y] = col
       endif
     endfor

     ;mwrfits, array, 'ofinal.fits', /create

    return, array

end

