function recolorlris, photin

   phot = photin

;  phot[l].Uj +=  0.00 * (photin[l].Uj-photin[l].Bj)
   phot[l].Bj +=  0.00 * (photin[l].Bj-photin[l].Vj)
   phot[l].Vj +=  0.00 * (photin[l].Vj-photin[l].Rc)
   phot[l].Rc +=  0.00 * (photin[l].Vj-photin[l].Rc)
   phot[l].Ic +=  0.30 * (photin[l].Rc-photin[l].Ic) ; based on a single star pair...

   phot[l].us +=  0.00 * (photin[l].Uj-photin[l].Bj)
   phot[l].gs +=  0.00 * (photin[l].Bj-0.14)-(photin[l].Vj-0.02)) ; zpt
;  phot[r].rs += 
;  phot[i].is += 
   phot[l].zs +=  0.00 * ((photin[l].Rc+0.17)-(photin[l].Ic+0.44))


end

