function connectspec, spec1in, spec2in, attachwav=attachwav

   ; connect two spectra with different wavelength ranges at the attachment point.
    
   spec1 = spec1in
   spec2 = spec2in
   if median(spec1.wav) gt median(spec2.wav) then begin
      spec2hold = spec2
      spec2 = spec1
      spec1 = spec2hold
   endif

   w1below = where(spec1.wav le attachwav, ct1)
   w2above = where(spec2.wav gt attachwav, ct2)
   n = ct1 + ct2

   spec = replicate({wav:0.0d, flux:0.0, sky:0.0, unc:0.0, pix:0L, resp:0.0},n)

   wav = [spec1[w1below].wav , spec2[w2above].wav]
   wsbelow = where(wav le attachwav)
   wsabove = where(wav gt attachwav)

   spec[wsbelow] = spec1[w1below]
   spec[wsabove] = spec2[w2above]

   return, spec

end

