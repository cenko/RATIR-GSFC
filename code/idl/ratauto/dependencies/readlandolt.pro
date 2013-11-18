function flt, inpval, defval
   c = strtrim(inpval,2)
   if c eq '' then return, defval
   return, float(inpval)
end


function getstarmag, landolt, filter, unc=unc, nickel=nickel
   unc = -999
   if filter eq 'U' then unc = landolt.eUj
   if filter eq 'B' then unc = landolt.eBj
   if filter eq 'V' then unc = landolt.eVj
   if filter eq 'R' then unc = landolt.eRc
   if filter eq 'I' then unc = landolt.eIc
   if filter eq 'u' then unc = landolt.eus
   if filter eq 'g' then unc = landolt.egs
   if filter eq 'r' then unc = landolt.ers
   if filter eq 'i' then unc = landolt.eis
   if filter eq 'z' then unc = landolt.ezs
   if keyword_set(nickel) then begin
      if filter eq 'B' then return, landolt.Bj - 0.06*(landolt.Bj - landolt.Vj)
      if filter eq 'g' then return, landolt.gs - 0.23*(landolt.Bj - landolt.Vj  - 0.741)
      if filter eq 'V' then return, landolt.Vj + 0.09*(landolt.Bj - landolt.Vj  - 0.741)
      if filter eq 'R' then return, landolt.Rc + 0.10*(landolt.Vj - landolt.Rc  - 0.426)
      if filter eq 'I' then return, landolt.Ic - 0.04*(landolt.Vj - landolt.Ic  - 0.833)
   endif
   if filter eq 'U' then return, landolt.Uj
   if filter eq 'B' then return, landolt.Bj
   if filter eq 'V' then return, landolt.Vj
   if filter eq 'R' then return, landolt.Rc
   if filter eq 'I' then return, landolt.Ic
   if filter eq 'u' then return, landolt.us
   if filter eq 'g' then return, landolt.gs
   if filter eq 'r' then return, landolt.rs
   if filter eq 'i' then return, landolt.is
   if filter eq 'z' then return, landolt.zs
   return, -999
end


function readlandolt, filename=filename

;#   1- 11 A11    ---      Star   Star name
;#  13- 14 I2     h        RAh    Hour of Right Ascension (J2000) 
;#  16- 17 I2     min      RAm    Minute of Right Ascension (J2000) 
;#  19- 24 F6.3   s        RAs    Second of Right Ascension (J2000) 
;#      26 A1     ---      DE-    Sign of the Declination (J2000)
;#  27- 28 I2     deg      DEd    Degree of Declination (J2000) 
;#  30- 31 I2     arcmin   DEm    Arcminute of Declination (J2000) 
;#  33- 37 F5.2   arcsec   DEs    Arcsecond of Declination (J2000) 
;#  39- 44 F6.3   mag      Vmag   The V band magnitude
;#  46- 51 F6.3   mag      B-V    The (B-V) color
;#  53- 58 F6.3   mag      U-B    The (U-B) color
;#  60- 65 F6.3   mag      V-R    The (V-R) color
;#  67- 72 F6.3   mag      R-I    The (R-I) color
;#  74- 79 F6.3   mag      V-I    The (V-I) color
;#  81- 83 I3     ---      n      Number of times star was observed
;#  85- 87 I3     ---      m      Number of nights star was observed
;#  89- 94 F6.4   mag    e_Vmag   ? Mean Error of the Mean of Vmag 
;#  96-101 F6.4   mag    e_B-V    ? Mean Error of the Mean of B-V
;# 103-108 F6.4   mag    e_U-B    ? Mean Error of the Mean of U-B
;# 110-115 F6.4   mag    e_V-R    ? Mean Error of the Mean of V-R
;# 117-122 F6.4   mag    e_R-I    ? Mean Error of the Mean of R-I
;# 124-129 F6.4   mag    e_V-I    ? Mean Error of the Mean of V-I

;# PG0029+024  00 31 42.20  +02 37 44.3  15.268 +0.362 -0.184 +0.251 +0.337 +0.593   5   2 0.0094 0.0174 0.0112 0.0161 0.0125 0.0067 

n = 595
landolt = replicate({name:'', ra:0.D, dec:0.D, rastr:'', decstr:'', Uj:0., Bj:0., Vj:0., Rc:0., Ic:0., us:0., gs:0., rs:0., is:0., zs:0., $
                                                                 euJ:0.,eBj:0.,eVj:0.,eRc:0.,eIc:0.,eus:0.,egs:0.,ers:0.,eis:0.,ezs:0., n:0, m:0}, n)

if n_elements(filename) eq 0 then filename = '~/progs/idl/phot/landolt2009.txt'
openr, 1, filename

iline = ''
l = 0l
while not eof(1) do begin
   readf, 1, iline
   if strmid(iline,0,1) eq '#' then continue
   if strlen(iline) lt 120 then continue

   landolt[l].name = strmid(iline,0,11)
   landolt[l].rastr = strmid(iline,12,10)
   landolt[l].decstr = strmid(iline,25,10)
   landolt[l].ra = 15*ten(landolt[l].rastr)
   landolt[l].dec = ten(landolt[l].decstr)
   landolt[l].Vj = flt(strmid(iline,38,6))
   landolt[l].Bj = landolt[l].Vj + flt(strmid(iline,45,6))
   landolt[l].Uj = landolt[l].Bj + flt(strmid(iline,52,6))
   landolt[l].Rc = landolt[l].Vj - flt(strmid(iline,59,6))
   landolt[l].Ic = landolt[l].Vj - flt(strmid(iline,73,6))
   landolt[l].n  = fix(strmid(iline,80,3))
   landolt[l].m  = fix(strmid(iline,84,3))
   landolt[l].eVj = flt(strmid(iline,88,6),0.1)
   landolt[l].eBj = sqrt(landolt[l].eVj^2 + flt(strmid(iline,95,6),0.173205))
   landolt[l].eUj = sqrt(landolt[l].eBj^2 + flt(strmid(iline,102,6),0.282843))  
   landolt[l].eRc = sqrt(landolt[l].eVj^2 + flt(strmid(iline,109,6),0.173205))
   landolt[l].eIc = sqrt(landolt[l].eVj^2 + flt(strmid(iline,123,6),0.173205))

    ; Jester 2005 transform

   landolt[l].gs = landolt[l].Vj + 0.60*(landolt[l].Bj-landolt[l].Vj) - 0.12  ;  0.02
   landolt[l].rs = landolt[l].Vj - 0.42*(landolt[l].Bj-landolt[l].Vj) + 0.11  ;  0.03
   landolt[l].us = landolt[l].gs + 1.28*(landolt[l].Uj-landolt[l].Bj) + 1.13  ;  0.06
   landolt[l].is = landolt[l].rs - 0.91*(landolt[l].Rc-landolt[l].Ic) + 0.20  ;  0.03
   landolt[l].zs = landolt[l].rs - 1.72*(landolt[l].Rc-landolt[l].Ic) + 0.41  ;  0.03

   landolt[l].egs = sqrt( landolt[l].eVj^2 + 0.60*flt(strmid(iline,95,6),0.173205)^2 + 0.02^2 )
   landolt[l].ers = sqrt( landolt[l].eVj^2 + 0.42*flt(strmid(iline,95,6),0.173205)^2 + 0.03^2 )
   landolt[l].eus = sqrt( landolt[l].egs^2 + 1.28*flt(strmid(iline,102,6),0.282843)^2 + 0.06^2 )
   landolt[l].eis = sqrt( landolt[l].ers^2 + 0.91*flt(strmid(iline,116,6),0.173205)^2 + 0.03^2 )
   landolt[l].ezs = sqrt( landolt[l].ers^2 + 1.72*flt(strmid(iline,116,6),0.173205)^2 + 0.03^2 )

   l = l + 1

endwhile

close, 1
return, landolt


end

