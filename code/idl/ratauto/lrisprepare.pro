;+
; NAME: lrisprepare.pro
;
; PURPOSE: 
;       All-purpose program to perform basic operations on an LRIS
;       image (or 2D spectrum) including bias subtraction, flagging
;       of bad pixels, header standardization, provisional astrometry,
;       etc.
;       Auto-detects whether it is red or blue, which LRIS version
;       it is, and whether it is imaging or spectroscopy.
;       This the latest, final, catch-all version.
;       It does not yet work with old-format images, but should work
;       with all modern LRIS imaging/spectroscopy regardless of amp
;       status, binning, etc.
;
; CALLING SEQUENCE: 
;       lrisprepare, filename
;
; INPUTS:
;	filename = name of the multi-HDU FITS file to prepare,
;                  or array of filenames, or file w/list of filenames
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       split       - write out left and right chips separately (default)
;       merge       - keep left and right chips together
;       rightonly 
;       leftonly
;       crops       - crop instruction
;       xcrop       - crop locations
;       lcrop       -
;       rcrop       -
;       ycrop       -
;       outname     - specify output file to write to disk
;       loutname
;       routname
;       prefix      - specify output prefix
;       lprefix     -
;       rprefix     -
;       suffix      - specify output suffix
;       lsuffix     -
;       rsuffix     -
;       nobias
;       scalarbias
;       gain        - gain correction for all four (two if LRISR1) amps
;       flag
;       verbose
;       namefix     - correct names using a catalog
;       timer       - print out processing times
;
; OPTIONAL OUTPUT KEYWORD PARAMETERS:
;       header (?)
; 
; OUTPUTS:
;       No formal outputs.
;       The prepared images are written to disk.
;
; REQUIREMENTS:
;       - Requires the IDL Astronomy User's Library routines 
;       (http://idlastro.gsfc.nasa.gov/) 
;       - Requires readmhdu.pro to process LRIS images from 2009/later.

; Note - now attempts to add physical coordinates based on binned data.
; However, there is a genuine physical shift for binned observations - when binned 2x2,
; an extra row of dead pixels is read out along the edge of the gap (on the r chip and 
; looks like on the l chip too.)  This applies to R3; I don't know if it applies to R2.

;----------------------------------------------------------------------------------------
pro lrisprepare, filename, split=split, merge=merge, rightonly=rightonly, leftonly=leftonly, $
        crops=crops, xcrop=xcrop, lcrop=lcrop, rcrop=rcrop, ycrop=ycrop, notranspose=notranspose, $
        outname=outname, loutname=loutname, routname=routname, $
        prefix=prefix, lprefix=lprefix, rprefix=rprefix, $
        suffix=suffix, lsuffix=lsuffix, rsuffix=rsuffix, $
        nobias=nobias, scalarbias=scalarbias, gain=gain, verbose=verbose, $
        fixdefect=fixdefect, flag=flag, namefixfiles=namefixfiles, header=header, timer=timer
   ; arguments repeated once

common lrisconfig, lrisautopath, refdatapath, defaultspath
common lrisbpmblock, bpm2, bpm3 
common lrisprepareprev, mask, prevcamera, prevversion, prevpane, prevxbin, prevybin
common lrisbadarc, badarcdir, badarclist

if n_elements(refdatapath) eq 0 then bpmpath = '/home/dperley/progs/idl/lris/refdata/' else bpmpath = refdatapath
if strlen(bpmpath) gt 0 then if strmid(bpmpath,strlen(bpmpath)-1,1) ne '/' then bpmpath += '/'

if keyword_set(merge) and keyword_set(split) then begin
   print, 'Cannot both split and merge.'
   return
endif
if keyword_set(merge) eq 0 and keyword_set(split) eq 0 then begin
   split = 1   
endif
merge = keyword_set(merge)
split = keyword_set(split)

if n_elements(crops) eq 0 and $
   n_elements(xcrop) eq 0 and n_elements(lcrop) eq 0 and n_elements(rcrop) eq 0 and n_elements(ycrop) eq 0 then crops = 'auto'


; ---------- Process input filename(s), call recursively if necessary----------

; Check for empty filename
if n_elements(filename) eq 0 then begin
   print, 'No filename specified.'
   return
endif

; Check for array input
if n_elements(filename) gt 1 then begin
   files = filename
endif

; Check for list input
if n_elements(filename) eq 1 then begin
   filearr = strsplit(filename,'.', /extract)
   fileroot = filearr[0]
   fileext = filearr[1]
   if fileext eq 'cat' or fileext eq 'lis' or fileext eq 'list' or fileext eq 'txt' then begin
     ;files = grabcolumn(filename,0,/str)
     readcol, filename, files, format='a'
   endif
endif


; Check for wildcard input
if n_elements(filename) eq 1 then begin
   starpos = strpos(filename,'*')
   qpos = strpos(filename,'?')
   if starpos ge 0 or qpos ge 0 then begin
      files = findfile(filename)
      if n_elements(files) eq 1 then begin
        if files[0] eq '' then print, 'Cannot find any files matching ', filename
        return
      endif
   endif
endif

; Check for array filename, loop recursively if necessary
if n_elements(files) gt 1 then begin
  for f = 0, n_elements(files)-1 do begin
    lrisprepare, files[f], split=split, merge=merge, rightonly=rightonly, leftonly=leftonly, $
        crops=crops, xcrop=xcrop, lcrop=lcrop, rcrop=rcrop, ycrop=ycrop, notranspose=notranspose, $
        outname=outname, loutname=loutname, routname=routname, $
        prefix=prefix, lprefix=lprefix, rprefix=rprefix, $
        suffix=suffix, lsuffix=lsuffix, rsuffix=rsuffix, $
        nobias=nobias, scalarbias=scalarbias, gain=gain, verbose=verbose, $
        fixdefect=fixdefect, flag=flag, namefixfiles=namefixfiles, timer=timer
  endfor
  return
endif


print, filename, format='($,A)'


; ---------- Read data and process header information ----------

; Read in the data and header

starttime = systime(/seconds) 

linebias = keyword_set(scalarbias) eq 0
;array = readmhdufitsdap(filename, header=header, nobias=nobias, linebias=linebias, verbose=verbose) ; do gain later
array = readfits(filename, header)
sin = size(array) ; input size
nxin = sin[1]
nyin = sin[2]

;; for old-style LRIS, at this point we will have to indicate where the amps are manually.
;
;if nxin eq 4096 and nyin eq 512 and sxpar(header,'INSTRUME') eq 'LRISBLUE' and clip(sxpar(header,'SLITNAME')) eq 'direct' then begin
;   print, ' - MIRA file; skipping.'
;   return
;endif


readtime = systime(/seconds)  ; mark time after file read from disk with readmhdu

; Get essential info about the data structure: version, camera, binning, windowing, lmode

;mjd = float(sxpar(header,'MJD-OBS'))
;if mjd gt 0     and mjd lt 54952 then lrisversion = 1
;if mjd ge 54952 and mjd lt 55562 then lrisversion = 2
;if mjd ge 55562 then lrisversion = 3
;if n_elements(lrisversion) eq 0 then begin
;   print, '[No MJD entry in header: assuming LRIS3]', format='($,A)'
;   lrisversion = 3
;endif

;instrume = clip(sxpar(header,'INSTRUME'))
camera = 'b'                    ;placeholder
;if instrume eq 'LRIS'     then camera = 'r'
;if instrume eq 'LRISBLUE' then camera = 'b'
;if camera eq 'r' then sxaddpar, header, 'CAMERA', 'red'
;if camera eq 'b' then sxaddpar, header, 'CAMERA', 'blue'
;if camera eq '' then begin
;   print, ': unrecognized camera (INSTRUME = ', instrume, ')'
;   return
;endif

sxaddpar, header, 'CAMERA', 'blue' ;placeholder
lrisversion = 3                    ;placeholder
;sxaddpar, header, 'LRISVER', lrisversion ;placeholder

xbin = 1                        ;placeholder
ybin = 1                        ;placeholder
;panestr = clip(sxpar(header,'PANE'))
;pane = fix(strsplit(panestr,',',/extract))
;binning = clip(sxpar(header,'BINNING'))
;bin = fix(strsplit(binning,',',/extract))
;xbin = bin[0]
;ybin = bin[1]
sxaddpar, header, 'XBIN', xbin, after='BINNING'  ; actually Y, but we rotate it latter for spectroscopy...
sxaddpar, header, 'YBIN', ybin, after='BINNING'
;;sxaddpar, header, 'CRPIX1', 0
;;sxaddpar, header, 'CRPIX2', 0
;;sxaddpar, header, 'CRVAL1', 0          ; doesn't seem to work, look at old Nickel data?
;;sxaddpar, header, 'CRVAL2', 0
;;sxaddpar, header, 'CDELT1', fix(xbin)
;;sxaddpar, header, 'CDELT2', fix(ybin)

;lmode = ''
;if camera eq 'b' and clip(sxpar(header,'GRISNAME'))  eq 'clear'  then lmode = 'img'
;if camera eq 'r' and clip(sxpar(header,'GRANAME'))  eq 'mirror' then lmode = 'img'
;if camera eq 'b' and strpos(sxpar(header,'GRISNAME'),'/') gt 0   then lmode = 'spec'
;if camera eq 'r' and strpos(sxpar(header,'GRANAME'),'/') gt 0    then
;lmode = 'spec'
lmode = 'img'                   ;placeholder

;smode = ''
;slit =  strlowcase(clip(sxpar(header, 'SLITNAME')))
;if slit eq 'direct' then smode = 'img'
;if strpos(slit,'_') gt 0 then smode = 'spec'
smode = 'img'                   ;placeholder

;if lmode eq 'img' and smode eq 'spec' or smode eq 'img' and lmode eq 'spec' then begin
;   print
;   print, '     WARNING: Mixed-up grating/slit settings!  Ignoring file...'
;   return
;endif
;if lmode eq '' and smode eq '' then begin
;   print
;   print, '     WARNING: grating and slit keywords corrupt, absent, or unusual.'
;  if camera eq 'b' then print, '      GRISNAME:', clip(sxpar(header,'GRISNAME')), '    SLITNAME:', clip(sxpar(header,'SLITNAME')), '  FILTER:', clip(sxpar(header,'BLUFILT')), '  :', lmode, '?'
;  if camera eq 'r' then print, '       GRANAME:', clip(sxpar(header,'GRANAME')), '  SLITNAME:', clip(sxpar(header,'SLITNAME')), '  FILTER:', clip(sxpar(header,'REDFILT')), '  :', lmode, '?'
;   print, '              Skipping...'
;   return
;endif
;if lmode eq '' and smode ne '' then lmode = smode

; Check for missing keywords
targname = sxpar(header, 'TARGNAME', count=ct)
if ct eq 0 then begin
   prpslid = sxpar(header,'PRPSLID',count =ct)
   vstid = sxpar(header, 'VSTID', count=ct)
   IF ct gt 0 then begin
      targname=strcompress(prpslid+'-'+'vis'+string(vstid),/remove_all)
      print,'Target Name Now ',targname
   ENDIF ELSE BEGIN
      print, '[TARGNAME keyword missing from header!]', format='($,A)'
      sxaddpar, header, 'TARGNAME', 'unknown'
   endelse
endif

; Tidy up the header by removing unhelpful information and 
; reorganizing some keywords.
; If anyone recognizes keywords here that do contain useful information
; for an observer, please let me know.

sxdelpar, header, ''           ; remove blank lines
sxdelpar, header, 'COMMENT' 
sxdelpar, header, 'VSTNM'
sxdelpar, header, 'VSTID'
sxdelpar, header, 'ECRD4TR'   
sxdelpar, header, 'ECRD4' 
sxdelpar, header, 'ECRD3TR'   
sxdelpar, header, 'ECRD3' 
sxdelpar, header, 'ECRD2TR'   
sxdelpar, header, 'ECRD2' 
sxdelpar, header, 'ECRD1TR'   
sxdelpar, header, 'ECRD1'
sxdelpar, header, 'ECRC4TR'   
sxdelpar, header, 'ECRC4' 
sxdelpar, header, 'ECRC3TR'   
sxdelpar, header, 'ECRC3' 
sxdelpar, header, 'ECRC2TR'   
sxdelpar, header, 'ECRC2' 
sxdelpar, header, 'ECRC1TR'   
sxdelpar, header, 'ECRC1'
sxdelpar, header, 'ECRBTR'
sxdelpar, header, 'ECRB'
sxdelpar, header, 'ECRATR'
sxdelpar, header, 'ECRA'
sxdelpar, header, 'ECRTMT'
sxdelpar, header, 'ECRPRTR'
sxdelpar, header, 'ECRPR'
sxdelpar, header, 'ECRPRT'
sxdelpar, header, 'ECRT'
sxdelpar, header, 'ECRACT'
sxdelpar, header, 'ECRAC'
sxdelpar, header, 'ECRRQACT'
sxdelpar, header, 'ECRRQAC'
sxdelpar, header, 'ECRSTT'
sxdelpar, header, 'ECRST'

;sxrenamepar, header, 'GRISM', 'GRISNUM'    ; these cause confusion.
;sxrenamepar, header, 'DICHROIC', 'DICHNUM'
;sxrenamepar, header, 'GRATING', 'GRANUM'
;sxrenamepar, header, 'OBJECT', 'TEMPOBJ' ; remove duplicate keyword
;sxrenamepar, header, 'TEMPOBJ', 'OBJECT'

;sxmovepar, header, 'TEMPDET'
;sxmovepar, header, 'TEMPSET'
;
;sxmovepar, header, 'DEC', after='DATE-OBS'
;sxmovepar, header, 'RA', before='DEC'
;sxmovepar, header, 'EQUINOX', after='DEC'
;sxmovepar, header, 'HA', before='AIRMASS'
;sxmovepar, header, 'OBJECT', after='TARGNAME'

if clip(sxpar(header, 'UTC')) eq '0' then sxaddpar, header, 'UTC', sxpar(header, 'UT'), after='UT'
if clip(sxpar(header, 'UT')) eq '0' then sxaddpar, header, 'UT', sxpar(header, 'UTC'), after='UTC'

;sxaddpar, header, 'BZERO', 0 
;sxaddpar, header, 'BSCALE', 1
;sxaddpar, header, 'NAXS1', 0, after='NAXIS' ; just to save the location
;sxaddpar, header, 'NAXIS2', 0, after='NAXIS1'

if n_elements(namefixfiles) gt 0 then begin
   namefix, header, namefixfiles
endif

if sxpar(header, 'EXPTIME') eq 0 then begin ; swarp needs this
    if sxpar(header, 'ELAPTIME') gt 0 then begin
       sxaddpar, header, 'EXPTIME', sxpar(header, 'ELAPTIME'), before='ELAPTIME'
    endif else begin
       if sxpar(header, 'TTIME') gt 0 then $
          sxaddpar, header, 'EXPTIME', sxpar(header, 'TTIME'), before='TTIME'
    endelse
endif


; Saturation on the help page keeps changing.  Currently claims both are 60k.
if camera eq 'b' then sxaddpar, header, 'SATURATE', 55000, before='BINNING'
;if camera eq 'r' and lrisversion eq 1 then sxaddpar, header, 'SATURATE', 55000, before='BINNING'
;if camera eq 'r' and lrisversion ge 2 then sxaddpar, header, 'SATURATE', 55000, before='BINNING'

;if camera eq 'b' then sxaddpar, header,'FILTER', sxpar(header,'BLUFILT'), before='BLUFILT'
;if camera eq 'r' then sxaddpar, header,'FILTER', sxpar(header,'REDFILT'), before='BLUFILT'

;if camera eq 'b' then sxdelpar, header, ['REDFILT','REDFNUM','REDFOCUS']
;if camera eq 'r' then sxdelpar, header, ['BLUFILT','BLUFNUM','BLUFTRAN','BTRAN','BLUFOCUS']

;if camera eq 'r' then begin
;   if lrisversion eq 1 then pxscale = 0.210*(xbin) $
;                       else pxscale = 0.135*(xbin)
;endif
;if camera eq 'b' then pxscale = 0.1357*(xbin)
pxscale = 0.3*(xbin)

sxaddpar, header, 'PXSCALE', pxscale, 'arcsec/pix', after='WAVELEN'
sxaddpar, header, 'OBSERVAT',"SPM", before='TELESCOP'
;sxaddpar, header, 'SLITPA', (float(sxpar(header, 'ROTPOSN'))+90) mod 360., 'Slit or right chip CCW from North', before='ROTPOSN'
sxaddpar, header, 'RDNOISE', 4, before='BINNING' ; default as placeholder

; --- Sun/moon ephemeris ---

mjd = sxpar(header, 'MJD-OBS')
jd = double(mjd) + 2400000.5
sxaddpar, header, 'JD', jd, after='MJD-OBS'

;;sxaddpar, header, 'CRVAL1', 0          ; doesn't seem to work, look at old Nickel data?
;;sxaddpar, header, 'CRVAL2', 0
;;sxaddpar, header, 'CDELT1', fix(xbin)
;;sxaddpar, header, 'CDELT2', fix(ybin)

lst = sxpar(header, 'ST')
zenithra = 15.*ten(lst)
zenithdec = 30.7500

; these are using the start of the exposure.  The middle (or "worst") would be more useful.
moonpos, jd, moonra, moondec  ; needs isarray, which seems to have vanished from the GSFC library.
sunpos, jd, sunra, sundec

;radeg = sxpar(header,'TELRA')
;decdeg = sxpar(header,'TELDEC')
;radeg = sxpar(header,'STRCURA')
;decdeg = sxpar(header,'STRCHDE')
;raoff = sxpar(header,'STRAPRAO')
;decoff = sxpar(header,'STRAPDEO')
;radeg = radeg-raoff
;decdeg = decdeg-decoff
;IF strcmp(filter,'H') OR strcmp(filter,'Y') THEN radeg = radeg-450.*pxscale_deg  ;account for the fact that we're dividing the image into two different planes
;IF strcmp(filter,'J') OR strcmp(filter,'Z') THEN radeg = radeg+450.*pxscale_deg
;decdeg = decdeg + 30./3600.     ;artificial offset
;pxscale = 0.3*(xbin)
;pxscale_deg = 0.3/3600.

radeg = sxpar(header,'ETRRQRA')
decdeg = sxpar(header,'ETRRQDE')
;radeg=120.53095833333333        ;06jd 
;decdeg=0.80875                  ;06jd
;radeg = 145.722208              ;10jl
;decdeg = +09.494944             ;10jl

;(2012A-0002,vis1 = 259.2807917) (2012A-0000-vis4 = 276.8063333) (2012A-0000-vis8 = 316.0614583) (2012A-0000-vis9 = 316.0614583) (2012A-0000-vis12 = 348.09)
;(2012A-0002,vis1 = 43.13594444) (2012A-0000-vis4 = 4.0526111) (2012A-0000-vis8 = 313.1970833) (2012A-0000-vis9 = 6.66805556) (2012A-0000-vis12 = 10.784472)

;radeg =35.665;2:22:39.6
;decdeg = 42.99205;43.03553;43:02:07.9

; 326.78725;
; -8.235444;

tmp=strsplit(filename,'.',/extract)
tmp=strsplit(tmp[0],'_',/extract)
filter=tmp[2]

gcirc, 2, zenithra, zenithdec, moonra, moondec, moonzdist
gcirc, 2, zenithra, zenithdec, sunra, sundec, sunzdist
gcirc, 2, radeg, decdeg, moonra, moondec, moondist

sunelev = 90. - sunzdist/3600.
moonelev = 90. - moonzdist/3600.
moonsep = moondist/3600

sxaddpar, header, 'SUNELEV', sunelev, (sunelev>0.?'Above horizon':'Below horizon'), before='ROTMODE'
sxaddpar, header, 'MOONELEV', moonelev, (moonelev>0?'Above horizon':'Below horizon'), before='ROTMODE'
sxaddpar, header, 'MOONDIST', moonsep, before='ROTMODE'

;if lmode eq 'spec' then begin
;  dispersion = 0.
;  if camera eq 'b' then begin
;    grism = clip(sxpar(header, 'GRISNAME'))
;    if grism eq  '300/5000' then dispersion = 1.43
;    if grism eq  '400/3400' then dispersion = 1.09
;    if grism eq  '600/4000' then dispersion = 0.63
;    if grism eq '1200/3400' then dispersion = 0.24
;    if dispersion eq 0. then print, '   WARNING: Unrecognized grism GRISNAME=', grism $
;    else sxaddpar, header, 'DISPERSN', dispersion*ybin, 'Ang/pix'  ; ybin (not rotated yet)
;    sxaddpar, header, 'GRNAME', grism
;  endif
;  if camera eq 'r' then begin
;   grating = clip(sxpar(header,'GRANAME'))
;   if grating eq  '150/7500' then dispersion = 3.0
;   if grating eq  '300/5000' then dispersion = 1.59
;   if grating eq  '400/8500' then dispersion = 1.16
;   if grating eq  '600/5000' then dispersion = 0.80
;   if grating eq  '600/7500' then dispersion = 0.80
;   if grating eq  '600/10500' then dispersion = 0.80
;   if grating eq  '831/8200' then dispersion = 0.58
;   if grating eq  '900/5500' then dispersion = 0.53
;   if grating eq '1200/7500' then dispersion = 0.40
;   if dispersion eq 0. then print, '   WARNING: Unrecognized grating GRANAME=', grating $
;   else sxaddpar, header, 'DISPERSN', dispersion*ybin, 'Ang/pix'  ; ybin (not rotated yet)
;   sxaddpar, header, 'GRNAME', grating
;  endif
;endif
if lmode eq 'img' then begin

;   radeg = sxpar(header,'TELRA')
;   decdeg = sxpar(header,'TELDEC')
;   radeg = sxpar(header,'STRCURA')
;   decdeg = sxpar(header,'STRCHDE')
;   raoff = sxpar(header,'STRAPRAO')
;   decoff = sxpar(header,'STRAPDEO')
;   radeg = radeg-raoff
;   decdeg = decdeg-decoff   
;   pxscale = 0.3*(xbin)
;   pxscale_deg = 0.3/3600.
;   IF strcmp(filter,'H') OR strcmp(filter,'Y') THEN radeg = radeg-450.*pxscale_deg ;account for the fact that we're dividing the image into two different planes
;   IF strcmp(filter,'J') OR strcmp(filter,'Z') THEN radeg = radeg+450.*pxscale_deg
;   decdeg = decdeg + 30./3600.  ;artificial offset

;   radeg = sxpar(header,'ETROBRA')
;   decdeg = sxpar(header,'ETROBDE')

   tmp=strsplit(filename,'.',/extract)
   tmp=strsplit(tmp[0],'_',/extract)
   filter=tmp[2]

   IF strcmp(filter,'J') OR strcmp (filter,'Z') THEN BEGIN
;      sxaddpar,header,'CRPIX1',899 
;      sxaddpar,header,'CRPIX2',850.
      sxaddpar,header,'CRPIX1', 700. 
      sxaddpar,header,'CRPIX2',900.
      sxaddpar, header, 'CD1_1',  -8.17218901106E-05
      sxaddpar, header, 'CD1_2',  3.7651099887E-06
      sxaddpar, header, 'CD2_1',  4.20827225712E-06
      sxaddpar, header, 'CD2_2',  8.26704009041E-05
   ENDIF ELSE IF strcmp(filter,'H') OR strcmp (filter,'Y') THEN BEGIN
      sxaddpar, header, 'CRPIX1', 200.
      sxaddpar, header, 'CRPIX2', 900.
      sxaddpar, header, 'CD1_1',  -8.17218901106E-05
      sxaddpar, header, 'CD1_2',  3.7651099887E-06
      sxaddpar, header, 'CD2_1',  4.20827225712E-06
      sxaddpar, header, 'CD2_2',  8.26704009041E-05
   ENDIF ELSE BEGIN
      sxaddpar, header, 'CD1_1', -8.80977078079E-05 
      sxaddpar, header, 'CD1_2',  1.86753101419E-06 
      sxaddpar, header, 'CD2_1',  1.86716671065E-06
      sxaddpar, header, 'CD2_2',  8.81208878042E-05 
      sxaddpar, header, 'CRPIX1', 512.
      sxaddpar, header, 'CRPIX2', 512.
   ENDELSE
   sxaddpar,header,'AIRMASS',1.
   
;   padeg = sxpar(header, 'ROTPOSN')
   padeg = 0                    ;placeholder (since everything has been rotated appropriately)
   parad = padeg*3.1415927/180.
   pxscaledeg = pxscale/3600.
   parity = 1
;   sxaddpar, header, 'CD1_1',  pxscaledeg * cos(parad)*parity
;   sxaddpar, header, 'CD1_2',  pxscaledeg * sin(parad)
;   sxaddpar, header, 'CD2_1', -pxscaledeg * sin(parad)*parity
;   sxaddpar, header, 'CD2_2',  pxscaledeg * cos(parad)
;   sxaddpar,header,'CD1_1', -4.30761532329e-06
;   sxaddpar,header,'CD1_2', -8.29197751839e-05
;   sxaddpar,header,'CD2_1', -8.29197751839e-05
;   sxaddpar,header,'CD2_2', 4.30761532329e-06

   sxaddpar, header, 'CRVAL1', float(radeg)
   sxaddpar, header, 'CRVAL2', float(decdeg)
   sxaddpar, header, 'CTYPE1', 'RA---TAN'
   sxaddpar, header, 'CTYPE2', 'DEC--TAN'
;;sxaddpar, header, 'CRPIX1', 0
;;sxaddpar, header, 'CRPIX2', 0
;;sxaddpar, header, 'CRVAL1', 0          ; doesn't seem to work, look at old Nickel data?
;;sxaddpar, header, 'CRVAL2', 0
;;sxaddpar, header, 'CDELT1', fix(xbin)
;;sxaddpar, header, 'CDELT2', fix(ybin)
   ; crpix will wait until post-cropping

   exptime = sxpar(header,'EXPTIME')
   sxaddpar,header,'ELAPTIME',exptime
endif

headtime = systime(/seconds)  ; mark time after header modifications, including epheremis


;; -------- Fix CCD bug (?) -----------
;
;; correct for the different pixel readout pattern in binned data (some bug on the CCD?)
;if lrisversion eq 3 and camera eq 'r' and xbin eq 2 then begin
;   ampl2 = fix(strsplit(sxpar(header,'AMP1'),',',/extract))
;   ampr1 = fix(strsplit(sxpar(header,'AMP4'),',',/extract))
;   hold = array[nxin-2:nxin-1,*]
;   array[ampr1[0]-1:nxin-3,*] = array[ampr1[0]+1:nxin-1,*]
;   ; array[nxin-2:nxin-1,*] = hold ; not actually sure if this is the best thing to do -
;                                   ; will mess up the edge of an imaging frame windowed and
;                                   ; binned for spectroscopy.
;                           
;endif
;
;; Should check again to verify that this is a real effect of binning and is necessary.
;
;; New exciting bug discovered - diminished columns (5-10% lower response) seem to move around
;; even within the same setting.
;
;
; -------------- Figure out the amp/chip boundaries ------------
;
;
;if camera eq 'b' then begin
;   ampl1 = fix(strsplit(sxpar(header,'AMP1'),',',/extract))
;   ampl2 = fix(strsplit(sxpar(header,'AMP2'),',',/extract))
;   ampr1 = fix(strsplit(sxpar(header,'AMP3'),',',/extract))
;   ampr2 = fix(strsplit(sxpar(header,'AMP4'),',',/extract))
;endif else begin
;   ampl1 = fix(strsplit(sxpar(header,'AMP2'),',',/extract))
;   ampl2 = fix(strsplit(sxpar(header,'AMP1'),',',/extract))  ; red amps are not in sequence
;   ampr1 = fix(strsplit(sxpar(header,'AMP4'),',',/extract))
;   ampr2 = fix(strsplit(sxpar(header,'AMP3'),',',/extract))
;endelse
;
;chipl = [min([ampl1,ampl2]), max([ampl1,ampl2])]
;chipr = [min(ampr1),         max([ampr1,ampr2])]
;
;chipboundaryl = chipl[1]
;chipboundaryr = chipr[0] ; should be separated by one...
;
;
; ---------------- Flag/fix bad pixels -------------------

;; Load bad pixel mask and rebin/crop appropriately
;
;if camera eq 'r' then begin
;
;   reusemask = 0
;   if n_elements(mask) gt 0 and n_elements(prevversion) gt 0 and n_elements(prevpane) eq n_elements(pane) then begin
;     if lrisversion eq prevversion and xbin eq prevxbin and ybin eq prevybin and total(pane eq prevpane) eq n_elements(pane) then $
;        reusemask = 1      
;   endif
;   ; reuse the cropped mask if we can, although this saves only trivial processing time (aside from the time saved by not rereading the file)
;
;   if reusemask eq 0 then begin ; load mask from disk
;      ; should confirm that this works as planned....
;      bpmfile = ''
;      if camera eq 'r' and lrisversion eq 2 and n_elements(bpm2) lt 2 then bpmfile = bpmpath+'bpm_red2.fits'
;      if camera eq 'r' and lrisversion eq 3 and n_elements(bpm3) lt 2 then bpmfile = bpmpath+'bpm_red3.fits'
;      if bpmfile ne '' then if file_test(bpmfile) eq 0 then begin
;         print
;         print, 'Bad pixel mask '+bpmfile+': check refdata directory'
;      endif else begin
;         if camera eq 'r' and lrisversion eq 2 and n_elements(bpm2) lt 2 then bpm2 = mrdfits(bpmfile,0,hbpm, /silent)
;         if camera eq 'r' and lrisversion eq 3 and n_elements(bpm3) lt 2 then bpm3 = mrdfits(bpmfile,0,hbpm, /silent)
;         if camera eq 'r' and lrisversion eq 2 then mask = bpm2[pane[0]:pane[0]+pane[2]-1,pane[1]:pane[1]+pane[3]-1]  ; pane is before binning
;         if camera eq 'r' and lrisversion eq 3 then mask = bpm3[pane[0]:pane[0]+pane[2]-1,pane[1]:pane[1]+pane[3]-1]
;         sbpm = size(mask)
;         if xbin gt 1 or ybin gt 1 then mask = rebin(mask, sbpm[1]/xbin, sbpm[2]/ybin)
;         ;if xbin gt 1 or ybin gt 1 then mask = shift(mask, 1, 0) ; mysterious why this is necessary.  should be revisited with flats.
;         ; should no longer be necessary with the problem now identified as a CCD bug, but revisit this to check someday
;      endelse
;   endif else begin
;      ; just reuse the old mask that's already in there.
;      sbpm = size(mask)
;   endelse
;   
;   if n_elements(flag) gt 0 then array[where(mask gt 0)] = !values.f_nan
;   if n_elements(fixdefect) gt 0 then begin
;     array[where(mask gt 0)] = !values.f_nan
;     ct = 1
;     while ct gt 0 do begin 
;       bad = where(finite(array) eq 0, ct)
;       for m = 0, ct-1 do begin
;         mm = bad[m]
;         array[mm] = median(array[mm-2:mm+2])
;       endfor
;     endwhile
;   endif
;
;   ; Recognize and replace hyper-saturated pixels (so saturated they register as dead)
;   ; from neon or argon lamps
;   ;  (could extend this to all files at the expense of processing time if it seems
;   ;   necessary)
;   lamps = sxpar(header, 'LAMPS')
;   lamparr = fix(strsplit(lamps,',',/extract))
;   if lamparr[1] or lamparr[2] then begin
;     truesaturation = max(array) ; both chips should saturate at the same place?
;     arrayl = array[chipl[0]:chipl[1]-4,*] 
;     arrayr = array[chipr[0]+4:chipr[1],*] 
;     arrayl = fixhypersaturation(arrayl, deadthresh=1500)
;     arrayr = fixhypersaturation(arrayr, deadthresh=500)
;     array[chipl[0]:chipl[1]-4,*] = arrayl
;     array[chipr[0]+4:chipr[1],*] = arrayr
;
;   endif
;
;   prevversion = lrisversion
;   prevxbin = xbin
;   prevybin = ybin
;   prevpane = pane
;   
;endif
;
;prevcamera = camera
;
;flagtime = systime(/seconds)  ; mark time after flagging
;

;; ------------- Gain correction -------------------
;
;if n_elements(gain) eq 0 then gain = [1.,1,1,1]
;
;if camera eq 'r' then sxaddpar, header, 'GAIN', 0.95/gain[2], after='RDNOISE' ; is this actually true?  It is close to 1, at any rate.
;if camera eq 'b' then sxaddpar, header, 'GAIN', 1.63/gain[2], after='RDNOISE'
;if max(ampl1) gt 0 then array[ampl1[0]:ampl1[1],*] = array[ampl1[0]:ampl1[1],*] * gain[0]
;if max(ampl2) gt 0 then array[ampl2[0]:ampl2[1],*] = array[ampl2[0]:ampl2[1],*] * gain[1]
;if max(ampr1) gt 0 then array[ampr1[0]:ampr1[1],*] = array[ampr1[0]:ampr1[1],*] * gain[2]
;if max(ampr2) gt 0 then array[ampr2[0]:ampr2[1],*] = array[ampr2[0]:ampr2[1],*] * gain[3]
;if max(ampl1) gt 0 then sxaddpar, header, 'INGAINL1', gain[0]
;if max(ampl2) gt 0 then sxaddpar, header, 'INGAINL2', gain[1]
;if max(ampr1) gt 0 then sxaddpar, header, 'INGAINR1', gain[2]
;if max(ampr2) gt 0 then sxaddpar, header, 'INGAINR2', gain[3]
;

; ---------- Figure out crop boundaries ----------

;if n_elements(ycrop) eq 0 then ycrop = [0,nyin-1]
;
;if n_elements(crops) eq 1 then begin
;  if crops eq 'none' then begin
;     lcrop = chipl
;     rcrop = chipr
;     ycrop = [0, nyin-1]
;     xcrop = [chipl[0], chipr[1]]
;  endif
;  if crops eq 'auto' then begin
;     ycrop = [0, nyin-1] ; default
;     if camera eq 'b' then begin
;        if lmode eq 'img' then begin
;           lcrop = [(360-pane[0])/xbin > chipl[0], chipl[1]]
;           rcrop = [chipr[0], (3790-pane[0])/xbin-1 < chipr[1]]
;           ycrop = [(760-pane[1])/ybin > 0, (3200-pane[1])/ybin-1 < nyin-1]
;           xcrop = [lcrop[0], rcrop[1]]
;        endif
;        if lmode eq 'spec' then begin
;           lcrop = [(1440-pane[0])/xbin > chipl[0],      chipl[1]]
;           rcrop = [chipr[0],      (2670-pane[0])/xbin-1 < chipr[1]]
;           xcrop = [lcrop[0], rcrop[1]]
;        endif
;     endif
;     if camera eq 'r' then begin
;        if (lrisversion eq 2 or lrisversion eq 3) and lmode eq 'img' then begin
;           lcrop = [(430-pane[0])/xbin > chipl[0], chipl[1]]
;           rcrop = [chipr[0],  (3778-pane[0])/xbin-1 < chipr[1]]
;           ycrop = [(638-pane[1])/ybin > 0, (3100-pane[1])/ybin-1 < nyin-1]
;           xcrop = [lcrop[0], rcrop[1]]
;        endif
;        if (lrisversion eq 2 or lrisversion eq 3) and lmode eq 'spec' then begin
;           lcrop = [(1520-pane[0])/xbin > chipl[0], chipl[1]]
;           rcrop = [chipr[0],  (2660-pane[0])/xbin-1 < chipr[1]] ; 2660
;           xcrop = [lcrop[0], rcrop[1]]
;        endif
;     endif
;  endif
;endif

dotpos = strpos(filename,'.',/reverse_search)
slashpos = strpos(filename,'/',/reverse_search)
fileroot = strmid(filename,slashpos+1,dotpos-slashpos-1)
fileext = strmid(filename,dotpos+1)

get_date, now
hist = 'Processed by lrisprepare '+now

;if ycrop[1] le ycrop[0] or xcrop[1] le xcrop[0] then begin
;   print
;   print, '  ERROR:  Invalid crop region.  Cannot prepare file.'
;   stop
;   return
;   ; This seems to occur when processing MIRA files.
;   ; investigate why, someday.  2012-04-27 run.
;endif


; ------ Do the cropping and make final edits to the header, then write to disk ------

sxaddpar, header, 'LTM1_1', 1./xbin
sxaddpar, header, 'LTM2_2', 1./ybin
header2 = header

n=n_elements(header2)
header = [header2[0:30],header2[n-25:n-1]]
filter = clip(sxpar(header2, 'FILTER'))
IF strcmp(filter,'BB') THEN filter='r'
sxaddpar, header, 'FILTER', filter

outfilename = outname
mwrfits, array, outfilename, header, /create
print, ' -> ', outfilename
;if merge then begin
;   array  = array[xcrop[0]:xcrop[1],  ycrop[0]:ycrop[1]]
;   croptime = systime(/seconds)  ; mark time after cropping
;
;   sxaddpar, header, 'LTV1', -pane[0]*1.0/xbin-xcrop[0] ; LTV1 and LTV2 will be reversed in transpose function
;   sxaddpar, header, 'LTV2', -pane[1]*1.0/ybin-ycrop[0]
;
;   if lmode eq 'img' then begin
;     if camera eq 'r' then begin
;        sxaddpar, header, 'CRPIX1', 2060 +  (29-xcrop[0]), after='CRVAL2'
;        sxaddpar, header, 'CRPIX2', 1150 +  (37-ycrop[0]), after='CRPIX1'
;     endif
;     if camera eq 'b' then begin
;        sxaddpar, header, 'CRPIX1', 180  + (2050-xcrop[0]), after='CRVAL2'
;        sxaddpar, header, 'CRPIX2', 1180 +  (768-ycrop[0]), after='CRPIX1'
;     endif
;   endif
;
;   ampl1 -= xcrop[0]
;   ampl2 -= xcrop[0]
;   ampr1 -= xcrop[0]
;   ampr2 -= xcrop[0]
;   if n_elements(ampl1) ge 2 then sxaddpar, header, 'AMPL1', clip(ampl1[0])+','+clip(ampl1[1]), ' (image coordinates)'
;   if n_elements(ampl2) ge 2 then sxaddpar, header, 'AMPL2', clip(ampl2[0])+','+clip(ampl2[1])
;   if n_elements(ampr1) ge 2 then sxaddpar, header, 'AMPR1', clip(ampr1[0])+','+clip(ampr1[1])
;   if n_elements(ampr2) ge 2 then sxaddpar, header, 'AMPR2', clip(ampr2[0])+','+clip(ampr2[1])
;   sxdelpar, header, 'AMP1'
;   sxdelpar, header, 'AMP2'
;   sxdelpar, header, 'AMP3'
;   sxdelpar, header, 'AMP4'
;
;   sxaddpar, header, 'COUNTS', median(array)
;
;   if camera eq 'r' then sxaddpar, header, 'RDNOISE', 4
;
;   if n_elements(outname) gt 0 then outfilename = outname
;   if n_elements(outfilename) eq 0 then begin
;      if n_elements(prefix) eq 0 and n_elements(suffix) eq 0 then begin
;         prefix = 'p'
;      endif
;      if n_elements(prefix) gt 0 then outfilename = prefix+fileroot+'.'+fileext
;      if n_elements(suffix) gt 0 then outfilename = fileroot+suffix+'.fits'
;   endif
;
;   if keyword_set(notranspose) eq 0 and lmode eq 'spec' then begin
;      ltranspose, array, header
;   endif
;
;   print, ' -> ', outfilename
;
;   sxaddhist, hist, header
;
;
;   mwrfits, array, outfilename, header, /create
;endif
;if split then begin
;   rheader = header
;   lheader = header
;
;   if keyword_set(leftonly)  and n_elements(lcrop) eq 0 and n_elements(xcrop) eq 2 then lcrop = xcrop
;   if keyword_set(rightonly) and n_elements(rcrop) eq 0 and n_elements(xcrop) eq 2 then rcrop = xcrop
;   doleft  = keyword_set(rightonly) eq 0
;   doright = keyword_set(leftonly)  eq 0
; 
;   if doleft  then arrayl = array[lcrop[0]:lcrop[1],ycrop[0]:ycrop[1]] 
;   if doright then arrayr = array[rcrop[0]:rcrop[1],ycrop[0]:ycrop[1]] 
;
;   croptime = systime(/seconds)  ; mark time after cropping
;
;
;   if lmode eq 'img' then begin
;     if camera eq 'r' then begin
;        sxaddpar, rheader, 'CRPIX1', 213  + (1657-rcrop[0]), after='CRVAL2'
;        sxaddpar, rheader, 'CRPIX2', 1150 +   (37-ycrop[0]), after='CRPIX1'
;         sxaddpar, lheader, 'CRPIX1', 2060 +   (29-lcrop[0]), after='CRVAL2'
;        sxaddpar, lheader, 'CRPIX2', 1150 +   (37-ycrop[0]), after='CRPIX1'
;     endif
;     if camera eq 'b' then begin
;        sxaddpar, rheader, 'CRPIX1', 180  + (2050-rcrop[0]), after='CRVAL2'
;        sxaddpar, rheader, 'CRPIX2', 1180 +  (768-ycrop[0]), after='CRPIX1'
;
;        sxaddpar, lheader, 'CRPIX1', 1970  + (374-lcrop[0]), after='CRVAL2'
;        sxaddpar, lheader, 'CRPIX2', 1180  + (768-ycrop[0]), after='CRPIX1'
;
;        leftrot = 0.5
;        sxaddpar, lheader, 'CD1_1',  pxscaledeg * cos(parad+leftrot)*parity
;        sxaddpar, lheader, 'CD1_2',  pxscaledeg * sin(parad+leftrot)
;        sxaddpar, lheader, 'CD2_1', -pxscaledeg * sin(parad+leftrot)*parity
;        sxaddpar, lheader, 'CD2_2',  pxscaledeg * cos(parad+leftrot)
;     endif
;   endif
;
;   if doleft  then sxaddpar, lheader, 'LTV1', -pane[0]/xbin-lcrop[0]  ; LTV1 and LTV2 will be reversed in transpose function
;   if doright then sxaddpar, rheader, 'LTV1', -pane[0]/xbin-rcrop[0]
;
;   if doleft  then sxaddpar, lheader, 'LTV2', -pane[1]/ybin-ycrop[0]
;   if doright then sxaddpar, rheader, 'LTV2', -pane[1]/ybin-ycrop[0]
;
;   if camera eq 'r' then begin
;      if doleft  then sxaddpar, lheader, 'RDNOISE', 4
;      if doright then sxaddpar, rheader, 'RDNOISE', 3.79 ; confusing, since the fits header says this is 1-13.
;   endif
;
;   if doleft  then sxaddpar, lheader, 'CHIP', 'left'
;   if doright then sxaddpar, rheader, 'CHIP', 'right'
;
;   if doleft  then sxaddpar, lheader, 'COUNTS', median(arrayl)
;   if doright then sxaddpar, rheader, 'COUNTS', median(arrayr)
;
;   if keyword_set(leftonly) and n_elements(loutfilename) eq 0 and n_elements(outfilename) gt 0 then loutfilename=outfilename
;   if keyword_set(rightonly) and n_elements(routfilename) eq 0 and n_elements(outfilename) gt 0 then routfilename=outfilename
;
;   if n_elements(prefix) gt 0 and n_elements(rprefix) eq 0 and keyword_set(rightonly) then rprefix = prefix
;   if n_elements(prefix) gt 0 and n_elements(lprefix) eq 0 and keyword_set(rightonly) then lprefix = prefix
;
;   if n_elements(loutname) gt 0 then loutfilename = loutname
;   if n_elements(loutfilename) eq 0 then begin
;      if n_elements(lprefix) eq 0 and n_elements(lsuffix) eq 0 then begin
;        lprefix = 'p'
;        lsuffix = 'l'
;      endif
;      if n_elements(lprefix) gt 0 then loutfilename = lprefix+fileroot+'.'+fileext
;      if n_elements(lsuffix) gt 0 then loutfilename = fileroot+lsuffix+'.fits'
;   endif
;   if n_elements(routname) gt 0 then routfilename = routname
;   if n_elements(routfilename) eq 0 then begin
;      if n_elements(rprefix) eq 0 and n_elements(rsuffix) eq 0 then begin
;         rprefix = 'p'
;         rsuffix = 'r'
;      endif
;      if n_elements(rprefix) gt 0 then routfilename = rprefix+fileroot+'.'+fileext
;      if n_elements(rsuffix) gt 0 then routfilename = fileroot+rsuffix+'.fits'
;   endif
;
;   if  keyword_set(notranspose) eq 0 and lmode eq 'spec' then  begin
;      if doleft  then ltranspose, arrayl, lheader 
;      if doright then ltranspose, arrayr, rheader 
;   endif
;
;   if doleft and doright        then print, ' -> ', loutfilename, ',', routfilename
;   if doleft and (doright eq 0) then print, ' -> ', loutfilename
;   if doright and (doleft eq 0) then print, ' -> ', routfilename
;
;   ;if doleft  and file_test(loutfilename) then print, 'Overwriting ', loutfilename
;   ;if doright and file_test(routfilename) then print, 'Overwriting ', routfilename
;
;   sxaddhist, hist, lheader
;   sxaddhist, hist, rheader
;
;   if doleft  then mwrfits, arrayl, loutfilename, lheader, /create
;   if doright then mwrfits, arrayr, routfilename, rheader, /create
;
;endif

finaltime = systime(/seconds)  ; mark time at end of file processing

if keyword_set(timer) then begin
print, 'Total duration:  ', finaltime-starttime
print, '  Read duration: ', readtime-starttime
print, '  Head duration: ', headtime-readtime
print, '  Flag duration: ', flagtime-headtime
print, '  Crop duration: ', croptime-flagtime
print, '  Write duration:', finaltime-croptime
endif

end

