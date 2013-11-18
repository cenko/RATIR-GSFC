function countlines, infile

if file_test(infile) eq 0 then return, 0l

openr, unit, infile, /get_lun

iline=''
nlines = 0l
while not eof(unit) do begin
  readf, unit, iline
  nlines = nlines + 1
endwhile
close, unit
free_lun, unit

return, nlines

end
