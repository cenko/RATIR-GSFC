; Count how many files match a wildcard string.

function countfiles, arg
 
  if file_test(arg) eq 0  then return, 0 $
                          else return, n_elements(where(findfile(arg) ne ''))
   
end

