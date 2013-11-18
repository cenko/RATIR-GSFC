pro hswap, header, hkey1, hkey2
   hold1 = sxpar(header, hkey1, count=check1)
   hold2 = sxpar(header, hkey2, count=check2)
   if check1 eq 0 and check2 eq 0 then return
   sxaddpar, header, hkey1, hold2
   sxaddpar, header, hkey2, hold1
end

pro ltranspose, array, header

   array = transpose(array)

   hswap, header, 'LTV1', 'LTV2'
   hswap, header, 'CDELT1', 'CDELT2'
   hswap, header, 'XBIN', 'YBIN'

   hswap, header, 'LTM2_1', 'LTM2_2'
   hswap, header, 'LTM1_2', 'LTM1_1'

end

