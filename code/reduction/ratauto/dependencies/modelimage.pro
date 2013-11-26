function modelsky, img, scale=scale
    ; Quickly remove the background.  Sort of like a median filter but faster.
    
      if n_elements(scale) eq 0 then scale = 9

      dims = size(img, /dimens)
      nx = dims[0]
      ny = dims[1]

      imgs0 = shift(img, scale, scale)
      imgs1 = shift(img, scale, 0)
      imgs2 = shift(img, scale, -scale)
      imgs3 = shift(img, 0, scale)
      imgs4 = shift(img, 0, -scale)
      imgs5 = shift(img, -scale, scale)
      imgs6 = shift(img, -scale, 0)
      imgs7 = shift(img, -scale, -scale)

      block = fltarr(dims[0], dims[1], 8)
      block[*,*,0] = imgs0
      block[*,*,1] = imgs1
      block[*,*,2] = imgs2
      block[*,*,3] = imgs3
      block[*,*,4] = imgs4
      block[*,*,5] = imgs5
      block[*,*,6] = imgs6
      block[*,*,7] = imgs7
   
      return, median(block, dimension=3, /even)

end


function modelsky2, img, scale=scale
    ; Quickly remove the background.  Sort of like a median filter but faster.

      if n_elements(scale) eq 0 then scale = 9

      dims = size(img, /dimens)
      block = fltarr(dims[0], dims[1], 4)
      block[*,*,0] = (shift(img, scale, 0) + shift(img, -scale, 0) ) / 2.
      block[*,*,1] = (shift(img, 0, scale) + shift(img, 0, -scale) ) / 2.
      block[*,*,2] = (shift(img, scale, scale) + shift(img, -scale, -scale) ) / 2.
      block[*,*,3] = (shift(img, scale, -scale) + shift(img, -scale, scale) ) / 2.
   
      return, median(block, dimension=3, /even)

end


function model2image, imgin

     idims = size(imgin, /dimens)
     inx = idims[0]
     iny = idims[1]
     nx = inx + 6
     ny = iny + 6

      ; buffer the image by making the buffer zones reflections over the edges

     img = fltarr(nx,ny)
     img[3:nx-4,3:ny-4] = imgin
     img[*,2] = img[*,3]
     img[*,1] = img[*,4]
     img[*,0] = img[*,5]
     img[*,ny-3] = img[*,ny-4]
     img[*,ny-2] = img[*,ny-5]
     img[*,ny-1] = img[*,ny-6]
     img[2,*] = img[3,*]
     img[1,*] = img[4,*]
     img[0,*] = img[5,*]
     img[nx-3] = img[nx-4]
     img[nx-2] = img[nx-5]
     img[nx-1] = img[nx-6]

     imgx = shift(img, 1, 0)
     ximg = shift(img, -1, 0)
     imgy = shift(img, 0, 1)
     yimg = shift(img, 0, -1)
     imgxy = shift(img, 1, 1)
     imgyx = shift(img, 1, -1)
     xyimg = shift(img, -1, -1)
     yximg = shift(img, -1, 1)
     imgxx  = shift(img, 2, 0)
     imgxxx = shift(img, 3, 0)
     xximg = shift(img, -2, 0)
     xxximg = shift(img, -3, 0)
     imgyy  = shift(img, 0, 2)
     imgyyy = shift(img, 0, 3)
     yyimg = shift(img, 0, -2)
     yyyimg = shift(img, 0, -3)
     imgxyxy = shift(img, 2, 2)
     imgxyxyxy = shift(img, 3, 3)
     xyxyimg =  shift(img, -2, -2)
     xyxyxyimg =  shift(img, -3, -3)
     imgyxyx = shift(img, 2, -2)
     imgyxyxyx =  shift(img, 3, -3)
     yxyximg = shift(img, -2, 2)
     yxyxyximg = shift(img, -3, 3)

     block = fltarr(nx, ny, 16)
     block[*,*,0] = (imgx + (imgx-imgxx))/2  + (ximg + (ximg-xximg))/2
     block[*,*,1] = (imgx+ximg)/2
     block[*,*,2] = (imgy + (imgy-imgyy))/2  + (yimg + (yimg-yyimg))/2
     block[*,*,3] = (imgy+yimg)/2
     block[*,*,4] = (imgxy + (imgxy-imgxyxy))/2  + (xyimg + (xyimg-xyxyimg))/2
     block[*,*,5] = (imgxy+xyimg)/2
     block[*,*,6] = (imgyx + (imgyx-imgyxyx))/2  + (yximg + (yximg-yxyximg))/2
     block[*,*,7] = (imgyx+yximg)/2  
     block[*,*,8] = (imgxx + (imgxx-imgxxx)*2)/2  + (xximg + (xximg-xxximg)*2)/2
     block[*,*,9] = (imgxx+xximg)/2
     block[*,*,10] = (imgyy + (imgyy-imgyyy)*2)/2  + (yyimg + (yyimg-yyyimg)*2)/2
     block[*,*,11] = (imgyy+yyimg)/2
     block[*,*,12] = (imgxyxy + (imgxyxy-imgxyxyxy)*2)/2  + (xyxyimg + (xyxyimg-xyxyxyimg)*2)/2
     block[*,*,13] = (imgxyxy+xyxyimg)/2
     block[*,*,14] = (imgyxyx + (imgyxyx-imgyxyxyx)*2)/2  + (yxyximg + (yxyximg-yxyxyximg)*2)/2
     block[*,*,15] = (imgyxyx+yxyximg)/2  


     ;return, ((img-imgx) > 0) + ((img-ximg) > 0) + ((img-yimg)>0) + ((img-imgy)>0) + $
     ;        ((img-imgxy)>0) + ((img-yimgx)>0) + ((img-ximgy)>0) + ((img-xyimg)>0)

     ;return, (img-(imgx+ximg)/2.)>0 + (img-(imgy+yimg)/2.)>0 + (img-(imgxy+xyimg)/2.)>0 + (img-(ximgy+yimgx)/2.)>0

     img = median(block, dimension=3, /even)

     return, img[3:nx-4,3:ny-4]

end

function modelimage, img

      ; Model the image using an interpolation scheme
      ; This is an average of a couple different methods and is not mathematically rigorous

      dims = size(img, /dimens)
      nx = dims[0]
      ny = dims[1]

      imgxx  = shift(img, 2, 0)
      imgxx[0,*] = img[0,*]
      imgxx[1,*] = img[0,*]
      imgxxx = shift(img, 3, 0)
      imgxxx[0,*] = img[0,*]
      imgxxx[1,*] = img[0,*]
      imgxxx[2,*] = img[0,*]

      xximg = shift(img, -2, 0)
      xximg[nx-1,*] = img[nx-1,*]
      xximg[nx-2,*] = img[nx-1,*]
      xxximg = shift(img, -3, 0)
      xxximg[nx-1,*] = img[nx-1,*]
      xxximg[nx-2,*] = img[nx-1,*]
      xxximg[nx-3,*] = img[nx-1,*]

      imgyy  = shift(img, 0, 2)
      imgyyy = shift(img, 0, 3)

      yyimg = shift(img, 0, -2)
      yyyimg = shift(img, 0, -3)

      imgxyxy = shift(img, 2, 2)
      imgxyxy[0,*] = img[0,*]
      imgxyxy[1,*] = img[0,*]
      imgxyxyxy = shift(img, 3, 3)
      imgxyxyxy[0,*] = img[0,*]
      imgxyxyxy[1,*] = img[0,*]
      imgxyxyxy[2,*] = img[0,*]

      xyxyimg =  shift(img, -2, -2)
      xximg[nx-1,*] = img[nx-1,*]
      xximg[nx-2,*] = img[nx-1,*]
      xyxyxyimg =  shift(img, -3, -3)
      xyxyxyimg[nx-1,*] = img[nx-1,*]
      xyxyxyimg[nx-2,*] = img[nx-1,*]
      xyxyxyimg[nx-3,*] = img[nx-1,*]

      imgyxyx = shift(img, 2, -2)
      imgyxyx[0,*] = img[0,*]
      imgyxyx[1,*] = img[0,*]
      imgyxyxyx =  shift(img, 3, -3)
      imgyxyxyx[0,*] = img[0,*]
      imgyxyxyx[1,*] = img[0,*]
      imgyxyxyx[2,*] = img[0,*]

      yxyximg = shift(img, -2, 2)
      yxyximg[nx-1,*] = img[nx-1,*]
      yxyximg[nx-2,*] = img[nx-1,*]
      yxyxyximg = shift(img, -3, 3)
      yxyxyximg[nx-1,*] = img[nx-1,*]
      yxyxyximg[nx-2,*] = img[nx-1,*]
      yxyxyximg[nx-3,*] = img[nx-1,*]

      block = fltarr(dims[0], dims[1], 8)
      block[*,*,0] = (imgxx + (imgxx-imgxxx)*2)/2  + (xximg + (xximg-xxximg)*2)/2
      block[*,*,1] = (imgxx+xximg)/2
      block[*,*,2] = (imgyy + (imgyy-imgyyy)*2)/2  + (yyimg + (yyimg-yyyimg)*2)/2
      block[*,*,3] = (imgyy+yyimg)/2
      block[*,*,4] = (imgxyxy + (imgxyxy-imgxyxyxy)*2)/2  + (xyxyimg + (xyxyimg-xyxyxyimg)*2)/2
      block[*,*,5] = (imgxyxy+xyxyimg)/2
      block[*,*,6] = (imgyxyx + (imgyxyx-imgyxyxyx)*2)/2  + (yxyximg + (yxyximg-yxyxyximg)*2)/2
      block[*,*,7] = (imgyxyx+yxyximg)/2  
      ;block = fltarr(dims[0], dims[1], 4)
      ;block[*,*,8] = (0.3*xximg + 1.5*imgxx - 0.8*imgxxx)/2. + $
      ;               (0.3*imgxx + 1.5*xximg - 0.8*xxximg)/2.
      ;block[*,*,9] = (0.3*yyimg + 1.5*imgyy - 0.8*imgyyy)/2. + $
      ;               (0.3*imgyy + 0.5*yyimg - 0.8*yyyimg)/2.
      ;block[*,*,10] = (0.3*xyxyimg + 1.5*imgxyxy - 0.8*imgxyxyxy)/2. + $
      ;               (0.3*imgxyxy + 1.5*xyxyimg - 0.8*xyxyxyimg)/2.
      ;block[*,*,11] = (0.3*yxyximg + 1.5*imgyxyx - 0.8*imgyxyxyx)/2. + $
      ;               (0.3*imgyxyx + 1.5*yxyximg - 0.8*yxyxyximg)/2.


      model = median(block, dimension=3, /even)

      ;model = fltarr(dims[0], dims[1])
      ;for x = 3, dims[0]-4 do begin
      ;  for y = 3, dims[1]-4 do begin
      ;      arr = block[x,y,*]
      ;      arr = arr[sort(arr)]
      ;      med = median(arr)
      ;      min = arr[0]
      ;      arr = arr[where(arr lt med+1.5*(med-min))]
      ;      model[x,y] = median(arr)
      ;  endfor
      ;endfor


      return, model

end

function usampmodelimage, img
      dims = size(img, /dimens)
      inblock = fltarr(dims[0], dims[1], 16)
      inblock[*,*,0] =  shift(img, 1, 0)
      inblock[*,*,1] =  shift(img,-1, 0)
      inblock[*,*,2] =  shift(img, 0, 1)
      inblock[*,*,3] =  shift(img, 0,-1)
      inblock[*,*,4] =  shift(img, 1, 1)
      inblock[*,*,5] =  shift(img,-1,-1)
      inblock[*,*,6] =  shift(img,-1, 1)
      inblock[*,*,7] =  shift(img, 1,-1)
      inblock[*,*,8] =  shift(img, 2, 0)
      inblock[*,*,9] =  shift(img,-2, 0)
      inblock[*,*,10] = shift(img, 0, 2)
      inblock[*,*,11] = shift(img, 0,-2)
      inblock[*,*,12] = shift(img, 2, 2)
      inblock[*,*,13] = shift(img,-2,-2)
      inblock[*,*,14] = shift(img, 2,-2)
      inblock[*,*,15] = shift(img,-2, 2)

      ;outblock = fltarr(dims[0], dims[1], 4)
      ;outblock[*,*,0] = 0.7*inblock[*,*,0] + 0.7*inblock[*,*,1] - 0.2*inblock[*,*,8]  - 0.2*inblock[*,*,9]
      ;outblock[*,*,1] = 0.7*inblock[*,*,2] + 0.7*inblock[*,*,3] - 0.2*inblock[*,*,10] - 0.2*inblock[*,*,11]
      ;outblock[*,*,2] = 0.7*inblock[*,*,4] + 0.7*inblock[*,*,5] - 0.2*inblock[*,*,12] - 0.2*inblock[*,*,13]
      ;outblock[*,*,3] = 0.7*inblock[*,*,6] + 0.7*inblock[*,*,7] - 0.2*inblock[*,*,14] - 0.2*inblock[*,*,15]

      outblock = fltarr(dims[0], dims[1], 8)
      outblock[*,*,0] = 1.4*inblock[*,*,0] - 0.4*inblock[*,*,8]
      outblock[*,*,1] = 1.4*inblock[*,*,1] - 0.4*inblock[*,*,9]
      outblock[*,*,2] = 1.4*inblock[*,*,2] - 0.4*inblock[*,*,10]   ; this does seem to do better than the above. 
      outblock[*,*,3] = 1.4*inblock[*,*,3] - 0.4*inblock[*,*,11]
      outblock[*,*,4] = 1.4*inblock[*,*,4] - 0.4*inblock[*,*,12] 
      outblock[*,*,5] = 1.4*inblock[*,*,5] - 0.4*inblock[*,*,13]
      outblock[*,*,6] = 1.4*inblock[*,*,6] - 0.4*inblock[*,*,14] 
      outblock[*,*,7] = 1.4*inblock[*,*,7] - 0.4*inblock[*,*,15]

      return, median(outblock,dimension=3,/even)
 
end



