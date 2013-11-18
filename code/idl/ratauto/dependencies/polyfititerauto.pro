function polyfititerauto, xin, yin, xout, yout, targetdev=targetdev, initrejectdev=initrejectdev, finalrejectdev=finalrejectdev, alwaysrejectdev=alwaysrejectdev, startorder=startorder, maxorder=maxorder, verbosein=verbosein

   ; Iterative polynomial fit with automatic order determination

   ; initrejectdev - absolute deviation threshold to remove outliers.  With each outlier removed this drops towards finalthreshold.
   ; finalrejectdev - absolute deviation threshold to stop removing outliers.
   ; alwaysrejectdev - never increase the order if there is an outlier worse than this in the current order.
   ; targetdev - stop fitting when there are no outliers and the median absolute deviation is less than this.

   ; right now used only within solvearc.  Please use caution in modifying it to avoid changing behavior of that code.

   rejectdev = initrejectdev

   if n_elements(verbosein) eq 0 then verbose = 0 else verbose = verbosein

   outlierndev = 3.5

   x = xin
   y = yin

   if n_elements(startorder) eq 0 then startorder = 1
   if n_elements(maxorder) eq 0 then maxorder = n_elements(xin)-1

   if verbose then $
      print, 'Starting order ', startorder

   absdev = [999]
   nreject = 0
   order = startorder
   iter = 0
   skiprefit = 0
   while 1 do begin  ; median(absdev) gt targetdev and max(absdev) gt finalrejectdev
      iter += 1
 
      if skiprefit eq 0 then begin
         fitpar = polyfitg(x, y, order, yfit=yfit)
      endif
      n = n_elements(y)
      absdev = abs(yfit - y)
      meddev = median(absdev, /even)
      maxdev = max(absdev)
      sdev = stdev(yfit-y)
      dof = n-order
      outlier = where(absdev gt outlierndev*meddev, noutlier)
      bad = where(absdev gt rejectdev, nbad)

      ; Characterize outliers and investigate what would happen if order were increased.
      ; Use a voting system taking various measurements into account.

      trypar = polyfitg(x, y, order+1, yfit=tryyfit, guess=fitpar)
      tryabsdev = abs(tryyfit - y)
      trysdev   = stdev(tryyfit - y)
      trydof = n-order-1
      trymeddev = median(tryabsdev, /even)
      trymaxdev = max(tryabsdev)
      tryoutlier = where(tryabsdev gt outlierndev*trymeddev, trynoutlier)
      trybad = where(absdev gt rejectdev, trynbad)
      if verbose gt 0 then begin
         print, ' order: ', rclip(order,5)+'   ',  rclip(order+1,5), '     order+1:  reject:  stop:'
         print, 'meddev: ', fpr(meddev, 3.3), ' ', fpr(trymeddev,3.3), ' / ', clip(targetdev,7),    '          ',clip(targetdev)
         print, 'maxdev: ', fpr(maxdev, 3.3), ' ', fpr(trymaxdev,3.3), ' / ', clip(alwaysrejectdev,7), '  ', clip(rejectdev),'  ',clip(finalrejectdev)
         print, 'stddev: ', fpr(sdev,   3.3), ' ', fpr(trysdev,3.3)
         print, '#outlr: ', rclip(noutlier,5)+'   ', rclip(trynoutlier,5)+'   ', ' / '
         print, '#bad:   ', rclip(nbad,5)+'   ',     rclip(trynbad,5)+'   ', ' / '
      endif
      voteincorder = 0.3
      votereject = 0.2          ; decimals are tiebreakers
      votestop = 0.1
      if order eq 2 or order eq 3 then votereject += 0.5
      if order ge 3 then votestop += 0.5
      ;                                                       Vote to reject if...
      if maxdev gt rejectdev then votereject += 1                    ; an point exceeds current threshold
      if maxdev gt alwaysrejectdev then votereject += 2              ; an point exceeds always threshold
      if noutlier ge 1 and trynoutlier ge 1 then votereject += 2     ; an outlier exists now and if order increases
      if noutlier ge 2 then votereject += 1                          ; there are many outliers
                                                            ; Vote to increase order if..,
      if meddev lt targetdev then voteincorder += 1                  ; we have not reached target deviation
      if trymeddev lt meddev then voteincorder += 1                  ; it would improve the median deviation
      if trysdev/trydof lt sdev/dof then voteincorder += 1           ; it would improve chi^2/dof
                                                            ; Vote to stop if...
      if order ge maxorder then votestop += 10
      if maxdev lt finalrejectdev then votestop += 1                 ; there are no points worse than final tolerance
      if meddev lt targetdev then votestop += 1.0                    ; the median deviation is within tolerance 
      if meddev lt targetdev*0.5 then votestop += 0.5                ; the median deviation is well within tolerance
      if sdev lt targetdev then votestop += 0.5                      ; the standard deviation is within tolerance
      if noutlier eq 0 then votestop += 1                            ; there are no outliers at default sigma 
      if maxdev/meddev lt 2.5 then votestop += 1                     ; there are no outliers even at 2.5 sigma
      if order gt 8 then votestop += 1                               ; very high order solution
      if n lt order + 5 then votestop += 1                           ; we are approaching few degrees of freedom
      if n lt order + 2 then votestop += 1                           ; we are approaching few degrees of freedom
      if n lt order + 1 then votestop += 1                           ; we are approaching few degrees of freedom
      maxvote = max([voteincorder, votereject, votestop])
      if verbose gt 0 then $
         print, 'Votes:  reject ', clipzero(votereject), '  increase ', clipzero(voteincorder), '  stop ', clipzero(votestop)

      if votestop eq maxvote then begin
         if verbose gt 0 then begin
           print, 'End condition reached.'
           print, 'Final order ', clip(order)
         endif
         break
      endif
      skiprefit = 0
      if voteincorder eq maxvote then begin
         ; if so - raise the order officially
         order += 1
         if verbose gt 0 then print, 'Raising solution order to ', clip(order)
         fitpar = trypar
         yfit = tryyfit
         skiprefit = 1
      endif 
      if votereject eq maxvote then begin
         ; Remove the most discrepant point
         notworstoutlier = where(absdev lt max(absdev), ctok, complement=worstoutlier, ncomplement=nworst)

         if verbose gt 0 then begin
              print, '          Rejecting ', clip(x[worstoutlier[0]],7), '=', clip(y[worstoutlier[0]],8), $
                              ' (dev = '+clip(absdev[worstoutlier[0]])+' / '+clip(rejectdev)+')'
         endif
         x = x[notworstoutlier]
         y = y[notworstoutlier]
         rejectdev *= 0.9
         nreject += nworst
      endif

   endwhile

    xout = x
    yout = y

    return, fitpar

end

