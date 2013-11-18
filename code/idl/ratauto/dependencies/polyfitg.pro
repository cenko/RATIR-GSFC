function polyfitg, x, y, degree, guess=guessin, yfit=yfit, _EXTRA=ex

; polyfit with initial guess.
; fits to the residuals.

if n_elements(guess) eq 0 then guess = fltarr(degree+1)
if (size(guess)) [0] gt 1 and (size(guess)) [2] gt 1 then guess = transpose(guess)
while n_elements(guess) lt degree+1 do guess = [guess, 0.]

model = poly(x, guess)
residual = y - model

result = poly_fit(x, residual, degree, yfit=ryfit, _EXTRA=ex)
if (size(result)) [0] gt 1 and (size(result)) [2] gt 1 then result = transpose(result)

yfit = model + ryfit

return, guess + result

end

