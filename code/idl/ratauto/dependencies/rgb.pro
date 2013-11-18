function rgb, r, g, b
  ; returns the color integer for X window plotting.

  return, b*256*256L+g*256L+r
end

