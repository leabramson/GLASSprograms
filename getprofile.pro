function getProfile, di

  if 0 then begin
;  data = mrdfits(filename, 'DSCI')
     s    = size(di, /dim)
     xran = s[0]
     x    = findgen(xran)
     profile = total(di[20:40,20:40], 1)
     profile /= total(profile)
     
     ;; Centroid
     g = gaussfit(x, profile, nterms = 3, out)
     ctr = out[1]
     
     ;; Fold the light curve
     ctr = value_locate(x, ctr)
     gc = profile[ctr:*] + reverse(profile[0:ctr-1])
     gc = total(gc, /cum) / total(gc)
     
     ;; Get "half light radius"
     re = value_locate(gc, 0.5)
  endif
  
;  peak = mpfit2dpeak(di[20:40,20:40], /moffat, /tilt, out)
;  re   = out[3] * cos(out[-2])  ;; Output is in rotated coordinate system

  ;; Or just fit a moffat profile
  ;; SPECTRAL DIRECTION
  profile = total(di[20:40,20:40], 2)
  s       = size(di, /dim)
  lsf     = mpfitpeak(findgen(n_elements(profile)), profile, /moffat, out)
  tx      = findgen(s[0]) - (out[1] + 20)
  u       = tx / out[2]
  lsf     = out[0] / (u^2 + 1)^out[3]
  lsf /= total(lsf)
  
  ;; SPATIAL DIRECTION
  spatial = total(di[20:40,20:40], 1)
  prof    = mpfitpeak(findgen(n_elements(spatial)), spatial, /moffat, out)
  tx      = findgen(s[1]) - (out[1] + 20)
  u       = tx/out[2]
  profile = out[0] / (u^2 + 1)^out[3]
  tot     = total(profile,/cum) / total(profile)
  q1 = value_locate(tot, 0.25)
  q2 = value_locate(tot, 0.75)
  re = 0.5 * (q2 - q1)  ;;  Half the light's in here

  out[1]  = 0
  dat = {SPATIAL_PARAMS: out, RE: re, LSF: LSF}
  
  RETURN, dat; profile
end
