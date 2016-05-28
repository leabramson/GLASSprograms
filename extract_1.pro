function extract_1, infile, z, $
                    PA = pa, $
;                    MASK = mask, $
                    REFIT = refit

  if N_ELEMENTS(z) eq 0 then z = -99
  if NOT keyword_SET(PA) then pa = -99
  if NOT keyword_set(MASK) then mask = 0 else mask = 1 ;; mask rather than subtract contamination
  if NOT keyword_set(REFIT) then refit = 0 else refit = 1  ;; refit continuum centerline
  
  ;; Read crap
  
  null   = mrdfits(infile, 0, thead, /silent)
  spec   = mrdfits(infile, 'sci', head, /silent)
  contam = mrdfits(infile, 'contam', /silent)
  sens   = mrdfits(infile, 'sens', /silent)
  trace  = mrdfits(infile, 'ytrace', /silent)
  err    = mrdfits(infile, 'wht', /silent)
  di     = mrdfits(infile, 'dsci', /silent)
  
  ;; Basic Info
  exptime = sxpar(thead, 'EXPTIME')   ;; SECONDS
  mag     = sxpar(thead, 'MAG')       ;; Probably F140W?

  ;; TO POISSON UNITS
  spec   *= exptime
  contam *= exptime
  var     = (err * exptime)^2
  
  ;; GET SPATIAL PROFILE/LSF
  prof    = getProfile(di)
  profile = prof.LSF       
  re      = ceil(prof.RE)        ;; Returns a pseudo "half-light radius"
  params  = prof.SPATIAL_PARAMS
  
  ;; MAKE A WAVELENGTH AXIS
  l0  = sxpar(head, 'CRVAL1')
  dl  = sxpar(head, 'CD1_1')
  run = sxpar(head, 'NAXIS1')
  lambda = findgen(run) * dl + l0

  yctr = sxpar(head, 'CRPIX2') - 1
  tx = findgen(sxpar(head, 'NAXIS2'))
  height = n_elements(tx)

  ;; MAKE A MASK IMAGE
  if mask then begin
     mask = contam
     mask[where(contam ge 2 * sqrt(var), compl = good)] = 1
     mask[good] = 0
  endif
  
  ;; MAKE A SCIENCE IMAGE
  if NOT mask then $
     sci   = (spec - contam) $
  else $
     sci = spec[where(~mask)]


  if NOT refit then begin
     ctroid = fltarr(2,run)
     ctroid[0,*] = trace
     ctroid[1,*] = re
  endif else begin
     ctroid = fltarr(4,run)
     for jj = 0, run - 1 do begin
        if 20 * mean(reform(sci[jj,*]) / sqrt(reform(var[jj,*]))) gt 10 then begin  ;;  "S/N in central region"
           g = gaussfit(tx, reform(sci[jj,*]), $
                        nterms = 3, out, $
                        sig = sig)
           ctroid[*,jj] = [out[1:2], sig[1:2]]
        endif
     endfor

     use = where(abs(ctroid[0,*] - height/2) le 10, nuse)
     tofit = use[where(abs(ctroid[0,use] - median(ctroid[0,use])) $
                       le 2 * stddev(ctroid[0,use]))]
     cline = poly_fit(lambda[tofit], ctroid[0,tofit], $
                      measure_err = ctroid[2, tofit], 1)
     trace = cline[0] + lambda * cline[1]
  endelse

  modIm = sci  ;; Will contain the profile for optimal extraction
  for jj = 0, n_elements(trace) - 1 do begin
     u = (tx - trace[jj]) / params[2]
     modIm[jj,*] = 1. / (u^2 + 1)^params[3]
     modIm[jj,*] /= total(modIm[jj,*])
  endfor
  extrIm = sci
  
  ;;  DEFINE EXTRACTION RADIAL REGIMES
  innerup = trace   + 0.25 * re
  innerdn = trace   - 0.25 * re
  interup = innerup + 0.50 * re
  interdn = innerdn - 0.50 * re
  outerup = interup + 0.75 * re
  outerdn = interdn - 0.75 * re

  skies  = dblarr(run)
  integ  = dblarr(run)
  outer  = dblarr(run)
  inter  = dblarr(run)
  inner  = dblarr(run)
  
  ivar   = dblarr(run)
  ovar   = dblarr(run)
  intvar = dblarr(run)
  innvar = dblarr(run)

  nAlls = intarr(run)
  nInns = intarr(run)
  nInts = intarr(run)
  nOuts = intarr(run)
  
  ;;  DO THE EXTRACTIONS
  for jj = 0, run - 1 do begin
     ;; Do the optimal extraction
     goods    = where(~mask[jj,*])
     opti     = total(modim[jj,*]*sci[jj,*]/var[jj,*]) $
                / total(modim[jj,*]^2/var[jj,*])
     opti_var = 1. / total(modim[jj,*]/var[jj,*])
;     optiMask = total(modim[jj,goods]*spec[jj,goods]/var[jj,goods]) $
;                / total(modim[jj,goods]^2/var[jj,goods])
;     optiMask_Var = 1. / total(modim[jj,goods]/var[jj,goods])

     ;; Do the binned extractions
     all = where(tx ge outerdn[jj] AND tx le outerup[jj], nall)
     out = where((tx ge interup[jj] AND tx lt outerup[jj]) $
                 OR $
                 (tx gt outerdn[jj] AND tx le interdn[jj]), nout)
     int = where((tx ge innerup[jj] AND tx lt interup[jj]) $
                 OR $
                 (tx ge interdn[jj] AND tx lt innerdn[jj]), nint)
     in  = where(tx gt innerdn[jj] AND tx lt innerup[jj], nin)
     
     sky = mean([sci[jj,0:10], sci[jj,-11:*]])

     skies[jj]  = sky
     integ[jj]  = total(sci[jj,all]) / nall
     outer[jj]  = total(sci[jj,out]) / nout
     inter[jj]  = total(sci[jj,int]) / nint
     inner[jj]  = total(sci[jj,in ]) / nin 
     
     ivar[jj]   = total(var[jj,all]) / nall
     ovar[jj]   = total(var[jj,out]) / nout
     intvar[jj] = total(var[jj,int]) / nint
     innvar[jj] = total(var[jj,in ]) / nin 

     nAlls[jj] = nall
     nInns[jj] = nin
     nInts[jj] = nint
     nOuts[jj] = nout
     
     extrIm[jj,out] = 10
     extrIm[jj,int] = 20
     extrIm[jj,in ] = 30
     
  endfor       

  ;;  WRITE OUTPUT
  savedata = {LAMBDA:      lambda, $
              Z:           z, $
              RA:          sxpar(thead, 'RA'), $
              DEC:         sxpar(thead, 'DEC'), $
              MAG:         mag, $
              F_TOT:       integ, $
              F_INNER:     inner, $
              F_INTER:     inter, $
              F_OUTER:     outer, $
              VAR_TOT:     ivar, $
              VAR_INNER:   innvar, $
              VAR_INTER:   intvar, $
              VAR_OUTER:   ovar, $
              NPIX_TOT:    nAlls, $
              NPIX_INNER:  nInns, $
              NPIX_INTER:  nInts, $
              NPIX_OUTER:  nOuts, $
              SENSITIVITY: sens, $
              EXPTIME:     exptime, $
              FILENAME:    infile, $
              MODEL_TRACE: modim, $
              EXTRACT_ZONES: extrIM, $
;              EXTRACT_MIDPTS: , $
              LSF: prof.LSF, $
              RE: prof.RE, $
              SPEC2D: sci, $
              DIRECT_IM: di, $
              TRACE: trace, $
              INNERUP: innerup, $
              INNERDN: innerdn, $
              INTERUP: interup, $
              INTERDN: interdn, $
              OUTERUP: outerup, $
              OUTERDN: outerdn, $
              PA: pa, $
              MASK: mask}
;  mwrfits, savedata, pa+'_extractions/'+string(ii, f = '(I04)')+'_radStuff.fits', /create
  
  RETURN, savedata
end
