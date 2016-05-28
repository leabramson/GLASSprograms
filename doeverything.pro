;; Cull and get the files for 1 field

pro doEverything, field, $
                  ZMIN = zmin, $
                  ZMAX = zmax, $
                  MINMAG = minmag, $
                  MAXMAG = maxmag, $
                  DOCONTAM = docontam, $
                  CLEVEL = clevel, $
                  BOTHPA = bothPa, $
                  OUTPREFIX = outprefix
  
  if NOT keyword_set(ZMIN) then zmin = 1.0
  if NOT keyword_set(ZMAX) then zmax = 1.8
  if NOT keyword_set(MINMAG) then minmag = 14
  if NOT keyword_set(MAXMAG) then maxmag = 21.5
  if NOT keyword_set(DOCONTAM) then docontam = 0 else docontam = 1
  if docontam then begin
     if NOT keyword_set(CLEVEL) then clevel = 1.5
     if NOT keyword_set(BOTHPA) then bothPA = 0 else bothPA = 1
  endif
  
  data_dir   = '../MasterCats_V1'
  cat_suffix = '_allAvailableData.fits' ;'_MASTER_z_GiG_Kpz.fits'

  case field of
     'a370'    : field = 'ABEL0370'
     'a0370'   : field = 'ABEL0370'
     'abel370' : field = 'ABEL0370'
     'abel0370': field = 'ABEL0370'
     'ABEL370' : field = 'ABEL0370'
     'ABEL0370': field = 'ABEL0370'

     'a2744'   : field = 'ABEL2744'
     'abel2744': field = 'ABEL2744'
     'ABEL2744': field = 'ABEL2744'
     
     'm0416'   : field = 'MACS0416'
     'macs0416': field = 'MACS0416'
     'MACS0416': field = 'MACS0416'
     
     'm0717'   : field = 'MACS0717'
     'macs0717': field = 'MACS0717'
     'MACS0717': field = 'MACS0717'
     
     'm0744'   : field = 'MACS0744'
     'macs0744': field = 'MACS0744'
     'MACS0744': field = 'MACS0744'

     'm1149'   : field = 'MACS1149'
     'macs1149': field = 'MACS1149'
     'MACS1149': field = 'MACS1149'

     'm1423'   : field = 'MACS1423'
     'macs1423': field = 'MACS1423'
     'MACS1423': field = 'MACS1423'
     
     'm2129'   : field = 'MACS2129'
     'macs2129': field = 'MACS2129'
     'MACS2129': field = 'MACS2129'

     'r1347'   : field = 'RXJC1347'
     'rxj1347' : field = 'RXJC1347'
     'rxjc1347': field = 'RXJC1347'
     'RXJ1347' : field = 'RXJC1347'
     'RXJC1347' : field = 'RXJC1347'

     'r2248'   : field = 'RXJC2248'
     'rxj2248' : field = 'RXJC2248'
     'rxjc2248': field = 'RXJC2248'
     'RXJ2248' : field = 'RXJC2248'
     'RXJC2248' : field = 'RXJC2248'
  endcase

  cat = data_dir+'/'+field+cat_suffix

  if NOT keyword_set(OUTPRIFIX) then begin
     zzmin = string(zmin  , f = '(F3.1)')
     zzmax = string(zmax  , f = '(F3.1)')
     lomag = string(minmag, f = '(F4.1)')
     himag = string(maxmag, f = '(F4.1)')

     if docontam then begin
        cclevel = 'c'+string(clevel, f = '(F3.1)')
        bbpa    = 'pacut'+string(bothPa, f = '(I1)')
     endif else begin
        cclevel = 'c99.9'
        bbpa    = 'pacut0'
     endelse

     outprefix = field+'_'+zzmin+'_'+zzmax+'_'$
                 +lomag+'_'+himag+'_'$
                 +cclevel+'_'+bbpa
     
  endif
  
  cutOnZ, cat, zmin = zmin, zmax = zmax, $
          /photoz, $
          docontam = docontam, clevel = clevel, bothPa = bothpa, $
          output = 'tmpZ.list'
  cutOnMag, 'tmpZ.list', minMag = minmag, maxMag = maxmag, $
            output = outprefix+'_FILELIST.list'
  
  toGet = dump1field(outprefix+'_FILELIST.list', field)

  stop
  
end
