pro makemaster_roman, FIELD = field, $
                      POINTING = pointing, $
                      OUTDIR = outdir

  catDir = '$DATA_DIR/GLASS_official_products/'+field+'/Catalogs'

  master    = outdir+'/'+field+'_glass_srclist.fits'    
  redshifts = outdir+'/'+field+'_glass_redshifts.fits'  
  gigCat    = outdir+'/'+field+'_glass_GiG_results.fits'
  photozOut = outdir+'/'+field+'_Kuang_photozs.fits'
  
  readglasslist, catDir+'/*glassmaster.cat'    , master
  readredshift , catDir+'/*redshiftcatalog.txt', redshifts
  readgig      , catDir+'/*glassgigcatalog.txt', gigCat

  photozs   = catDir+'/*_photoz.fits'
  translatekuang, photozs, photozOut

  print, ''
  print, 'APPENDING GiG RESULTS TO MASTER GLASS SOURCE CATALOG ... '
  print, ''
  join2cats, master, gigCat, 'tmp1.fits', $
             striptag = 'GLASS_ID'
  print, ''
  print, 'APPENDING GiGZ RESULTS TO GLASS SOURCE + GiG CATALOG ... '
  print, ''
  join2cats, 'tmp1.fits', redshifts, 'tmp2.fits', $
             striptag = ['GLASS_ID', 'RA', 'DEC']
  print, ''
  print, 'APPENDING KUANG/ROMAN PHOTO-ZS TO GLASS SOURCE + GiG + GiGz CATALOG ... '
  print, ''
  join2cats, 'tmp2.fits', photozOut, outdir+'/'+field+'_allAvailableData.fits', $
             striptag = ['RA', 'DEC'], $
             dposname = 'ROMAN2GLASS_DPOS', $
             POINTING = pointing
  
  print, ''
  print, ' >>> ALL DONE! <<< '
  
end
