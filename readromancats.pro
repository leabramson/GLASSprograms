;; Read in both the roman photmoetry file AND the derived photo-zs,
;; merge them, and spit-out a FITS binary table

pro readRomanCats, photometryCat, redshiftCat, $
                   OUTFITS = outfits

  ;; Read the photometry catalog
  readcol, photometryCat, $
           ID, RA, DEC, $
           B435, V606, I814, $
           Y105, J125, JH140,$
           H160, Ks, CH1,$
           CH2 ,$
           errB435, errV606, errI814,$
           errY105, errJ125, errJH140,$
           errH160, errKs, errCH1, errCH2, $
           f = 'L,D,D,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,F', comment = '#'
  ngals = n_elements(ID)
  
  ;; Read the photo-z catalog
  readcol, redshiftCat, $
           zID, ZBEST, $
           ZSPEC_FLAG, ZSPEC_ID, $
           PHOTOZ_FLAG, $
           f = 'L,'+$
           'F,F,'+$
           'F,F,'+$
           'F', comment = '#', /fast

  ;; Check that they have the same number of objects.
  ;; If so, dump to output.
  if total(ID - zID) ne 0 then begin
     print, ''
     print, ' !!! CATALOGS DO NOT ALIGN --> ABORTING !!!'
     stop
  endif else begin
     print, ''
     print, ' >>> Catalogs look aligned; dumping to output ... '

     savedata = {$
                ID:               0L , $
                RA:               0.d, $
                DEC:              0.d, $
                ZBEST_ROME:       0. , $
                ZSPEC_FLAG_ROME:  0  , $
                ZSPEC_ID_ROME:    0L , $
                PHOTOZ_FLAG_ROME: 0B , $
                B435R:            0. , $
                V606R:            0. , $
                I814R:            0. , $
                Y105R:            0. , $
                J125R:            0. , $
                JH140R:           0. , $
                H160R:            0. , $
                KsR:              0. , $
                CH1R:             0. , $
                CH2R:             0. , $
                errB435R:         0. , $
                errV606R:         0. , $
                errI814R:         0. , $
                errY105R:         0. , $
                errJ125R:         0. , $
                errJH140R:        0. , $
                errH160R:         0. , $
                errKsR:           0. , $
                errCH1R:          0. , $
                errCH2R:          0.   $
                }
     
     savedata = replicate(savedata, ngals)
     for ii = 0, ngals - 1 do begin
        savedata[ii].ID               = ID[ii]
        savedata[ii].RA               = RA[ii]              
        savedata[ii].DEC              = DEC[ii]                           
        savedata[ii].ZBEST_ROME       = ZBEST_ROME[ii]             
        savedata[ii].ZSPEC_FLAG_ROME  = ZSPEC_FLAG_ROME[ii]   
        savedata[ii].ZSPEC_ID_ROME    = ZSPEC_ID_ROME[ii]       
        savedata[ii].PHOTOZ_FLAG_ROME = PHOTOZ_FLAG_ROME[ii] 
        savedata[ii].B435R            = B435R[ii]                       
        savedata[ii].V606R            = V606R[ii]                       
        savedata[ii].I814R            = I814R[ii]                       
        savedata[ii].Y105R            = Y105R[ii]                       
        savedata[ii].J125R            = J125R[ii]                       
        savedata[ii].JH140R           = JH140R[ii]                     
        savedata[ii].H160R            = H160R[ii]                       
        savedata[ii].KsR              = KsR[ii]                           
        savedata[ii].CH1R             = CH1R[ii]                         
        savedata[ii].CH2R             = CH2R[ii]                         
        savedata[ii].errB435R         = errB435R[ii]                 
        savedata[ii].errV606R         = errV606R[ii]                 
        savedata[ii].errI814R         = errI814R[ii]                 
        savedata[ii].errY105R         = errY105R[ii]                 
        savedata[ii].errJ125R         = errJ125R[ii]                 
        savedata[ii].errJH140R        = errJH140R[ii]               
        savedata[ii].errH160R         = errH160R[ii]                 
        savedata[ii].errKsR           = errKsR[ii]                     
        savedata[ii].errCH1R          = errCH1R[ii]                   
        savedata[ii].errCH2R          = errCH2R[ii]         
     endfor
     
     mwrfits, savedata, outfits, /create
     print, ''
     print, ' >>> Catalog output to : '+outfits

  endelse

end
