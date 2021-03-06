;; Convert a GLASS Master ASCII Catalog to a FITS table

pro readGlassList, incat, outname

;NUMBER
;X_IMAGE Y_IMAGE
;X_WORLD Y_WORLD
;A_IMAGE B_IMAGE THETA_IMAGE
;A_WORLD B_WORLD THETA_WORLD
;FLUX_APER FLUX_APER2 FLUX_APER3
;FLUXERR_APER FLUXERR_APER2 FLUXERR_APER3
;MAG_APER MAG_APER2 MAG_APER3
;MAGERR_APER MAGERR_APER2 MAGERR_APER3
;FLUX_AUTO FLUXERR_AUTO
;MAG_AUTO MAGERR_AUTO

;KRON_RADIUS PETRO_RADIUS
;BACKGROUND THRESHOLD
;XWIN_IMAGE YWIN_IMAGE
;AWIN_IMAGE BWIN_IMAGE THETAWIN_IMAGE
;MU_THRESHOLD FLAGS FWHM_IMAGE
;FLUX_RADIUS FLUX_RADIUS2 CLASS_STAR X_FLT Y_FLT
  
  ;; Do fluxes (mainly)
  print, ''
  print, 'Reading MASTER GLASS catalog ... '
  print, ''
  readcol, incat, $
           NUMBER, $
           X_IMAGE, Y_IMAGE, $
           X_WORLD, Y_WORLD, $
           A_IMAGE, B_IMAGE, THETA_IMAGE, $
           A_WORLD, B_WORLD, THETA_WORLD, $
           FLUX_APER, FLUX_APER2, FLUX_APER3, $
           FLUXERR_APER, FLUXERR_APER2, FLUXERR_APER3, $
           MAG_APER, MAG_APER2, MAG_APER3, $
           MAGERR_APER, MAGERR_APER2, MAGERR_APER3, $
           FLUX_AUTO, FLUXERR_AUTO, $
           MAG_AUTO, MAGERR_AUTO, $
           f = $
           'F,'+$
           'F,F,'+$
           'D,D,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,'+$
           'F,F'
  number = long(number)
  
  ;; Do the rest
  readcol, incat, $
           KRON_RADIUS, PETRO_RADIUS, $
           BACKGROUND, THRESHOLD, $
           XWIN_IMAGE, YWIN_IMAGE, $
           AWIN_IMAGE, BWIN_IMAGE, THETAWIN_IMAGE, $
           MU_THRESHOLD, FLAGS, FWHM_IMAGE, $
           FLUX_RADIUS, FLUX_RADIUS2, CLASS_STAR, X_FLT, Y_FLT, $
           f = 'X,'+$
           'X,X,'+$
           'X,X,'+$
           'X,X,X,'+$
           'X,X,X,'+$
           'X,X,X,'+$
           'X,X,X,'+$
           'X,X,X,'+$
           'X,X,X,'+$
           'X,X,'+$
           'X,X'+$  
           'F,F,'+$
           'F,F,'+$
           'F,F,'+$
           'F,F,F,'+$
           'F,L64,F,'+$
           'F,F,F,F,F'

  savedata = {$
             NUMBER:        0L  ,$
             X_IMAGE:       0.  ,$
             Y_IMAGE:       0.  ,$
             X_WORLD:       0.d ,$
             Y_WORLD:       0.d ,$
             A_IMAGE:       0.  ,$
             B_IMAGE:       0.  ,$
             THETA_IMAGE:   0.  ,$
             A_WORLD:       0.  ,$
             B_WORLD:       0.  ,$
             THETA_WORLD:   0.  ,$
             FLUX_APER:     0.  ,$
             FLUX_APER2:    0.  ,$
             FLUX_APER3:    0.  ,$
             FLUXERR_APER:  0.  ,$
             FLUXERR_APER2: 0.  ,$
             FLUXERR_APER3: 0.  ,$
             MAG_APER:      0.  ,$
             MAG_APER2:     0.  ,$
             MAG_APER3:     0.  ,$
             MAGERR_APER:   0.  ,$
             MAGERR_APER2:  0.  ,$
             MAGERR_APER3:  0.  ,$
             FLUX_AUTO:     0.  ,$
             FLUXERR_AUTO:  0.  ,$
             MAG_AUTO:      0.  ,$
             MAGERR_AUTO:   0.  ,$
             KRON_RADIUS:   0.  ,$
             PETRO_RADIUS:  0.  ,$
             BACKGROUND:    0.  ,$
             THRESHOLD:     0.  ,$
             XWIN_IMAGE:    0.  ,$
             YWIN_IMAGE:    0.  ,$
             AWIN_IMAGE:    0.  ,$
             BWIN_IMAGE:    0.  ,$
             THETAWIN_IMAGE:0.  ,$
             MU_THRESHOLD:  0.  ,$
             FLAGS:  long64(0)  ,$
             FWHM_IMAGE:    0.  ,$
             FLUX_RADIUS:   0.  ,$
             FLUX_RADIUS2:  0.  ,$
             CLASS_STAR:    0.  ,$
             X_FLT:         0.  ,$
             Y_FLT:         0.   $
             }
  savedata = replicate(savedata, n_elements(X_WORLD))

  for ii = 0, n_elements(savedata) - 1 do begin
     savedata[ii].NUMBER        = NUMBER[ii]
     savedata[ii].X_IMAGE       = X_IMAGE[ii]      
     savedata[ii].Y_IMAGE       = Y_IMAGE[ii]  
     savedata[ii].X_WORLD       = X_WORLD[ii]            
     savedata[ii].Y_WORLD       = Y_WORLD[ii]
     savedata[ii].A_IMAGE       = A_IMAGE[ii]           
     savedata[ii].B_IMAGE       = B_IMAGE[ii]           
     savedata[ii].THETA_IMAGE   = THETA_IMAGE[ii]   
     savedata[ii].A_WORLD       = A_WORLD[ii]             
     savedata[ii].B_WORLD       = B_WORLD[ii]             
     savedata[ii].THETA_WORLD   = THETA_WORLD[ii]    
     savedata[ii].FLUX_APER     = FLUX_APER[ii]         
     savedata[ii].FLUX_APER2    = FLUX_APER2[ii]    
     savedata[ii].FLUX_APER3    = FLUX_APER3[ii]       
     savedata[ii].FLUXERR_APER  = FLUXERR_APER[ii]  
     savedata[ii].FLUXERR_APER2 = FLUXERR_APER2[ii]
     savedata[ii].FLUXERR_APER3 = FLUXERR_APER3[ii] 
     savedata[ii].MAG_APER      = MAG_APER[ii]           
     savedata[ii].MAG_APER2     = MAG_APER2[ii]         
     savedata[ii].MAG_APER3     = MAG_APER3[ii]         
     savedata[ii].MAGERR_APER   = MAGERR_APER[ii]     
     savedata[ii].MAGERR_APER2  = MAGERR_APER2[ii]   
     savedata[ii].MAGERR_APER3  = MAGERR_APER3[ii]  
     savedata[ii].FLUX_AUTO     = FLUX_AUTO[ii]         
     savedata[ii].FLUXERR_AUTO  = FLUXERR_AUTO[ii]   
     savedata[ii].MAG_AUTO      = MAG_AUTO[ii]           
     savedata[ii].MAGERR_AUTO   = MAGERR_AUTO[ii]     
     savedata[ii].KRON_RADIUS   = KRON_RADIUS[ii]     
     savedata[ii].PETRO_RADIUS  = PETRO_RADIUS[ii]   
     savedata[ii].BACKGROUND    = BACKGROUND[ii]    
     savedata[ii].THRESHOLD     = THRESHOLD[ii]       
     savedata[ii].XWIN_IMAGE    = XWIN_IMAGE[ii]      
     savedata[ii].YWIN_IMAGE    = YWIN_IMAGE[ii]    
     savedata[ii].AWIN_IMAGE    = AWIN_IMAGE[ii]        
     savedata[ii].BWIN_IMAGE    = BWIN_IMAGE[ii]        
     savedata[ii].THETAWIN_IMAGE= THETAWIN_IMAGE[ii]
     savedata[ii].MU_THRESHOLD  = MU_THRESHOLD[ii] 
     savedata[ii].FLAGS         = FLAGS[ii]  
     savedata[ii].FWHM_IMAGE    = FWHM_IMAGE[ii]        
     savedata[ii].FLUX_RADIUS   = FLUX_RADIUS[ii]
     savedata[ii].FLUX_RADIUS2  = FLUX_RADIUS2[ii]    
     savedata[ii].CLASS_STAR    = CLASS_STAR[ii]        
     savedata[ii].X_FLT         = X_FLT[ii]                  
     savedata[ii].Y_FLT         = Y_FLT[ii] 
  endfor

  print, 'Saving to FITS ... '
  print, ''
  mwrfits, savedata, outname, /create
  print, ' >>> New FITS catalog output to : '+outname
  print, ''
  
end

  
