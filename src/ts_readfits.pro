FUNCTION TS_READFITS, filename, extension=extension, $
            t=t, dt=dt, start_t=start_t, scale=scale, cols=cols

; ----------------------------------------------------------
;+
; NAME:
;       TS_READFITS 
;
; PURPOSE:
;       Read two columns data (e.g. time, flux) from a FITS file
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       data = TS_READFITS('mydata.fits')
;
; INPUTS:
;       filename  - (string) file name  
;
; OPTIONAL INPUTS:
;       extension - (integer) FITS extension to read (default 1)
;       scale     - (logical) if set then scale times to units of dt
;                             starting at zero.
;       cols      - (vector) two numbers specifying columns to read 
;                            [time,flux]. default is [1,2]
;
; OUTPUTS:
;       x         - 1 dimensional time series (vector)
;
; OPTIONAL OUTPUTS:
;       t         - (vector) sample times 
;       dt        - (scalar) sampling interval
;       start_t   - (scalar) start time
;
; PROCEDURES CALLED:
;       READFITS, TBGET, SXPAR
;
; DETAILS:
;       Read two columns of data from a FITS file. Use for reading
;       time series (time,flux) from FITS files. Based around 
;       READFITS which can work on gzip or Unix compressed FITS files.
;
; Example calls:
; IDL> x = ts_readfits('mydata.fits',t=t,dt=dt)
; IDL> y = ts_readfits('mydata.FTZ',t=t,dt=dt,/scale)
;
; HISTORY:
;       25/04/07  - v1.0 - first working version
;
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments

; is the file name defined?
  if (n_elements(filename) eq 0) then begin
      filename=''
      read,'-- Enter file name (ENTER to list current directory): ',filename
      if (filename eq '') then begin
          list = findfile()
          print, list
          read,'-- Enter file name: ',filename
      endif
  endif

; is the extension parameter defined?
  if (n_elements(extension) eq 0) then extension=1

; is the cols parameter defined?
  if (n_elements(cols) eq 0) then cols=[1,2]

; ----------------------------------------------------------
; Extract the data from file into arrays

; call READFITS to read the data from FITS file
  data=readfits(filename,htab,EXTEN_NO=extension,/SILENT)

; check for errors from READFITS
  if (n_elements(data) eq 1) then begin
      if (data eq -1) then begin
          print,'** Error calling READFITS in TS_READFITS'
          print,'** '+!ERROR_STATE.MSG
          return,0
      endif
  endif

; now convert first column from binary to floating point
  t=tbget(htab,data,cols[0])

; likewise for second column
  x=tbget(htab,data,cols[1])

; record, then subtract the start time 
  start_t=t[0]

; record the sampling period dT
  dt=sxpar(htab,'TIMEDEL',count=count)
  if (count ne 1) then begin
      print,'** No TIMEDEL parameter in file. Setting dt=1'
      dt=1.0
  endif

; if SCALE is set then convert time into units of dT from t0
  if keyword_set(scale) then $
    t=temporary((t-start_t)/dt)

; ----------------------------------------------------------
; Return the data array to the user

  return,x

END
