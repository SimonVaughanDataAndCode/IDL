FUNCTION SPEC_READFITS, filename, extension=extension, cols=cols, $
                        double=double, texp=texp

; ----------------------------------------------------------
;+
; NAME:
;       SPEC_READFITS 
;
; PURPOSE:
;       Read specified columns of data (e.g. energy, counts) from a FITS file
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       data = SPEC_READFITS('mydata.fits')
;
; INPUTS:
;       filename  - (string) file name  
;
; OPTIONAL INPUTS:
;       extension - (integer) FITS extension to read (default 1)
;       cols      - (vector) two numbers specifying columns to read 
;                            Default is [1,2]
;       double    - (logical) use double precision?
;
; OUTPUTS:
;       x         - Two dimensional data array (N*M)
;
; OPTIONAL OUTPUTS:
;       texp      - Expsure time (from EXPOSURE keyword)
;
; PROCEDURES CALLED:
;       READFITS, TBGET, SXPAR
;
; DETAILS:
;       Read columns of data from a FITS file. Output is a (floating?)
;       array of dimensions N*M where N = number of "rows" and M =
;       number of "columns" in the FITS file.
;       Based around READFITS by W.B. Landsman which can work on gzip
;       or Unix compressed FITS files. 
;
; Example calls:
;
; IDL> data = spes_readfits('mydata.fits', extension=1, cols=[2,3,5])
;
; HISTORY:
;       26/11/09  - v1.0 - first working version
;       15/12/09  - v1.1 - added TEXP output
;
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)

  COMPILE_OPT idl2  

; watch out for errors

  on_error,2

; ----------------------------------------------------------
; Check the arguments

; is the file name defined?

  if (n_elements(filename) eq 0) then begin
      filename = ''
      read,'-- Enter file name (ENTER to list current directory): ',filename
      if (filename eq '') then begin
          list = findfile()
          print, list
          read,'-- Enter file name: ',filename
      endif
  endif

; is the extension parameter defined?

  if (N_ELEMENTS(extension) eq 0) then extension = 1

; is the cols parameter defined?

  m = N_ELEMENTS(cols)
  if (m eq 0) then begin
      cols = [1, 2]
      m = 2
  endif

; ----------------------------------------------------------
; Extract the data from file into arrays

; call READFITS to read the data from FITS file

  filedata = READFITS(filename, htab, EXTEN_NO=extension, /SILENT)

; check for errors from READFITS

  if (N_ELEMENTS(filedata) le 1) then begin
      if (filedata eq -1) then begin
          print,'** Error calling READFITS in SPEC_READFITS'
          print,'** '+!ERROR_STATE.MSG
          return, 0
      endif
  endif

; how many rows are there in the file?

  s = SIZE(filedata, /DIMENSIONS)
  n = s[1]

; prepare data array for output

  data = MAKE_ARRAY(n, m, DOUBLE=KEYWORD_SET(double))

; extract the M columns from the binary table

  for i = 0, m-1 do begin

; now convert columns from binary to floating point

      x = TBGET(htab, filedata, cols[i])
      data[0,i] = x[*]

; extract the exposure time

  texp = SXPAR(htab,"EXPOSURE")

  endfor

; ----------------------------------------------------------
; Return the data array to the user

  return, data

END
