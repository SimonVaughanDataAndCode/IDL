PRO write_table, data, filename, WIDTH=width, DOUBLE=double, $
                                 _EXTRA=EXTRA_KEYWORDS, $
                                 INTEGER=integer

; ----------------------------------------------------------
;+
; NAME:
;       WRITE_TABLE
;
; PURPOSE:
;       Write a data array to an ASCII table file.
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       WRITE_TABLE, data, filename='file.dat'
;
; INPUTS:
;       data    - (array) 1 or 2 dimensional data array
;       file    - (string) file name
;
; OPTIONAL INPUTS:
;       width   - (integer) width of lines in characters
;       double  - (logical) whether to use double prec. format
;       integer - (logical) whether to use integer format
;
; OUTPUTS:
;       ASCII file
;
; DETAILS:
;       The 1d or 2d 'data' array is written to an ASCII file.
;       If 1-dimensional it is written as a column vectors,
;       otherwise the array is written 'as is'.
;       If there are many columns the width parameter may
;       need to be increased to fit the columns on one line.
;       Note, data are simple type - ie: all floating point,
;       all integers, etc... no mixing is allowed (yet).
;
;       Using the DOUBLE keyowrd will force all columns to
;       be printed in double precision format.
;
;       If more specific formats are needed, they may be
;       passed using the _EXTRA keyword (with DOUBLE not set).
;
; EXAMPLE CALL:
;  IDL> WRITE_TABLE, double(x), 'temp.txt', /DOUBLE
;
; HISTORY:
;       01/02/2007 -  v1.0  - first working version
;       27/04/2007 -  v1.1  - added EXTRA and DOUBLE options
;       07/07/2013 -  v1.2 - changed FORMAT to Gw.d (from Dw.d)
;       10/03/2014 -  v1.3 - added INTEGER format output
;
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)

  COMPILE_OPT idl2

; watch out for errors

  ON_ERROR, 0

; ----------------------------------------------------------
; Check the arguments

  if (N_ELEMENTS(data) eq 0) THEN BEGIN
      PRINT, '** No data in WRITE_TABLE.'
      RETURN
  ENDIF

; supply default filename if none supplied

  if (N_ELEMENTS(filename) eq 0) then filename='idl.out'

; supply default width if none supplied

  if (N_ELEMENTS(width) eq 0) then width=150

; use SINGLE format if DOUBLE and INTEGER are not being used

  if (KEYWORD_SET(double) eq 0 and KEYWORD_SET(integer) eq 0) THEN single=1
  
; ----------------------------------------------------------
; Prepare for writing data array

; determine dimensions of data array
  s = SIZE(data)
  d = s[0]

  if (d eq 2) THEN BEGIN
      m=s[1]                    ; columns if 2d
      n=s[2]                    ; rows 
      type=s[3]                 ; variable type
  ENDIF ELSE BEGIN
      m=1
      n=s[1]                    ; rows if 1d
      type=s[2]                 ; variable type
  ENDELSE

  if (d gt 2) THEN BEGIN
      PRINT, '** Data array has >2 dimensions in WRITE_TABLE.'
      RETURN
  ENDIF

; ----------------------------------------------------------
; Open ASCII file, write data, then close file

  OPENW, lun, filename, /GET_LUN, WIDTH=width
  if KEYWORD_SET(double) THEN BEGIN

      if (d eq 2) THEN BEGIN
          PRINTF, lun, data, FORMAT='('+strtrim(m,2)+'(G30.16))',_EXTRA=extra_keywords
      ENDIF
      if (d eq 1) THEN BEGIN
          PRINTF, lun, TRANSPOSE(data), FORMAT='(G)', _EXTRA=extra_keywords
      ENDIF

  ENDIF 
  
  if KEYWORD_SET(single) THEN BEGIN

      if (d eq 2) THEN BEGIN
          PRINTF, lun, data, FORMAT='('+strtrim(m,2)+'(G18.7))', _EXTRA=extra_keywords
      ENDIF
      if (d eq 1) THEN BEGIN
          PRINTF, lun, TRANSPOSE(data), FORMAT='(G)', _EXTRA=extra_keywords
      ENDIF

  ENDIF

  if KEYWORD_SET(integer) THEN BEGIN

      if (d eq 2) THEN BEGIN
          PRINTF, lun, data, FORMAT='('+strtrim(m,2)+'(I10))', _EXTRA=extra_keywords
      ENDIF
      if (d eq 1) THEN BEGIN
          PRINTF, lun, TRANSPOSE(data), FORMAT='(I10)', _EXTRA=extra_keywords
      ENDIF

  ENDIF

  FREE_LUN, lun

; ----------------------------------------------------------
; Return to user

END
