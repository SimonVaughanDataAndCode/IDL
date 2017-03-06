FUNCTION READ_TABLE, filename, $
                     columns=columns, nrows=nrows, $
                     nmax=nmax, double=double, $
                     text=text, head=head

; ----------------------------------------------------------
;+
; NAME:
;       READ_TABLE
;
; PURPOSE:
;       Read an ASCII table into a data array.
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       data = READ_TABLE('file.dat')
;
; INPUTS:
;       file - (string) file name 
;
; OPTIONAL INPUTS:  
;       columns - (integer vector) which columns of table to retain 
;       nrows   - (integer) number of lines to read 
;       nmax    - (integer) minimum size of file, default is 100,000
;       double  - (logical) whether to use double or single prec.
;       text    - (logical) whether to load text from file       
;       head    - (integer) number of lines to skip in header
;
; OUTPUTS:
;       2-dimensional data array (floating point values)
;
; DETAILS:
;       Data are assumed to be floating point numbers, by
;       default single precision, seperated by spaces in a
;       rectangular grid (ncols*nrow)
;
; Example calls:
;
; IDL> data = read_table('table.txt')
; IDL> data = read_table('table.txt', columns=[2,3])
; IDL> data = read_table('table.txt', n=100, HEAD=1, /TEXT)
;
; HISTORY:
;       Based losely on DDREAD.PRO by Frank Knight (MIT)
;       11/01/2007 -  v1.0 - first working version
;       27/04/2007 -  v1.1 - added TEXT option
;       22/11/2007 -  v1.2 - added HEAD option
;       22/09/2008 -  v1.3 - HEAD now takes integer input
;       16/06/2010 -  v1.4 - fixed handling of text input
;                             using STRSPLIT function.
;       20/07/2012 -  v1.5 - If error encountered then RETURN
;                             a value -1, for consistency with
;                             other routines.
;       28/07/2012 - v1.6  - Replace obsolete FINDFILE with
;                             FILE_SEARCH. Upon error return
;                             !NULL. (New to IDL 8). Use 
;                             QUERY_ASCII to check the file is
;                             indeed ASCII, and to define the
;                             number of rows in the file.
;       20/08/2014 - v1.7 - Fixed minor bug occuring when 
;                             MAX(columns) = ncols-1. Thanks 
;                             to Tatjana Koukal for spotting this.  
;                              
; NOTES:
;       This is similar to the IDL READ_ASCII function. But
;       it allow you to force the data type. This is useful
;       if the input file comprises columns of different types,
;       e.g. 1 column of strings, 3 columns of floats. In this
;       case use READ_TABLE(..., /TEXT) which will force all
;       columns to be read as strings, and then convert the
;       numerical columns from string to float as needed.
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)

  COMPILE_OPT idl2, HIDDEN

; watch out for errors

  ON_ERROR, 0

; ----------------------------------------------------------
; Check the arguments

; is the file name defined?

  IF (N_ELEMENTS(filename) eq 0) THEN BEGIN
      filename = ''
      READ, '-- Enter file name (ENTER to list current directory): ', filename
  ENDIF

  IF (filename eq '') THEN BEGIN
      list = FILE_SEARCH(/TEST_READ, /TEST_REGULAR)
      PRINT, list
      READ,'-- Enter file name: ', filename
  ENDIF

; are we reading in single (=4) or double precision (=5)?
 
  type=4
  IF (KEYWORD_SET(double)) THEN type=5

; are we reading numbers or text?
 
  IF (KEYWORD_SET(text)) THEN type=7
 
; ----------------------------------------------------------
; Checks of the file existance and shape

; check the file exists

  file = (FILE_SEARCH(filename, /TEST_READ, /TEST_REGULAR))
  IF (file[0] eq '') THEN BEGIN
      PRINT,'** File not found.'
      RETURN, !NULL
  ENDIF
  IF (N_ELEMENTS(file) ne 1) THEN BEGIN
      PRINT,'** File not found.'
      RETURN, !NULL
  ENDIF

; check the file is ASCII and count the total number of lines

  check = QUERY_ASCII(file, file_info)
  IF (check eq 0) THEN BEGIN
      PRINT, '** File is not ASCII.'
      RETURN, !NULL
  ENDIF
  file_row = file_info.lines
  IF KEYWORD_SET(head) THEN file_row = file_row - head
  IF NOT KEYWORD_SET(nrows) THEN nrows = file_row
  nrows = (nrows < file_row)
  
; find the number of columns in the file by reading first line
; into a string (tmp)

  ncols = 0
  tmp = ''
  OPENR, lun, file, /GET_LUN
  if KEYWORD_SET(head) THEN BEGIN
      FOR i=0, head-1 DO READF, lun, tmp ; skip header
  ENDIF
  READF, lun, tmp
  FREE_LUN, lun

; remove whitespace

  tmp = ' ' + STRCOMPRESS(STRTRIM(tmp, 2))

; count the spaces (there is one per column)

  FOR i=0, STRLEN(tmp)-1 DO BEGIN
      ncols = ncols + (STRPOS(tmp, ' ', i) eq i)
  END

; ----------------------------------------------------------
; load the data into an array

; define the data array ready to receive data

  data = MAKE_ARRAY(size=[2, ncols, nrows, type, ncols*nrows])

; define a single line (row) array for reading each line
; except for text which is loaded a whole line at a time 

  IF NOT KEYWORD_SET(text) THEN BEGIN
      record = MAKE_ARRAY(size=[1, ncols, type, ncols])
  ENDIF ELSE BEGIN
      record = ''
  ENDELSE

; Open the file ready to read

  OPENR, lun, file, /GET_LUN

; skip header line if HEAD keyword is set

  if KEYWORD_SET(head) THEN BEGIN
      FOR i=0, head-1 DO READF, lun, tmp 
  ENDIF

; Read each line one at a time, until either end-of-file
; or we reach nrows.

  n = 0L
  WHILE (eof(lun) ne 1) DO BEGIN
      ON_IOERROR, IOERR
      error = 1
      READF, lun, record
      error = 0
      IF KEYWORD_SET(text) THEN BEGIN
          data[*, n] = STRSPLIT(record, ' ', /EXTRACT)
      ENDIF ELSE BEGIN
          data[*,n] = record
      ENDELSE
      n = n + 1
      IF (n eq nrows) THEN BREAK

      IOERR:
      IF (error eq 1) THEN BEGIN
          PRINT, '** Error reading line', n, ' in READ_TABLE'
          FREE_LUN, lun
          RETURN, !NULL
      ENDIF

  ENDWHILE

; Close the file

  FREE_LUN, lun

; ----------------------------------------------------------
; Return the data array to the user

; if no column selection, RETURN entire array 

  if (N_ELEMENTS(columns) eq 0) THEN RETURN, data

; otherwise remove unwanted columns before RETURNing

  indx = WHERE((columns ge ncols), count)
  IF (count eq 0) THEN BEGIN
      data = data[columns, *]
  ENDIF ELSE BEGIN
      PRINT, '** Requested columns outside allowed range'
      PRINT, '** Returning all columns from READ_TABLE'
  ENDELSE

  RETURN, data

END
; ----------------------------------------------------------
