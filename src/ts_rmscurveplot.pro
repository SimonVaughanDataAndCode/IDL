PRO TS_RMSCURVE_PLOT, filename

; ----------------------------------------------------------
;+
; NAME:
;       TS_RMSCURVE_PLOT
;
; PURPOSE:
;       Plot the output from TS_RMSCURVE
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       rms = TS_RMSCURVE_PLOT('filename.txt')
;
; INPUTS:
;       filename  - (string) name of file listing FITS time series
;
; OUTPUTS:
;       <files>   - series of ASCII files 
;
; DETAILS:
;       The filename given on input should point to an
;       ASCII file listing the names of a series of files
;       produced by TS_RMSCURVE. These will be read into
;       memory, concatonated into a single dataset and the
;       two time series (one for mean flux, one for rms)
;       will be plotted.
;
; Example call:
; IDL> rms = TS_RMSCURVE('filelist.txt',n_seg=64,/pois)
;
; PROCEDURES CALLED:
;       TS_READFITS, TS_SEGMENT, TS_RMSFLUX, WRITE_TABLE
;
; HISTORY:
;       27/04/07  - v1.0 - first working version
;
; NOTES:
;       + Add deadtime correction (to TS_RMSFLUX)
;
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments

; check filename has been set
  if (n_elements(filename) eq 0) then begin
      print,'** FILENAME not set in TS_RMSCURVE_PLOT'
      return
  endif

; ----------------------------------------------------------
; Main part of procedure

; Load list of FITS file names from the file
  filelist = read_table(filename,/text)

; Determine how many FITS files are to be processed
  N = n_elements(filelist)
  print,'-- Files to process',N

; Loop over each file to process
  for i=0,N-1 do begin

      print,'-- Loading file',i+1

; Read the file
      x = read_table(filelist[0,i],/double)

; add to the data arrays

      if (i eq 0) then begin
          time = x[0,*]
          flux = x[1,*]
          rms  = x[2,*]
      end else begin
          time = [time,x[0,*]]
          flux = [flux,x[1,*]]
          rms  = [rms,x[2,*]]
      endelse

; End loop over files

  endfor

; subtract the start time

;  time = time - time[0]

; ----------------------------------------------------------
; Return control to the user

  return

END
