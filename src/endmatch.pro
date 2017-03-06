FUNCTION endmatch, y, x=x

; ----------------------------------------------------------
;+
; NAME:
;       EndMatch
;
; PURPOSE:
;       Perform 'end matching' of a dataset
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       result = ENDMATCH(y, x=x)
;
; INPUTS:
;       y      - (vector) the data to be end-matched
;
; OPTIONAL INPUTS:
;       x      - (vector) positions of the data y
;
; OUTPUTS:
;       result - (vector) y positions after end matching
;
; OPTIONAL OUTPUTS:
;       none
;
; DETAILS:
;       Perform a simple 'end matching' - i.e. subtract
;       a linear function such that the first and last data points
;       have equal values: y[0] = y[N-1].
;
;       This can be useful in spectral analysis of time
;       series. Fougere (1985) showed how end matching can reduced red
;       noise 'leakage' and help recover intrinsically steep power
;       spectra. See
;
;         Fougere, P. F., 1985 JOURNAL OF GEOPHYSICAL RESEARCH,
;          VOL. 90, NO. A5, PP. 4355-4366, 1985 
;
; EXAMPLE USAGE:
;       
;        N = 100
;        x = INDGEN(N)
;        r = RANDOMN(seed, N)
;        y = TOTAL(r, /CUMULATIVE)
;        PLOT, x, y
;        yout = ENDMATCH(y, x=x)
;        OPLOT, x, yout, LINESTYLE=1
;
; HISTORY:
;        13/02/2012 - v1.0 - first working version
;
; NOTES:
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)
  COMPILE_OPT idl2

; watch out for errors
  on_error,2

; ---------------------------------------------------------
; Check arguments

; check the shape of the input array

  N = N_ELEMENTS(y)
  if (N le 2) then begin
    PRINT, '** Array Y is too small in ENDMATCH.'
    RETURN, -1
  endif

  s=size(y)
  if (s[0] ne 1) then begin
     PRINT, '** Array Y has ',s[0],' dimensions in ENDMATCH'
     RETURN, -1
  endif

; check to see if X positions are supplied

  Nx = N_ELEMENTS(x)
  if (Nx gt 0) then begin
     if (Nx ne N) then begin
         PRINT, '** Arrays X and Y are different sizes in ENDMATCH.'
         RETURN, -1
     endif 
  endif else begin
     x = INDGEN(N)
  endelse

; ---------------------------------------------------------
; Main routine

; compute the overall mean of Y

  ymean = MEAN(y)

; compute the linear function Y = mX + C joining y[0] and y[N-1]

  m = (y[N-1] - y[0]) / (x[N-1] - x[0])
  c = y[0] - m*x[0]
  ymod = m*x + c

; subtract the linear function

  yout = y - ymod

; make sure the mean value is returned to original value

  yout = yout - MEAN(yout) + ymean

; ---------------------------------------------------------
; Return to user

  RETURN, yout

END

