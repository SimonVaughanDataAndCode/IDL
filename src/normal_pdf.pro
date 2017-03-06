FUNCTION normal_pdf, x, MEAN=mean, SD=sd

; ----------------------------------------------------------
;+
; NAME:
;       NORMAL_PDF
;
; PURPOSE:
;       Compute density for chi-sq distribution
;
; AUTHOR:
;       Simon Vaughan
;
; CALLING SEQUENCE:
;       d = NORMAL_PDF(x, mean=0, sd=1.0)
;
; INPUTS:
;       x       - (float vector) N values at which to compute density
;       mean    - (float) value or N-vector of means
;       sd      - (float) value or N-vector of std devs.
;
; OUTPUTS:
;       d       - (float vector) densities pdf(x)
;       
; OPTIONAL OUTPUTS:
;       none
;       
; DETAILS:
;       Compute the probability density function (PDF) for the 
;       normal distribution.
;       
;       Example usage
;       IDL> x = INDGEN(100, /FLOAT) * 0.1
;       IDL> d = NORMAL_PDF(x, mean=4.3, sd=1.3)
;       IDL> PLOT, x, d
;
; HISTORY:
;       18/07/13  - v1.0 - first working version
;
; USES:
;       none
;
; NOTES: 
;       none
;
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)

  COMPILE_OPT idl2

; watch out for errors

  ON_ERROR, 3

; ----------------------------------------------------------
; check arguments

  n = N_ELEMENTS(x)

; is the sampling period dt supplied? (default = 1.0)

  if (n eq 0) then begin
    PRINT, '-- Missing X values in NORMAL_PDF'
    RETURN, !NULL
  ENDIF

  if NOT KEYWORD_SET(mean) THEN mean=0.0D
  if NOT KEYWORD_SET(sd) THEN mean=1.0D
  
; ----------------------------------------------------------
; main routine

  z = (x-mean) / sd 
  d = EXP(-0.5*z*z) / SQRT(2.0D*!pi) / sd

; ----------------------------------------------------------
; return values and end

  RETURN, d

END
  