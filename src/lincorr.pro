FUNCTION LINCORR, X, Y, DOUBLE=double, T_STAT=t_stat

; ----------------------------------------------------------
;+
; NAME:
;       LINCORR
;
; PURPOSE:
;       Linear correlation coefficient with p-value
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       result = LINCORR(x, y)
;
; INPUTS:
;       x         - (vector) list of x values
;       y         - (vector) list of y values
;
; OPTIONAL INPUTS:
;       double    - (logical) use double prec for r calculation?
;
; OUTPUTS:
;       result    - list containing r and p
;
; OPTIONAL OUTPUTS:
;       t_stat    - (float) the t statistic
;
; DETAILS:
;       Computes the linear (Pearson) correlation coefficient for Y
;       upon X using the CORRELATE command. Then calculates the t
;       statistic and applies a two-sided t-test to the result and
;       returns the p-value. Small p-value --> "significant" result.
;
; EXAMPLE USAGE:
;       x = RANDOMN(seed, 20)
;       y = -0.75 * x + RANDOMN(seed, 20) 
;       result = LINCORR(x, y, T_STAT=t1)
;       PRINT, result
;       PRINT, t1
;
;        Try the regression route and show that t=slope/error
;        is the same t statistic (and therefore gives the same
;        p-value) 
;
;      slope = REGRESS(x, y, sigma=sigma)
;      t2 = slope/sigma
;      PRINT, t2
;
; HISTORY:
;       15/12/2010 - v1.0 - first working version
;       24/02/2011 - v1.1 - added a catch to cope with 
;                           div by zero error when r=1.
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

  n = N_ELEMENTS(x)
  if (n le 3) then begin
      print,'** Insufficient data in LINCORR'
      return, -1
  endif

  ny = N_ELEMENTS(y)
  if (ny ne n) then begin
      print,'** X and Y arrays are not the same size in LINCORR'
      return, -1
  endif

; ---------------------------------------------------------
; Main routine

; compute the correlation coefficient

  r = CORRELATE(x, y, DOUBLE=KEYWORD_SET(double))

; compute the degrees of freedom

  df = n - 2

  if (1.0 - r le 1.0E-7) then begin

       p = 0.0

  endif else begin

; compute the t statistic

        t = r/sqrt((1.0-r*r)/FLOAT(df))

; calculate the two-side 'tail area' probability 

        p = 2.0 * (1.0 - T_PDF(ABS(t), df))

  endelse

; ---------------------------------------------------------
; Return to user

  t_stat = t
  return, [r, p]

END
