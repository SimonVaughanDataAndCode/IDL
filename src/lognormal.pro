FUNCTION LOGNORMAL,x,mu,sigma,thresh

; ----------------------------------------------------------
;+
; NAME:
;       LOGNORMAL
;
; PURPOSE:
;       Compute density of lognormal function
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       y = LOGNORMAL(x,0.0,1.0,1.0)
;
; INPUTS:
;       x      - (float/array) position to evaulate density
;       mu     - (scalar) location parameter
;       sigma  - (scalar) shape parameter
;       thresh - (scalar) threshold parameter
;
; OUTPUTS:
;       y      - (float/array) densities evaluated at x
;
; DETAILS:
;       Compute the lognormal density function at position(s) x.
;       
;         y(x) = exp[ -(log[x-thresh] - mu)^2 / 2*sigma^2 ]       
;                -------------------------------------
;                      (x-thresh)*sigma*sqrt(2*pi)
;
;
; HISTORY:
;        09/05/2007 -- v1.0 -- first working version
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments

  if (n_elements(mu) eq 0) then mu = 0.0
  if (n_elements(sigma) eq 0) then sigma = 1.0
  if (n_elements(thresh) eq 0) then thresh = 0.0

; ----------------------------------------------------------
; evaluate density

  y = exp( -(alog(x-thresh) - mu)^2 / (2.0*sigma^2) ) / $ 
           ((x-thresh)*sigma*sqrt(2*!pi))

; ----------------------------------------------------------
; return to caller
  return,y

END
