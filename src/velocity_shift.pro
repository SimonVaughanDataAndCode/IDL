FUNCTION VELOCITY_SHIFT, w, w0, VELOCITY=velocity

; ----------------------------------------------------------
;+
; NAME:
;      VELOCITY_SHIFT
;
; PURPOSE:
;      Convert from observed and rest frame wavelength to velocity
;
; AUTHOR:
;      Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;      result = VELOCITY_SHIFT(21.11, 21.601)
;
; INPUTS:
;      W        - (float/vector) observed wavelength (velocities)
;      W0       - (float) rest-frame wavelength 
;      VELOCITY - (logical) reverse transformation?
;
; OPTIONAL INPUTS:
;      none
;
; OUTPUTS:
;      result   - output vector of velocities
;
; OPTIONAL OUTPUTS:
;      none
;
; DETAILS:
;      Computes the line-of-sight Doppler shift. 
;      Given an observed wavelength (w)
;      and a rest-frame wavelength (w0) the fractional change in
;      velocity due to the Doppler shift is: 
;      
;        w/w0 = sqrt((1+beta)/(1-beta))
;        
;      where beta = v/c (c = speed of light)
;      In input w can be a vector of values and the velocity
;      is computed for each one with respect to the rest-frame
;      wavelength w0. Negative velocity means blue-shift - 
;      motion towards the observed.
;      
;      If the VELOCITY keyword is set then the reverse transformation
;      from velocity to wavelength is applied, again with respect
;      to the rest-frame wavelength w0. 
;
; EXAMPLE USAGE:
;      result = VELOCITY_SHIFT(21.11, 21.601)
;      result = VELOCITY_SHIFT(-6891.8, 21.601, /VELOCITY)
;
; PROCEDURES CALLED:
;      none
;
; HISTORY:
;      22/12/2009 - v1.0 - first working version
;      12/08/2014 - v1.1 - changed loop to begin at i=0
;
; NOTES:
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)

  COMPILE_OPT idl2

; watch out for errors

  ON_ERROR, 0

  ; speed of light (km/s)

  c = 299792.458D
  
  delta = w / w0
  delta2 = delta * delta
  v = (delta2 - 1.0D) / (delta2 + 1.0D) * c
  
  IF KEYWORD_SET(velocity) THEN BEGIN
    beta = w/c
    v = w0 * SQRT((1+beta)/(1-beta))
  ENDIF
  
  RETURN, v

END