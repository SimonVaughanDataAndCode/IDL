function int_powerlaw,indx,n0,minx,maxx

; ----------------------------------------------------------
;+
; NAME:
;       INT_POWERLAW
;
; PURPOSE:
;       Find definite integral of power law function
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       z = INT_POWERLAW(-2.5,100.0,25,100)
;
; INPUTS:
;       indx  - (float) index of power law
;       n0    - (float) normalisation
;       minx  - (float) lower limit of integration
;
; OPTIONAL INPUTS:
;       maxx  - (float) upper limit of integration
;
; OUTPUTS:
;       z     - (float) definite integral
;
; DETAILS:  
;       Calculates integral of power law function
;   
;              f(x) = N0 * x^INDX
;
;       between limits MINX and MAXX.
;
;       If MAXX is not given at input then it is assumed
;       to be infinite, the integral is carried out over
;       [MINX,infinity] where possible. If the integral
;       diverges (at x=0 when indx <= -1, or at x=infinity
;       when indx >= -1) then the result is returned as
;       'Inf' (without warning).
;
; HISTORY:
;       30/07/2007  -- v1.0 -- first working version
;
; NOTES:
;       Makes use of the special floating point values
;       !VALUES.F_NAN (when f(x) not real-valued) and 
;       !VALUES.F_INFINITY (when integral diverges)
;       See on-line help regarding the use of these values.
;
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments and input

; check minimum number of arguments given
  n = n_params()
  if (n lt 3) then begin
      print,'** Usage: z=INT_POWERLAW(indx,norm,min,max)'
      return,0
  endif

; if MAXX is not set then default to infinity
  if (n eq 3) then maxx = !VALUES.F_INFINITY

; check that limits are in correct order: MINX < MAXX
  if (minx ge maxx) then begin
      print,'** INT_POWERLAW: MIN/MAX values are faulty.'
      return,0
  endif

; ----------------------------------------------------------
; Main section

; check that function f(x) is real-valued
  if (minx lt 0.0) then begin
      if ((indx mod 1) ne 0) then begin
          print,'** INT_POWERLAW: f(x) not real-valued in [MIN,MAX]'
          return,!VALUES.F_NAN
      endif
  endif

; check for divergence at x = 0
  if ((indx lt -1) && (minx le 0.0) && (maxx ge 0.0)) then begin
      return,!VALUES.F_INFINITY
  endif
  if ((indx eq -1) && (minx le 0.0)) then begin
      return,!VALUES.F_INFINITY
  endif
  
; check for divergence at x = infinity
  if ((indx ge -1) && (maxx eq !VALUES.F_INFINITY)) then begin
      return,!VALUES.F_INFINITY
  endif

; calculate the integral
  if (indx ne -1) then begin
      indx1 = indx + 1.0
      z = n0 / indx1 * (maxx^indx1 - minx^indx1)
  endif else begin
      z = n0 * (alog(maxx) - alog(minx))
  endelse

; ----------------------------------------------------------
; Return the result to the user
  return, z

END
