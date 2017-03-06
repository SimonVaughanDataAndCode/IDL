function seq,x0,xn,dx,n=n,endpoint=endpoint,double=double

; ----------------------------------------------------------
;+
; NAME:
;       SEQ
;
; PURPOSE:
;       Generate sequence of evenly spaced real numbers
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       x = seq(0,1)
;
; INPUTS:
;       x0     - (float) start value of sequence
;       xn     - (float) end value of sequence
;
; OPTIONAL INPUTS:
;       dx     - (float) step size for sequence
;       n      - (length) of sequence (default = 101)
;       endpoint - (logical) flag to exclude endpoint
;       double - (logical) double prec output 
;
; OUTPUTS:
;       x     - (vector) sequence of numbers over (x0,xn)
;
; DETAILS:  
;       Generate a vector containing a sequence of 
;       floating point numbers evenly spaced along a user
;       defined interval (X_0,X_N). The start and end 
;       points of the interval must be specified and then
;       either the length of the sequence N or the step
;       size DX must be specified. If both are specified
;       then DX takes priority over N.
;
;       If DX is specified then the length of the sequence
;       is given by 
;                   N = (X_N - N_0)/DX + 1
;       The extra one is to make sure the sequence is inclusive
;       of the end points of the interval (X_0,X_N).
;       This can be switch off useing the ENDPOINT flag, in
;       which case the interval is non-inclusive (X_0,X_1].
;       ENDPOINT is only used if DX but not N is not given.
;
; EXAMPLE USAGE:
;       y = seq(0.0,64,0,dx=0.1)
;
; HISTORY:
;       19/10/2007  -- v1.0 -- first working version
;
; NOTES:
;       + Treat double prec calculations better
;
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments and input

; check interval (X0,XN) is given
  if (n_elements(x0) eq 0) OR (n_elements(xn) eq 0) then begin
      print,'** SEQ: START and/or STOP values are not set.'
      return,0
  endif

; check interval makes sense
  if (x0 ge xn) then begin
      print,'** SEQ: STOP value > START value.'
      return,0
  endif

; check if N is defined
  if (n_elements(n) gt 0) then old_n = n

; determine N from DX value if given
  if (n_elements(dx) eq 1) then begin
      n = (xn - x0)/dx + 1
      if keyword_set(endpoint) then n = n -1
  endif

; make sure N exists (either from input N or DX)
; if not, use default
  if (n_elements(n) eq 0) then n=101
  
; make sure N is a worthwhile size
  if (n lt 2) then begin
      print,'** SEQ: Sequence is too short.'
      return,0
  endif

; make sure N is (long) integer form
  n = long(n)

; ----------------------------------------------------------
; Main section

; determine correct DX
  if (n_elements(dx) eq 0) then begin
      dx = (xn - x0)/float( n - 1 )
  endif

; define sequence of length N
  if keyword_set(double) then begin
      y = dindgen(n)
  endif else begin
      y = findgen(n)
  endelse

; scale to to cover correct interval (X0,X1)
  y = (y * dx) + x0

; ----------------------------------------------------------
; replace value of N if not used, before returning to user
  if (n_elements(old_n) gt 0) then n = old_n

; Return the result to the user
  return, y

END
