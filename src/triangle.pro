FUNCTION TRIANGLE, x, r=r

; ----------------------------------------------------------
;+
; NAME:
;       TRIANGLE
;
; PURPOSE:
;       Generate a simple triangular function
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       y = triangle(x,r=10)
;
; INPUTS:
;       x         - the dependent variable (vector or scalar)
;
; OPTIONAL INPUTS:
;       r         - the half-width of the triangle
;
; DETAILS:
;       The function is 
;
;              f(x) = 1-|x/r|   where |x| < r
;              f(x) = 0         elsewhere
;       
;       The output is therefore always non-negative.
;       The triangle is centered on x=0.
;
;      ^
;  y=1 |           ^  
;      |          / \ 
;      |         /   \
;      |        /     \
;      |       /       \
;      |      /         \
;  y=0 +------+----+----+-----> x
;          x=-r   x=0   x=r
;
; HISTORY:
;       19/07/07  - v1.0 - first working version
;
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

  y = ( (1-abs(float(x)/r)) > 0 )

  return, y

END
