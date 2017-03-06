PRO ARROW_PLOT, GRAD, X, Y, $
                SPREAD_x=spread_x, SCALE=scale, NOPLOT=NOPLOT, $
                _EXTRA=EXTRA_KEYWORDS, FUDGE=fudge

; ----------------------------------------------------------
;+
; NAME:
;       ARROW_PLOT
;
; PURPOSE:
;       Plot scatter diagram using arrows of given gradients
;
; AUTHOR:
;       Simon Vaughan (U.Leicester) 
;
; CALLING SEQUENCE:
;       ARROW_PLOT, grad, mean_x, mean_y
;
; INPUTS:
;       grad          - (vector) gradients (dy/dx) of each arrow
;       x             - (vector) x-position of arrow centre
;       y             - (vector) y-position of arrow centre
;
; OPTIONAL INPUTS:  
;       spread_x      - (vector) x-length of arrows
;       scale         - (float) re-scaling factor for arrows
;       noplot        - (logical) do not make a new plot?
;       _extra        - (any) any other keywords to pass to ARROW
;       fudge         - (flat) extend axis ranges by this factor
;
; OUTPUTS:
;       none
;
; OPTIONAL OUTPUTS:  
;       none
;
; DETAILS:
;       Produce a scatter plot using arrows instead of points. 
;       The directions of the are given by the GRAD input and 
;       the (centre) positions are given by X and Y.
;       The length of the arrow is scaled to SPREAD_X. An optional
;       input SCALE specifies a re-scaling factor for all lengths.
;
; EXAMPLE USAGE:
;
;       x = INDGEN(10)+1
;       y = 0.7*x + RANDOMN(seed, N_ELEMENTS(x))
;       grad = RANDOMU(seed, N_ELEMENTS(x))-0.5 
;       spread_x = sqrt(x+1)
;       
;       ARROW_PLOT, grad, x, y, spread_x=spread_x, scale = 0.5 
;
;         ...or using the NOPLOT option...
;
;       PLOT, x, y, xrange=[0, 12], yrange=[-3, 10], $
;             /XSTYLE, /YSTYLE, PSYM=7, XTITLE="x", YTITLE="y"
;       ARROW_PLOT, grad, x, y, spread_x=spread_x, scale = 0.5, $
;                   /SOLID, /NOPLOT
;
;          ...and arrow heads that scale with arrow length...
;
;       ARROW_PLOT, grad, x, y, spread_x=spread_x, scale = 0.5, $
;                   /SOLID, /NOPLOT, HSIZE=-0.1
;
;
; HISTORY:
;        23/06/2010 - v1.0 - first working version
;
; NOTES:
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)

  COMPILE_OPT idl2  
  
; watch out for errors

  on_error, 2

; ---------------------------------------------------------
; Check arguments

; is input array large enough to bother with?

  N = N_ELEMENTS(x)
  if (N le 2) then begin
      PRINT,'** Too few data in ARROW_PLOT'
      RETURN
  endif

; check the other arrays match in size

  if (N_ELEMENTS(y) ne N) then begin
      PRINT, '** Arrays X and Y do not match in ARROW_PLOT'
      RETURN
  endif

  if (N_ELEMENTS(grad) ne N) then begin
      PRINT, '** Arrays X and GRAD do not match in ARROW_PLOT'
      RETURN
  endif

  if KEYWORD_SET(spread_x) then begin
    if (N_ELEMENTS(spread_x) ne N) then begin
        PRINT, '** Arrays X and SPREAD_X do not match in ARROW_PLOT'
        RETURN
    endif
  endif else begin
     spread_x = MAKE_ARRAY(N, VALUE=1.0)
  endelse

  if NOT KEYWORD_SET(scale) then scale = 1.0

  if NOT KEYWORD_SET(fudge) then fudge = 0.1

; ---------------------------------------------------------
; Main routine

; compute the scaled DX and DY

  dx = spread_x * scale    
  dy = grad * dx

; make the plot

  x0 = x - dx/2.0
  x1 = x + dx/2.0
  y0 = y - dy/2.0
  y1 = y + dy/2.0

  if NOT KEYWORD_SET(noplot) then begin

      xspan = MAX(x) - MIN(x) 
      xf = xspan * fudge
      xrange = [MIN(x)-xf, MAX(x)+xf]
      
      yspan = MAX(y) - MIN(y)
      yf = yspan * fudge
      yrange = [MIN(y)-yf, MAX(y)+yf]
  
      PLOT, x, y, XRANGE=xrange, YRANGE=yrange, /NODATA

  endif

  ARROW, x0, y0, x1, y1, /DATA, _EXTRA=extra_keywords

; ---------------------------------------------------------

END
