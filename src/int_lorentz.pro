FUNCTION int_lorentz, x0, dx, N, minx=minx, maxx=maxx, Q=Q


; ----------------------------------------------------------
;+
; NAME:
;       INT_LORENTZ
;
; PURPOSE:
;       Calculate the definite integral of a Lorentzian function
;
; AUTHOR:
;       Simon Vaughan (U.Leicester) 
;
; CALLING SEQUENCE:
;       result = INT_LORENTZ(x0, dx, N)
;
; INPUTS:
;       x0       - Line centre
;       dx       - Line FWHM (Q = x0/dx)
;       N        - Line normalisation
;
; OPTIONAL INPUTS:  
;       minx     - Lower limit of integration
;       maxx     - Upper limit of integration
;
; OUTPUTS:
;       result   - Definite integral from x_lo to x_hi
;
; OPTIONAL OUTPUTS:  
;       Q        - Quality factor = x0/(dx)
;
; DETAILS:
;       Integrates the Lorentzian line profile
;
;                    N * dx/2pi
;         f(x) = -------------------
;                (x-x0)^2 + (dx/2)^2
; 
;       over the range MINX to MAXX. If these limits are not specified
;       we integrate over the positive real line (0, +infty). The
;       above can be re-written in terms of x0 and the quality factor
;       Q = x0/dx as
;
;                  2N * Q * x0 / pi
;         f(x) = -------------------
;                x0^2 + 4Q^2*(x-x0)^2
; 
;       Its analytical integral is:
;
;         F(x) = N/pi * arctan([x-x0] / [dx/2])
;
;       The values at x=0 and x=+infty are
;
;         F(0) = N/pi * arctan(-2Q)
;
;         F(+infty) = N/2
;
;       Over the full real range (-infty,+infty) this is equal to
;       N. (As we expect since it is equivalent to the Cauchy
;       distribution which is normalised such.) And as x0 or Q tend to
;       +infty we also get integral --> N.
;
; EXAMPLE USAGE:
;       result = INT_LORENTZ(1.0, 0.1, 1.0)
;
;       print, INT_LORENTZ(1.0, 0.1, 1.0, minx=0.1, maxx=10)
;
; HISTORY:
;       14/06/2010 - v1.0 - first working version
;       17/06/2010 - v1.1 - changed definition of Lorentzian to that
;                    used in Pottschmidt et al. (2003; AA, 407, 1039)
;                    which is more standard. This means Q=x0/dx. This
;                    makes no difference to the integral results. 
;        3/08/2010 - v1.2 - explained df parameter better in header
;
; NOTES:
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)
  COMPILE_OPT idl2  
  
; watch out for errors
  ON_ERROR,2

; ---------------------------------------------------------
; Check arguments

  N_par = N_PARAMS()
  if (N_par lt 3) then BEGIN
      PRINT,'** Usage: z=INT_LORENTZ(x0, dx, N)'
      RETURN, -1
  endif

; Set lower/upper bounds to 0, +infty if not set

  if not KEYWORD_SET(minx) then minx = 0.0
  if not KEYWORD_SET(maxx) then maxx = !VALUES.F_INFINITY

;  check that limits are in correct order: MINX < MAXX

  if (minx ge maxx) then begin
      PRINT,'** INT_LORENTZ: MIN/MAX values are faulty.'
      RETURN, -1
  endif

; ---------------------------------------------------------
; Main routine

; calcualate 'quality factor'

  Q = x0 / dx

; evaluate integral at lower/upper bounds

  y0 = ATAN((minx-x0)/(dx*0.5))
  y1 = ATAN((maxx-x0)/(dx*0.5))

; integral over lower-upper range is difference

  result = N / !pi * (y1 - y0)

; ---------------------------------------------------------
; RETURN to user

  RETURN, result

END

; ---------------------------------------------------------
;
;
;
;
;
;
; ---------------------------------------------------------

FUNCTION int_lorentz_model, x, parm, dx=dx 

; ----------------------------------------------------------
;+
; NAME:
;       INT_LORENTZ_MODEL
;
; PURPOSE:
;       Calculate the integrated Lorentzian function at several points
;
; AUTHOR:
;       Simon Vaughan (U.Leicester) 
;
; CALLING SEQUENCE:
;       result = INT_LORENTZ_MODEL(x, parm, dx=dx)
;
; INPUTS:
;       x        - bin centre, at which to compute model
;       dx       - bin width, over which to integrate model
;       parm     - (vector) parameters of Lorentzian (centre, width, norm)
;
; OPTIONAL INPUTS:  
;       none
;
; OUTPUTS:
;       result   - Definite integral over [x-dx/2, x+dx/2] for each x
;
; OPTIONAL OUTPUTS:  
;       none
;
; DETAILS:
;       The parm input vector specifies the three Lorentzian function
;       parameters: [centre, width, norm]. Where centre in the
;       location, width gives the FWHM (in the same units) and norm
;       gives the overal normalisation (i.e. integral from -infty to
;       +infty). 
;;
;       Integrates the Lorentzian line profile
;
;                    N * FWHM/2pi
;         f(x) = ---------------------
;                (x-x0)^2 + (FWHM/2)^2
; 
;       over each bin [x-dx/2, x+dx/2]. 
;
;       Its analytical integral is:
;
;         F(x) = N/pi * arctan([x-x0] / [FWHM/2])
;
;
; EXAMPLE USAGE:
;       x = INDGEN(30)+800
;       parm = [815.4, 3, 1]
;       dx = 0.5
;       result = INT_LORENTZ_MODEL(x, parm, dx=dx)
;       PLOT, x, result, PSYM=10
;       x = INDGEN(300)/10.0+800
;       lorz = parm[2]*parm[1]/(2*!pi)/((0.5*parm[1])^2 + (x-parm[0])^2)
;       OPLOT, x, lorz
;
;
; HISTORY:
;       15/06/2011 - v1.0 - first working version
;
; NOTES:
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)
  COMPILE_OPT idl2  
  
; watch out for errors
  ON_ERROR,2

; ---------------------------------------------------------
; Check arguments

  Nx = N_ELEMENTS(x)
  if (Nx lt 1) then begin
      PRINT,'** INT_LORENTZ_MODEL: not enough data points in X'
      RETURN, -1
  endif

  Np = N_ELEMENTS(parm)
  if (Np ne 3) then begin
      PRINT,'** INT_LORENTZ_MODEL: not enough parameters in PARM'
      RETURN, -1
  endif

  if (N_ELEMENTS(dx) eq 0) then dx = x[1]-x[0]

; ---------------------------------------------------------
; Main routine

  x0 = parm[0]
  fwhm = parm[1]
  N = parm[2]

  y = MAKE_ARRAY(Nx)

  for i = 0, Nx-1 do begin
      minx = x[i] - dx*0.5
      maxx = x[i] + dx*0.5
      y[i] = INT_LORENTZ(x0, fwhm, N, MINX=minx, MAXX=maxx)
  endfor

  y = y/dx

; ---------------------------------------------------------
; RETURN to user

  RETURN, y

END
