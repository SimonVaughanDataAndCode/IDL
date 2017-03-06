
; ----------------------------------------------------------
;+
;
; Define model: broken power law plus constant
;
;  y(x) = N*(x/x_b)^-a1 + C       where x <  x_b
;         N*(x/x_b)^-a2 + C       where x >= x_b
;
; where x_b is the break point.
;
; The parameters are stored in the one dimensional array P
; P = [N, a1, a1, fb, C]
;      N   = normalisation (input as log[N])
;      a1  = power law index at x < x_b
;      a2  = power law index at x > x_b
;      x_b = break point
;      C   = constant
;
; x_b is passed to the function in log_10 units.
;
; EXAMPLE:
;
;   f = indgen(1000)+1.0
;   p = [0.0, 1.0, 2.0, 1.7, 0.1]
;   plot, f, p_model(f, p), /xlog, /ylog
;
;-
; ----------------------------------------------------------

FUNCTION P_MODEL1, x, p

  Norm = 10.0^float(p[0])
  a1   = float(p[1])
  a2   = float(p[2])
  x_b  = 10.0^float(p[3])
  C    = float(p[4])

  if (C lt 0.0) then C = 0.0

; is input X a scalar or vector?

  dim = size(x,/DIMENSION)
 
; Normalisation above (N2) and below (N1) the break x_b

  N1 = Norm*x_b^a1
  N2 = Norm*x_b^a2
 
; if F is a scalar comput scalar output PER

  if (dim eq 1) then begin

     if (x lt x_b) then begin

         y = N1*x^(-a1) 
  
     endif else begin

         y = N2*x^(-a2) 

     endelse

 endif else begin

; create output PER same size as input F

     n = n_elements(x)
     y = make_array(n,value=0.0)

     low_x = where(x lt x_b, COUNT, COMPLEMENT=high_x)

     if (COUNT gt 0) then y[low_x]  = N1*x[low_x]^(-a1)
     if (COUNT lt n) then y[high_x] = N2*x[high_x]^(-a2)

 endelse

; add the constant, force to be non-negative

  y = abs(y) + C

; return to calling function

  return, y

END

; ----------------------------------------------------------
;+
;
;
; Define model: single power law plus constant 
;
;  y(x) = Nx^-a + C
;
; the parameters are stored in the one dimensional array P
;
; P = [N, a, C]
;      N   = normalisation (input in log[N])
;      a   = power law index 
;      C   = constant
;
; EXAMPLE:
;
;   f = indgen(1000)+1.0
;   p = [1.0,1.0,0.1]
;   plot, f, p_model1(f, p), /xlog, /ylog
;
;-
; ----------------------------------------------------------

FUNCTION P_MODEL0, x, p

  Norm = 10.0^float(p[0])
  a    = float(p[1])
  C    = float(p[2])

  if (C lt 0.0) then C = 0.0

; create output PER same size as input F

  y = Norm * x^(-a)

; add the constant

  y = y + C

; force to be non-negative

  y = abs(y)

; return to calling function

  return, y

END

; ----------------------------------------------------------
;+
;
; Define log likelihood function for periodogram
;
; The probability density function (PDF) for a standard
; periodogram ordinate is
;
;   p( I_j | M_j ) = (1/M_j) exp(-I_j / M_j)
;
; where I_j is the periodgram ordinate at frequency f_j
; and M_j is the "true" value. This can be used to formulate
; the likelihood function of the model M_j given the data I_j
;
;  L = prod_{j} [ (1/M_j) exp(-I_j / M_j) ]
;
; which is better expressed as S = -2 ln[L]
;
;  S = 2 sum_{j} ln[M_j] + I_j/M_j
;
; For more details see the Appendix of Vaughan (2005; A&A, 431, 391)
;
; INPUTS:
;      p   - array of parameters values for model
;
;
; VIA COMMON BLOCK:
;      x - array of frequencies
;      y - periodogram value (data) at each frequency
; 
;
; -
; ----------------------------------------------------------

FUNCTION MLOGL1, p

; access the common block containing the data points

  common per_data, x, y

; compute model periodgram at frequencies f

  model = p_model1(x, p)

; check for zeros

  mask = where(model le 0.0, count)
  if (count gt 0) then begin
     print,"** Warning: detected zero or -ve power in model"
     print,"** at ",count," frequencies: ",x[mask]
     model[mask] = abs(model[mask]) + 1.0e-12
 endif

; calculate l_j = ln[L_j] given model and data

  l = (alog(model) + y/model)

; calculate S = -2 sum_j ln[l_j] given model and data

  S = 2.0 * total(l,/DOUBLE)

 return, S

END

; ----------------------------------------------------------

FUNCTION MLOGL0, p

; compute model periodgram at frequencies f
; ++ Alternative version for simple power law model ++

; access the common block containing the data points

  common per_data, x, y

 model = p_model0(x, p)

; check for zeros

  mask = where(model le 0.0, count)
  if (count gt 0) then begin
     print,"** Warning: detected zero or -ve power in model"
     print,"** at ",count," frequencies: ",x[mask]
     model[mask] = abs(model[mask]) + 1.0e-12
  endif

; calculate l_j = ln[L_j] given model and data

  l = (alog(model) + y/model)

; calculate S = -2 sum_j ln[l_j] given model and data

  S = 2.0 * total(l,/DOUBLE)

 return, S

END
