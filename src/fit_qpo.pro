

; ----------------------------------------------------------
;+
; NAME:
;       P_MODEL0
;
; PURPOSE:
;       Define model: power law plus constant
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; USED BY:
;       QPO_FIT
;
; DETAILS:
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
; N and C are passed to the function in log_10 units.
;
; To help alleviate "underflow" errors the computions
; are done in log-space, and converted to linear space
; at the end.
;
; EXAMPLE:
;
;   f = indgen(1000)+1.0
;   p = [1.0,1.0,0.1]
;   plot, f, p_model0(f, p), /xlog, /ylog
;
;-
; ----------------------------------------------------------

FUNCTION P_MODEL0, x, p

  Norm = p[0]
  a    = p[1]
  C    = p[2]

; Minimum allowed value - 10^logmin - used to prevent underflow

  logmin = -30.0

; create output array of cirrect size

  n = N_ELEMENTS(x)
  y = MAKE_ARRAY(n, /DOUBLE)

; convert to log_10 units

  logx = ALOG10(x)

; create output PER same size as input F

  logy = Norm - a*logx

; watch out for very low/high values (might cause underflow/overflow)

  lo_y = WHERE(logy lt logmin, n_l)
  hi_y = WHERE(logy gt -logmin, n_h)
  me_y = WHERE(logy gt logmin and logy lt -logmin, n_m)

  IF (n_m gt 0) THEN y[me_y] = DOUBLE(10.0^logy[me_y])
  IF (n_h gt 0) THEN y[hi_y] = DOUBLE(10.0^(-logmin))

; add the constant

  IF (C gt logmin) THEN y = y + DOUBLE(10.0^C)

; return to calling function

  return, y

END

;
;
;
;
;
;
;
;
;
;
;
;



; ----------------------------------------------------------
;+
; NAME:
;       P_MODEL1
;
; PURPOSE:
;       Define model: bending power law plus constant
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; USED BY:
;       QPO_FIT
;
; DETAILS:
;
; Define model: bending power law plus constant
;
; y(x) = N*x^(-a1) / [1 + (x/x_b)^(a2-a1) ]
;
; This is a power law of slope -a1 at x << x_b
; and a power law of slope -a2 at x >> x_b.
; Between these two extremes is smoothly changes
; from one slope to the other.
; (x_b is the break point.)
;
; The parameters are stored in the one dimensional array P
; P = [N, a1, a1, fb, C]
;      N   = normalisation (input as log[N])
;      a1  = power law index at x << x_b
;      a2  = power law index at x >> x_b
;      x_b = break point
;      C   = constant
;
; x_b, N and C are passed to the function in log_10 units.
;
; To help alleviate "underflow" errors the computions
; are done in log-space, and converted to linear space
; at the end.
;
; EXAMPLE:
;
;   f = indgen(1000)+1.0
;   p = [0.0, 1.0, 2.0, 1.7, -1.0]
;   plot, f, p_model1(f, p), /xlog, /ylog
;
;-
; ----------------------------------------------------------

FUNCTION P_MODEL1, x, p

  Norm = p[0]
  a1   = 1.0
  a2   = p[1]
  x_b  = p[2]
  C    = p[3]

  IF (a1 lt 0.0) THEN a1 = 0.0
  IF (a1 gt 10.0) THEN a1 = 10.0
  IF (a2 lt 0.0) THEN a2 = 0.0
  IF (a2 gt 10.0) THEN a2 = 10.0

; Minimum allowed value - 10^logmin - used to prevent underflow

  logmin = -30.0

; is input X a scalar or vector?

  n = N_ELEMENTS(x)
 
; convert to log_10 units

  logx = ALOG10(x)

; create output array Y same size as input array X

  y = MAKE_ARRAY(n, /DOUBLE)
  z = MAKE_ARRAY(n, /DOUBLE)
  logq = MAKE_ARRAY(n, /DOUBLE)
  
; Calculate bending factor q=1+z with z = (x/x_b)^(a2-a1).
; But to avoid overflow/underflow we calculate
;   log[z] = (a2 - a1) * (log[x] - log[x_b])

  logz = (a2-a1)*(logx - x_b)

; split into regions of low, medium and high z.
; i.e. where x << x_b, x ~ x_b and x >> x_b
; and treat these seperately to avoid underflow/
; overflow errors

  lo_z = WHERE(logz lt -7.0, n_l)
  hi_z = WHERE(logz gt +7.0, n_h)
  me_z = WHERE(logz ge -7.0 and logz le 7.0, n_m)

  IF (n_l gt 0) THEN logq[lo_z] = ALOG10(1.0)
  IF (n_h gt 0) THEN logq[hi_z] = logz[hi_z]
  IF (n_m gt 0) THEN BEGIN
      logq[me_z] = ALOG10( 1.0 + DOUBLE(10.0^(logz[me_z])) )
  ENDIF

; calculate log[y]

  logy = Norm - a1*logx - logq

; watch out for very low/high values (might cause underflow/overflow)

  lo_y = WHERE(logy lt logmin, n_l)
  hi_y = WHERE(logy gt -logmin, n_h)
  me_y = WHERE(logy gt logmin and logy lt -logmin, n_m)

  IF (n_m gt 0) THEN y[me_y] = DOUBLE(10.0^logy[me_y])
  IF (n_h gt 0) THEN y[hi_y] = DOUBLE(10.0^(-logmin))

; add the constant

  IF (C gt logmin) THEN y = y + DOUBLE(10.0^C)

; return to calling function

  RETURN, y

END

;
;
;
;
;
;
;
;
;
;
;
;
; ----------------------------------------------------------
;+
; NAME:
;       P_MODEL
;
; PURPOSE:
;       Define the model
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; USED BY:
;       QPO_FIT
;
; DETAILS:
;
;-
; ----------------------------------------------------------

FUNCTION P_MODEL, x, p, model

  IF (model eq 0) then y = P_MODEL0(x, p)
  IF (model eq 1) then y = P_MODEL1(x, p)

  RETURN, y

END

;
;
;
;
;
;
;
;
;
;
;
;
; ----------------------------------------------------------
;+
; NAME:
;       MLOGL
;
; PURPOSE:
;       Define log likelihood function for periodogram 
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; USED BY:
;       QPO_FIT
;
; DETAILS:
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
;      q - flag for dataset: 1, 2, 3, ..., N_files
;  model - which model to fit (choice of 0, 1, 2, ...) [unused]
;
; -
; ----------------------------------------------------------

FUNCTION MLOGL, p

; access the common block containing the data points

  COMMON per_data, x, y, q, model

; New lines to cope with multiple spectra         xxx

  N_parm = N_ELEMENTS(p)
  N_spec = MAX(q)

; compute a new parameter vector with only one C parameter

  N_nonC = N_parm - N_spec
  C_i = p[N_nonC:(N_nonC+N_spec-1)]
  parm = MAKE_ARRAY(N_nonC+1)
  parm[0:N_nonC-1] = p[0:N_nonC-1]
  S_i = MAKE_ARRAY(N_spec)

  for i = 0, N_spec-1 do begin

; compute model periodgram at frequencies f

      mask = WHERE(q eq i+1)
      xx = x[mask]
      yy = y[mask]

      parm[N_nonC] = C_i[i]
      model_y = p_model(xx, parm, model) 

; check for zeros

      mask = WHERE(model_y le 0.0, count)
      IF (count gt 0) THEN BEGIN
          PRINT,"** Warning: detected zero or -ve power in model"
          PRINT,"** at ",count," frequencies: ",x[mask]
          model[mask] = ABS(model_y[mask]) + 1.0e-12
      ENDIF

; calculate l_j = ln[L_j] given model and data

      l = DOUBLE(ALOG(model_y) + yy/model_y)

; calculate S = -2 sum_j ln[l_j] given model and data

      S_i[i] = 2.d * TOTAL(l, /DOUBLE)

  endfor

  S = TOTAL(S_i, /DOUBLE)

;  print,'--',S,S_i,'--',p

 RETURN, S

END

;
;
;
;
;
;
;
;
;
;
;
;
;
; ----------------------------------------------------------
;+
; NAME:
;       QPO_FIT_MODEL
;
; PURPOSE:
;       Perform the model fitting for QPO_FIT
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; INPUTS:
;       f      - (vector) independent variable (frequency)
;       per    - (vector) response variable (periodogram ordinate)
;       p0     - (vector) starting values of model parameter(s)
;
; OPTIONAL INPUTS:
;       scale0 - (vector) likely range for model parameter(s)
;       graph  - (integer) display optional graphics output
;       pvar   - (vector) stand. dev. for randomising starting values
;       N_fit  - (integer) number of trial starting values to use
;       model  - (integer) code for which model is being fitted
;       nodisp - (logical) switch off on-screen updates
;
; OUTPUTS:
;       nu     - (integer) degrees of freedom in fit
;       Smin   - (real) best-fit S = -2*ln[likelihood]
;       res    - (vector) data/model residuals from best fit
;       prob   - (real) p-value from KS-test
;
; USED BY:
;       QPO_FIT
;
; DETAILS:
;
;
;-
; ----------------------------------------------------------

FUNCTION QPO_FIT_MODEL, p0, nu=nu, Smin=Smin, res=res, graph=graph, $
                        pvar=pvar, N_fit=N_fit, prob=prob, scale0=scale0, $
                        nodisp=nodisp, norenorm=norenorm

; ----------------------------------------------------------
; access the common block containing the data points

  COMMON per_data, x, y, q, model

; ----------------------------------------------------------
;
; options for compilation (recommended by RSI)

  COMPILE_OPT idl2

; watch out for errors

  ON_ERROR, 2

; ----------------------------------------------------------
; Check the arguments

; number of frequencies

  N_spec = MAX(q)                           ; xxx

  n_f = N_ELEMENTS(x)

  IF (n_f lt 5) THEN BEGIN
      PRINT,"** Not enough data points to fit."
      RETURN,0
  ENDIF

  IF (N_ELEMENTS(y) ne n_f) THEN BEGIN
      PRINT, "** Array size mismatch in QPO_FIT_MODEL."
      RETURN, 0
  ENDIF

; check input parameter values

  IF (NOT KEYWORD_SET(p0)) THEN BEGIN
      PRINT, "** Parameters P0 not set in QPO_FIT_MODEL."
      RETURN, 0
  ENDIF

  N_par = N_ELEMENTS(p0)

  IF (N_ELEMENTS(pvar) eq 0) THEN pvar=make_array(N_par, value=1.0)

  IF (N_ELEMENTS(pvar) ne N_par) THEN BEGIN
      PRINT, "** PVAR size does not match P0 in QPO_FIT_MODEL:", $
        N_par, N_ELEMENTS(pvar)
      RETURN, 0
  ENDIF

; number of trial values to use of the starting position of fit

  IF (NOT KEYWORD_SET(N_fit)) THEN N_fit=20

; check that SCALE0 is set if we want to use AMEOBA

  IF (N_fit eq 1) THEN BEGIN
      IF (N_ELEMENTS(scale0) ne N_par) THEN BEGIN
          PRINT,"** SCALE0 is incorrect in QPO_FIT_MODEL."
          RETURN, 0
      ENDIF
  ENDIF

; ----------------------------------------------------------

  IF (KEYWORD_SET(nodisp) eq 0) THEN BEGIN
      PRINT,"-- [Iteration, S_min, parameters, ...]"
  ENDIF

  FOR k = 1, N_fit DO BEGIN

      p_k = p0

      IF (k gt 1) THEN p_k = p_k + RANDOMN(seed, N_ELEMENTS(pvar)) * pvar

; re-normalise before starting interative minimisation by stepping
; through values of NORM and selecting the one that minimises the fit
; function (with all other parameters fixed at their starting values).

      IF KEYWORD_SET(norenorm) THEN BEGIN

          p_k = p_k

      ENDIF ELSE BEGIN

          nn=200
          norm = (INDGEN(nn) - nn/2.0) * 0.1 + p0[0]
          ll = MAKE_ARRAY(nn)
          p_i = p0
          FOR i = 0, nn-1 DO BEGIN
              p_i[0] = norm[i]
              ll[i] = mlogl(p_i)
          ENDFOR
          result = MIN(ll, j)

; use the best normalisation that was found

          p_k[0] = norm[j]

      ENDELSE

; now, to fit the model - by minimising S
; using the POWELL routine

      minl = mlogl(p_k)
      IF (KEYWORD_SET(nodisp) eq 0) THEN BEGIN
          PRINT, FORMAT='("-- [", i5, f10.2, 8(f8.2), "]")', k, minl, p_k
      ENDIF

      IF (N_fit gt 1) THEN BEGIN

          Xi = MAKE_ARRAY(N_par, N_par, VALUE=0.0)
          ii = INDGEN(N_par)
          Xi[ii, ii] = 1.0
          POWELL, p_k, Xi, 1.0e-6, minl, "MLOGL", /DOUBLE, ITMAX=10000L

      ENDIF ELSE BEGIN

          p_k = AMOEBA(1.0e-7, FUNCTION_NAME="MLOGL", FUNCTION_VALUE=minl, $ 
                       P0=p_k, SCALE=scale0, NMAX=10000L)

      ENDELSE

      IF (KEYWORD_SET(nodisp) eq 0) THEN BEGIN
          PRINT, FORMAT='("-- [", 5x, f10.2, 4(f8.2), "]")', minl[0], p_k
      ENDIF

; check for errors

      IF (N_ELEMENTS(p_k) eq 1) THEN BEGIN
          IF (p_k eq -1) THEN BEGIN
              PRINT,"** Fit failed in QPO_FIT_MODEL."
              RETURN, 0
          ENDIF
      ENDIF

; If this is the first fit attempt then store the results

      IF (k eq 1) THEN BEGIN
          Smin = minl[0]
          par = p_k
      ENDIF 

; check whether current fit is the best one so far.

      IF (minl[0] lt Smin) THEN BEGIN
          Smin = minl
          par = p_k
      ENDIF
  
  ENDFOR

; ----------------------------------------------------------
; calculate degrees of freedom

  nu = n_f - N_par

; calculate the Poisson noise level
;                                                 xxx

  N_nonC = N_par-N_spec
  C = par[N_nonC:(N_par-1)]
  IF (MIN(C) gt -30) THEN BEGIN 
      C_n = DOUBLE(10.0^C)
  ENDIF ELSE BEGIN
      C_n = REPLICATE(0.0, N_spec)
  ENDELSE

; calculate the best-fit model (incl. constant)

  par0 = MAKE_ARRAY(N_nonC+1)
  par0[0:(N_nonC-1)] = par[0:(N_nonC-1)]

  prob = MAKE_ARRAY(N_spec)

  for i = 0, N_spec-1 do begin

      mask = WHERE(q eq i+1)

      par0[N_nonC] = C[i]

      model_y = p_model(x[mask], par0, model) 

; calculate the residuals res = 2*data/model if the model is ok these
; should follow a chi-square distribution (with 2 dof) see:
; http://en.wikipedia.org/wiki/Chi-square

      res = y[mask]/model_y

; plot the data and model in one panel
; plus the data/model residuals in a second panel

      IF KEYWORD_SET(graph) THEN BEGIN

; plot the periodogram data

          mask = WHERE(q eq i+1)

          str = STRTRIM(STRING(i+1),2)
          IF (model eq 0) THEN title="Periodogram "+str+" with MODEL_0"
          IF (model eq 1) THEN title="Periodogram "+str+" with MODEL_1"
          
          xrange = [MIN(x)*0.8, MAX(x)/0.8]
          yrange = [MIN(y)*0.8, MAX(y)/0.8]
          PLOT, x[mask], y[mask], /XLOG, /YLOG, PSYM=10, POSITION=[0.15,0.35,0.95,0.95], $
            YTITLE="Power ([rms/mean]!U2!N Hz!U-1!N)", $
            XRANGE=xrange, xtickname=REPLICATE(' ', 30), /XSTYLE, $
            YRANGE=yrange, /YSTYLE, TITLE=title, /NODATA
          
          OPLOT, x[mask], y[mask], COLOR=200, PSYM=10
          
          OPLOT, x[mask], model_y

; overlay the Poisson noise level (if not too small)
          
          IF (MIN(C) gt -30) THEN BEGIN 
              OPLOT, x[mask], model_y - C_n[i], LINESTYLE=2
              OPLOT, xrange, [C_n[i], C_n[i]], LINESTYLE=2
          ENDIF

; *** OPTION [added 23/06/2009] ***
; overlay the 90% confidence band on the power model

          pr = 1.00 - (1.00 - 0.1)^(1.0/n_f)
          cl = CHISQR_CVF(pr, 2)
          OPLOT, x[mask], model_y*(cl/2), LINESTYLE=3

; plot the residuals in second panel

          yrange = [MIN(res)/1.5, MAX(res)*1.5]
          PLOT, x[mask], res, /XLOG, /YLOG, PSYM=10, POSITION=[0.15,0.15,0.95,0.35], $
            XTITLE="Frequency (Hz)", YTITLE="!7c!X = I!Dj!N/M!Dj!N", $
            XRANGE=xrange, /noerase, /XSTYLE, YRANGE=yrange, /YSTYLE
          OPLOT, xrange, [1,1], LINESTYLE=1
;          tmp=""
;          READ, tmp, PROMPT="-- Press ENTER to continue> "
      ENDIF

; compare the distribution of residuals with the chi-square
; distribution. Use KSONE.PRO from IDL Astronomy Library to compare
; data and model - the model comes from the CHISQ2_CDF.PRO function
; which is the cumulative distribution function for chi-square (with 2
; dof). see: http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test 

      KSONE, res*2.0, 'chisq2_cdf', D, pval
      prob[i] = pval

  ENDFOR

; ----------------------------------------------------------
; return the best-fitting parameters to the calling routine

  RETURN, par

END

;
;
;
;
;
;
;
;
;
;
;
;
; ----------------------------------------------------------
;+
; NAME:
;       LRT
;
; PURPOSE:
;       Perform a Likelihood Ratio Test (LRT)
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; INPUTS:
;       Smin_0   - (real) value of -2ln[L] for model_0
;       Smin_1   - (real) value of -2ln[L] for model_1
;       nu_0     - (real) degrees of freedom (dof) for model_0
;       nu_1     - (real) degrees of freedom (dof) for model_1
;
; OPTIONAL INPUTS:
;       disp     - (logical) whether to print results to terminal
;
; OUTPUTS:
;      p_test    - (real) p-value for test
;
; USED BY:
;       QPO_FIT
;
; DETAILS:
;       The two models are then compared using a Likelihood Ratio Test
;       (see Protassov et al. 2002). The LRT statistic is T =
;       -2*ln[L_0/L_1] which is equal to S_0 - S_1 (where L_0 is
;       max(Likelihood) for MODEL_0 and S_0 = -2ln[L_0]).Under certain
;       regularity conditions this should be distributed as a
;       chi-square with nu_0 - nu_1 degrees of freedom (where nu_0 = dof
;       for MODEL_0). We therefore compare T to CHISQR_PDF with
;       nu_0-nu_1 dof and calculate the (upper) tail area
;       probability. This gives a p-value for the test: if p is small
;       (p < ALPHA) then we favour the more complex model (MODEL_1)
;       over the null model. Otherwise we favour the null (MODEL_0). 
;
;       NOTE: As discussed by Freeman et al. (1999) the LRT may be
;       sub-optimal when the null values of additional parameters are
;       not well defined.
;
; Refs:
;       Freeman P., et al., 1999, ApJ, v524, p753
;       Protassov, R., et al., 2002, ApJ, v571, p545
;
;
;-
; ----------------------------------------------------------

FUNCTION LRT, Smin_0, Smin_1, nu_0, nu_1, disp=disp

; check arguments

  IF (N_PARAMS() lt 4) THEN BEGIN
      PRINT,"** Missing arguments in call to LRT."
      RETURN, -1
  ENDIF

; change in 'fit statistic' = LRT

  dS = Smin_0 - Smin_1

; change in degrees of freedom

  dnu = nu_0 - nu_1

; compare likelihood ratio to appropriate chi-square
; distribution, and calculate (upper) tail area

  p_test = 1.0 - CHISQR_PDF(dS, dnu)

; display results

  IF KEYWORD_SET(disp) THEN BEGIN
      PRINT,"-------------------------------------"
      PRINT,"-- Results of Likelihood Ratio Test"
      PRINT, FORMAT='("-- L-ratio:           ",f9.2)',dS
      PRINT, FORMAT='("-- Change in DOF:     ",i9)',dnu
      PRINT, FORMAT='("-- p-value:           ",e9.2)',p_test
      PRINT,"-------------------------------------"
  ENDIF

  RETURN, p_test

END
   
;
;
;
;
;
;
;
;
;
;
;
; ----------------------------------------------------------
;+
; NAME:
;       FIT_RESULTS
;
; PURPOSE:
;       Display the results of model fitting
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; INPUTS:
;       par      - (vector) list of parameter values
;       Smin     - (real) value of -2ln[L] for current model
;       nu       - (real) degrees of freedom (dof) for current model
;       prob     - (real) p-value from KS-test
;       model    - (integer) which model was used
;
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; USED BY:
;       QPO_FIT
;
; DETAILS:
;
;
;-
; ----------------------------------------------------------

PRO FIT_RESULTS, par, err, Smin, nu, prob, model

; check arguments

  IF (N_PARAMS() lt 6) THEN BEGIN
      PRINT,"** Missing arguments in call to FIT_RESULTS."
      RETURN
  ENDIF

  N_par = N_ELEMENTS(par)

  N_spec = N_par - 3
  if (model eq 0) then N_spec = N_par - 2

  N_nonC = N_par-N_spec
  C = par[N_nonC:(N_par-1)]
  IF (MIN(C) gt -30) THEN BEGIN 
      C_n = DOUBLE(10.0^C)
  ENDIF ELSE BEGIN
      C_n = REPLICATE(0.0, N_spec)
  ENDELSE

  PRINT,"-------------------------------------"
  PRINT,"-- Best fitting parameters for model"

  IF (model eq 0) THEN BEGIN
      PRINT, "-- MODEL 0: single power law + const"
      PRINT, "-- Best fitting parameters for model"
      PRINT, FORMAT='("-- ln[Norm]:          ",f9.5," [",2(f9.5),"]")',par[0], err[*, 0]
      PRINT, FORMAT='("-- slope:             ",f9.5," [",2(f9.5),"]")',par[1], err[*, 1]
      FOR i = 0, N_spec-1 do begin
          PRINT, FORMAT='("-- log[C_N]           ",f9.5," [",2(f9.5),"]")',par[2+i], err[*, 2+i]
      ENDFOR
  ENDIF

  IF (model eq 1) THEN BEGIN
      PRINT,"-- MODEL 1: bending power law + const"
      PRINT, FORMAT='("-- log[Norm]:         ",f9.5," [",2(f9.5),"]")',par[0], err[*, 0]
      PRINT, FORMAT='("-- slope_1:           ",f9.5," [",2(f9.5),"]")',1.0,[0.0,0.0]
      PRINT, FORMAT='("-- slope_2:           ",f9.5," [",2(f9.5),"]")',par[1], err[*, 1]
      PRINT, FORMAT='("-- log[f_break]       ",f9.5," [",2(f9.5),"]")',par[2], err[*, 2]
      FOR i = 0, N_spec-1 do begin
          PRINT, FORMAT='("-- log[C_N]           ",f9.5," [",2(f9.5),"]")',par[3+i], err[*, 3+i]
      ENDFOR
  ENDIF

  PRINT, FORMAT='("-- Noise level(s) C_P:",f15.8," [fitted]")',C_n
  PRINT, FORMAT='("-- Fit statistic S:   ",f9.2)',Smin
  PRINT, FORMAT='("-- Degrees of freedom:",i9)',nu
  PRINT,"-- Result of KS test:"
  PRINT, FORMAT='("--   p =              ",f9.5)',prob
;  tmp=""
;  READ, tmp, PROMPT="-- Press ENTER to continue> "

END

;
;
;
;
;
;
;
;
;
;
;
;
; ----------------------------------------------------------
;+
; NAME:
;       LOCATE
;
; PURPOSE:
;       Search a vector for first occurence of specific value
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; INPUTS:
;       x      - (vector) vector to search
;       value  - (float) value to search for
;
; OPTIONAL INPUTS:
;       none
;
; OUTPUTS:
;       pos  - (integer) location of best-match
;
; USED BY:
;       QPO_FIT_ERRORS
;
; DETAILS:
;       A replacement for the built-in routine VALUE_LOCATE. This
;       version simply finds the lowest j for which VALUE is
;       between X[j] and X[j+1]. For small vectors (X) only. If X is
;       large then VALUE_LOCATE is faster.
;
;-
; ----------------------------------------------------------

FUNCTION LOCATE, x, value

; ----------------------------------------------------------

; check arguments

  IF (N_PARAMS() lt 2) THEN BEGIN
      PRINT,"** Missing arguments in call to LOCATE."
      RETURN, -1
  ENDIF

  N = N_ELEMENTS(x)

  IF (N le 2) THEN BEGIN
      PRINT,"** Not enough data to run LOCATE."
      RETURN, -1
  ENDIF

; ----------------------------------------------------------

  pos = -1

  FOR i = 0, N-1 DO BEGIN

      IF (x[i] ge value) THEN BREAK

  ENDFOR

  pos = i - 1

; ----------------------------------------------------------
; return to calling procedure

  RETURN, pos

END

; ----------------------------------------------------------
;
;
;
;
;
;
;
;
;
;
;
; ----------------------------------------------------------
;+
; NAME:
;       QPO_FIT_ERRORS
;
; PURPOSE:
;       Estimate the parameter uncertainties from QPO_FIT
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; INPUTS:
;       par    - (vector) starting values of model parameter(s)
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;      err     - (vector) "errors" from covariance matrix
;
; USED BY:
;       QPO_FIT
;
; DETAILS:
;
;
;-
; ----------------------------------------------------------

FUNCTION QPO_FIT_ERRORS, par, delta=delta, diagn=diagn
 
; ----------------------------------------------------------
; access the common block containing the data points

  COMMON per_data, x, y, q, model

; check arguments

  IF (N_PARAMS() lt 1) THEN BEGIN
      PRINT,"** Missing arguments in call to QPO_FIT_ERRORS."
      RETURN, -1
  ENDIF

  N_p = N_ELEMENTS(par)

  IF (NOT(KEYWORD_SET(delta))) THEN delta=1.0


; ----------------------------------------------------------
; calculate the errors by "brute force" using dS = 1.0

; calculate change in S = -2ln[L] as function of one parameter
; while allowing all other parameters to vary in the fit
; NOTE: do not bother with normalisation (since each fit
; is renormalised anyway)

  err_del = MAKE_ARRAY(2, N_p)

  S0 = MLOGL(par)
  err_0 = MAKE_ARRAY(N_p)

; loop over each parameter in turn

  FOR i = 0, N_p-1 DO BEGIN

      xtitle="Parameter"+STRING(i+1,format='(i1)')
      title="MODEL_"+STRING(model,format='(i1)')
      xtitle=STRCOMPRESS(xtitle)
      title=STRCOMPRESS(title, /REMOVE_ALL)

; explore increasing parameter values

      p = par
      scale = ABS(p) * 0.15
      scale[i] = 0.0
      nn = 15
      dp = INDGEN(nn) / FLOAT(nn)
      dS = MAKE_ARRAY(nn, /DOUBLE)

      FOR j = 0, nn-1 DO BEGIN
          p = par
          p[i] = p[i] + dp[j]
          temp = qpo_fit_model(p, Smin=Smin, N_fit=1, scale=scale, /NODISP, /NORENORM)  
          dS[j] = DOUBLE(Smin - S0)
          IF KEYWORD_SET(diagn) THEN PRINT,"-- ", j, dp[j], dS[j]
      ENDFOR
 
      pos = LOCATE(dS, delta*2.0)

      IF (pos ge nn-1) THEN BEGIN
          c_u = dp[pos] * 5.0
          xrange=[-5,5]
      ENDIF ELSE BEGIN
          c_u = dp[pos+1]
          xrange=[-1,1]
      ENDELSE

      xrange = xrange + par[i]
      yrange = [-0.5, delta*4]
      PLOT, [0,0], [0,0], YRANGE=yrange, /YSTYLE, XTITLE=xtitle, $
        YTITLE="dS = S[x]-S[x0]", TITLE=title, XRANGE=xrange, /NODATA, /XSTYLE
      OPLOT, xrange, [delta, delta], LINESTYLE=2
      OPLOT, dp+par[i], dS, LINESTYLE=2

      nn = 50
      dp = INDGEN(nn) / FLOAT(nn) * c_u
      dS = MAKE_ARRAY(nn, /DOUBLE)

      FOR j = 0, nn-1 DO BEGIN
          p = par
          p[i] = p[i] + dp[j]
          temp = qpo_fit_model(p, Smin=Smin, N_fit=1, scale=scale, /NODISP, /NORENORM)  
          dS[j] = DOUBLE(Smin - S0)
          IF KEYWORD_SET(diagn) THEN PRINT,"-- ", j, dp[j], dS[j]
      ENDFOR

      OPLOT, dp+par[i], dS

      pos = LOCATE(dS, delta)
      IF (pos lt nn-1) THEN BEGIN
          c_u = dp[pos+1]
          S_crit = dS[pos+1]
      ENDIF ELSE BEGIN
          c_u = !VALUES.F_Infinity
          S_crit = dS[pos]
      ENDELSE

      OPLOT, [c_u]+par[i], [S_crit], PSYM=6

; explore decreasing parameter values

      nn = 15
      dp = INDGEN(nn) / FLOAT(nn)
      dS = MAKE_ARRAY(nn, /DOUBLE)
      dp = -dp

      FOR j = 0, nn-1 DO BEGIN
          p = par
          p[i] = p[i] + dp[j]
          temp = qpo_fit_model(p, Smin=Smin, N_fit=1, scale=scale, /NODISP, /NORENORM)  
          dS[j] = DOUBLE(Smin - S0)
          IF KEYWORD_SET(diagn) THEN PRINT,"-- ", j, dp[j], dS[j]
      ENDFOR
      OPLOT, dp+par[i], dS, LINESTYLE=2

      pos = LOCATE(dS, delta*2.0)
      IF (pos ge nn-1) THEN BEGIN
          c_l = dp[pos] * 5.0
      ENDIF ELSE BEGIN
          c_l = dp[pos+1]
      ENDELSE

      nn = 50
      dS = MAKE_ARRAY(nn, /DOUBLE)
      dp = INDGEN(nn) / FLOAT(nn) * c_l

      FOR j = 0, nn-1 DO BEGIN
          p = par
          p[i] = p[i] + dp[j]
          temp = qpo_fit_model(p, Smin=Smin, N_fit=1, scale=scale, /NODISP, /NORENORM)  
          dS[j] = DOUBLE(Smin - S0)
          IF KEYWORD_SET(diagn) THEN PRINT,"-- ", j, dp[j], dS[j]
      ENDFOR

      OPLOT, dp+par[i], dS

      pos = LOCATE(dS, delta)
      IF (pos lt nn-1) THEN BEGIN
          c_l = dp[pos+1]
          S_crit = dS[pos+1]
      ENDIF ELSE BEGIN
          c_l = -!VALUES.F_Infinity
          S_crit = dS[pos]
      ENDELSE

      OPLOT, [c_l]+par[i], [S_crit], PSYM=6

      tmp=""
 ;     READ, tmp, PROMPT="-- Press ENTER to continue> "

      err_del[0, i] = c_l
      err_del[1, i] = c_u

  ENDFOR

; ----------------------------------------------------------
; return the errors to the calling function

  RETURN, err_del

END
        
; ----------------------------------------------------------
;
;
;
;
;
;
;
;
;
;
;
;


PRO QPO_FIT, filename, parm0=p0, parm1=p1, mask=mask, alpha=alpha, $   ; xxx
             beta=beta, plot=plot, ps=ps, N_fit=N_fit, $
             diagn=diagn, delta=delta, p_test=p_test

; ----------------------------------------------------------
;+
; NAME:
;       QPO_FIT
;
; PURPOSE:
;       Search for QPOs in time series
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       qpo_fit,"full.lc"
;
; INPUTS:
;       filename - List of N filenames, each containing data
;
; OPTIONAL INPUTS:
;       parm0 - (vector) initial parameters for MODEL_0
;       parm0 - (vector) initial parameters for MODEL_1
;       mask  - (vector) 2*N-component list of [start, stop] 
;               times of a subset of data to use.
;       alpha - (scalar) defines the significance level of the
;                 test between the two models. Must be in range [0,1]
;       beta  - (scalar) defines the significance level of the
;                 test between for QPO peaks. Must be in range [0,1]
;       plot  - (logical) show graphical output of results?
;       ps    - (logical) output graphics to PostScript files
;       diagn - (logical) display diagnostics while running
;       delta - (float) delta fit-statistic value to use for
;                 calculating confidence intervals.
;                = 1.00 --> 68.3% [default] ("1 sigma")
;                = 2.71 --> 90.0%
;                = 4.00 --> 95.4% ("2 sigma") 
;                = 6.63 --> 99.0%
;
; OUTPUTS:
;       none
;
; OPTIONAL OUTPUTS:
;       p_test  - (float) p-value for LRT (model0 vs. model1)
;
; DETAILS:
;       This procedure performs a search for QPOs in an input time
;       series. The input data are read from anASCII file, with time
;       in the first column and the time series variable in the second
;       column.
;
;       In outline the steps are:
;        1. Calculate periodogram of data
;        2. Fit periodogram with simple continuum model
;              (power law plus constant: MODEL_0)
;        3. Fit periodogram with complex continuum model
;              (broken power law plus constant: MODEL_1)
;        4. Select between the two models using a Likelihood Ratio Test
;        5. Examine the residuals from the selected model
;        6. Find the most extreme outlier of the residuals
;        7. Compute where this outlier is in the expected distribution
;              of residuals
;        8. Apply a correction of the number of data points examined 
;              (Bonferroni correction)
;        9. Repeat the examination of residuals but after smoothing 
;              the data with a 2-bin width filter
;
;       The periodogram of the data is computed. This is fitted with a
;       model using the method of Maximum Likelihood to optimise the
;       parameters. In fact it is the minus log-likelihood function S
;       = -2ln[L] that is minimised - see the appendix of Vaughan
;       (2005) for details and Anderson et al. (1990). The
;       minimisation is done using the POWELL routine.
;
;       Two models are used as standard:
;         * MODEL_0 - power law + constant (3 params)
;         * MODEL_1 - broken power law + constant (5 params)
;       For each model we minimise S to find 
;         - the best-fitting parameters  
;            (i.e. the Maximum Likelihood Estimates; MLEs)
;         - the miminim of S, i.e. S_min
;         - the degrees of freedom: nu = N - M
;         - result of a KS-test 
;       For each model we perform a KS-test comparing the data/model
;       residuals to their expected distribution if the model is 
;       correct (a chi-square distribution). This is a goodness-of-fit 
;       test -  a simple check of how good/bad the model fit is.
;
;       NOTE: In order to make the model parameters "well behaved"
;       - i.e. easier and more stable to fit - the fitting is done
;       using the logarithm of the normalisation and the break
;       frequency.
;
;       The two models are then compared using a Likelihood Ratio Test
;       (see Protassov et al. 2002). The LRT statistic is T =
;       -2*ln[L_0/L_1] which is equal to S_0 - S_1 (where L_0 is
;       max(Likelihood) for MODEL_0 and S_0 = -2ln[L_0]).Under certain
;       regularity conditions this should be distributed as a
;       chi-square with nu_0 - nu_1 degrees of freedom (where nu_0 = dof
;       for MODEL_0). We therefore compare T to CHISQR_PDF with
;       nu_0-nu_1 dof and calculate the (upper) tail area
;       probability. This gives a p-value for the test: if p is small
;       (p < ALPHA) then we favour the more complex model (MODEL_1)
;       over the null model. Otherwise we favour the null (MODEL_0). 
;
;       NOTE: As discussed by Freeman et al. (1999) the LRT may be
;       sub-optimal when the null values of additional parameters are
;       not well defined.In this case the break frequency of MODEL_1
;       has no defined null value in MODEL_0. 
;
;       With a model selected we now examine the data/model residuals
;       to search for peaks that could be QPOs. We compute RES =
;       DATA/MODEL. If the model is correctthen 2*RES should be
;       distributed as a chi-square variablewith 2 dof. We can use
;       this to convert from RES to Prob(RES > RES_obs) and find out
;       where the observed values lie in the expected
;       distribution. The next step is to correct for the fact that we
;       are interested in the most extreme of N_F frequency points,so
;       we apply a Bonferroni correction to account forsearching
;       through so much data to find one extreme point.
;       
;       The process of examining the residuals and comparing to the
;       expected chi-square distribution is then repeated after 
;       smoothing the residuals using a 3-bin filter. Neighbouring
;       pairs of periodogram data points are averagedto improve the
;       sensitivity to QPOs that lie between two Fourier frequencies
;       (that have their power splitbetween two frequencies rather 
;       than concentrated in one. See section 6 of van der Klis,
;       1989). The smoothed periodogram is then compared to the best
;       fitting model (from above) to produce smoothed residuals. As
;       the sum of 3 chi-square (dof=2) variables is one chi-square
;       (dof=6) variable we can compare the residuals to the
;       chi-square (dof=6) distribution, and compute p-values for the
;       extreme points as before.
;
; Refs:
;       Anderson E. R., et al. 1990, ApJ, v364, p699
;       Freeman P., et al., 1999, ApJ, v524, p753
;       Protassov, R., et al., 2002, ApJ, v571, p545
;       van der Klis, M. 1989, in "Timing Neutron Stars", p27
;       Vaughan S., 2005, A&A, v431, p391
;
; EXAMPLE USAGE:
;       qpo_fit, "full.lc", N_fit=50, /PLOT
;       qpo_fit, "full.lc", mask=[24500, 91400]
;
; EXTERNAL PROCEDURES CALLED:
;       KSONE.PRO
;       PERIODGRAM.PRO
;       PLOT_ERR.PRO
;       READ_TABLE.PRO
;       PS_OPEN.PRO
;       PS_CLOSE.PRO
;
; INTERNAL PROCEDURES CALLED:
;       P_MODEL0
;       P_MODEL1
;       P_MODEL
;       MLOGL
;       QPO_FIT_MODEL
;       LRT
;       FIT_RESULTS
;       LOCATE
;       QPO_FIT_ERRORS
;
;
; HISTORY:
;         13/02/2009 - v1.0 - first version
;                             (based on old script)
;         17/02/2009 - v1.1 - added PS keyword
;                             Modularised the fitting.
;                             Added multiple initial positions
;                              to check for global minima.
;         19/02/2009 - v1.2 - Replaced models with single function
;                              to complute models. Now use a P_MASK
;                              array to fix/free parameters.
;         20/02/2009 - v1.3 - Added function to calculate numerical
;                              derivatives of fit function.
;         23/02/2009 - v1.4 - Added step where the largest outlying
;                              data point is excluded from the
;                              continuum fit prior to calculating its
;                              p-value.
;                             Put in a seperate procedure the display
;                              of model parameters after fitting.
;                             Put in a seperate procedure the LRT.
;         02/03/2009 - v1.5 - Use POWELL for minimisation (replaced
;                             AMOEBA). Use smoothly bending power law
;                             model in place of sharply broken
;                             model. Included steps to avoid
;                             overflow/underflow errors in model
;                             calculation. Added routine to compute
;                             derivative of log[L] w.r.t. parameters. 
;         09/03/2009 - v1.6 - Added routine to compute confidence.
;                             intervals on the fitted parameters.
;                             Uses the 'graphical method' on finding
;                             where S = -2ln[L] changes by some
;                             amount. (Similar to the delta chi^2
;                             method.) The confidence level is set
;                             using the DELTA keyword.
;         29/04/2010 - v1.7 - Changed MODEL1 to have a fixed lower
;                             frequency slope. This allows for much
;                             simpler fitting of the model and less
;                             degeneracy in the parameters (and their
;                             errors). 
;         23/06/2011 - v1.8 - Removed HESSIAN function and covariance
;                             calculation (for speed). Changed C_n to
;                             a free parameter in the model(s) rather
;                             than a (fixed) constant. Replaced 2-bin
;                             smoothing with 3-bin smoothing in the
;                             final stage of the analysis. Changed
;                             residual plots to log[data/model]
;         25/08/2011 - v1.9 - Modified several routines to allow
;                             simultaneous fitting of multiple
;                             datasets using the same model (but
;                             different Poisson noise levels). Changed 
;                             FINDFILE for FILE_SEARCH in IDL 8+.
;
; NOTES:
;
;-
; ----------------------------------------------------------
;
; options for compilation (recommended by RSI)

  COMPILE_OPT idl2

; watch out for errors

  ON_ERROR, 0

; ----------------------------------------------------------
; Check the arguments

; is the file name defined?

  N_spec = N_ELEMENTS(filename)              ; xxx
  IF (N_spec eq 0) THEN BEGIN
      filename=''
      READ,'-- Enter file name (ENTER to list current directory): ',filename
      IF (filename eq '') THEN BEGIN
          list = FILE_SEARCH()
          PRINT, list
          READ,'-- Enter file name: ',filename
      ENDIF
  ENDIF

; initialise the parameters for MODEL_0 and MODEL_1

  N0 = N_ELEMENTS(parm0)                     ; xxx
  N1 = N_ELEMENTS(parm1)
  IF (N0 ne 3) THEN parm0 = DOUBLE([-3.0, 1.5, 0])
  IF (N1 ne 4) THEN parm1 = DOUBLE([-3.0, 3.5, -3.5, 0])
  N0 = 3
  N1 = 4

; select a "significance level" for selecting
; between MODEL_0 and MODEL_1
; if LRT p-value < alpha then select MODEL_1
; E.g. alpha = 0.05 <--> "95% significance"

  alpha0 = 0.05
  IF (NOT KEYWORD_SET(alpha)) THEN alpha = alpha0
  IF (alpha ge 1.0 OR alpha le 0.0) THEN BEGIN
      PRINT,"** ALPHA is not within [0,1]; reset to default ", alpha0
      alpha = alpha0
  ENDIF

; select a "significance level" for QPO search
; if global p-value < beta then detect QPO

  beta0 = 0.01
  IF (NOT KEYWORD_SET(beta)) THEN beta = beta0
  IF (beta ge 1.0 OR beta le 0.0) THEN BEGIN
      PRINT,"** BETA is not within [0,1]; reset to default ", beta0
      beta = beta0
  ENDIF

; define graphics setting.
; 0 = no graphics
; 1 = graphics on screen 
; 2 = graphics to PS file

  IF KEYWORD_SET(plot) THEN graph=1
  IF KEYWORD_SET(ps) THEN graph=2
  IF (NOT KEYWORD_SET(graph)) THEN graph=0

; open PS device if needed

  IF (graph eq 2) THEN PS_OPEN

; number of trial values to use of the starting position of fit

  IF (NOT KEYWORD_SET(N_fit)) THEN N_fit=20

; number of trial values to use of the starting position of fit

  IF (NOT KEYWORD_SET(delta)) THEN delta=1.0
  IF (delta le 0.0) THEN BEGIN
       PRINT,"** DELTA <= zero in QPO_FIT: ", delta
       RETURN
  ENDIF


; ----------------------------------------------------------
;
; load data

  for i_file = 0, N_spec-1 do begin                 ; xxx

      data = READ_TABLE(filename[i_file], /head, /double)

; if no data loaded then finish

      IF (N_ELEMENTS(data) eq 1) THEN RETURN

; are error supplied?

      ndim = (SIZE(data, /DIM))[0]

; sort into vectors for time (t), count rate (r) and its error (dr)

      t = REFORM(data[0,*])
      r = REFORM(data[1,*])
      dr = r-r
      if (ndim ge 3) then dr = REFORM(data[2,*])
      dt = t[1] - t[0]
      n_t = N_ELEMENTS(t)

; select a subset of data if requested (XXX this needes fixing for
; N_spec > 1 case XXX)

      IF (N_ELEMENTS(mask) eq 2) THEN BEGIN
          subset = WHERE((t gt mask[0] and t lt mask[1]), count)
          IF (count eq 0) THEN PRINT,"** No data in selected time interval"
          t = t[subset]
          r = r[subset]
          dr = dr[subset]
      ENDIF

; estimate the Poisson noise level in the data

      C_n = MEAN(dr*dr)*2.0*dt/MEAN(r)^2
      if KEYWORD_SET(pois) then C_n = 2.0/MEAN(r)
      C_n = (C_n > 1e-8)

; give some feedback to the user

      PRINT,"-------------------------------------"
      PRINT,"-- Time series loaded      ",filename[i_file]   ; xxx
      PRINT,"-- Number of data points:  ",n_t
      PRINT,"-- Sampling rate:          ",dt        
      PRINT,"-- Mean count rate         ",MEAN(r)
      PRINT,"-- Estimated Poisson level ",C_n

; ----------------------------------------------------------
; plot the time series data
; with errors (using PLOT_ERR.PRO)

      IF KEYWORD_SET(graph) THEN BEGIN
          xrange=[MIN(t), MAX(t)]
          PLOT, t, r, XRANGE=xrange, /XSTYLE, TITLE="Time series", $
            XTITLE="Time (s)", YTITLE="Rate (ct/s)"
          PLOT_ERR, t, r, dr
          tmp=""
;          READ, tmp, PROMPT="-- Press ENTER to continue> "
          IF (tmp eq "q" OR tmp eq "Q") THEN RETURN
      ENDIF

; ----------------------------------------------------------
; calculate periodogram (in RMS units) using PERIODOGRAM.PRO
; see: http://en.wikipedia.org/wiki/Discrete_Fourier_transform

      per = PERIODOGRAM(r,dt=dt,f=f,df=df,/rms)
      n_f = N_ELEMENTS(f)
      flag = REPLICATE(i_file+1, n_f)          ; xxx

      if (i_file eq 0) then begin
          merged_per = per
          merged_nf = n_f
          merged_f = f
          merged_flag = flag
      endif else begin
          merged_per = [merged_per, per]
          merged_nf = merged_nf + n_f
          merged_f = [merged_f, f]
          merged_flag = [merged_flag, flag]
      endelse

  endfor

  per = merged_per
  flag = merged_flag
  f = merged_f
  n_f = merged_nf

; ----------------------------------------------------------
; place the data points into a common block so they are accessible to 
; the (likelihood) function

                                         ; xxx
  COMMON per_data, x, y, q, model            
  x = f
  y = per
  q = flag

; ----------------------------------------------------------
; start with simple power law model; model_0
; parameters that are NOT to be varied have p_mask[m] = 0

  model = 0

; convert Poisson noise level to log units

  C = ALOG10(C_n)

; set-up the parameter vector                 xxx
; (with one C value for each spectrum if they are multiple)

  p0 = MAKE_ARRAY(N0-1+N_spec)
  p0[0:(N0-2)] = parm0[0:(N0-2)]
  p0[(N0-1):(N0-2+N_spec)] = REPLICATE(C, N_spec)

; define the parameter settings
; P_MASK shows which parameters are to be varied (0 = fixed)
; P0 is the initial setting of the parameter values prior to fitting
; PVAR gives the standard deviation used to distribute starting value
;  when the fit is repeated for randomly spaced starting positions

  scale0 = ABS(p0)*0.5                                 ; xxx

; perform the fitting 

  par_0 = qpo_fit_model(p0, nu=nu_0, Smin=Smin_0, res=res_0, graph=graph, $ 
                        N_fit=N_fit, prob=prob, scale=scale0)


; calculate errors

  err_0 = QPO_FIT_ERRORS(par_0, DELTA=delta, DIAGN=KEYWORD_SET(diagn))

; display results

  FIT_RESULTS, par_0, err_0, Smin_0, nu_0, prob, model

; ----------------------------------------------------------
; fit broken power law model: model_1

  model = 1

; set-up the parameter vector                 xxx
; (with one C value for each spectrum if they are multiple)

  p1 = MAKE_ARRAY(N1-1+N_spec)
  p1[0:(N1-2)] = parm1[0:(N1-2)]
  p1[(N1-1):(N1-2+N_spec)] = REPLICATE(C, N_spec)

; update parameter settings

  scale1 = ABS(p1)*0.1
 
; perform the fitting 

  par_1 = qpo_fit_model(p1, nu=nu_1, Smin=Smin_1, res=res_1, $
                        graph=graph, N_fit=N_fit, prob=prob, $
                        scale=scale1)  


; calculate errors

  err_1 = QPO_FIT_ERRORS(par_1, DELTA=delta, DIAGN=KEYWORD_SET(diagn))

; display results

  FIT_RESULTS, par_1, err_1, Smin_1, nu_1, prob, model

; ----------------------------------------------------------
; perform Likelihood Ratio Test (LRT) to
; compare model_0 to model_1.

  p_test = LRT(Smin_0, Smin_1, nu_0, nu_1, /DISP)

; make a decision: select MODEL_1 (more complex model) only if 
; p_test < alpha, where alpha is a small number e.g. 0.05 
; (= "95% significance"). Otherwise, chose MODEL_0 by default.

  m_flag = 0
  IF (p_test le alpha) THEN m_flag=1

  IF (m_flag eq 0) THEN PRINT,"-- selected:        MODEL_0"
  IF (m_flag eq 1) THEN PRINT,"-- selected:        MODEL_1"

; pick residuals and parameters for the SELECTED model

  res = res_0
  par = par_0
  scale=scale0
  model = 0
  IF (m_flag eq 1) THEN begin
      res = res_1
      par = par_1
      scale = scale1
      model = 1
  endif

; ----------------------------------------------------------
; plot and calculate probs
 
  N_par = N_ELEMENTS(par)
  N_nonC = N_par-N_spec
  C = par[N_nonC:(N_par-1)]
  IF (MIN(C) gt -30) THEN BEGIN 
      C_n = DOUBLE(10.0^C)
  ENDIF ELSE BEGIN
      C_n = REPLICATE(0.0, N_spec)
  ENDELSE

; calculate the best-fit model (incl. constant)

  par_ref = MAKE_ARRAY(N_nonC+1)
  par_ref[0:(N_nonC-1)] = par[0:(N_nonC-1)]

  for k = 0, N_spec-1 do begin                     ; xxx

      mask = WHERE(q eq k+1)
      xx = f[mask]
      yy = per[mask]

 
      mask = WHERE(q eq k+1)

      par_ref[N_nonC] = C[k]

      model_y = p_model(xx, par_ref, model) 
      res = yy/model_y

; ----------------------------------------------------------
; plot the residuals as probabilities (convert from x to Prob(>x)
; using chi-square) NOTE: this is the probability per frequency, so it
; does not account for the number of independent frequencies
; observed. NOTE: there is no account of the uncertainty in the the
; model. Therefore this analysis is valid assuming the model is
; exactly right.

      prob = 1.0 - CHISQR_PDF(2*res, 2)

; make the correction for the number of independent frequencies
; observed using 
;   p[multiple] = 1 - ( 1 - p[single] )^N
; see: http://en.wikipedia.org/wiki/Multiple_testing
;
; To avoid arithmetic underflow errors, do the calculation using
; logarithms:
;  log( 1 - p[multiple] ) = N*log( 1 - p[single] )

      nprob = MAKE_ARRAY(SIZE(prob,/DIM))

      log_pp = n_f*ALOG(DOUBLE(1.0-prob))

      indx = WHERE(log_pp gt -30.0, count, COMPLEMENT=indx_c)
      IF (count gt 0) THEN nprob[indx] = exp(log_pp[indx])
      IF (count lt n_f) THEN nprob[indx_c] = 0.0
      
      nprob = DOUBLE(1.0 - nprob)
          
; find the minimum p-value - i.e. most "significant" point

      p_min = MIN(nprob, pos)
      f_min = xx[pos]

      PRINT,"-------------------------------------"
      PRINT,"-- Most significant outlier in spectrum", k+1
      PRINT,"-- Frequency (mHz):   ",f_min*1.0e3
      PRINT,"-- Peak residual:     ",res[pos]
      PRINT,"-- Single-trail p:    ",prob[pos]
      PRINT,"-- Global p-value:    ",p_min
      IF (p_min lt beta) THEN BEGIN
          PRINT,"-- POSSIBLE QPO DETECTION ******"

; and plot the global p-values

          IF KEYWORD_SET(graph) THEN BEGIN
              ymin = MIN([beta, nprob])*0.5
              IF (ymin le 1.0e-10) then ymin = 1.0e-10
              yrange=[ymin, 1.0]
              PLOT, xx, nprob, /XLOG, /YLOG, PSYM=10, YRANGE=yrange, $
                XTITLE="Frequency (Hz)",YTITLE="Prob(>x)", $
                TITLE="Global p-values per frequency bin"

; plot the significance theshold (0.01 = "99% signif.")

              OPLOT,[1e-8, 1.0], [beta, beta], LINESTYLE=2
;              READ, tmp, PROMPT="-- Press ENTER to continue> "
              IF (tmp eq "q" OR tmp eq "Q") THEN RETURN
          ENDIF
      ENDIF ELSE BEGIN
          PRINT,"-- No QPO detected"
      ENDELSE


; ----------------------------------------------------------
; Now repeat the analysis of residuals to the selected model,
; but using the average of two neighbouring data points.
; I.e. the data are smoothed with a running 3-bin wide filter.
; This improves sensitive to peaks that lie between two bins,
; which have power spread between the bins rather than concentrated 
; all in a single bin.

; calculate the running three-bin mean of the periodogram
;                                                         xxx

      i = INDGEN(N_ELEMENTS(xx)-2)
      bin_per = (yy[i] + yy[i+1] + yy[i+2]) / 3.0
      bin_f = xx[i+1]

; plot the smoothed data against the maximum likelihood
; model - obtained (above) by fitting the unbinned data

      model_y = P_MODEL(bin_f, par_ref, model) 

; residuals: data/model

      res = bin_per/model_y

; plot the binned periodogram

      str = STRTRIM(STRING(k+1),2)
      title = "Periodogram "+str+" after 3-bin smoothing"
      IF KEYWORD_SET(graph) THEN BEGIN
          xrange = [MIN(f)*0.8, MAX(f)/0.8]
          yrange = [MIN(per)*0.8, MAX(per)/0.8]
          PLOT, bin_f, bin_per, /XLOG, /YLOG, PSYM=10, POSITION=[0.15,0.35,0.95,0.95], $
            YTITLE="Power ([rms/mean]!U2!N Hz!U-1!N)", $
            XRANGE=xrange, xtickname=REPLICATE(' ', 30), /XSTYLE, $
            YRANGE=yrange, /YSTYLE, TITLE=title, /NODATA
          OPLOT, bin_f, bin_per, COLOR=200, PSYM=10

; overlay the Poisson noise level

;      C_n = par[N_ELEMENTS(par)-1]
          IF (min(C) gt -30) THEN BEGIN 
              OPLOT, bin_f, model_y - C_n[k], LINESTYLE=2
              OPLOT, xrange, [C_n[k], C_n[k]], LINESTYLE=2
          ENDIF

; plot the model (previously fitted to the unbinned data)

          OPLOT, bin_f, model_y
 
; plot the residuals in second panel

          yrange = [MIN(res)/1.5, MAX(res)*1.5]
          PLOT, bin_f, res, /XLOG, /YLOG, PSYM=10, POSITION=[0.15,0.15,0.95,0.35], $
            XTITLE="Frequency (Hz)", YTITLE="data/model", $
            XRANGE=xrange, /noerase, /XSTYLE, YRANGE=yrange, /YSTYLE
          OPLOT, xrange, [1,1], LINESTYLE=1
;          READ, tmp, PROMPT="-- Press ENTER to continue> "
          IF (tmp eq "q" OR tmp eq "Q") THEN RETURN
      ENDIF

; plot the residuals as probabilities 
; NOTE: all the caveats given above hold. In addition: as we have
; averaged together two bins with degrees-of-freedom is now
; 2+2=4. Furthermore: the number of 'independent trials' is no longer
; well defined (because the smoothed bins are not independent, they
; are overlapping) 

; NOTE: The factor of two in "res*2" used so that we are comparing the
; sum of two (dof=2) variables (not their average) with a dof=4
; distribution. 

      df = 6
      prob = 1.0 - CHISQR_PDF(res*df, df)

;  PLOT, bin_f, prob, /XLOG, /YLOG, PSYM=10, $
;    XTITLE="Frequency (Hz)", YTITLE="Prob(>x)", $
;    TITLE="Smoothed (2 bin) residuals - single trial"
;
; pause
;
;  READ, tmp, PROMPT="-- Press ENTER to continue> "

; make the correction for the number of independent
; frequencies observed using 
;   p[multiple] = 1 - ( 1 - p[single] )^N
; see: http://en.wikipedia.org/wiki/Multiple_testing
;
; where p is very small  p[multiple] ~ p[single]*N

      n_f = n_elements(bin_f)
      nprob = MAKE_ARRAY(SIZE(prob,/DIM))
      
      log_pp = n_f*ALOG(DOUBLE(1.0-prob))
      
      indx = WHERE(log_pp gt -30.0, count, COMPLEMENT=indx_c)
      IF (count gt 0) THEN nprob[indx] = exp(log_pp[indx])
      IF (count lt n_f) THEN nprob[indx_c] = 0.0
      
      nprob = DOUBLE(1.0 - nprob)
      
; find the minimum p-value - i.e. most "significant" point

      p_min = MIN(nprob,pos)
      f_min = f[pos]
      
      PRINT,"-------------------------------------"
      PRINT,"-- Most significant outlier (after binning)"
      PRINT,"-- Frequency (mHz):   ",f_min*1.0e3
      PRINT,"-- Single-trail p:    ",prob[pos]
      PRINT,"-- Global p-value:    ",p_min
      IF (p_min lt beta) THEN BEGIN
          PRINT,"-- POSSIBLE QPO DETECTION ******"
          
          IF KEYWORD_SET(graph) THEN BEGIN
              ymin = MIN([beta,nprob])*0.5
              IF (ymin le 1.0e-10) then ymin = 1.0e-10
              yrange=[ymin, 1.0]

; and plot the global p-values

              PLOT, bin_f, nprob, /XLOG, /YLOG, PSYM=10, YRANGE=yrange, $
                XTITLE="Frequency (Hz)",YTITLE="Prob(>x)", $
                TITLE="Global p-values per frequency bin (2-bin)"

; plot the significance theshold (0.01 = "99% signif.")

              OPLOT,[1e-8,1.0],[beta, beta],LINESTYLE=2
          ENDIF
      ENDIF ELSE BEGIN
          PRINT,"-- No QPO detected"
      ENDELSE
      PRINT,"-------------------------------------"
      
  ENDFOR

; close PS device if opened

  IF (graph eq 2) THEN PS_CLOSE

; ---------------------------------------------------------

END
