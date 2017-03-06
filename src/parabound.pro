
FUNCTION parabound, x, par, covar, conf=conf, plotq=plotq, $
                    nsim=nsim

; ----------------------------------------------------------
;+
; NAME:
;       PARABOUND
;
; PURPOSE:
;       Calculate confidence bounds on quadratic model
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       result = PARABOUND(x, par, covar)
;
; INPUTS:
;       x       - (array) x-values at which to evaluate model
;       par     - (array) best fitting parameter values
;       covar   - (array) covariance matrix for parameters 
;
; OPTIONAL INPUTS:
;       conf    - (float) confidence level to use (default = 0.9)
;       plotq   - (logical) plot results? (detault = FALSE)
;       nsim    - (integer) number of simulations (default = 500)
;
; OUTPUTS:
;       result  - (array) upper and lower limits to model
;
; OPTIONAL OUTPUTS:
;       none
;
; DETAILS:
;       Calculate the lower and upper (alpha%) confidence bounds
;       around a quadratic model. The best fitting parameters of the
;       model are input as PAR, and the covariance matrix of the
;       parameters are input as COVAR. The best fitting model is then
;       evalated at the values of X. The confidence intervals at each
;       X value are calculated by Monte Carlo as follows.
;
;       We use the covariance matrix to drawn random numbers from a
;       multi-variate Normal distributuion, and add the best fitting
;       values. This then gives a set of random numbers drawn from the
;       joint distribution of the parameter estimates. The random
;       numbers are generated using the NORMAL function. We then
;       calculate a model from each of these sets of parameters,
;       giving a distribution of Y at each X, and find
;       the alpha*100 and (1-alpha)*100 percentiles at each of the X
;       values. The lower and upper limits of the interval at each X
;       value are returned in the output array. 
;       
;       The number of simulations used is set using the NSIM keyword
;       (detault = 500). 
;
;       Note that the confidence level is related to alpha by alpha =
;       (1-conf)/2, so that with conf=0.90 we calculate the lower 5%
;       and upper 95% limits, between which lies the 90% confidence
;       interval. 
;
; EXAMPLE USAGE:
;
;  n = 6                         ; define some X data
;  x = INDGEN(n)
;  
;  par = [1.0, -1.0, 0.5]        ; define a Y model
;  y = par[2]*x*x + par[1]*x + par[0]
;  err = MAKE_ARRAY(n, value=0.4)
;  y = y + RANDOMN(seed, n)*err  ; randomise Y data
;  
;  PLOT, x, y, PSYM=1, XRANGE=[-0.5,8], yrange=[-2,20], /XSTYLE, /YSTYLE
;  PLOT_ERR, x, y, err
;  result = POLY_FIT(x, y, 2, COVAR=covar, MEASURE_ERRORS=err)
;  par = REFORM(result)                     ; reformat output
;  
;  xmod = INDGEN(100)/10.0 - 2.0            ; define model X values
;  cbounds = PARABOUND(xmod, par, covar)    ; find intervals
;  OPLOT, xmod, cbounds[0,*], LINESTYLE=1   ; lower limit
;  OPLOT, xmod, cbounds[1,*], LINESTYLE=1   ; upper limit
;
; PROCEDURES CALLED:
;          NORMAL
;
; HISTORY:
;       05/05/2010 -  v1.0  - first working version
;;
; NOTE:
;       none
;
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)
  COMPILE_OPT idl2

; watch out for errors
  on_error,2
 
; ----------------------------------------------------------
; Check the arguments

; is the input data array well-defined?

  n = n_elements(x)
  if (n eq 0) then begin
      print,'** No data in PARABOUND.'
      return,0
  endif

; are the parameters and their covariance defined?

  d = N_ELEMENTS(par)
  if (d ne 3) then begin
      print,'** ',d,'parameters found. Expecting 3 in PARABOUND.'
      return,0
  endif

  dd = N_ELEMENTS(covar)
  if (dd ne 9) then begin
      print,'** ',dd,'covariance elements found. Expecting 9 in PARABOUND.'
      return,0
  endif

; has the user selected the confidence level?
  if not KEYWORD_SET(conf) then conf=0.90

; has the user selected the number of simulations?
  if not KEYWORD_SET(nsim) then nsim=500

; ----------------------------------------------------------
; Main part of routine

; draw NSIM random numbers from a multi-variate Normal distribution
; with covariance matrix COVAR

  params = NORMAL(seed, covar=covar, n=nsim)

; add the best fitting values from PAR to 'center' the distribution
; over the best fit 

  params = params + REBIN(par, SIZE(params, /DIMENSIONS))

; now loop over every set of parameters. 
; For each one calculate a model Y = a + bX + cX^2 

  ysims = MAKE_ARRAY(nsim, n)

  for i = 0, nsim-1 do begin
      a = params[0,i]
      b = params[1,i]
      c = params[2,i]
      ysims[i,*] = c*x*x + b*x + a
;      OPLOT, x, ysims[i,*]
  endfor

; find the upper (alpha*100%) and lower ([1-alpha]*100%) percintiles
; of the distribution of Y's at each X value. Do this by sorting the Y
; values at each X into accending order, then finding the element
; closest to the required percentile.

  alpha = (1.0-conf)/2.0

  n_hi = ROUND((1.0-alpha)*nsim)
  n_lo = ROUND(alpha*nsim)

  ylo = MAKE_ARRAY(n)
  yhi = MAKE_ARRAY(n)

;  plot, params[0,*], params[1,*], psym=1

  for j = 0, n-1 do begin
      yslice = ysims[*,j]
      dist = yslice[SORT(yslice)]
      ylo[j] = dist[n_lo]
      yhi[j] = dist[n_hi]
  endfor

; ----------------------------------------------------------
; Optional plot

  if KEYWORD_SET(plotq) then begin

      ymod = par[2]*x*x + par[1]*x + par[0]
      plot, x, ymod
      oplot, x, ylo, LINESTYLE=4
      oplot, x, yhi, LINESTYLE=4

  endif

; ----------------------------------------------------------
; finished. Return the result to the user

  result = TRANSPOSE([[ylo], [yhi]])
  RETURN, result

END

; ----------------------------------------------------------
