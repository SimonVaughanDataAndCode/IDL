PRO LOG_PERIODOGRAM_FIT, filename, FREQ_RANGE=freq_range, DT=dt

; ----------------------------------------------------------
;+
; NAME:
;       LOG_PERIODGRAM_FIT
;
; PURPOSE:
;       Simple power law fit to periodogram data.
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       qpo_fit,"full.lc"
;
; INPUTS:
;       filename - Name of file containing data
;
; OPTIONAL INPUTS:
;       FREQ_RANGE - (array) Lower, upper frequencies to fit
;       DT         - (float) sampling interval
;
; OUTPUTS:
;       none
;
; OPTIONAL OUTPUTS:
;       none
;
; DETAILS:
;       The periodogram of the data is computed. This is fitted with a
;       simple power law model using the method of Maximum Likelihood
;       to optimise the 
;       parameters. In fact it is the minus log-likelihood function S
;       = -2ln[L] that is minimised - see the appendix of Vaughan
;       (2005), Vaughan (2010) and Anderson et al. (1990). The
;       minimisation is done using the POWELL routine.
;
;       For each model we perform a simple goodness-of-fit test
;       comparing the data/model 
;       residuals to their expected distribution if the model is 
;       correct (a chi-square distribution). 
;
;       With the model fitted we examine the data/model residuals
;       to search for peaks that could be QPOs. We compute RES =
;       2*DATA/MODEL. If the model is correct then RES should be
;       distributed as a chi-square variable with 2 dof. We can use
;       this to convert from RES to prob(RES > RES_obs) and find out
;       where the observed values lie in the expected
;       distribution. The next step is to correct for the fact that we
;       are interested in the most extreme of N_F frequency points,so
;       we apply a Bonferroni correction to account forsearching
;       through so much data to find one extreme point.
;       
;       The process of examining the residuals and comparing to the
;       expected chi-square distribution is then repeated after 
;       smoothing the residuals using a 2-bin filter. Neighbouring
;       pairs of periodogram data points are averaged to improve the
;       sensitivity to QPOs that lie between two Fourier frequencies
;       (that have their power splitbetween two frequencies rather 
;       than concentrated in one. See section 6 of van der Klis,
;       1989). The smoothed periodogram is then compared to the best
;       fitting model (from above) to produce smoothed residuals. As
;       the sun of two chi-square (dof=2) variables is one chi-square
;       (dof=4) variable we can compare the residuals to the
;       chi-square (dof=4) distribution, and compute p-values for the
;       extreme points as before.
;
; Refs:
;       Anderson E. R., et al. 1990, ApJ, v364, p699
;       van der Klis, M. 1989, in "Timing Neutron Stars", p27
;       Vaughan S., 2005, A&A, v431, p391
;
; EXAMPLE USAGE:
;       LOG_PERIDOGRAM_FIT, "mydata.dat"
;
; EXTERNAL PROCEDURES CALLED:
;       PERIODGRAM.PRO
;       READ_TABLE.PRO
;
; INTERNAL PROCEDURES CALLED:
;       none
;
; HISTORY:
;         30/10/2010 - v1.0 - first version
;
;-
; ----------------------------------------------------------
;
; options for compilation (recommended by RSI)

  COMPILE_OPT idl2

; watch out for errors

  ON_ERROR, 2

; ----------------------------------------------------------
; Check the arguments

; is the file name defined?

  IF (N_ELEMENTS(filename) eq 0) THEN BEGIN
      filename=''
      READ,'-- Enter file name (ENTER to list current directory): ',filename
      IF (filename eq '') THEN BEGIN
          list = FINDFILE()
          PRINT, list
          READ,'-- Enter file name: ',filename
      ENDIF
  ENDIF

; ----------------------------------------------------------
; load the time series data

  data = READ_TABLE(filename,/head,/double)

; if no data loaded then finish

  IF (N_ELEMENTS(data) eq 1) THEN RETURN

; sort into vectors for time (t), count rate (r)

  t = REFORM(data[0,*])
  x = REFORM(data[1,*])
  N = N_ELEMENTS(x)

; check the sampling rate is regular

  if KEYWORD_SET(dt) then begin
      print,'-- Even sampling rate assumed = ',dt
  endif else begin
      dt = t[1:(N-1)]-t[0:(N-2)]
      dt_frac = STDDEV(dt) / MEAN(dt)
      dt = MEAN(dt)
      if (dt_frac gt 1.0E-3) then begin
          print,'-------------------------------------------------'
          print,'** Uneven sampling found in time series'
          print,'** Mean sampling interval = ',dt
          print,'** standard deviation / mean of sampling interval = ',dt_frac
          print,'** Analysis will proceed assuming even sampling at mean rate'
          print,'** This may cause spurious results. Please check carefully.'
          print,'-------------------------------------------------'
      endif else begin
          print,'-- Even sampling rate assumed = ',dt
      endelse
  endelse

; ----------------------------------------------------------
; calculate periodogram 

  per = PERIODOGRAM(x, DT=dt, F=f, DF=df, /RMS)
  nf = N_ELEMENTS(f)

; plot on log scale

  yrange = [MIN(per)/1.5, MAX(per)*2]
  xrange = [MIN(f)/1.5, MAX(f)*1.5]

  PLOT, f, per, PSYM=10, /XLOG, /YLOG, XTITLE="Frequency (Hz)", $
    YTITLE="Power ([rms/mean]!U2!N Hz!U-1!N)", YRANGE=yrange, $
    XRANGE=xrange, /XSTYLE, /YSTYLE

; ----------------------------------------------------------
; Check the frequency range to fit over and define it in dimensionless
; units [0,1,...,N/2-1] 

  frange = INTARR(2)

; if FREQ_RANGE=[f1, f2] not defined then set frange to max/min
; possible values, i.e. [f_1, f_{N/2}] 

  if (N_ELEMENTS(freq_range) ne 2) then begin
      frange = [0, nf-1]
  endif else begin

; check the frequency range is sane (f1 < f2).
; if so, define FRANGE in dimensionless units

      if (freq_range[0] lt freq_range[1]) then begin
          frange[0] = FLOOR(freq_range[0] * dt * N) - 1
          frange[1] = FLOOR(freq_range[1] * dt * N) - 1
      endif else begin
          frange = [0, nf-1]
      endelse
  endelse

; also clip f1,f2 at min(f),max(f) if they overrun

  frange[0] = (frange[0] > 0)
  frange[1] = (frange[1] < (nf-1))

; number of frequencies fitted

  n_fband = frange[1] - frange[0] + 1

; ----------------------------------------------------------
; linear regression of LOG(pow)-LOG(freq)

  mask = WHERE(f GE f[frange[0]] AND f LE f[frange[1]], count)
  err = MAKE_ARRAY(count, VALUE=1.0)
  sigma2 = (!pi/6.0/ALOG(10))^2
  per_err = MAKE_ARRAY(count, VALUE=SQRT(sigma2))
  result = LINFIT(ALOG10(f[mask]), ALOG10(per[mask]), COVAR=covar, $
                  MEASURE_ERROR=per_err)

; define the model

  alpha = result[1]
  norm = 10.0^(result[0] + 0.25068)
  model = norm * f^alpha
  OPLOT, f[mask], model[mask], THICK=3
  print,'-- Log regression index     = ',alpha,' +/-',SQRT(covar[1,1])
  print,'-- Log regression log[norm] = ',ALOG10(norm),' +/-',SQRT(covar[0,0])

; calculate uncertainty on the model

  err2 = covar[1,1] * ALOG10(f)^2 + covar[0,0] + 2.0*covar[0,1] * ALOG10(f)
  err = SQRT(err2)
  mod_up = model * 10.0^(err)
  mod_do = model / 10.0^(err)
  OPLOT, f[mask], mod_up[mask], LINESTYLE=0, COLOR=130, THICK=3
  OPLOT, f[mask], mod_do[mask], LINESTYLE=0, COLOR=130, THICK=3

; ----------------------------------------------------------
; find the 3-sigma confidence level above the best fitting model 
; using no model uncertainty correction

  p0 = 1-GAUSS_PDF(3.0)
  pN = p0 / count
  conf_factor = CHISQR_CVF(pN, 2)/2.0

  OPLOT, f[mask], model[mask]*conf_factor, LINESTYLE=2

; ----------------------------------------------------------
; check the fit

  resid = 2.0 * per[mask] / model[mask]
  sum_rat = TOTAL(resid)
  dof = 2*count
  p = CHISQR_PDF(sum_rat, dof)
  print,'-- Goodness-of-fit test p = ',p

; ----------------------------------------------------------
  RETURN

END
