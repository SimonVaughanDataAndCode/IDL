
FUNCTION rms_mean_fit, x, DT=dt, ERROR=error, NBIN=nbin, N_SEG=n_seg, $
                       SEG_LIST=seg_list, FREQ_RANGE=freq_range, QPS=qps, $
                       COVAR=covar, QPLOT=qplot, POIS=pois, $
                       PSD_I=psd_i, ERR_I=err_i, FREQ_i=freq_i, NOISE_I=noise_i, $
                       MEAN_I=mean_i, SILENT=silent, BINFACTOR=binfactor, TITLE=title, $
                       NUM_SEG=Num_seg, RMS=rms, ERMS_I=erms_i, ERMSERR_I=ermserr_i

; ----------------------------------------------------------
;+
; NAME:
;       rms_mean_fit
;
; PURPOSE:
;       Calculate, plot and fit the rms-mean relation given input time series
;
; AUTHOR:
;       Simon Vaughan (U.Leicester) 
;
; CALLING SEQUENCE:
;       result = RMS_MEAN_FIT(x, ...)
;
; INPUTS:
;       X          - (vector) one or more time series (concatinated)
;
; OPTIONAL INPUTS:  
;       DT         - (float) time bin size
;       ERROR      - (vector) errors on the time series (concatinated)
;       SEG_LIST    - (array) [n_seg,2] array listing first and last
;                     element numbers of the combined time series
;       NBIN       - (integer) Bin the periodograms into NBIN spectra
;       N_SEG      - (integer) Length of time series segments to use
;       FREQ_RANGE - (array) Lower, upper frequencies to inegrate over
;       QPS        - (logical) Make a PostScript plot?
;       COVAR      - (array) Covariance matrix of the model fit
;       QPLOT      - (logical) Plot the rms-mean relation?
;       POIS       - (logical) subtract off Poisson noise?
;       SILENT     - (logical) suppress on-screen output?
;       BINFACTOR  - (float) fractional binwidth (default=1.2) for
;                            plotting PSDs
;       RMS        - (logical) perform periodogram averaging in
;                            fractional (i.e. rms/mean^2) units?
;
; OUTPUTS:
;       RESULT     - (array) a list of the fit results 
;
; OPTIONAL OUTPUTS:  
;       PSD_I      - (array) [NBIN, N_F] PSDs for each of NBIN flux bins
;                    at each of N_F frequencies
;       ERR_I      - (array) error bars on the flux-binned PSDs
;       ERMS_I      - (vector) rms estimates for each mean_i bin
;       ERMSERR_I   - (vector) error on rms estimates
;       FREQ_I     - (vector) list of N_F frequencies of PSDs
;       NOISE_I    - (array) estimated Poisson noise levels per PSD
;       MEAN_I     - (vector) mean flux for each PSD
;       NUM_SEG    - (integer) number of intervals used in analysis
;
; DETAILS:
;       Takes as input one (or more) time series X and computes the rms
;       as a function of mean. The time series is broken into
;       consecutive non-overlapping intervals of length N_SEG. For
;       each of these we calculate the mean and the periodogram using
;       DYNAMIC_PSD. The individual periodograms are then sorted into
;       order of the mean <X>, and averaged together into NBIN
;       bins. These show how the power spectrum varies as a function
;       of the mean level. The errors on the averaged periodograms are
;       computed by combining the variances of the individual
;       chi-square distributed periodogram values.
;
;       The averaged periodograms are then integrated over a specified
;       frequency range FREQ_RANGE to give the variance, the square
;       root of which is the rms. A background noise level can be
;       removed first, if needed. If the keyword POIS is TRUE then the
;       data are assumed to be count rates (in ct/s with DT in
;       units of s), and the Poisson noise level is 2.<X>. If the
;       errors on X are supplied (as ERROR) then this is used to
;       estimate the power density due to measurement errors from the
;       mean square error of each interval - as 2.DT.<ERROR^2>.
;
;       The result is an estimate of the rms, with error, for
;       different average levels <X>. These data are then fitted with
;       a linear function, using LINFIT. The resulting parameters, and
;       the quality of the fit as determined using a chi-square test,
;       are given in the output vector:
;
;         result[0] = slope
;         result[1] = slope error
;         result[2] = y-intercept
;         result[3] = y-intercept error
;         result[4] = x-intercept
;         result[5] = x-intercept error
;         result[6] = fractional x-intercept
;         result[7] = chi-square of the linear fit
;         result[8] = degrees of freedom of the linear fit
;         result[9] = p-value from chi-square test
;        
;       The NBIN means levels output as MEAN_I, the
;       frequencies as FREQ_I, and the averaged periodogram powers as
;       PSD_I. The errors on the averaged peridograms are given in
;       ERR_I, and the estimated noise levels (see above) in NOISE_I. 
;
;       The input time series may be a concatination of several time
;       series of the same variable. If this is the case then the
;       start and stop positions of the individual
;       continuous/contiguous time series within the input array
;       should be specified using the SEG_LIST keyword. For example if
;       X is the result of combining three time series of lengths 1000,
;       500 and 2000, we would have SEG_LIST like this:
;
;       SEG_LIST = [[  0, 1000, 1500],
;                   [999, 1499, 3499]
;
;       Each segment is treated individually - meaning that each one
;       is divided into intervals of length NSEG, and no intervals
;       cross the boudary between segments of the input data - but are
;       then merged prior to binning according to the mean <X>.
;
; EXAMPLE USAGE:
;
;
; HISTORY:
;       12/07/2010 - v1.0 - first working version  
;       14/07/2010 - v1.1 - added BINFACTOR keyword and log-rebin of
;                            PSD plots  
;       16/07/2010 - v1.2 - added quick exit if NBIN=1
;       02/09/2010 - v1.3 - added NUM_SEG output
;       16/09/2010 - v1.4 - added RMS input
;       25/09/2010 - v1.5 - replaced LOG_REBIN with GROUP
;       05/10/2010 - v1.6 - added a catch for when
;                            DYNAMIC_PSD returns no useful data
;       30/10/2010 - v1.7 - changed error calculation. fixed it.
;       10/03/2011 - v1.8 - replaced GROUP with REGROUP
;
; PROCEDURES CALLED:
;       PLOTSYM, PLOT_ERR, DYNAMIC_PSD, PS_OPEN, PS_CLOSE, REGROUP
;
; NOTES:
;-
; ----------------------------------------------------------

; ----------------------------------------------------------
; Check the arguments

; check the shape of the input array

  s = SIZE(x)
  if (s[0] gt 1) then begin
      PRINT, '** Array x has ',s[0],' dimensions in RMS_MEAN_FIT'
      RETURN, -1
  endif

; do we have enough data to make this worthwhile?

  N = N_ELEMENTS(x)
  if (N le 256) then begin
      PRINT, '** Array X is too small in RMS_MEAN_FIT'
      RETURN, -1
  endif

; if number of flux bins not defined, set to 5

  if (N_ELEMENTS(NBIN) eq 0) then nbin=5

; if N_SEG not defined, set default

  if (N_ELEMENTS(n_seg) eq 0) then n_seg=256

; if DT not defined, set default

  if (N_ELEMENTS(dt) eq 0) then dt = 1.0

; if SEG_LIST not defined, set default

  if (N_ELEMENTS(seg_list) lt 2) then begin
      seg_list = MAKE_ARRAY(1, 2, /LONG)    
      seg_list[0,0] = 0
      seg_list[0,1] = N-1
  endif

; if BINFACTOR not defined, set default

 if (N_ELEMENTS(binfactor) eq 0) then binfactor = 1.2

; ----------------------------------------------------------
; Main routine

; Determine how many time series are to be processed

  N = (SIZE(seg_list, /DIMENSIONS))[0]
  PRINT, '-- Time series to process',N
  PRINT, '-- Sampling rate',dt

  first = 0
  for i = 0, N-1 do begin

; load a time series

       y = x[seg_list[i,0]:seg_list[i,1]]
      if (N_ELEMENTS(error) ne 0) then begin
          dy = error[seg_list[i,0]:seg_list[i,1]]
      endif else begin
          dy = MAKE_ARRAY(N_ELEMENTS(y))
      endelse

; compute the average periodogram, and the flux and rms in 
; segments of length N_SEG*DT

       psd = DYNAMIC_PSD(y, n_seg, error=dy, dt=dt, f=f, $
                        df=df, mean=flux, POIS=KEYWORD_SET(pois), $
                        /NOSUB, P_N=P_N, RMS=KEYWORD_SET(rms))

; skip ahead if insufficient data

      print,seg_list[i,0],seg_list[i,1],n,i,N_ELEMENTS(flux),N_ELEMENTS(psd)
      if (N_ELEMENTS(psd) eq 1) then CONTINUE

; concatinate the arrays into one big array
; npsd = noise PSDs, mpsd = main PSDs

      if (first eq 0) then begin
          first = 1
          mpsd = psd
          mm = flux
          npsd = P_N
      endif else begin
          mpsd = [mpsd, psd]
          mm = [mm, flux]
          npsd = [npsd, P_N]
      endelse

  endfor

; ----------------------------------------------------------
; now sort the periodograms into accending flux order

  indx = sort(mm)
  mf = mm[indx]
  psd = mpsd[indx,*]
  nspsd = npsd[indx,*]

; now average the periodograms in flux bins
; and calculate the error by propagating the
; chi-square variances of each periodogram ordinate

  nf = N_ELEMENTS(f)
  N = N_ELEMENTS(mf)
  npb = N/nbin
  PRINT, '-- Number of intervals used',N
  PRINT, '-- Summed duration of intervals',N*N_seg*dt

  psd_i = MAKE_ARRAY(nf, nbin)
  err_i = MAKE_ARRAY(nf, nbin)
  noise_i = MAKE_ARRAY(nf, nbin)
  mean_i = MAKE_ARRAY(nbin)

  for i = 0, nbin-1 do begin

      j = i*npb
      k = (i+1) * npb - 1
      dpsd = psd[j:k, *]
      noise = nspsd[j:k, *]
;      print,i,j,k,N,nbin,npb
;      print,size(dpsd),size(psd_i)
      psd_i[*,i] = TOTAL(dpsd, 1) 
      mean_i[i] = MEAN(mf[j:k])
      noise_i[*,i] = TOTAL(noise, 1) 
 
  endfor

  psd_i = psd_i/npb
  err_i = psd_i/SQRT(npb > 1)
  noise_i = noise_i/npb

; ----------------------------------------------------------
; Check the frequency range and define it in dimensionless
; units [0,1,...,N/2-1]

  frange = INTARR(2)

; if FREQ_RANGE=[f1,f2] not defined then set frange to max/min
; possible values, i.e. [f_1, f_{N/2}] 

  if (N_ELEMENTS(freq_range) ne 2) then begin
      frange = [0, nf-1]
  endif else begin

; check the frequency range is sane (f1 < f2).
; if so, define FRANGE in dimensionless units

      if (freq_range[0] lt freq_range[1]) then begin
          frange[0] = FLOOR(freq_range[0] * dt * n_seg) - 1
          frange[1] = FLOOR(freq_range[1] * dt * n_seg) - 1
      endif else begin
          frange = [0, nf-1]
      endelse
  endelse

; also clip f1,f2 at min(f),max(f) if they overrun

  frange[0] = (frange[0] > 0)
  frange[1] = (frange[1] < (nf-1))

; number of frequencies integrated over

  n_fband = frange[1] - frange[0] + 1

; ----------------------------------------------------------
; subtract the noise levels

  psd_i = psd_i - noise_i

; ----------------------------------------------------------
; loop over the flux-binned periodograms and calculate
; the rms by integrating over the frequency range freq_range[0]-freq_range[1]
; Also, plot the periodograms

  rms = MAKE_ARRAY(nbin)
  err = MAKE_ARRAY(nbin)
  pow = MAKE_ARRAY(nbin)
  dpow = MAKE_ARRAY(nbin)
  col = 200 - INDGEN(nbin)*150.0/nbin

  for i = 0, nbin-1 do begin

      yrange = [min(ABS(psd_i))*0.9, max(psd_i)*1.1]
      xrange = [min(f)*0.9, max(f)*1.1]

      bin_y = REGROUP(psd_i[*,i], f, dy=err_i[*, i], binwidth=-binfactor, $
                        bin_x=bin_x, bin_dy=bin_dy, bin_l=bin_l, bin_u=bin_u)

      if (i eq 0) then begin
          PLOT, bin_x, bin_y, /xlog, /ylog, yrange=yrange, xrange=xrange, psym=10, $
            /xstyle, /ystyle, xtitle="Frequency (Hz)", ytitle="Power ([ct s!U-1!N]!U2!N Hz!U-1!N)", /NODATA
      endif 

      for j = 0, N_ELEMENTS(bin_x)-1 do begin
          clip = [xrange[0],yrange[0],xrange[1],yrange[1]]
          plots, [bin_l[j],bin_u[j]], [bin_y[j],bin_y[j]], CLIP=clip, NOCLIP=0, /DATA, COLOR=col[i]
      endfor
      for j = 0, N_ELEMENTS(bin_x)-2 do begin
          plots, [bin_u[j],bin_u[j]], [bin_y[j],bin_y[j+1]], CLIP=clip, NOCLIP=0, /DATA, COLOR=col[i]
      endfor

;      OPLOT, bin_x, bin_y, psym=10, COLOR=col[i]
      PLOT_ERR, bin_x, bin_y, (bin_dy < bin_y*0.99999), COLOR=col[i]

  
; sum the powers, and their squares

      pow[i]  = TOTAL( psd_i[frange[0]:frange[1],i] ) * df
      pow2 = TOTAL( err_i[frange[0]:frange[1],i]^2 ) * (df*df)

; form the error on the summed power (by adding in quadrature)

      dpow[i] = SQRT(pow2)

; rms is sqrt(power)

      rms[i] = SQRT(MAX([pow[i], 0.0]))

; error on rms, i.e. error on sqrt(power)

      err[i] = 0.5 * dpow[i]/SQRT(pow[i])
      
  endfor

; ----------------------------------------------------------
; remove any dodgy (i.e. -ve power) data

  mask = WHERE(rms gt 0.0, count)
  if (count gt 0) then begin
      rms = rms[mask]
      err = err[mask]
      mean_i = mean_i[mask]
  endif

; ----------------------------------------------------------
; if only one flux bin, finish now...

  if (nbin eq 1) then begin
      freq_i = f
      RETURN, -1
  endif

; ----------------------------------------------------------
; fit rms-flux with linear function

  result = LINFIT(mean_i, rms, MEASURE_ERRORS=err, chisq=chisq, $
                  prob=prob, sigma=sigma, yfit=yfit, covar=covar)

  dof = N_ELEMENTS(rms) - 2
  off = -result[0]/result[1]
  doff = ABS(off) * sqrt( (sigma[0]/result[0])^2 + (sigma[1]/result[1])^2 )

  if NOT KEYWORD_SET(silent) then begin

      PRINT, "Fit results:"
      PRINT, "Slope:      ", result[1], " +/-", sigma[1]
      PRINT, "y Intercept:", result[0], " +/-", sigma[0]
      PRINT, "x Intercept:", off, " +/-", doff
      PRINT, "x Intercept:", off/MEAN(mean_i), " (fractional)"
      PRINT, "Chi-square =", chisq 
      PRINT, "dof        =", dof
      PRINT, "p-value    =", prob

  endif

  outp = [result[1], sigma[1], result[0], sigma[0], off, doff, $
              off/MEAN(mean_i), chisq, dof, prob]
      

; plot RMS-FLUX relation

  xrange = [0, MAX(mean_i)*1.1]
  yrange = [0, MAX(rms+err)*1.1]

  if KEYWORD_SET(qplot) then begin
      if KEYWORD_SET(qps) then ps_open, THICK=5
      
      PLOTSYM, 0, 1, /fill
      PLOT, mean_i, rms, psym=8, xtitle="Flux (ct s!U-1!N)", ytitle="rms (ct s!U-1!N)", $
        xrange=xrange, yrange=yrange, /xstyle, /ystyle, title=title, POSITION=[0.15,0.15,0.95,0.95]
      PLOT_ERR, mean_i, rms, err
      
      OPLOT, mean_i, yfit, COLOR=170
      OPLOT, mean_i, rms, psym=8
      if KEYWORD_SET(qps) then ps_close
  endif
      

; ----------------------------------------------------------
; Return the data array to the user

  Num_seg = N
  freq_i = f
  erms_i = rms
  ermserr_i = err
  RETURN, outp

; ----------------------------------------------------------
END
