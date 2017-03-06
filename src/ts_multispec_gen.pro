FUNCTION TS_MULTISPEC_GEN, dt=dt, n_seg=n_seg, n_sim=n_sim, $
          simlen=simlen, mfreq=mfreq, mpow=mpow, seed=seed, $
          spile=spline, double=double, rebin=rebin, $
          rms=rms, dp=dp, freq=freq, f0=f0, f1=f1
          
; ----------------------------------------------------------
;+
; NAME:
;       TS_MULTISPEC_GEN
;
; PURPOSE:
;       Calculate average periodogram for artificial data
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       p = TS_MULTISPEC_GEN(dt=1.0,n_seg=1024,mpow=mpow,mfreq=mfreq,freq=freq)
;
; INPUTS:
;
; OPTIONAL INPUTS:
;  For simulated data:
;       dt        - (scalar) Sampling period (default = 1.0)
;       n_seg     - (integer) number of points in each segment
;       n_sim     - (integer) number of simulations to perform
;       simlen    - (integer) length of each simulation (>= n_seg) 
;       mfreq     - (vector) frequencies at which spectrum is known
;       mpow      - (vector) spectral density at frequencies FREQ
;       seed      - (long integer) seed for random number generator
;       spline    - (logical) use cubic spline interpolation
;       double    - (logical) perform FFT in double prec.
;  For output spectrum:
;       rebin     - (scalar) frequency rebinning factor
;       rms       - (logical) use [rms/mean]^2 norm for spectrum?
;
; OUTPUTS:
;       p         - averaged power density spectrum
;
; OPTIONAL OUTPUTS:
;       dp        - (vector) error bar on averaged spectrum
;       freq      - (vector) frequencies (centre of bin)
;       f0        - (vector) lower bound of frequency bin
;       f1        - (vector) upper bound of frequency bin
;
; DETAILS:
;       Estimate the 'observed' power spectrum of an ensamble of 
;       time series generated from an input 'model' power spectrum.
;       The input power spectrum is specified by MPOW and MFREQ,
;       and the sampling interval and length of each simulated
;       time series are given by DT and SIMLEN. The number of 
;       simulations is given by N_SIM.
;
;       For each generated time series a segment of length N_SEG
;       is extracted and a periodogram is 
;       computed. These are ensamble averaged, and then 
;       rebinned in frequency to give a consistent estimate
;       of the power spectrum (with error bars).
; 
;       The amount of rebinning is controlled with the REBIN
;       parameter. If a positive integer, we rebin with REBIN
;       points per bin. If a negative floating point number, we
;       rebin logarimically every factor of (1-REBIN). For example
;       REBIN = -0.1 will bin over ranges f -> 1.1f.
;       If REBIN=1 then no frequency binning is applied (only 
;       ensamble averaging). The default value is REBIN = 50.
;
;       The error bars on the power densities are computed from
;       the chi-square distribution of powers, i.e. error is
;       1/sqrt{M} where M is the number of periodogram ordinates
;       averaged in a given frequency bin (either ensamble or
;       frequency averaging). In the limit of large M (>50)
;       the averaged power densities are normally distributed.
;
;       The time series generation uses TS_GEN and the
;       average periodogram works as does TS_MULTISPEC
;
;
; PROCEDURES CALLED:
;       TS_GEN, DYNAMIC_PSD, LIN_REBIN, LOG_REBIN
;
; HISTORY:
;       15/05/07  - v1.0 - first working version
;
; NOTES:
;       + Add deadtime correction (to DYNAMIC_PSD)
;
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments

; if SIMLEN not defined, set default
  if (n_elements(simlen) eq 0) then simlen = 65536

; make sure SIMLEN is even
  if ((simlen mod 2) ne 0) then begin
      print,'** Please make SIMLEN even in TS_MULTISPEC_GEN'
      return,0
  endif

; check N_SIM is defined, or assume default
  if (n_elements(n_sim) eq 0) then n_sim = 128

; check the shape of the input array
  nf = n_elements(mfreq)
  np = n_elements(mpow)
  if (nf ne np) then begin
      print,'** MFREQ and MPOW of differing sizes in TS_MULTISPEC_GEN.'
      return,0
  endif

; if MFREQ is not defined, set-up default (flat) spectrum
  if (nf eq 0) then begin
      mfreq = [0.0,0.5/dt]
      mpow = [1.0,1.0]
  endif

; if sampling period DT not defined, set to 1.0
  if (n_elements(dt) eq 0) then dt = 1.0

; if N_SEG not defined, set default
  if (n_elements(n_seg) eq 0) then n_seg=128

; if REBIN not defined, set default
  if (n_elements(rebin) eq 0) then rebin=50

; ----------------------------------------------------------
; Main part of procedure

  ndata = 0

; Loop over each FITS file to process
  for i=0,n_sim-1 do begin

      print,'-- Processing time series',i+1

; Generate artificial time series
      data = ts_gen(simlen,dt=dt,freq=mfreq,pow=mpow,seed=seed, $
             spline=keyword_set(spline), double=keyword_set(double))

; make exponential transform
      data = exp(data)

; Compute set of periodograms
      dpsd = DYNAMIC_PSD(data,n_seg,dt=dt,t=t,f=f,df=df,mean=ave, $
                         pois=keyword_set(pois),rms=keyword_set(rms))

; Add all periodograms together
      if (i eq 0) then begin
          pow = total(dpsd,1)
      endif else begin
          pow = temporary(pow) + total(dpsd,1)
      endelse

; calculate number of periodograms (NDATA)
      ndata = ndata + (size(dpsd))[1]

; End loop over simulationd
  endfor

; normalise by number of data segments (to get ensamble average power)
  pow = temporary(pow) / float(ndata)

; if REBIN is positive then use linear frequency bins
  if (rebin gt 1) then begin
      binwidth = floor(rebin)*df
      bin_p = lin_rebin(f, pow, binwidth=binwidth, bin_x=freq, $ 
                        bin_l=f0, bin_u=f1, bin_n=bin_n)
  endif

; if REBIN is negative use logarithmic binning
  if (rebin lt 0) then begin
      bin_p = log_rebin(f, pow, binfactor=1.0-rebin, bin_x=freq, $ 
                        bin_l=f0, bin_u=f1, bin_n=bin_n)
  endif

; otherwise do not rebin
  if ((rebin ge 0) and (rebin le 1)) then begin
      bin_n = intarr(n_elements(f)) + 1
      bin_p = pow
      freq = f
      f0 = f - 0.5*df
      f1 = f + 0.5*df
  endif

; calculate total number of periodogram ordinates in each frequency
; bin
  bin_n = bin_n * ndata

; calculate error on average power in each frequency bin
  dp = bin_p / sqrt(bin_n)

; ----------------------------------------------------------
; Return the data array to the user

  return,bin_p

END
