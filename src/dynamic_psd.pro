FUNCTION DYNAMIC_PSD, X, N_SEG, ERROR=error,  $ 
                      DT=dt, T=t, F=f, DF=df, RMS=rms, MEAN=mean, $
                      POIS=pois, LEAHY=leahy, P_N=P_N, NOSUB=nosub, $
                      ERR2=err2

; ----------------------------------------------------------
;+
; NAME:
;       DYNAMIC_PSD
;
; PURPOSE:
;       Calculate the dynamic (time-resolved) power spectrum.
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       img = DYNAMIC_PSD(x,2048)
;
; INPUTS:
;       x       - (vector) evenly sampled time series 
;       n_seg   - (int) number of points in each segment 
;
; OPTIONAL INPUTS:
;       dt     - (float) the sampling rate (default = 1.0)
;       rms    - (logical) whether to use rms/mean normalisation
;       pois   - (logical) subtract Poisson noise level?
;       leahy  - (logical) use Leahy normalisation
;       error  - (vector) errors on the time series
;       nosub  - (logical) do not subtract noise level? (TRUE)
;
; OUTPUTS:
;       dper   - 2-d array listing power at [frequency,time]
;
; OPTIONAL OUTPUTS:
;       df     - frequency resolution
;       f      - array of frequencies 
;       t      - array of times
;       mean   - array of mean(x) in each segment
;       P_N    - array of Poisson noise levels for each periodogram
;       err2   - array of mean square errors
;
; DETAILS:
;       Calculates the dynamic periodogram of a univariate time series
;       from the modulus square of the Fourier transform.
;       The time series is split into M segments of length N_seg
;       and the periodogram for each is computed.
;       The output is a 2-d array listing power density at
;       frequencies f_j = j / (N_seg*dt), and times 
;       t_i= i * (N_seg*dt), where dt is the sampling interval.
;       NB: Only non-zero positive frequencies are returned.

;       If dt is supplied and the rms keyword
;       is set them the 'fractional' normalisation 2*dT / N*mean(x)^2
;       is used. This yields power density in fractional units of
;       [(rms/mean)^2 Hz^-1], popular in X-ray astronomy.
;
;       If LEAHY keyword set then use normalisation from 
;       Leahy et al. (1983; ApJ, 266, 160) which is defined
;       so that, if the input time series is nothing but Poisson
;       noise, the output periodogram has a mean power of 2
;       and follows a chi_2^2 distribution. The RMS keyword
;       overrules the LEAHY keyword.
;
;       If keyword POIS is set then expected Poisson noise level
;       is subtracted from the powers. In this case the data are 
;       assumed to be in count/sec units and the Poisson noise level 
;       is expected to be 2<x> for absolute normalisation or 2/<x> 
;       for fractional rms normalisation. 
;
;       If ERROR is supplied then the noise level will be assumed to
;       be 2.dT.<err2> for absolute rms normalisation, where
;       <err2> is the mean square error (i.e. mean[ERROR^2]). 
;
;       If NOSUB keyword it TRUE then we do not actually subtract the
;       Poisson noise level, just return the expected levels as P_N.
;
;
; EXAMPLE USAGE:
;        the following shows how to create a simple
;        dynamic (time - frequency) periodogram.         
;
;          y = randomn(seed,8192L)           
;          p = dynamic_psd(y,64,f=f,t=t)        
;          shade_surf,p,t,f,xtitle="Time",ytitle="Frequency", $
;               ztitle="Power",charsize=2.0
;
;        The following uses DYNAMIC_PSD to calculate
;        periodograms from several segments of one longer
;        time series, then averages them together to make
;        a single periodogram with errors.
;
;          y = randomn(seed,32768L)           
;          p = dynamic_psd(y,256,f=f)        
;          m = (size(p))[1]                   
;          mean_p = total(p,1)/m             
;          err = mean_p/sqrt(m)                
;          plot,f,mean_p, xrange=[0,max(f)*1.02],yrange=[0,4], $
;               /ystyle, /xstyle, xtitle="Frequency",ytitle="Power",psym=10
;          plot_err,f,mean_p,err
;          oplot,f,rebin([2.0],m),linestyle=2
;
; HISTORY:
;       11/01/2007 -  v1.0  - first working version
;       17/01/2007 -  v1.1  - added POIS keyword
;       26/04/2007 -  v1.2  - added check for n_seg integer  
;                              fixed bug for POIS keyword check
;       03/12/2007 -  v1.3  - added LEAHY keyword
;       06/07/2010 -  v1.4  - added ERROR keyword and new P_N
;                              calculation
;       05/10/2010 -  v1.5  - added catches for when the number of
;                              good intervals, m, is small (0 or 1)
;
; NOTE:
;       If memory space is an issue try 
;       freeing arrays once they are no longer needed.
;
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)

  COMPILE_OPT idl2

; watch out for errors

  on_error, 2

; ----------------------------------------------------------
; Check the arguments and input

; is the time series data array well-defined?

  n = N_ELEMENTS(x)
  if (n eq 0) then begin
      PRINT, '** No data in DYNAMIC_PSD.'
      RETURN, -1
  endif

; is the segment length supplied?

  if (N_ELEMENTS(n_seg) eq 0) then begin
      PRINT, '** Must supply segment length to DYNAMIC_PSD.'
      RETURN, -1
  endif

; check n_seg is an integer

  if ((SIZE(n_seg))[1] ne 2) then n_seg=floor(n_seg)

; and check it is large enough to work

  if (n_seg lt 4) then begin
      PRINT, '** N_SEG too small in DYNAMIC_PSD'
      RETURN, -1
  endif

; is the sampling period dt supplied?

  if (N_ELEMENTS(dt) eq 0) then dt = 1.0

; ----------------------------------------------------------
; Break the time series into M segments of N_seg length each

  m = n/n_seg

  if (m eq 0) then begin
      print,"** No complete intervals in DYNAMIC_PSD.", m, n, n_seg
      RETURN, -1
  endif

  n_cut = m*n_seg

  data = REFORM(x[0:n_cut-1], n_seg, m)

; array of times for each segment

  t = FINDGEN(m) * (n_seg*dt) + 0.5*dt

; no. +ve frequencies

  nf = n_seg/2                      

; frequency resolution

  df = 1.0/(dt*n_seg)

; frequency array

  f = (FINDGEN(nf)+1) * df 

; Nyquist frequency

  f_n = 0.5/dt

; ----------------------------------------------------------
; Calculate the dynamic periodogram using FFT command

; calculate mean value (DC component) for each segment

  mean = TOTAL(data, 1)/FLOAT(n_seg)

; if ERROR supplied then calculate mean square error for each segments

  if KEYWORD_SET(error) then begin
      data_err = REFORM(error[0:n_cut-1], n_seg, m)
      err2 = TOTAL(data_err^2, 1)/FLOAT(n_seg)
  endif

; subtract mean (DC component)

  data = TEMPORARY(data) - TRANSPOSE(REBIN(mean, m, n_seg, /sample))

; Calculate the Fourier transform of each segment (row)

  dft = FFT(data, 1, dimension=1, /double)

; extract only positive (non-zero) frequency components
; and square to get power

  dper = ABS(dft[1:nf,*])
  dper = dper * dper

; apply the absolute normalisation

  norm = (2.0 * dt / FLOAT(n_seg))
  dper = TEMPORARY(dper)*norm

; calculate the Poisson noise component if requested
; use either 2.dT.<err2> or 2.<x> formulation depending on inputs 

  if KEYWORD_SET(error) then begin
      P_i = 2.0 * dt * err2 
  endif else begin
      if KEYWORD_SET(pois) then begin
          P_i = 2.0 * mean
      endif
  endelse

  if (N_ELEMENTS(P_i) gt 0) then begin
;      if (m gt 1) then begin
          P_N = TRANSPOSE(REBIN(REFORM(P_i), m, nf, /sample))
;      endif else begin
;          P_N = P_i
;      endelse
  endif

; subtract off the Poisson noise level unless NOSUB=TRUE

  if (NOT KEYWORD_SET(nosub) and N_ELEMENTS(P_N) gt 0) then $
    dper = TEMPORARY(dper) -  P_N 


; apply the fractional rms normalisation if requested
; also re-normalise P_N (Poisson noise level)

;  if (m gt 1) then begin
      mean_i = TRANSPOSE(REBIN(REFORM(mean), m, nf, /sample))
;    endif else begin
;      mean_i = mean[0]
;  endelse
  
  if KEYWORD_SET(rms) then begin
      dper = TEMPORARY(dper) / mean_i^2
      P_N = TEMPORARY(P_N) / mean_i^2
  endif else begin
      if KEYWORD_SET(leahy) then dper = TEMPORARY(dper) / mean_i
      if KEYWORD_SET(leahy) then P_N = TEMPORARY(P_N) / mean_i
  endelse

; ----------------------------------------------------------
; return the data array to the user
; after rotating them to M*N order

  if (N_ELEMENTS(P_N) gt 0) then begin
      P_N = TRANSPOSE(P_N)
  endif else begin
      P_N = [0]
  endelse

  dper = TRANSPOSE(dper)
  
  RETURN, dper

END
