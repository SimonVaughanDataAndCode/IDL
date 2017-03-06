FUNCTION PERIODOGRAM, x, dt=dt, f=f, df=df, rms=rms, leahy=leahy

; ----------------------------------------------------------
;+
; NAME:
;       PERIODOGRAM
;
; PURPOSE:
;       Calculate periodogram (mod-square Fourier transform)
;       of a univariate time series (x)
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       per = PERIODOGRAM(x)
;
; INPUTS:
;       x      - (array) evenly sample time series
;
; OPTIONAL INPUTS:  
;       dt     - (float) the sampling rate (default = 1.0)
;       rms    - (logical) whether to use rms/mean normalisation
;       leahy  - (logical) use Leahy normalisation
;
; OUTPUTS:
;       per    - 1-d array listing frequencies
;
; OPTIONAL OUTPUTS:
;       df     - frequency resolution 
;       f      - array of frequencies 
;
; DETAILS:
;       Calculates the periodogram of a univariate time series
;       from the modulus square of the Fourier transform using
;       the FFT command. The output is a 1-d array listing power 
;       density at frequencies f_j = j / N.dt, where N is the
;       number of data points in the time series, dt is the sample
;       interval of the time series and j=1,2,...,N/2.
;       NB: Only non-zero positive frequencies are returned.
;       There are three possible normalisations. The simplest
;       (default) is 2/N. If dt is supplied then instead the
;       normalisation is 2*dT/N which yields power density in 
;       absolute units [X^2 Hz^-1], where X is the units of
;       the time series, x. 
;
;       If dt is supplied and the rms keyowrd
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
; EXAMPLE USAGE:
;       IDL> dt = 0.1
;       IDL> t = indgen(1024)*dt
;       IDL> x = 0.4*sin(4.0*t) + randomn(seed,1024)
;       IDL> p = periodogram(x,dt=dt,/rms,f=f)
;       IDL> plot,f,p,xtitle="Frequency (Hz)", $
;                   ytitle="Power ([rms/mean]!U2!N Hz!U-1!N)"
;
;       will return the periodogram in fractional rms^2 units.
;
; HISTORY:
;       11/01/2007 -  v1.0  - first working version
;       03/12/2007 -  v1.1  - added LEAHY keyword
;
; NOTES:
;
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)
  COMPILE_OPT idl2

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments and input

; is the time series data array well-defined?
  n = n_elements(x)
  if (n lt 4) then begin
      print,'** Not enough data in PERIODOGRAM.'
      return,0
  endif

; is the sampling period dt supplied?
  if (n_elements(dt) eq 0) then dt=1.0

; make sure input array is 1 dimensional
  s=size(x)
  if (s[0] gt 1) then x=reform(x,n,/overwrite)

; ----------------------------------------------------------
; Calculate the periodogram

; no. +ve frequencies
  nf = n/2                      

; frequency resolution
  df = 1.0/(dt*n)

; frequency array
  f = (findgen(nf)+1) * df 

; calculate mean (DC component)
  mean_x = mean(x)

; Calculate the Fourier transform of mean-subtracted data
  dft = fft(x-mean_x, 1, dimension=1)  ; (NB: 1/N factor when for forward DFT)

; extract only positive (non-zero) frequency components
; and square to get power

  per = abs(dft[1:nf])^2

; apply the appropriate normalisation

  norm = (2.0 * dt / float(n))
  if keyword_set(rms) then norm=norm/(mean_x^2)

  per = temporary(per)*norm

; ----------------------------------------------------------
; Return the data array to the user

  return,per

END
