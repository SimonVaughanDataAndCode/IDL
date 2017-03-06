FUNCTION CDYNAMIC_PSD, x, n_seg, step, $ 
          dt=dt, t=t, f=f, df=df, rms=rms, meanx=meanx, pois=pois

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
;       x     - (float) evenly sampled time series 
;       n_seg - (int) number of points in each segment 
;
; OPTIONAL INPUTS:
;       dt   - (float) the sampling rate (default = 1.0)
;       rms  - (logical) whether to use rms/mean normalisation
;       pois - (logical) subtract Poisson noise level?
;
; OUTPUTS:
;       dper - 2-d array listing power at [frequency,time]
;
; OPTIONAL OUTPUTS:
;       df   - frequency resolution
;       f    - array of frequencies 
;       t    - array of times
;       mean - array of mean(x) in each segment
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
;
;       Example, array x sampled every 0.0625 sec 
;
;       IDL> img = DYNAMIC_PSD(x,2048,dt=0.0625,t=t_img,f=f_img,/rms)
;
;       this returns an image array (img) showing the periodogram
;       calculated every 128 sec, in fractional rms^2 units.
;       The time and frequencies scales are in t_img, f_img.
;       If keyword POIS is set then expected Poisson noise level
;       is subtracted from the powers. In this case the data are 
;       assumed to be in count/sec units and the Poisson noise level 
;       is expected to be 2<x> for abs normalisation or 2/<x> 
;       for fractional rms normalisation
;
; HISTORY:
;       11/01/2007 -  v1.0  - first working version
;       17/01/2007 -  v1.1  - added POIS keyword
;
; NOTE:
;       If memory space is an issue try using TEMPORARY()
;       and freeing arrays once they are no longer needed.
;
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments and input

; is the time series data array well-defined?
  n = n_elements(x)
  if (n eq 0) then begin
      print,'** No data in DYNAMIC_PSD.'
      return,0
  endif

; is the segment length supplied?
  if (n_elements(n_seg) eq 0) then begin
      print,'** Must supply segment length to DYNAMIC_PSD.'
      return,0
  endif

; is the sampling period dt supplied?
  if (n_elements(dt) eq 0) then dt=1.0

; ----------------------------------------------------------
; Determine segment size and periodogram frequencies
  m = n/n_seg * step

; no. +ve frequencies
  nf = n_seg/2                      

; frequency resolution
  df = 1.0/(dt*n_seg)

; frequency array
  f = (findgen(nf)+1) * df 

; ----------------------------------------------------------
; Calculate the dynamic periodogram using FFT command

  t = fltarr(m-step+1)
  cdper = fltarr(nf,m-step+1)

  for i=0,m-step do begin

; Break the time series into segment of length N_seg 
      n1 = i*n_seg/step
      n2 = n1 + n_seg - 1
      data = x[n1:n2]
      t[i] = (n1*dt)

; subtract mean value (DC component)
      meanx = mean(data)
      data = data - meanx

; Calculate the Fourier transform of each segment (row)
      dft=fft(data,1,dimension=1)

; extract only positive (non-zero) frequency components
; and square to get power
      dper = abs(dft[1:nf])^2

; apply the absolute normalisation
      norm = (2.0 * dt / float(n_seg))
      dper=dper*norm

; subtract the Poisson noise component if requested
      if keyword_set(pois) then $
        dper = dper - 2.0*meanx

; apply the fractional rms normalisation if requested
      if keyword_set(rms) then $
        dper = dper / meanx^2

; transfer to output array
      cdper[*,i] = dper
  
  end

; ----------------------------------------------------------
; Return the data array to the user

  return,transpose(cdper)

END
