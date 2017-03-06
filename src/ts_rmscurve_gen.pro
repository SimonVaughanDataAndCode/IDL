PRO TS_RMSCURVE_GEN, dt=dt, n_seg=n_seg, n_sim=n_sim, $
                     simlen=simlen, mfreq=mfreq, mpow=mpow, seed=seed, $
                     spile=spline, double=double, freq=freq, $
                     p_n=p_n, pois=pois, var=var

; ----------------------------------------------------------
;+
; NAME:
;       TS_RMSCURVE_GEN
;
; PURPOSE:
;       Calculate average flux and rms time series for artificial data
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       p = TS_RMSCURVE_GEN(dt=1.0,n_seg=1024,mpow=mpow,mfreq=mfreq,freq=freq)
;
; INPUTS:
;       filename  - (string) name of file listing FITS time series
;       n_seg     - (integer) number of points in each rms segment
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
;       pois      - (logical) subtract Poisson noise level?
;       freq      - (vector) list of two frequencies defining
;                   the range over which to calculate the rms
;       var       - (logical) output variance rather than RMS?
;
; OUTPUTS:
;       <files>   - series of ASCII files 
;       rmscurve_gen.list - ASCII file listing the output filenames
;
; OPTIONAL OUTPUTS:
;       p_n       - noise level of the variance periodogram
;
; DETAILS:
;       Calculate mean flux and rms in segments of length n_seg
;       for several time series generated from an input power spectrum.
;       The input power spectrum is specified by MPOW and MFREQ,
;       and the sampling interval and length of each simulated
;       time series are given by DT and SIMLEN. The number of 
;       simulations is given by N_SIM.
;
;       Each series is broken into segments that are contiguous
;       and from each is computed a time series of the mean
;       flux and the rms. The mean is calculated in the
;       standard fashion, and the rms is calculated by integrating
;       the periodogram over the frequency range freq=[f1,f2].
;       If the time (dt) is in units of sec, then f1,f2 should
;       be in units of sec^-1 = Hz.
;
;       The periodograms may be optionally corrected for 
;       Poisson noise, in which case it is assumed the time series are 
;       in units of ct/s and dt is in units of s.
;
;       The time series generation uses TS_GEN and the
;       average periodogram works as does TS_RMSCURVE
;
; PROCEDURES CALLED:
;       TS_GEN, TS_RMSFLUX, WRITE_TABLE, STATUSLINE
;
; HISTORY:
;       16/05/07  - v1.0 - first working version
;       17/05/07  - v1.1 - used STATUSLINE for progress report
;
; NOTES:
;
;-
; ----------------------------------------------------------

; watch out for errors
;  on_error,2

; ----------------------------------------------------------
; Check the arguments

; if SIMLEN not defined, set default
  if (n_elements(simlen) eq 0) then simlen = 65536

; make sure SIMLEN is even
  if ((simlen mod 2) ne 0) then begin
      print,'** Please make SIMLEN even in TS_RMSCURVE_GEN'
      return
  endif

; check N_SIM is defined, or assume default
  if (n_elements(n_sim) eq 0) then n_sim = 128

; check the shape of the input array
  nf = n_elements(mfreq)
  np = n_elements(mpow)
  if (nf ne np) then begin
      print,'** MFREQ and MPOW of differing sizes in TS_RMSCURVE_GEN.'
      return
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

; if freq=[f1,f2] not defined then set to zero
; (signals to subsequent routines that it is not set)
  if (n_elements(freq) ne 2) then begin
      freq=[0,0]
  endif 

; ----------------------------------------------------------
; Main part of procedure

; Make array to store output filenames
  outlist = strarr(n_sim)

; initialise noise level
  p_n = 0.0

; Loop over each FITS file to process
  for i=0,n_sim-1 do begin

      string = '-- Processing time series '+strtrim(strupcase(i+1),2)
      string = string+' of '+strtrim(strupcase(n_sim),2)
      statusline,string

; Generate artificial time series
      data = ts_gen(simlen,dt=dt,freq=mfreq,pow=mpow,seed=seed, $
             spline=keyword_set(spline), double=keyword_set(double))

; make exponential transform
      data = exp(data)

; Compute time series for mean and rms
      rms = ts_rmsflux(data,dt=dt,n_seg=n_seg,freq=freq,mean=mean,t=t, $
                       pois=keyword_set(pois), var=keyword_set(var),pow=pow,df=df)

      if (i eq 0) then begin
          p = pow
      endif else begin
          p = temporary(p) + pow
      endelse

; Generate a filename for the output
      obsid = 'rmscurve_gen'
      flag_segm = strtrim(strupcase(i+1),2)
      filename = obsid+'_'+flag_segm+'.rms'

; put time (t), mean flux (mean) and rms (rms) in data array
      data = transpose([[t],[mean],[rms]])
      if (n_elements(t) eq 1) then begin
          data = [t[0],mean[0],rms[0]]
          data = reform(data,3,1)
      endif

; Save the output as an ASCII file
      write_table,data,filename,/double
      outlist[i] = filename

; End loop over simulation
  endfor

; blank line
  print

; normalise to get average power spectrum
  p = p / float(n_sim)

; variance-of-variance term : V[sigma] = \sum p(f)^2 df^2
  f = (indgen(n_elements(p))+1)*df
  mask = where((f ge freq[0]) and (f le freq[1]))
  E_sigma = total(p[mask]) * df
  V_sigma = total(p[mask]^2) * df^2 

; convert from variance to expected power density
; (formula is different depending on whether VAR or RMS is measured)
  if keyword_set(var) then begin
      p_n = V_sigma * (2.0*n_seg*dt)
  endif else begin
      p_n = 0.25 * V_sigma / E_sigma * (2.0*n_seg*dt)
  endelse

; output the list of filenames 
  write_table,outlist[0:n_sim-1],'rmscurve_gen.list'

; ----------------------------------------------------------
; Return the data array to the user

  return

END
