
FUNCTION TS_RMSFLUX, x, dt=dt, n_seg=n_seg, freq=freq, mean=mean, t=t, $
        pois=pois, deadtime=deadtime, Nseg=Nseg, pow=pow, df=df, var=var, $
        dpow=dpow, drms=drms

; ----------------------------------------------------------
;+
; NAME:
;       TS_RMSFLUX
;
; PURPOSE:
;       Calculate mean and rms at even intervals in a series
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       rms = TS_RMSFLUX(x,dx=1.0)
;
; INPUTS:
;       x         - (vector) list of samples, x(t)
;
; OPTIONAL INPUTS:
;       dt        - (scalar) Sampling period (default 1.0)       
;       n_seg     - (integer) number of points in each segment
;       freq      - (vector) list of two frequencies defining
;                   the range over which to calculate the rms
;       pois      - (scalar/vector) subtract Poisson noise level?
;       deadtime  - (logical) correct for detector deadtime?
;       var       - (logical) output variance rather than RMS
;
; OUTPUTS:
;       rms       - (vector) rms values for each segment
;
; OPTIONAL OUTPUTS:
;       mean      - (vector) mean values for each segment
;       Nseg      - (integer) number of segments
;       t         - times of outputs (units of dt)
;       v_err     - power density for the 'noise' in the variance itself
;       df        - frequency resolution of periodogram (p)
;       pow       - average periodogram
;
; DETAILS:
;       Calculate mean flux and rms in segments of length n_seg
;       from an input time series. The mean is calculated in the
;       standard fashion, and the rms is calculated by integrating
;       the periodogram over the frequency range freq=[f1,f2].
;       If the time (dt) is in units of sec, then freq should
;       be in units of sec^-1 = Hz.
;
;       The periodograms may be optionally corrected for 
;       Poisson noise. Setting POIS=-1 will assume x is in
;       units of ct/s and dt is in units of s, and subtract
;       the expected noise level (P_N = 2/<x> in rms normalisation).
;       Setting POIS=z will subtract the specified value P_N = z.
;
;       It is assumed that the input values are in time order
;       in intervals of dt, with no gaps.
;
; Example call:
; IDL> rms = TS_RMSFLUX(x, dt=0.25, n_seg=1024, mean=flux)
;
; PROCEDURES CALLED:
;       DYNAMIC_PSD
;
; HISTORY:
;       26/04/07  - v1.0 - first working version
;       27/04/07  - v1.1 - added test f1>=min(f) and f2<=max(f)
;       01/05/07  - v1.2 - removed bug in FREQ scaling
;                          added keyword_set(pois) check in call to
;                           DYNAMIC_PSD
;       03/05/07  - v1.3 - added DF and P optional outputs
;       16/10/09  - v2.0 - changed POIS input so it can be
;                           a user-defined scalar/vector
;                          Added "error" calculation to give
;                          the uncertainty on the outputs.
;
; NOTES:
;       + Add deadtime correction
;
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments

; if sampling period dt not defined, set to 1.0
  if (n_elements(dt) eq 0) then dt=1.0

; if n_seg not defined, set default
  if (n_elements(n_seg) eq 0) then n_seg=128

; check the shape of the input array
  s=size(x)
  if (s[0] gt 1) then print,'** Array x has ',s[0],' dimensions'

; do we have enough data to make this worthwhile?
  n=s[1]
  if (n le n_seg) then begin
      print,'** Array x too small in TS_RMSFLUX'
      return,0
  endif

; ----------------------------------------------------------
; Check the frequency range and define it in dimensionless
; units [0,1,...,N/2-1]

  frange=intarr(2)

; if freq=[f1,f2] not defined then set frange to max/min
; possible values, i.e. [f_1, f_{N/2}] 
  if (n_elements(freq) ne 2) then begin
      frange=[0,n_seg/2-1]
  endif else begin

; check the frequency range is sane (f1 < f2).
; if so, define FRANGE in dimensionless units
      if (freq[0] lt freq[1]) then begin
          frange[0]=floor(freq[0]*dt*n_seg)-1
          frange[1]=floor(freq[1]*dt*n_seg)-1
      endif else begin
          frange=[0,n_seg/2-1]
      endelse
  endelse

; also clip f1,f2 at min(f),max(f) if they overrun
  frange[0] = (frange[0] > 0)
  frange[1] = (frange[1] < (n_seg/2-1))

; number of frequencies integrated over
  n_fband = frange[1] - frange[0] + 1

; ----------------------------------------------------------
; Main part of procedure

; calculate dynamic power spectrum
; i.e. periodogram for each segment of n_seg points
  dpsd = DYNAMIC_PSD(x,n_seg,dt=dt,t=t,f=f,df=df,mean=mean)

; dimensions of the output array
  m = n_elements(t)
  nf = n_elements(f)

; calculate Poisson noise level
  if keyword_set(pois) then begin

      if (pois eq -1) then begin
          p_n = 2.0*transpose(rebin(mean,m,nf,/sample))
      endif else begin
          p_n = pois
      endelse

  endif else begin
      p_n = 0.0
  endelse

; integrate over frequency range to get variance in each time slice
  v = total( dpsd[*,frange[0]:frange[1]], 2 ) * df
  
; calculate the "error" (rms) on the variance. 
; (This comes from propagating the rms on the individual chi_2^2
; distributed powers.)
  dv = sqrt( total( (dpsd[*,frange[0]:frange[1]])^2, 2 ) ) * df

; subtract Poisson noise level (if requested)
  v = v - p_n * (n_fband * df)

; average over time to get average periodogram (without Poisson subtraction)
  pow = total(dpsd,1)/float(m)

; and calculate the error on the average periodogram
  dpow = sqrt( total((dpsd^2),1) )/float(m)

; check for any negative variances (after Poisson correction)
  mask = where(v lt 0.0,count)
  if (count gt 0) then begin
      print,'** Negative Var[x]-Var[Poisson] in TS_RMSFLUX',count
  endif

; if requested take square root of variance to give rms (zero if -ve variance)
; (and convert from error on V to error on rms)
; otherwise return the variances
  n_rms = n_elements(v)
  mask = where(v gt 0)
  if keyword_set(var) then begin
      rms = v
      drms = dv
  endif else begin
      rms = sqrt(v > 0)
      drms = make_array(n_rms)
      drms[mask] = dv[mask]/(2.0*rms[mask])
  endelse

; ----------------------------------------------------------
; Return the data array to the user

  return, rms

END
