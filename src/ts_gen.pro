FUNCTION TS_GEN, n, dt=dt, freq=freq, pow=pow, seed=seed, time=time, $
                 spline=spline, double=double, phase=phase, coh=coh, $
                 log=log, y=y

; ----------------------------------------------------------
;+
; NAME:
;       TS_GEN
;
; PURPOSE:
;       Generate a random time series from a power spectrum model 
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       x = TS_GEN(65536, dt=0.1)
;
; INPUTS:
;       n         - (scalar) length of time series (default = 65536) 
;
; OPTIONAL INPUTS:
;       dt        - (scalar) Sampling period (default = 1.0)
;       freq      - (vector) frequencies at which spectrum is known
;       pow       - (vector) spectral density at frequencies FREQ
;       phase     - (vector) phase shift at frequencies FREQ
;                   (default=0)
;       coh       - (vector) coherence of two signals (default = 0)
;       seed      - (long integer) seed for random number generator
;       spline    - (logical) use cubic spline interpolation
;       log       - (logical) interpolate linearly in log-log space
;       double    - (logical) perform FFT in double prec.
;
; OUTPUTS:
;       x         - (vector) output time series
;
; OPTIONAL OUTPUTS:
;       time      - (vector) sampling times [0, n-1]*dt
;       y         - (vector) output time series #2
;
; DETAILS:
;       Generate an evenly-sampled time series for a noise (random)
;       process with a power spectrum specified by POW and FREQ.  
;
;       The method comes from:
;         Davies R. B., Harte D. S., 1987, Biometrika, v74, pp95-101
;       and was introduced to astronomy by
;         Timmer J., Konig M., 1995, A&A, v300, pp707-710
;
;       The time series is generated using the following algorithm: 
;       
;        1. The "true" power spectrum is specified using POW(FREQ)
;      
;        2. Define the Fourier frequencies for the output time series
;        as f_j = j/(N*dT) with j=1,2,...,N/2. Use interpolation
;        to find the power spectrum at f_j from input POW(FREQ)
;
;        3. The "true" power spectrum is converted from power
;        (non-negative) to a "true" DFT for the process, using the
;        fact that POW = |DFT|^2, so we have a complex-valued
;        DFT = complex(sqrt(POW),sqrt(POW)) at each frequency f_j
;
;        4. Draw two sets of N/2 normal deviates (random numbers from
;        "normal" Gaussian distribution.
;
;        5. Multiply the real and imaginary parts of the DFT by the
;        deviates. This randomised the "true" DFT and gives it the
;        distribution expected for an observed or "sampled" DFT from a
;        single time series of a random process.
;         X(f_j) = DFT(f_j) * eps_j   where eps_j is a normal deviate
;
;        6. Use the inverse FT to convert from the frequency domain to
;        the time domain, i.e. from x(t_i) = FFT[X(f_j)] 
;
;        7. Fill-in the array of times t_i = i*dT for i=0,...,N-1
;
;       The randomisation step (5) is equivalent to drawing the square
;       amplitude of the DFT from a chi-squared distribution (with two
;       degrees of freedom), and the phase of the DFT from a uniform
;       distribution over the range [0,2*pi]. These are the expected
;       sampling distributions from a random time series.
;
;       Note that in reality the DFT is also defined for negative
;       Fourier frequencies j=-N/2,...,-1. In order for the resulting
;       time series to be real we require that the X(f_j) = X'(-f_j),
;       so the negative frequencies carry the complex conjugate of the
;       positive frequencies. Each side of the DFT is normalised by
;       1/2 so that the sum over all (-ve and +ve) frequencies is
;       equal to the total variace (the integral of the power
;       spectrum). 

;       Also, the DFT at the Nyquist frequency j=N/2 is always real
;       when N is even, so the imaginary part is set to zero.  
;       The DFT at zero frequency (j=0) determines the mean (DC
;       component) of the resulting time series. Here we generate
;       zero-mean data, so this is set to zero, i.e. X(f_j = 0) = 0.
;
;       The spectrum is specified by the vectors FREQ and POW, which
;       are interpolated as needed to populate the periodogram needed
;       for the generation (step 2). Interpolation is linear unless SPLINE
;       keyword is set (in which case it is cubic spline). If FREQ and
;       POW are not specified, the spectrum is assumed to be flat
;       (i.e. POW = constant). 
;
;       WARNING: The routine needs to know the power density at
;       frequencies f_j = j/(N*dT) with j=1,2,...,N/2. You need
;       to make sure your input spectrum spans this full range,
;       i.e. that MIN(FREQ) <= 1/(N*dT) and MAX(FREQ) >= 1/2dT.
;       This may involve simply adding more extra point to the
;       input power spectrum at a very low or high frequency.
;       If this is the case the program will return a WARNING
;       but carry on by extrapolating the data outside the
;       range of the input data.
;       
;       If the input power spectrum is a power law it may be best
;       to use the LOG keyword. This forces the interpolation to
;       be done using the log of the power spectrum. I.e. it
;       interpolates log(pow) - log(freq) data, and then converts
;       the result back to linear-space. In effect, it interpolates
;       between points using a power law model.       
; 
;       As the routine uses the FFT function, it works fastest if N is
;       a power of 2 (default = 2^16) 
;
;       There are addition optional inputs PHASE and COH. This allows 
;       two time series (x and y) to be generated with an optional
;       phase shift and coherence applied between them. Since the phase is randomly and
;       uniformly distibuted over the range [0,2*pi] this makes no
;       difference for a single time series. NOTE: coherence is not
;       currently in use, output is fully coherent.

;       It is possible to generate multiple phase lagged two time
;       series by calling the routine again using the same 
;       random number seed but different POW or PHASE values. The 
;       result will be time series that differ only in their power
;       spectrum (modulus square of DFT), phase (argument of DFT) and coherence. 
;       If X = A*exp(i*theta) then applying a phase shift phi we get
;       X' = A*exp(i*[theta + phi]) = X * exp(i*phi).
;
; EXAMPLE USAGE:
;
;   Generate time series with 1/f spectrum
;
;     IDL> freq = (INDGEN(512)+1)/1024.0
;     IDL> pow = freq^(-1)
;     IDL> x = TS_GEN(1024, dt=1.0, freq=freq, pow=pow,time=time, $
;                     seed=seed)
;     IDL> plot,time,x
;
;   Generate time series with 1/f spectrum  making use of LOG keyword
;
;     IDL> freq = [1e-6,100]
;     IDL> pow = 0.01*freq^(-1)
;     IDL> x = TS_GEN(1024, dt=1.0, freq=freq, pow=pow,time=time, $
;                     seed=seed,/log)
;     IDL> plot,time,x
;
;   Because the spectrum is a power law, we only need define two
;   end points at let the interpolation (in log-log space) do the rest.
;   (NB: try this without the LOG keyword to see how it goes wrong!)
;
;   Generate two time series with constant phase delay of pi/2 using a
;   1/f^2 spectrum 
;
;     IDL> freq = (INDGEN(512)+1)/1024.0
;     IDL> pow = freq^(-2)
;     IDL> s = 123L
;     IDL> x = TS_GEN(1024, dt=1.0, freq=freq, pow=pow,time=time, $
;                     seed=s,phase=0)
;     IDL> plot,time,x
;     IDL> phase = !pi/2
;     IDL> s = 123L
;     IDL> x = TS_GEN(1024, dt=1.0, freq=freq, pow=pow,time=time,$
;                     seed=s,phase=phase) 
;     IDL> oplot,time,x,color=170
;
; NB: A constant time delay of tau can be produced using 
; phase(j) = 2.0*!pi*tau*freq(j) 
;
; HISTORY:
;       14/05/07  - v1.0 - first working version
;       15/05/07  - v1.1 - bug fix: INDGEN now uses L64 keyword
;                           this is needed for N > 2^15
;       20/12/07  - v1.2 - added PHASE keyword
;       15/01/09  - v1.3 - added LOG keyword
;       19/01/09  - v.14 - added check that range of FREQ
;                          spans [f_min = 1/NdT, f_max = 1/2dT]
;       22/09/10  - v1.5 - added clauses to allow for integer DT
;                           values
;       28/11/11  - v2.0 - added COH keyword to simulate incoherent
;                           signals. Not yet implimented!
;
; NOTES:
;       + uses built in random number generator
;       + use COH to simulate partially coherent bivariate outputs
;
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)

  COMPILE_OPT idl2

; watch out for errors

  on_error, 2

; ----------------------------------------------------------
; Check the arguments

; if N not defined, set default

  if (N_ELEMENTS(n) eq 0) then n = 65536

; make sure N is even

  if ((n mod 2) ne 0) then begin
      PRINT,'** Please make N even in TS_GEN'
      RETURN, -1
  endif

; if sampling period dt not defined, set to 1.0

  if (N_ELEMENTS(dt) eq 0) then dt = 1.0

; check the shape of the input array

  nf = N_ELEMENTS(freq)
  np = N_ELEMENTS(pow)
  if (nf ne np) then begin
      PRINT,'** FREQ and POW of differing sizes in TS_GEN.'
      RETURN, -1
  endif
  if (nf lt 2) then begin
      PRINT,'** FREQ too small in TS_GEN.'
      RETURN, -1
  endif

; if FREQ is not defined, set-up default (flat) spectrum

  if (nf eq 0) then begin
      freq = [0.0,0.5/dt]
      pow = [1.0,1.0]
  endif

; if PHASE is not defined, set-up default (zero) phase shift

  np = N_ELEMENTS(phase)
  if (np ne nf and np gt 1) then begin
      PRINT,'** FREQ and PHASE of differing sizes in TS_GEN.'
      RETURN, -1        
  endif
  if (np eq 0) then phi = MAKE_ARRAY(nf,/DOUBLE,value=0.0D)
  if (np eq 1) then phi = MAKE_ARRAY(nf,/DOUBLE,value=phase[0])
  if (np eq nf) then phi = phase

; check that PHI is within range [0,2*pi]

  pi2 = 2.0*!pi
  phi = phi mod pi2

; if COH is not defined, set-up default (unit coherence)

  nc = N_ELEMENTS(coh)
  if (nc ne nf and nc gt 1) then begin
      PRINT,'** FREQ and COH of differing sizes in TS_GEN.'
      RETURN, -1        
  endif
  if (nc eq 0) then co = MAKE_ARRAY(nf, /DOUBLE, value=1.0D)
  if (nc eq 1) then co = MAKE_ARRAY(nf, /DOUBLE, value=coh[0])
  if (nc eq nf) then co = coh

; force CO to be within range [0, 1]

  co = (ABS(co) < 1)

; ----------------------------------------------------------
; check the range of input frequencies spans the range
; needed to generate the simulation

  f_min = 1.0/(n*dt)
  f_max = 1.0/(2.0*dt)

  if (MIN(freq) gt f_min) then begin
      PRINT,"-- WARNING. MIN(FREQ) > f_min in TS_GEN."
      PRINT,"-- MIN(FREQ) = ",MIN(freq)," f_min = ",f_min
      PRINT,"-- Data will be extrapolated. You may prefer to EXPand the range of FREQ."
  endif
  if (MAX(freq) lt f_max) then begin
     PRINT,"** MAX(FREQ) < f_max in TS_GEN. EXPand range of FREQ"
     PRINT,"-- MIN(FREQ) = ",MIN(freq)," f_min = ",f_min
     PRINT,"-- Data will be extrapolated. You may prefer to EXPand the range of FREQ."
  endif

; ----------------------------------------------------------
; Main part of procedure

; number of positive frequencies in Fourier Transform

  nf = n/2

; make array for Fourier Transform 
; (need room for positive and negative frequencies)
 
  x = MAKE_ARRAY(2*nf, /COMPLEX)
  y = MAKE_ARRAY(2*nf, /COMPLEX)

; make array for frequencies

  f = INDGEN(nf+1, /DOUBLE) / (n*dt)

; interpolation of input power spectrum, to fill in any gaps
; interpolate given spectrum POW(FREQ) onto required frequencies F

  if KEYWORD_SET(log) then begin

; convert to log-log space, interpolate there, and convert back

      lpow = ALOG(pow)
      lfreq = ALOG(freq)
      lf = ALOG(f[1:nf])
      lspec = INTERPOL(lpow, lfreq, lf, SPLINE=KEYWORD_SET(spline))
      spec = EXP(lspec)
      spec = [0.0, spec]
      lpow = 0
      lfreq = 0
      lf = 0

  endif else begin

; or just interpolate in lin-lin space as default

  spec = INTERPOL(pow, freq, f, SPLINE=KEYWORD_SET(spline))

  endelse

; set DC value ( psd(f=0) ) to zero

  spec[0] = 0.0

; interpolate phase shift spectrum (PHI)

  phi = INTERPOL(TEMPORARY(phi), freq, f, SPLINE=KEYWORD_SET(spline))
  phi = phi mod pi2

; interpolate coherence spectrum (CO)

  co = INTERPOL(TEMPORARY(co), freq, f, SPLINE=KEYWORD_SET(spline))
  co = (ABS(co) < 1)

; normalise spectrum

  spec = TEMPORARY(spec) * nf / (2.0 * dt)

; ----------------------------------------------------------
; now for the Davies-Harte algorithm

; take square root of spectrum

  spec = SQRT(TEMPORARY(spec))

; populate Fourier Transform of time series with SQRT(spectrum)
; multiplied by normal deviate 
; (independent for real and complex components)

; first positive frequencies

  x[1:nf-1] = COMPLEX(spec[1:nf-1] * RANDOMN(seed, nf-1), $ 
                      spec[1:nf-1] * RANDOMN(seed, nf-1))

; apply phase shift Y(f) = X(f) * EXP(i*phi)
  
  y[1:nf-1] = x[1:nf-1] * EXP(COMPLEX(0.0, phi[1:nf-1]))

; FT must be real at Nyquist frequency

  x[nf] = COMPLEX(spec[nf] * RANDOMN(seed, 1), 0.0)
  y[nf] = COMPLEX(spec[nf] * RANDOMN(seed, 1), 0.0)

; make sure FT at negative frequencies is conjugate: X(-f) = X*(f)

  x[nf+1:2*nf-1] = REVERSE(CONJ(x[1:nf-1]))
  y[nf+1:2*nf-1] = REVERSE(CONJ(y[1:nf-1]))

; then inverse Fourier Transform into time domain: X(f) -> x(t)

  x = FFT(x, -1, DOUBLE=KEYWORD_SET(double), /OVERWRITE) 
  y = FFT(y, -1, DOUBLE=KEYWORD_SET(double), /OVERWRITE) 

; drop imaginary part (which is zero)

  x = REAL_PART(TEMPORARY(x))
  y = REAL_PART(TEMPORARY(y))

; calculate TIME if needed

  time = INDGEN(n) * dt

; ----------------------------------------------------------
; Return the data array to the user

  RETURN, x

END
