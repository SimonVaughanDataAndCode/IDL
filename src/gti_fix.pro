FUNCTION gti_fix, t, x, gti, DT=dt, FRAC=frac, MINF=minf, POIS=pois, $
                  INTP_FRAC=intp_frac, SEED=seed, INTP_NUM=intp_num, $
                  WHICH_RESCALE=which_rescale, WHICH_INTP=which_intp

; ----------------------------------------------------------
;+
; NAME:
;       GTI_FIX
;
; PURPOSE:
;       Correct for gaps (GTI) in EPIC time series by interpolation
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       xfix = GTI_FIX(t, x, gti)
;
; INPUTS:
;       t      - (vector) times of the data
;       x      - (vector) values of time series
;       gti    - (array) list of start/stop times for GTIs
;
; OPTIONAL INPUTS:
;       dt     - (float) sampling interval
;       minf   - (float) minimum fraction exposure to allow
;       pois   - (logical) add Poisson noise to interpolated data
;       seed   - (integer) optional seed for random number generation
;
; OUTPUTS:
;       xfix   - (vector) the corrected time series
;
; OPTIONAL OUTPUTS:
;       frac          - (vector) the fractional exposure per bin
;       intp_frac     - (float) the fraction of interpolated data points
;       intp_num      - (integer) the number of interpolated data points
;       which_rescale - (vector) list of array elements that were rescaled
;       which_intp    - (vector) list of array elements that were interpolated
;
; DETAILS:
;       Take as input a time series (t, x) of binned EPIC event data
;       and a table listing the start/stop times of the Good Time
;       Intervals (GTIs). Use the GTIs to define Bad Time Intervals
;       (BTIs) where no data were taken. Correct points in the binned
;       time series for BTIs by either: 
;        (i) correcting for the less than 100% exposure in the time bin 
;         or 
;        (ii) interpolating between neightbouring (good) data points if
;         there is too little good exposure in a bin. 
;       The keyword MINF specifies the mimimum acceptable fractional exposure 
;       per bin - bins with less than this get the interpolation treatment.
;
;       The keyword POIS specifies that Poisson noise should be
;       'added' to the interpolated data to randomise it. This assumes
;       the input array x is in units of counts/bin.
;
;       This function is used by the MAKE_EPIC_LC function.
;
; EXAMPLE USAGE:
;       yfix = GTI_FIX(t, y, gti)
;
; PROCEDURES CALLED:
;       POIDEV
;
; HISTORY:
;       13/05/2010 - v1.0 - first working version
;       17/06/2010 - v1.1 - fixed bug: check for -1 
;                            values output from VALUE_LOCATE.
;       27/10/2010 - v1.2 - added INTP_FRAC output
;       21/12/2011 - v1.3 - change to interpolation step. Now uses
;                            nearest good data with frac > MINF.
;       26/04/2012 - v1.4 - minor fixes: replaced condition "frac gt
;                            minf" in MASK computation to "frac ge
;                            minf". Added SEED keyword. Added INTP_NUM
;                            output. Changed the 'inner' BTI loop to 
;                            subtract exposure lost within BTIs rather
;                            than add non-missing exposure. This allows 
;                            for multiple short BTIs within a single time bin.
;       23/07/2012 - v1.5 - Added XFIX = DOUBLE(x) to ensure output and
;                            intermediate calculations are in double prec.
;       29/01/2014 - v1.6 - Return !NULL if error encountered.
; NOTES:
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)

  COMPILE_OPT idl2

; watch out for errors

  ON_ERROR, 0

; ---------------------------------------------------------
; Check arguments

; check the shape of the input array

  s=SIZE(x)
  IF (s[0] gt 1) THEN BEGIN
      PRINT,'** Array x has ',s[0],' dimensions in GTI_FIX'
      RETURN, !NULL
  ENDIF

; do we have enough data to make this worth while?

  n=s[1]
  IF (n le 4) THEN BEGIN
      PRINT,'** Array x too small in GTI_FIX'
      RETURN, !NULL
  ENDIF

; do we have a GTI table to use?

  s = SIZE(gti)
  IF (s[0] ne 2) THEN BEGIN
      PRINT,'** Missing or incorrect GTI table supplied in GTI_FIX'
      RETURN, !NULL
  ENDIF

; how many GTIs in table?

  ngti = s[1]

; check that t and x have same number of elements

  IF (N_ELEMENTS(t) ne n) THEN BEGIN
      PRINT,'** Number of elements in T and X differ in GTI_FIX'
      RETURN, !NULL
  ENDIF

; have we been given dt?

  IF (N_ELEMENTS(dt) ne 1) THEN dt = t[1]-t[0]

; have we been given minf?

  IF (N_ELEMENTS(minf) ne 1) THEN minf = 0.3D

; ---------------------------------------------------------
; Main routine

; exposure in each time bin

  texp = MAKE_ARRAY(n, VALUE=dt, /DOUBLE)

; subtract exposure lost at start of observation

  gti_start = VALUE_LOCATE(t, gti[0,0])
  IF (gti_start ge 0) THEN BEGIN
      texp[0:gti_start] = 0.0
      texp[gti_start] =  dt - (gti[0,0] - t[gti_start])
  ENDIF

; subtract exposure lost at end of observation

  gti_end = VALUE_LOCATE(t+dt, gti[-1,1]) 
  IF (gti_end lt n-1) THEN BEGIN
      gti_end = gti_end + 1
      texp[gti_end:-1] = 0.0
      texp[gti_end] = gti[-1,1] - t[gti_end] 
  ENDIF

; make list of bad time intervals (BTIs)
; (if there are any...)

  IF (ngti gt 1) THEN BEGIN

      nbti = ngti - 1
      bti = MAKE_ARRAY(nbti, 2)
      FOR i = 0, nbti-1 DO BEGIN
          bti[i, 0] = gti[i, 1]
          bti[i, 1] = gti[i+1, 0]
      ENDFOR

; subtract exposure lost during each BTIs

      FOR i = 0, nbti-1 DO BEGIN
          
          bti_start = VALUE_LOCATE(t, bti[i,0])
          bti_end = VALUE_LOCATE(t, bti[i,1])
          bti_start = MAX([0, bti_start])
          IF (bti_end ge 0) THEN BEGIN
              size_of_gap = bti_end - bti_start
              IF (size_of_gap eq 0) THEN BEGIN
                  texp[bti_start] = texp[bti_start] - (bti[i,1]-bti[i,0])
              ENDIF 
              IF (size_of_gap gt 1) THEN BEGIN
                 texp[(bti_start+1):(bti_end-1)] = 0.0D
              ENDIF
              IF (size_of_gap gt 0) THEN BEGIN
                texp[bti_start] = texp[bti_start] - (t[bti_start] + dt - bti[i,0])
                texp[bti_end] = texp[bti_end] - (bti[i,1] - t[bti_end])
              ENDIF
          ENDIF

      ENDFOR

  ENDIF

; fractional exposure in each bin

  frac = DOUBLE(texp/dt)

; ---------------------------------------------------------
; correct for exposure loss

  xfix = DOUBLE(x)
  mask = WHERE(frac lt 1.0 and frac ge minf, count, /NULL)
  IF (count gt 0) THEN xfix[mask] = xfix[mask] / frac[mask]
  which_rescale = mask

; linearly interpolate bins with too little exposure
; between the closest neighbouring 'good' points

  gap = WHERE(frac lt minf, count, COMPLEMENT=full, /NULL)
  which_intp = gap
  IF (count gt 0) THEN BEGIN
      goodx = xfix[full]
      goodt = t[full]
      xfix = INTERPOL(goodx, goodt, t)

; add Poisson noise to interpolated values?

      IF KEYWORD_SET(pois) THEN BEGIN
          xfix[gap] = POIDEV(xfix[gap], SEED=seed)
      ENDIF
  ENDIF

  intp_frac = FLOAT(count) / n 
  intp_num = count
  
; ---------------------------------------------------------
; Return to user

  RETURN, xfix

END
