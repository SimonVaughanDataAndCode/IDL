FUNCTION LIN_REBIN, x, y, dy=dy, $
          binwidth=binwidth, minbin=minbin, varerr=varerr, $
          bin_x=bin_x, bin_l=bin_l, bin_u=bin_u, $
          bin_dy=bin_dy, bin_n=bin_n, hist=hist, binmap=binmap, $
                    sort=sort

; ----------------------------------------------------------
;+
; NAME:
;       LIN_REBIN
;
; PURPOSE:
;      Rebin one dimensional data into equal bins in X
;
; AUTHOR:
;      Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;      result = LIN_REBIN(x,y,binwidth=5.0)
;
; INPUTS:
;      x - vector containing abcissae.
;      y - vector containing ordinates.
;
; OPTIONAL INPUTS:
;      dy        - (vector) errors on ordinates
;      binwidth  - (float) bin width (default=5.0)
;      minbin    - (integer) minimum no. data points per bin
;      varerr    - (logical) whether to calculate error based
;                    on sqrt(variance) within bin (default=no)
;      hist      - (logical) produce histogram of x values
;      sort      - (logical) sort into accending order of x?
;
; OUTPUTS:
;      bin_y     - (vector)binned ordinates
;
; OPTIONAL OUTPUTS:
;      bin_x     - (vector) binned abcissae
;      bin_l     - (vector) lower edge of abcissae bin
;      bin_u     - (vector) upper edge of abcissae bin
;      bin_dy    - (vector) error on ordinate bin
;      bin_n     - (vector) no. points contained in bin
;      binmap    - (vector) array mapping input data points to output bins
;
; DETAILS:
;      Bin the X,Y using evenly spaced bins in X
;      i.e. each bin spans X -> X + binwidth.
;      The error on the bin is computed in various ways.
;      If errors on Y as supplied (dy) then the combined
;      error is given. If errors are not supplied they
;      are assumed to be unity. If the /varerr keyword
;      is set the error is calculated by the RMS scatter
;      within the bin (i.e. error = sqrt[variance]). 
;      The bin locations (along X) are returned as the
;      arithmetic mean X position or as the lower/upper
;      edges of each contiguous bin.
;
;      If BINWIDTH is positive then the bin width is in
;      the same units as X. If BINWIDTH is negative then
;      the bin width is in multiples of dX (where the
;      data are assumed evenly sampled in X). 
;
; HISTORY:
;      Based on LOG_REBIN.PRO
;      24/01/2007  -- v1.0 -- first version
;      01/05/2007  -- v1.1 -- use long integer counter for loop
;                             added check for zero good bins
;      08/05/2007  -- v1.2 -- bug fix: nbins determination fixed
;      09/05/2007  -- v1.3 -- added HIST option for histograms
;      27/05/2007  -- v1.4 -- added BINMAP output
;      01/07/2007  -- v1.5 -- added -ve option for BINWIDTH
;      08/06/2009  -- v1.6 -- bug fix. inside the crucial double loop
;                             there was a blunder. now fixed line is
;                             bin_l[bin] = bin_u[bin-1]. Without this
;                             it was not always true that bin_l[bin] =
;                             bin_u[bin-1]. 
;      01/02/2010  -- v1.7 -- bug fixes. Added SORT keyword upfront
;                             (was previously missing). Proper
;                             treatment of HIST option so that counts
;                             and sqrt(counts) error are returned as
;                             output bin_y and bin_dy.
;      09/05/2013  -- v1.8 -- minor bug fix, added BIN_U[mask] filter
;
; NOTES:
;      + Add lower/upper bin boundaries
;      + Option for variable width bins (using MINBIN)
;      + Option to open new plotting window
;
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)
  COMPILE_OPT idl2, HIDDEN

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments and input

; is the data array well-defined?
  n = n_elements(x)
  if (n lt 2) then begin
      print,'** Not enough data in LIN_REBIN.'
      return,0
  endif

; is y array not supplied, or is HIST keyword specified?
  if ((n_elements(y) eq 0) or (keyword_set(hist))) then begin
      y = intarr(n)
      y[*] = 1
  endif

; are the X, Y arrays the same size? 
  if (n_elements(y) ne n) then begin
      print,'** X and Y lengths differ in LIN_REBIN.'
      return,0
  endif

; are the errors on Y supplied?
  if (n_elements(dy) ne n) then dy=make_array(n, /float, value=1.0)

; is the bin size supplied
  if (n_elements(binwidth) eq 0) then binwidth=5.0

; is the bin size sensible?
  if (binwidth eq 0.0) then begin
      print,'** Rebinning factor 0 in LIN_REBIN.'
      return,0
  endif

; if BINWIDTH < 0 then treat as factor
  if (binwidth lt 0.0) then begin
      dx = (x[1]-x[0])*(-binwidth)
  endif else begin
      dx = binwidth
  endelse

; is the minimum number of points per bin supplied?
  if (n_elements(minbin) eq 0) then minbin=1

; make sure it's an integer
  minbin=round(minbin)

; make sure input arrays are 1 dimensional
  s=size(x)
  if (s[0] gt 1) then x=reform(x,n,/overwrite)
  s=size(y)
  if (s[0] gt 1) then y=reform(y,n,/overwrite)

; ----------------------------------------------------------
; determine number and location of bins

; if SORT is true then sort data into acending X order
  if keyword_set(sort) then begin
      indx = sort(x)
      xx = x[indx]
      yy = y[indx]
      dyy = dy[indx]
  endif else begin
      xx = x
      yy = y
      dyy = dy
  endelse

; find no. bins between min-max
  min_x = min(xx,max=max_x)
  nbins = ceil( (max_x-min_x)/dx )

; make arrays for binned data
  bin_x  = make_array(nbins,/float,value=0.0)
  bin_y  = make_array(nbins,/float,value=0.0)
  bin_dy = make_array(nbins,/float,value=0.0)
  bin_l  = make_array(nbins,/float,value=0.0)
  bin_u  = make_array(nbins,/float,value=0.0)
  bin_n  = make_array(nbins,/long,value=0)
  binmap = make_array(n,/long,value=0)
  
; loop over every data point
  bin=0
  bin_l[bin] = xx[0]
  bin_u[bin] = bin_l[bin] + dx
  for i=0L,n-1L do begin
      if (xx[i] gt bin_u[bin]) then begin          ; reached upper edge of bin?
          if (bin_n[bin] ge minbin) then begin     ; enough points in bin?
              bin = bin + 1                        ; if so, move onto the next bin
              if (bin gt nbins-1) then break
              bin_l[bin] = bin_u[bin-1]
              bin_u[bin] = bin_l[bin] + dx
          endif else begin
              bin_u[bin] = xx[i]                   ; if not, expand upper edge 
          endelse
      endif

      binmap[i] = bin
      bin_x[bin] = bin_x[bin] + xx[i]
      bin_y[bin] = bin_y[bin] + yy[i]
      bin_n[bin] = bin_n[bin] + 1
      if keyword_set(varerr) then begin
          bin_dy[bin] = bin_dy[bin] + yy[i]*yy[i]
      endif else begin
          bin_dy[bin] = bin_dy[bin] + dyy[i]*dyy[i]
      endelse

  endfor

; remove bins with no data
  mask = where(bin_n gt minbin, nbins)
  bin_x = bin_x[mask]
  bin_y = bin_y[mask]
  bin_dy = bin_dy[mask]
  bin_l = bin_l[mask]
  bin_u = bin_u[mask]
  bin_n = bin_n[mask]

; convert sums to averages
  bin_x = bin_x/float(bin_n)
  if (not KEYWORD_SET(hist)) then begin
      bin_y = bin_y/float(bin_n)
  endif

  if keyword_set(varerr) then begin
      bin_dy = bin_dy - float(bin_n)*bin_y*bin_y
      mask = where(bin_n gt 1,count)
      if (count gt 0) then begin
          bin_dy[mask] = bin_dy[mask] / float(bin_n[mask] - 1)
      endif
      if (count lt nbins) then begin
          mask = where(bin_n eq 1)
          bin_dy[mask] = bin_y[mask]*bin_y[mask]
      endif
      bin_dy = sqrt(bin_dy > 0.0) / sqrt(float(bin_n > 1))
  endif else begin
      bin_dy = sqrt(bin_dy) / float(bin_n)
  endelse

  if KEYWORD_SET(hist) then begin
      bin_dy = sqrt(bin_y)
  endif

; ----------------------------------------------------------
; return result 

  return,bin_y

END

