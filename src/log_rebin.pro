FUNCTION LOG_REBIN, x, y, dy=dy, $
          binfactor=binfactor, minbin=minbin, varerr=varerr, $
          sort=sort, bin_x=bin_x, bin_l=bin_l, bin_u=bin_u, $
          bin_dy=bin_dy, bin_n=bin_n, binmap=binmap

; ----------------------------------------------------------
;+
; NAME:
;       LOG_REBIN
;
; PURPOSE:
;      Rebin one dimensional data into equal logarithmic bins in X
;
; AUTHOR:
;      Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;      result = LOG_REBIN(x, y, binfactor=1.2)
;
; INPUTS:
;      x - vector containing abcissae.
;      y - vector containing ordinates.
;
; OPTIONAL INPUTS:
;      dy        - (float vector) errors on ordinates
;      binfactor - (float) fractional binwidth (default=1.2)
;      minbin    - (integer) minimum no. data points per bin
;      varerr    - (logical) whether to calculate error based
;                  on sqrt(variance) within bin (default=no)
;      sort      - (logical) if true, sort data into acending
;                  X order
;
; OUTPUTS:
;      bin_y - binned ordinates
;
; OPTIONAL OUTPUTS:
;      bin_x   - binned abcissae
;      bin_l   - lower edge of abcissae bin
;      bin_u   - upper edge of abcissae bin
;      bin_dy  - error on ordinate bin
;      bin_n   - no. points contained in bin
;      binmap  - array mapping input data points to output bins
;
; DETAILS:
;      Bin the X,Y using logarithmically spaced bins
;      in X, i.e. each bin spans X -> X*binfactor.
;      The error on the bin is computed in various ways.
;      If errors on Y as supplied (dy) then the combined
;      error is given. If errors are not supplied they
;      are assumed to be unity. If the /varerr keyword
;      is set the error is calculated by the RMS scatter
;      within the bin (i.e. error = sqrt[variance]). 
;      The bin locations (along X) are returned as the
;      geometric mean X position or as the lower/upper
;      edges of each contiguous bin.
;
;      The input data arrays are assumed to be in
;      sorted into increasing order of X. If not, use
;      the SORT optional keyword.
;
; HISTORY:
;      17/01/2007  -- v1.0 -- first version
;      09/02/2007  -- v2.0 -- rewritten as a 1D loop
;      01/05/2007  -- v2.1 -- use long-integer counter for loop
;      27/05/2007  -- v2.2 -- added BINMAP output
;      07/08/2007  -- v2.3 -- changed BIN_N to long integer and
;                             added SORT option
;      05/06/2009  -- v2.4 -- bug fix. binu[mask] wasn't filtered
;      08/06/2009  -- v2.5 -- bug fix. inside the crucial double loop
;                             there was a blunder. now fixed line is
;                             bin_l[bin] = bin_u[bin-1]. Without this
;                             it was not always true that bin_l[bin] =
;                             bin_u[bin-1]. 
; 
; Doesn't yet deal with identical x values...
;
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments and input

; is the data array well-defined?
  n = n_elements(x)
  if (n lt 2) then begin
      print,'** Not enough data in LOG_REBIN.'
      return,0
  endif

; is the Y array empty, if so, replicate X
  if (n_elements(y) eq 0) then y=x

; are the X, Y arrays the same size? 
  if (n_elements(y) ne n) then begin
      print,'** X and Y lengths differ in LOG_REBIN.'
      return,0
  endif

; are the errors on Y supplied?
  if (n_elements(dy) ne n) then dy=make_array(n, /float, value=1.0)

; is the fractional bin size supplied
  if (n_elements(binfactor) eq 0) then binfactor=1.2

; is the fractional bin size sensible?
  if ((binfactor) le 1.0) then begin
      print,'** Rebinning factor is less than 1 in LOG_REBIN.'
      return,0
  endif

; is the minimum number of points per bin supplied?
  if (n_elements(minbin) eq 0) then minbin=1

; make sure it's an integer
  minbin=round(minbin)

; ----------------------------------------------------------
; determine number and location of bins

; check that x > 0
  if (min(x) le 0.0) then begin
      print,'-- Negative or zero data in X in LOG_REBIN.'
      return,0
  endif

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

; convert X axis to log_10 units
  logx = alog10(xx)

; convert bin width dX to log_10 units
  lbinfactor = alog10(binfactor)

; find no. bins between min-max
  min_x = min(logx,max=max_x)
  nbins = ceil( (max_x-min_x)/lbinfactor )

; make arrays for binned data
  bin_x = make_array(nbins,/float,value=0.0)
  bin_y = make_array(nbins,/float,value=0.0)
  bin_dy = make_array(nbins,/float,value=0.0)
  bin_l = make_array(nbins,/float,value=0.0)
  bin_u = make_array(nbins,/float,value=0.0)
  bin_n = make_array(nbins,/long,value=0)
  binmap = make_array(n,/long,value=0)

; loop over every data point
  bin=0
  bin_l[bin] = logx[0]
  bin_u[bin] = bin_l[bin] + lbinfactor
  for i=0L,n-1L do begin
      if (logx[i] gt bin_u[bin]) then begin        ; reached upper edge of bin?
          if (bin_n[bin] ge minbin) then begin     ; enough points in bin?
              bin = bin + 1                        ; if so, move onto the next bin
              if (bin gt nbins-1) then break
              bin_l[bin] = bin_u[bin-1]
              bin_u[bin] = bin_l[bin] + lbinfactor
          endif else begin
              bin_u[bin] = logx[i]                 ; if not, expand upper edge 
          endelse
      endif

;      print,i,bin,bin_l[bin],bin_u[bin],logx[i]

      binmap[i] = bin
      bin_x[bin] = bin_x[bin] + logx[i]
      bin_y[bin] = bin_y[bin] + yy[i]
      bin_n[bin] = bin_n[bin] + 1
      if keyword_set(varerr) then begin
          bin_dy[bin] = bin_dy[bin] + yy[i]*yy[i]
      endif else begin
          bin_dy[bin] = bin_dy[bin] + dyy[i]*dyy[i]
      endelse

  endfor

; remove bins with no data
  mask = where(bin_n gt 0, nbins)
  bin_x = bin_x[mask]
  bin_y = bin_y[mask]
  bin_dy = bin_dy[mask]
  bin_l = bin_l[mask]
  bin_u = bin_u[mask]
  bin_n = bin_n[mask]

; convert sums to averages
  bin_x = bin_x/float(bin_n)
  bin_y = bin_y/float(bin_n)

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
      
; convert x-axis back from log to linear

  bin_x = 10.0^bin_x
  bin_l = 10.0^bin_l
  bin_u = 10.0^bin_u

; ----------------------------------------------------------
; return result 

  return,bin_y

END

