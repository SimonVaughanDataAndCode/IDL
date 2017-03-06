FUNCTION CLIP_RANGE, x, percent=percent

; ----------------------------------------------------------
;+
; NAME:
;       CLIP_RANGE
;
; PURPOSE:
;       Find the lower and upper per cent bounds on input data
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       img_range = CLIP_RANGE(img)
;
; INPUTS:
;       x - data to be examined (floating array)
;
; OPTIONAL INPUTS:
;       percent - (scalar) threshold in per cent (default is 5.0%)
;
; OUTPUTS:
;       two element vector giving the lower and upper
;       boundaries.
;
; DETAILS:
;       Called with
;
;       IDL> range = CLIP_RANGE(data,2.0)
;
;       this function will return the values of the data at
;       the 2nd and 98th percentile.
;
; HISTORY:
;       Based on a routine IMCLIP by Liam E. Gumley.
;       15/01/2007  -- v1.0 -- working version
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments and input

; is the data array well-defined?
  n = n_elements(x)
  if (n eq 0) then begin
      print,'** No data in CLIP_RANGE.'
      return,0
  endif

; is the percentile supplied? if not use default.
  if (n_elements(percent) eq 0) then percent = 5.0

; ----------------------------------------------------------

; find min and max of input data
  min_val = min(x,max=max_val)

; calculate the density distribution
  nbins=100
  if (percent le 1 && n ge 1000) then nbins=1000
  hist = histogram(float(x), nbins=nbins, locations=bin)

; calculate cumulative, normalised sum
  sum = total(hist,/cumulative) 
  sum = sum * (100.0/float(n))

; find the percentile ranges
  range = [min_val,max_val]
  index = where((sum ge percent) and (sum le (100.0-percent)),count) 
  if (count ge 2) then  range = [bin[index[0]],bin[index[count-1]]]

; ----------------------------------------------------------
; Return the data array to the user

  return,range

END
