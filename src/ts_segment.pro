FUNCTION TS_SEGMENT, x, dx=dx, Nseg=Nseg, minseg=minseg

; ----------------------------------------------------------
;+
; NAME:
;       TS_SEGMENT
;
; PURPOSE:
;       Define contiguous segments of an input series
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       seglist = TS_SEGMENT(x,dx=1.0)
;
; INPUTS:
;       x         - (vector) list of time points
;
; OPTIONAL INPUTS:
;       dx        - (scalar) Sampling period (default x[1]-x[0])       
;       minseg    - (scalar) minumum length of acceptable segment
;
; OUTPUTS:
;       seglist   - [n_seg,2] array listing first and last element
;                   numbers of n_seg contiguous segments
;
; OPTIONAL OUTPUTS:
;       Nseg      - number of segments
;
; DETAILS:
;       Takes as input a vector x and scans through marking
;       the first and last elements of contiguous segments.
;       x[i] and x[i+1] are said to be non-contiguous if
;       x[i+1]-x[i] > dx, where dx is the sampling time.
;       The first and last elements of the ith segment -
;       [a_i,b_i] - are collected in a 2*Nseg array that
;       forms the output. Nseg is the number of segments.
;       Any segment shorter than minseg is ignored.
;
;       It is assumed that the input values are in increases
;       order.
;
; Example call:
; IDL> seglist = TS_SEGMENT(x, dx=0.5, Nseg=Nseg, minseg=50)
;
; HISTORY:
;       25/04/07  - v1.0 - first working version
;       26/04/07  - v1.1 - made variable names consistent (t -> x)
;       01/05/07  - v1.2 - added Nseg=0 if input array too small
;
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments

; check the shape of the input array
  s=size(x)
  if (s[0] gt 1) then print,'** Array x has ',s[0],' dimensions'

; do we have enough data to make this worthwhile?
  n=s[1]
  if (n le 4) then begin
      print,'** Array x too small in TS_SEGMENT'
      Nseg=0
      return,0
  endif

; if minseg not defined, set default
  if (n_elements(minseg) eq 0) then minseg=2

; if sampling period dx not defined, set to 1.0
  if (n_elements(dx) eq 0) then dx=x[1]-x[0]

; ----------------------------------------------------------
; Main part of procedure

; define numerical tolerance
  tol = 1.0e-5

; find locations where x[i+1]-x[i] > sampling interval
  xstep = (1.0+tol)*dx
  gaps = where((x[1:n-1]-x[0:n-2] gt xstep),count)

; total number of contiguous segments
  Nseg = count + 1

; define seglist array to store start,stop locations
  seglist = [[0,gaps+1],[gaps,n-1]]

; check for segments with fewer than MINSEG elements
  seglength = seglist[*,1]-seglist[*,0]
  segmask = where(seglength ge minseg,count)
  Nseg = count

; ignore short segments by removing their locations from seglist
  if (Nseg ge 1) then begin
      seglist = seglist[segmask,*]
  endif else begin
      seglist = 0
  endelse

; ----------------------------------------------------------
; Return the data array to the user

  return,seglist

END
