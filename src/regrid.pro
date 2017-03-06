FUNCTION REGRID, X, V, U, SORT=sort

; ----------------------------------------------------------
;+
; NAME:
;      REGRID
;
; PURPOSE:
;      Redistribute one binned array (histogram) onto a different
;      grid, conserving the integral.
;
; AUTHOR:
;      Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;      z = REGRID(inbreak, y, outbreak)
;
; INPUTS:
;      X         - (vector) boundaries of input bins
;      V         - (vector) data value within input bins
;      U         - (vector) boundaries of output bins
;
; OPTIONAL INPUTS:
;      SORT      - (logical) do we need to sort the arrays?
;
; OUTPUTS:
;      Z         - (vector) the output array 
;
; OPTIONAL OUTPUTS:
;      none
;
; DETAILS:
;      Input X is renamed INBREAK,  U is renamed OUTBREAK,
;      and V is renamed Y to ensure they are not alterned upon output.
;      
;      The goal is to map an input array of data Y onto differently 
;      gridded output array. The input array is a vector of values 
;      Y[i] with i = 0...N-1, the bin boundaries are INBREAK which is an
;      array of N+1 elements. INBREAK[i]-INBREAK[i+1] for i=0,...,N-1
;      defines the range of the ith data bin.
;      This is to be mapped onto an array Z, with M elements, whose 
;      bins cover the ranges INBREAK[j]-INBREAK[j+1], with j=0,...,M-1. 
;      The arrays may be unevenly gridded, i.e. 
;      BINWIDTH[i] = INBREAK[i+1] - INBREAK[i] may vary with i.
;
;      Where a single bin of the output Z spans two or more bins of 
;      the input Y, the output is a simple weighted mean of the 
;      input values, with weight W given by the width of the input 
;      bin falling in the range of the output bin.
;      Weighted mean:
;      
;           Z = (sum_i=0,N-1: Y[i]*W[i]) / (sum_i=0,N-1: W[i])
;           
;                i =     5    6   7   8   9 10  11   12  
;      Input array  ...----+-----+-+----+--+--+----+----+--...
;      
;                j =     3       4     5       6     7
;      Output array ...------+------+------+------+-----+--...
;      
;      Here, the output bin 4 overlaps (partially or wholly) 
;      with input bins 6, 7, and 8. The output element Z[4] 
;      is a wieghted average of input data Y[6], Y[7] and Y[8].
;      All other data points Y[i] have weight zero.
;      Y[6] is weighted by the overlap between input bin 6 and
;      output bin 4, i.e. INBREAK[7]-OUTBREAK[4]. Y[7] is wholly
;      inside output bin 4 so has weight equal to its width,
;      i.e. INBREAK[8]-INBREAK[7]. Y[8] is partially inside
;      output bin 4 and has weight OUTBREAK[5]-INBREAK[8].
;      More generally output bin j is a weighted sum where
;      weights are assigned to be
;      
;       W[i] = upper - lower
;       upper = MIN(OUTBREAK[j+1], INBREAK[i+1]) 
;       lower = MAX(OUTBREAK[j], INBREAK[i]) 
;      
;      in cases where there is no overlap upper=lower=W=0.
;      
;      The result result is another histogram, with different
;      bin boundaries, which conserves the integral of the
;      original histogram. E.g. sum(y*dx) = sum(z*du), assuming
;      the output array (u) covers entirely the input array x.
;
;      We do assume both the input and output arrays are in acending 
;      order. If they are not use the SORT keyword.
;
;      Example usage (using PLOT_HIST.PRO):
;
;        x = INDGEN(10)
;        y = RANDOMU(seed,9)
;        nx = N_ELEMENTS(x)
;      
;        u = [-0.5, -0.3, TOTAL(RANDOMU(seed, 13), /CUMULATIVE)]
;        u = u / MAX(u) * MAX(x) + 0.1
;        nu = N_ELEMENTS(u)
;      
;        PLOT, x[0:(nx-2)], y, /NODATA, XRANGE=[-1,12], /XSTYLE
;        PLOT_HIST, y, x0=x[0:(nx-2)], X1=x[1:(nx-1)], /NOPLOT
;        z = REGRID(x, y, u)
;        indx = INDGEN(N_ELEMENTS(u)-1)
;        PLOT_HIST, z, X0=u[indx], X1=u[indx+1], COLOR='ff0000'x, /NOPLOT
;       
;        dx = DOUBLE(x[1:(nx-1)] - x[0:(nx-2)])
;        PRINT, 'Input integral ', TOTAL(y*dx, /DOUBLE)
;        du = DOUBLE(u[1:(nu-1)] - u[0:(nu-2)])
;        PRINT, 'Output integral ', TOTAL(z*du, /DOUBLE)
; 
; EXAMPLE USAGE:
;      z = REGRID(X, V, U, /SORT)
;
; PROCEDURES CALLED:
;      none
;
; HISTORY:
;      19/08/2014 - v0.1 - First working version.
;
; NOTES:
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)

  COMPILE_OPT idl2

; watch out for errors

  ON_ERROR, 0

; ----------------------------------------------------------
; Check arguments

  inbreak = x
  outbreak = u
  y = v

; check number of bin to process

  N = N_ELEMENTS(y)
  N1 = N_ELEMENTS(inbreak)
  M = N_ELEMENTS(outbreak) - 1

  if (M EQ 0 or N EQ 0) then begin
      PRINT,'** REGRID: requires more data for input/output arrays.'
      RETURN, !NULL
  endif

  if (N1 ne N+1) then begin
      PRINT,'** REGRID: INBREAK and Y arrays sizes do not match.'
      RETURN, !NULL
  endif
 
; sort arrays if requested, otherwise assume sorted 
  
  IF KEYWORD_SET(sort) then begin
    indx = SORT(inbreak)
    inbreak = inbreak[sort]
    y = y[sort]
    indx = SORT(outbreak)
    outbreak = outbreak[sort]
  ENDIF

; ----------------------------------------------------------
; pad the input array with zeroes at both ends

  y = [0, y, 0]
  inbreak = [MIN([outbreak, inbreak])-1, inbreak, MAX([outbreak, inbreak])+1]
  N = N_ELEMENTS(y)

; ----------------------------------------------------------
; main routine

  in_upper = INBREAK[1:N]
  in_lower = INBREAK[0:(N-1)]

  z = MAKE_ARRAY(m, /DOUBLE)

  FOR j = 0, M-1 do begin
    upper = (OUTBREAK[j+1] < in_upper) 
    lower = (OUTBREAK[j] > in_lower)
    w = DOUBLE(upper - lower) > 0
    z[j] = TOTAL(w*y, /DOUBLE)
    tw = TOTAL(w, /DOUBLE)
    if (tw gt 0) then z[j] = z[j] / tw
  ENDFOR

; ----------------------------------------------------------
; return to calling routine

  RETURN, z

END
