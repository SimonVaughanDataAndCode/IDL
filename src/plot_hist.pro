PRO plot_hist, Y, X0=x0, X1=x1, DX=dx, NOPLOT=noplot, $
               _EXTRA=EXTRA_KEYWORDS, FUDGE=fudge, FILL=fill, $
               COL_FILL=col_fill, COL_LINE=col_line, DY=dy, $
               NOLINES=nolines

; ----------------------------------------------------------
;+
; NAME:
;       PLOT_HIST
;
; PURPOSE:
;       Make a nice histogram plot
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       PLOT_HIST, data, bins
;
; INPUTS:
;       y             - (vector) bin values
;
; OPTIONAL INPUTS:
;       x0            - (vector) location of bin / lower edge
;       x1            - (vector) location of bin / upper edge
;       dx            - (float) bin width
;       noplot        - (logical) do not make a new plot?
;       fudge         - (flat) extend axis ranges by this factor
;       col_fill      - (integer) colour for filling histogram
;       col_line      - (integer) colour for line drawing
;       nolines       - (logical) set to turn off the vertical lines
;       dy            - (vector) errors on y values
;
; OUTPUTS:
;       none
;
; OPTIONAL OUTPUTS:
;       none
;
; DETAILS:
;       Provide a nicer alternative to the usual PLOT...PSYM=10 method
;       of plotting histograms in IDL. Input are the histogram values
;       Y and the bin locations. Bins are specified by the lower and
;       upper edges (BIN_L, BIN_U) and/or the bin width (DX). It is
;       assumed that the bins are in accending order (i.e. increasing
;       X value) and that the bins are contiguous, i.e. upper edge of
;       bin[j] is the lower edge of bin[j+1].
;
;       This differs from the usual (PLOT,x,y,PSYM=10) method in a few
;       respects. First the bins are drawn from the lower edge of the
;       lowest bin to the upper edge of the highest bin - the usual
;       method crops the bottom and top bins. Secondly vertical lines
;       are added to define each bin properly (this can be turned off
;       using the NOLINES keyword). The histogram can be filled with a
;       solid colour by setting the COL_FILL keyword.
;
;       By default a new graphics plot is opened with X and Y axis
;       ranges determined from the ranges of the data, with a 4 per
;       cent margin either side. This margin can be adjusted using the
;       FUDGE keyword. Setting the NOPLOT keyword means the existing
;       graphics plot is used, i.e. the histogram is overplotted on
;       the current graphics device.
;
; EXAMPLE USAGE:
;
;       n = 15
;       dx = 2.0
;       bin_l = INDGEN(n) * dx + 15.0
;       bin_u = bin_l + dx
;       bin_y = bin_l + RANDOMN(seed, n) * 5
;       PLOT_HIST, bin_y, x0=bin_l, x1=bin_u
;
;         And the same plot 'filled'
;
;       PLOT_HIST, bin_y, x0=bin_l, x1=bin_u, col_fill=200, THICK=2
;
;         Same plot without vertical lines around each bin
;
;       PLOT_HIST, bin_y, x0=bin_l, x1=bin_u, /nolines
;
;         The same data but using user defined plotting region/axes
;       
;       PLOT, [10,35],[0,50], /NODATA, POSITION=[0.1,0.3,0.9,0.9], XTICKNAME=REPLICATE(' ', 30)
;       PLOT_HIST, bin_y, x0=bin_l, x1=bin_u, col_fill=200, /NOPLOT
;
; HISTORY:
;       23/06/2011 - v1.0 - first working version
;       06/05/2013 - v1.1 - added DY input keyword
;
; NOTES:
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)

  COMPILE_OPT idl2, HIDDEN

; watch out for errors

  ON_ERROR, 2

; ---------------------------------------------------------
; Check arguments

  N = N_ELEMENTS(y)

; check there are enough data to plot

  if (N le 2) then begin
      PRINT,'** Too few data in PLOT_HIST'
      RETURN
  endif
 
  if NOT KEYWORD_SET(fudge) then fudge = 0.04

; make sure that the bin lower edges are defined

  n1 = N_ELEMENTS(x0)
  n2 = N_ELEMENTS(x1)
  if (n1 eq 0) then begin
      if KEYWORD_SET(dx) then begin
          x0 = INDGEN(n)
          bin_l = x0 - 0.5*dx
          bin_u = x0 + 0.5*dx
      endif else begin
          PRINT, '** No X positions specified in PLOT_HIST'
          RETURN
      endelse
  endif else begin
      if (n1 ne n) then begin
          PRINT, '** X and Y arrays of different size in PLOT_HIST'
          RETURN
      endif
      if (n2 eq n1) then begin
          bin_l = x0
          bin_u = x1

; make sure that the bin upper edges are defined

      endif else begin
          if (n2 eq 0) then begin
              if NOT KEYWORD_SET(dx) then dx = x0[1] - x0[0]
              bin_l = x0
              bin_u = [x0[1:(n-1)], x0[n-1]+dx]
          endif else begin
              PRINT, '** X0, X1 of different size in PLOT_HIST'
              RETURN
          endelse
      endelse
  endelse

; ---------------------------------------------------------
; Main routine

; first open a graphics window unless told otherwise (NOPLOT). With
; the ranges determined either automatically from the data ranges or
; adjusted using the FUDGE input.

  x = (bin_l + bin_u) * 0.5
  xrange = !X.RANGE
  yrange = !Y.RANGE

  if NOT KEYWORD_SET(noplot) then begin

      xspan = MAX(bin_u) - MIN(bin_l) 
      xf = xspan * fudge
      xrange = [MIN(bin_l)-xf, MAX(bin_u)+xf]
      
      yspan = MAX(y) - MIN(y)
      yf = yspan * fudge
      yrange = [MIN(y)-yf, MAX(y)+yf]
  
      PLOT, x, y, XRANGE=xrange, YRANGE=yrange, /NODATA, /YSTYLE, /XSTYLE
 
  endif

; next plot a filled/coloured histogram using the POLYFILL
; command, if requested - i.e., if COL_FILL is set. The filled
; histogram is plotted before the line histogram.

  miny = yrange[0]
  clip = [xrange[0], yrange[0], xrange[1], yrange[1]]

  if KEYWORD_SET(col_fill) then begin

      xx = MAKE_ARRAY(2*n+2)
      yy = MAKE_ARRAY(2*n+2)
      yy[0] = miny < MIN(y)
      for i = 0, N-1 do begin
          xx[2*i] = bin_l[i]
          xx[2*i+1] = bin_l[i]
          yy[2*i+1] = y[i]
          yy[2*i+2] = y[i]
      endfor
      xx[2*n] = bin_u[n-1]
      xx[2*n+1] = bin_u[n-1]
      yy[2*n+1] = miny
      
      POLYFILL, xx, yy, COLOR=col_fill, NOCLIP=0

; re-plot the axes - which may have been partially covered-up by the
; filled histogram.

      AXIS, XAXIS=0, XRANGE=xrange, /XSTYLE, xtickname=REPLICATE(' ', 30)
      AXIS, YAXIS=0, YRANGE=yrange, /YSTYLE, ytickname=REPLICATE(' ', 30)

  endif

; set the line colour if not already supplied.

  if NOT KEYWORD_SET(col_line) then col_line = !P.COLOR
      
; plot the line histogram over the pre-existing plot. 

  for i = 0, N-1 do begin
      PLOTS, [bin_l[i], bin_u[i]], [y[i], y[i]], color=col_line, NOCLIP=0, _EXTRA=extra_keywords
  endfor
  for i = 0,N-2 do begin
      PLOTS, [bin_u[i], bin_u[i]], [y[i], y[i+1]], color=col_line, NOCLIP=0, _EXTRA=extra_keywords
  endfor

; plot vertical lines around each bin - unless NOLINES is set

  if NOT KEYWORD_SET(nolines) then begin

      for i = 0, N-1 do begin
          PLOTS, [bin_l[i], bin_l[i]], [miny, y[i]], color=col_line, NOCLIP=0, _EXTRA=extra_keywords
      endfor
      PLOTS, [bin_u[N-1], bin_u[N-1]], [miny, y[N-1]], color=col_line, NOCLIP=0, _EXTRA=extra_keywords

  endif

  IF (N_ELEMENTS(dy) gt 0) THEN BEGIN
    PLOT_ERR, x, y, dy
  ENDIF

; ---------------------------------------------------------
; Return to user

END
