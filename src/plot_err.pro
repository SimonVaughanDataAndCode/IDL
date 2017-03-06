PRO PLOT_ERR, x, y, dy, dy2=dy2, dx1=dx1, dx2=dx2, width=width, $
              newplot=newplot,_EXTRA=EXTRA_KEYWORDS

; ----------------------------------------------------------
;+
; NAME:
;       PLOT_ERR
;
; PURPOSE:
;       Plot error bars on x,y plot (in y and x).
;
; AUTHOR:
;       Simon Vaughan
;
; CALLING SEQUENCE:
;       PLOT_ERR,x,y,dy
;
; INPUTS:
;       x       - x position (ordinate)
;       y       - y position (abscissa)
;       dy      - error bar on y (assumed symmetrical)
;
; OPTIONAL INPUTS:
;       width   - width of horizontal bars (default = 0)
;       dy2     - upper error bar on y (assumed asymmetrical if supplied)
;       dx1     - error bar on x (assumed symmetrical)
;       dy2     - upper error bar on y (assumed asymmetrical if supplied)
;       newplot - (logical) start a new graphics window?
;
; OUTPUTS:
;       none
;
; DETAILS:
;       Adds error bars on y positions of an x,y plot.
;       If only dy is supplied then errors are assumed
;       symmetrical, i.e. y +/- dy. If both dy and dy2 are
;       supplied then errors are from y-dy to y+dy2.
;       If dx1 is supplied then also plot x +/- dx1 bars.
;       (And asymmetrical x bars if dx2 also supplied.)
;       Width is in fraction of screen size (0.0 to 1.0).
;
; HISTORY:
;       Based on the routine ERR_PLOT by Liam E. Gumley.
;       27/04/2007  -- v1.1 -- added NEWPLOT option
;       09/05/2007  -- v1.2 -- set dy=1.0 if not set correctly
;
;-
; ---------------------------------------------------------

; watch out for errors
 on_error,2

; ---------------------------------------------------------
; Check arguments

  n = n_elements(x)
  if (n eq 0) then begin
      print,'** X is undefined in PLOT_ERR.'
      return
  endif
  if (n_elements(y) ne n) then begin
      print,'** Y and X different lengths in PLOT_ERR.'
      return
  endif
  if (n_elements(dy) ne n) then begin
      print,'** Arrays X, Y and DY of different length in PLOT_ERR. Continuing....'
      help,x, y, dy
      dy = make_array(n,value=1.0)
  endif

; Check for asymmetric error bars
  if (n_elements(dy2) eq 0) then dy2=dy

; Use zero width bars if no width supplied
  if (n_elements(width) eq 0) then width = 0.0

; Check for errors on x positions
  if (n_elements(dx1) ne 0) then begin
      if (n_elements(dx1) ne n) then begin
          print,'** X and DX1 of different lengths in PLOT_ERR.'
          return
      endif
      if (n_elements(dx2) eq 0) then dx2=dx1
      if (n_elements(dx2) ne n) then begin
          print,'** X and DX2 of different lengths in PLOT_ERR.'
          return
      endif
  endif 

; ---------------------------------------------------------
; Plot the error bars

  if keyword_set(newplot) then begin
      xspan=max(x)-min(x)
      yspan=max(y)-min(y)
      offset=0.1
      xrange=[min(x)-offset*xspan,max(x)+offset*xspan]
      yrange=[min(y-dy)-offset*yspan,max(y+dy)+offset*yspan]
      plot, x, y, /data, noclip=0, psym=3, xrange=xrange, $
         xstyle=1, yrange=yrange, ystyle=1, _extra=extra_keywords
  endif 

; Define upper/lower y points
  yhigh = y + dy2
  ylow  = y - dy

;- Plot the error bars
  for i = 0L,n-1L do begin

  ;- Plot vertical bar using data coordinates
      xdata = [x[i], x[i]]
      ydata = [ylow[i], yhigh[i]]
      plots, xdata, ydata, /data, noclip=0, _extra=extra_keywords


; Plot error bars on x if requested
      if (n_elements(dx1) eq n) then begin
          xdata = [x[i]-dx1[i], x[i]+dx2[i]]
          ydata = [y[i], y[i]]
          plots, xdata, ydata, /data, noclip=0, _extra=extra_keywords
      endif

      if (width gt 0.0) then begin
;- Compute horizontal bar width in normal coordinates
          normalwidth = (!x.window[1] - !x.window[0]) * width

;- Plot horizontal bar using normal coordinates
          lower = convert_coord(x[i], ylow[i], /data, /to_normal)
          upper = convert_coord(x[i], yhigh[i], /data, /to_normal)
          xdata = [lower[0] - 0.5 * width, lower[0] + 0.5 * width]
          ylower = [lower[1], lower[1]]
          yupper = [upper[1], upper[1]]
          plots, xdata, ylower, /normal, noclip=0, _extra=extra_keywords
          plots, xdata, yupper, /normal, noclip=0, _extra=extra_keywords
      endif


  endfor

; ---------------------------------------------------------

END
