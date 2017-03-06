PRO TS_SPEC_PLOT, freq, pow, dp=dp, frange=frange, $
                  fp=fp, ps=ps, log=log, col=col, $
                  _EXTRA=EXTRA_KEYWORDS

; ----------------------------------------------------------
;+
; NAME:
;       TS_SPEC_PLOT
;
; PURPOSE:
;       Plot a binned periodogram, from e.g. TS_MULTISPEC
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       TS_SPEC_PLOT,f,p,dp=dp,/fp
;
; INPUTS:
;       freq      - (vector) frequencies [abcissa]
;       pow       - (vector) powers [ordinates]
;
; OPTIONAL INPUTS:
;       dp        - (vector) error bar on powers
;       fp        - (logical) plot in power*frequency ?
;       frange    - (vector) lower/upper frequency to plot
;       ps        - (logical) produce PostScript output?
;       log       - (logical) plot on log-log scale
;       col       - (scalar) colour for data points
;
; OUTPUTS:
;       <plot>    - power spectral plot
;
; OPTIONAL OUTPUTS:
;       idl.ps    - (file) PS figure output
;
; DETAILS:
;       Plot a binned periodogram as produced by e.g. TS_MULTISPEC.
;
; PROCEDURES CALLED:
;       PLOT_ERR
;
; HISTORY:
;       02/05/07  - v1.0 - first working version
;       04/05/07  - v1.1 - added EXTRA_KEYWORDS option
;       17/05/07  - v1.2 - added COL option
;       05/06/09  - v1.3 - changed to histrogram style plot
;                          added check for -ve values when using
;                          log scale
;
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments

; is input array large enough to bother with?
  N = n_elements(freq)
  if (N le 2) then begin
      print,'** Too few data in TS_SPEC_PLOT'
      return
  endif

; if frange not defined, set default
  if (n_elements(frange) ne 2) then begin
      mask = indgen(N)
  endif else begin
      mask = where((freq ge frange[0]) and (freq le frange[1]))
  endelse

; ----------------------------------------------------------
; Main part of procedure

; select only good data to plot
  x=freq[mask]
  y=pow[mask]
  N = n_elements(x)
  if (n_elements(dp) eq N) then begin
      dy=dp[mask]
      if keyword_set(fp) then dy=x*dy
  endif

; if FP selected use freq*pow as ordinate
  if keyword_set(fp) then y=freq[mask]*pow[mask]

; define plotting symbol (filled circle)

  A = FINDGEN(17) * (!PI*2/16.)
  USERSYM, COS(A), SIN(A), /FILL

; open the PS device if requested
  if keyword_set(ps) then ps_open

; define X-axis range (FREQ)
  xoffset=0.1
  xspan=max(x)-min(x)
  xrange=[min(x)-xoffset*xspan,max(x)+xoffset*xspan]
  if keyword_set(log) then begin
       xspan=alog10(max(x)/min(x))
       xrange[0]=alog10(min(x))-xoffset*xspan
       xrange[1]=alog10(max(x))+xoffset*xspan
       xrange=10.0^xrange
  endif

; define Y axis range (FLUX)
  yoffset=0.1
  yspan=max(y)-min(y)
  yrange=[min(y)-yoffset*yspan,max(y)+yoffset*yspan]
  if keyword_set(log) then begin
       mask = WHERE(y gt 0.0)
       yspan=alog10(max(y)/min(y[mask]))
       yrange[0]=alog10(min(y[mask]))-yoffset*yspan
       yrange[1]=alog10(max(y))+yoffset*yspan
       yrange=10.0^yrange
  endif

  if keyword_set(fp) then begin
      ytitle = "freq * power density ([rms/mean]!U2!N)"
  endif else begin
      ytitle = "power density ([rms/mean]!U2!N Hz!U-1!N)"
  endelse

; Define plot region, axes and labels
  plot,x,y,psym=8,symsize=0.5,position=[0.15,0.15,0.95,0.95], $
    noclip=0, xrange=xrange, xstyle=1, yrange=yrange, ystyle=1, $
    xlog=keyword_set(log), ylog=keyword_set(log), /nodata, $
    xtitle="frequency (Hz)", ytitle=ytitle, _extra=extra_keywords

; plot data as histogram
    oplot,x,y,psym=10

; add data points in COLOR if needed
  if keyword_set(col) then begin
      plots,x,y,psym=8,symsize=0.5,color=col, _extra=extra_keywords
  endif else begin
      plots,x,y,psym=8,symsize=0.5, _extra=extra_keywords
  endelse

; Plot error bars if requested
  if (n_elements(dy) eq N) then begin
      if keyword_set(col) then begin
          plot_err,x,y,dy,color=color
      endif else begin
          plot_err,x,y,dy
      endelse

; mark the points where error bar goes -ve on log scale
      if keyword_set(log) then begin
          mask = WHERE( y le dy, count)
          if (count gt 0) then begin
              for i=0,count-1 do begin
                  j = mask[i]
                  xx = [x[j], x[j]]
                  yy = [yrange[0], y[j] + dy[j]]
                  plots, xx, yy, /data, noclip=0
              endfor
          endif
      endif

  endif

; close the PS device if requested
  if keyword_set(ps) then ps_close

; ----------------------------------------------------------
; Return the data array to the user

  return

END
