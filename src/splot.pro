PRO splot, data, POS=pos, PLINE=pline, LABELS=labels, EXPAND=expand, $
           SIZE=size, COLOR=color, PSYM=psym, SYMSIZE=symsize, PCORR=pcorr, $
           CORR_METH=corr_meth, CORR_INDEX=corr_index, PSOUT=psout, $
           AXIS_SIZE=axis_size, XLABOFFSET=xlaboffset, $
           YLABOFFSET=ylaboffset, XYZERO=xyzero

; ----------------------------------------------------------
;+
; NAME:
;       SPLOT
;
; PURPOSE:
;       Produce an NxN matrix of scatter plots
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       splot, data
;
; INPUTS:
;       data      - (array) an NxM array of data to be plotted
;
; OPTIONAL INPUTS:
;       POS       - (vector) specify plotting range
;       PLINE     - (logical) plot lines between points?
;       LABELS    - (string, vector) list of labels for the variables
;       EXPAND    - (float) expansion factor of panel axis ranges
;       SIZE      - (float) size of the text labels
;       COLOR     - (integer) colour for line drawing
;       PSYM      - (integer) choice of plotting symbol
;       SYMSIZE   - (float) size of plotting symbols
;       PCORR     - (logical/int) display correlation coefficients?
;       CORR_METH - (integer) 1=Spearman rank, 2=Kendall tau
;       CORR_INDEX - (float) scaling of text size with correlation
;       PSOUT     - (logical) produce a PostScript output?
;       AXIS_SIZE - (float) scaling for axis labels sizes
;       XLABOFFSET - (float) fractional x-offset for the text labels
;       YLABOFFSET - (float) fractional y-offset for the text labels
;       XYZERO    - (array) force which plots to include x=0 point
;
; OUTPUTS:
;       none
;
; OPTIONAL OUTPUTS:
;       none
;
; DETAILS:
;       Based on the PAIRS() command in R/SPLUS.
;
;       Given an NxM array of input, containing M observations of N
;       variables we generate an NxN matrix of scatterplots, with each
;       panel showing one of the pairings of variables. Along the
;       diagonal (i=j) we print the name of the variable. Axes are
;       drawn on the far top, bottom, left and right sides only, and
;       alternative between sides to give more room of the
;       characters. 
;
;       The PSYM and PSYMSIZE set the symbol type and size. These are
;       the same as used in PLOT, except by default we use the filled
;       circle symbol generated with the PLOTSYM routine.
;
;       The axes within each panel encompass the range of the data
;       values plus a constant value either size propotional to the
;       full range. This constant is set using the EXPAND keyword
;       The default value of 1.1 means the axes extend in both
;       directions by 10% of the full range. 
;
;       The LABELS option specifies the names of each variable, to
;       appear in the diagonal panels. The SIZE keyword sets the size
;       of the text labels in the diagonal panels.
;
;       The AXIS_SIZE keyword sets the character size of the axis
;       labels. 
;
;       The PLINE keyword determines whether to draw lines between the
;       points in the lower-left panels. The colour of the line can be
;       set with the COLOR keyword.
;
;       The PCORR keywords determines whether to use the upper-right
;       panels to display the correlation coefficient for each
;       pairwise compbination of variables. There is a choice of
;       correlation coefficient: Pearson's linear one (default), or
;       Spearman rank-order (CORR_METH=1), or Kendall's tau
;       (CORR_METH=2). The CORR_INDEX keyword sets a scaling index -
;       if this is set then the correlation coefficient is printed
;       using larger text for larger values, this parameter sets the
;       power law index of scaling between values and text
;       size. Default value is constant equal to SIZE (same as text in
;       the diagonal panels). Set CORR_INDEX = 1 to make the text size
;       scaling linearly with coefficient (I find 0.5 a good value to
;       use). Setting PCORR=1 (or TRUE) plots coefficients, setting 
;       PCORR=2 plots the p-values.
;
;       The PSOUT keyword sends the output to a PostScript file called
;       idl.ps. 
;         
; EXAMPLE USAGE:
;         N = 30
;         w = RANDOMN(seed, N)
;         x = RANDOMN(seed, N)
;         y = x*x + 1.0 + randomn(seed, N)*0.5
;         z = y + x + randomn(seed, N)*0.25
;         data = TRANSPOSE([[w],[x],[y],[z]])
;
;         splot, data
;
;            ...and now with more options...
;
;         splot, data, labels=["w","x","y","z"], /PLINE, /PCORR,  $
;                corr_index=0.5, expand=1.25, symsize=0.7
;
; HISTORY:
;       02/08/10 - v1.0 - first working version
;                  12/08/10 - v1.1 - bug fix. Added /[XY]STYLE to AXIS
;                             command to align axis labels
;                             properly. Added XYZERO and [XY]LABOFFSET
;                             keywords. 
;                  24/02/11 - v1.2 added p-value output (PCORR)
;
; PROCEDURES CALLED:
;       PLOTSYM, NUMBER_FORMATTER, PS_OPEN, PS_CLOSE
;
; NOTES:
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)

;  COMPILE_OPT idl2

; watch out for errors

;  on_error, 2

; ---------------------------------------------------------
; Check arguments

  s = SIZE(data)

  if (s[0] ne 2) then begin
     print,'** DATA not of dimension 2 in SPLOT.'
     return
  endif

; number of variables

  N = s[1]

; number of data points in N dimensional space

  M = s[2]

; if POS not set then use default

  if NOT KEYWORD_SET(pos) then pos = [0.1, 0.1, 0.9, 0.9]

; check that POS has the right shape

  sp = SIZE(pos)
  if (sp[0] ne 1 or sp[1] ne 4) then begin
      print, '** POS has the incorrect form in SPLOT.'
      return
  endif

; set EXPAND factor for axes if not set

  if NOT KEYWORD_SET(expand) then expand = 1.1

; set labels for axes if not set

  if NOT KEYWORD_SET(labels) then begin
      labels = INDGEN(N, /STRING)
      labels = STRCOMPRESS(labels, /REMOVE_ALL)
  endif

  if ( N_ELEMENTS(labels) ne N ) then begin
      print, '** Incorrect size of LABELS in SPLOT'
      return
  endif

; set text size of not set

  if NOT KEYWORD_SET(size) then size = 2.0

; set plotting symbol if not already set

  if NOT KEYWORD_SET(psym) then psym=8

; set SYMSIZE if not already set

  if NOT KEYWORD_SET(symsize) then symsize=1.0

; set CORR_INDEX if not already set

  if NOT KEYWORD_SET(corr_index) then corr_index=0

; set AXIS_SIZE if not already set

  if NOT KEYWORD_SET(axis_size) then axis_size=1

; set [XY]LABOFFSET if not already set

  if NOT KEYWORD_SET(xlaboffset) then xlaboffset=0.0
  if NOT KEYWORD_SET(ylaboffset) then ylaboffset=0.0

; set [XY]ZERO if not already set

  if NOT KEYWORD_SET(xyzero) then xyzero = MAKE_ARRAY(N)

; ---------------------------------------------------------
; Main routine

; open PS device?

  if KEYWORD_SET(psout) then PS_OPEN

; define additional plotting symbols

  plotsym, 0, symsize, /FILL, THICK=5

; define the POSITIONS of each of the NxN panels

  xmin = pos[0]
  ymin = pos[1]
  xmax = pos[2]
  ymax = pos[3]
  dx = (xmax - xmin)/N
  dy = (ymax - ymin)/N

  i = INDGEN(N+1)
  x0 = i*dx + xmin
  y0 = i*dy + ymin
  x1 = x0 + dx
  y1 = y0 + dy

  for i = 0, N-1 do begin

     for j = 0, N-1 do begin

; set a flag to open a new graphics window if this is the first panel
; to be plotted

         noe = 1
         if (i eq 0 and j eq 0) then noe = 0

; set a flag to open panel but not plot the data for the diagonal
; panels, i.e. when i = j

         nod = 0
         if (i eq j) then nod=1
         if KEYWORD_SET(pcorr) then begin
             if (i gt j) then nod = 1
         endif

; define the position in the window for this panel

         k = (N-1) - j
         position = [x0[i], y0[k], x1[i], y1[k]]

; extract the two variables to be plotted

         x_i = REFORM(data[i, *])
         y_j = REFORM(data[j, *])

; compute the correlation coefficient
; use SPEARMAN or KENDALL rank correlation methods
; if the CORR_METH keyword is set

         if NOT KEYWORD_SET(corr_meth) then begin
             r = CORRELATE(x_i, y_j, /DOUBLE)
             df = M - 2 ; compute the degrees of freedom
             if (1.0 - r le 1.0E-7) then begin
                  prob = 0.0
             endif else begin
                  t = r/sqrt((1.0-r*r)/FLOAT(df)) ; compute the t statistic
                  prob = 2.0 * (1.0 - T_PDF(ABS(t), df))  ; calculate the two-side 'tail area' probability 
             endelse
         endif

         if KEYWORD_SET(corr_meth) then begin
             if (corr_meth eq 1) then begin
                 cor = (R_CORRELATE(x_i, y_j))
             endif else begin
                 cor = (R_CORRELATE(x_i, y_j, /KENDALL))
             endelse
             r = cor[0]
             prob = cor[1]
         endif

; convert correlation to a string for adding to plot

         r_out = NUMBER_FORMATTER(r, DECIMAL=2)
         if KEYWORD_SET(PCORR) then begin
            if (PCORR eq 2) then r_out = NUMBER_FORMATTER(prob, DECIMAL=3)
         endif 

; define a suitable range for the plot

         xspan = (max(x_i) - min(x_i))
         yspan = (max(y_j) - min(y_j))

         dx = xspan * (expand-1.0)
         dy = yspan * (expand-1.0)

         xrange = [min(x_i)-dx, max(x_i)+dx]
         yrange = [min(y_j)-dy, max(y_j)+dy]

         if XYZERO[i] then xrange[0] = (xrange[0] < 0)
         if XYZERO[j] then yrange[0] = (yrange[0] < 0)

         xcent = MEAN(xrange)
         ycent = MEAN(yrange)

         xcentxt = xcent+xlaboffset*xspan
         ycentxt = ycent+ylaboffset*yspan

; plot the data in the panel, without axes

         plot, x_i, y_j, POSITION=position, NOERASE=noe, $
           XTICKFORMAT='(A1)', YTICKFORMAT='(A1)', XRANGE=xrange, $
           YRANGE=yrange, /XSTYLE, /YSTYLE, NODATA=nod, PSYM=psym, $
           SYMSIZE=symsize

; if a diagonal panel (i=j) then print the label instead of data

         if (nod eq 1 and i eq j) then begin           
             XYOUTS, xcentxt, ycentxt, labels[i], ALIGNMENT=0.5, $
               CHARSIZE=size
         endif

; if in the lower left panels plot optional lines between points

         if KEYWORD_SET(pline) then begin
             if (i lt j) then begin
                 OPLOT, x_i, y_j, COLOR=color
                 OPLOT, x_i, y_j, PSYM=psym, SYMSIZE=symsize
             endif
         endif

; if in the upper right panels display optional correlation values

         if KEYWORD_SET(pcorr) then begin
             if (i gt j) then begin
                 corr_size = size * (ABS(r)^corr_index) + size/2
                 XYOUTS, xcent, ycent, r_out, ALIGNMENT=0.5, $
                   CHARSIZE=corr_size
             endif
         endif

; add axes if on first/last row/column of the matrix plot
; but alternate either side to make sure they don't overlap

         if (i eq 0) then begin
             if (j mod 2 eq 0) then begin
                 if NOT (KEYWORD_SET(pcorr) and j eq 0) then begin
                     AXIS, YAXIS=0, CHARSIZE=axis_size, /YSTYLE
                 endif
             endif 
         endif
         if (i eq N-1) then begin
             if (j mod 2 eq 1) then begin
                 AXIS, YAXIS=1, CHARSIZE=axis_size, /YSTYLE
             endif 
         endif

         if (j eq 0) then begin
             if (i mod 2 eq 1) then begin
                 if NOT (KEYWORD_SET(pcorr) and i eq N-1) then begin
                     AXIS, XAXIS=1, CHARSIZE=axis_size, /XSTYLE
                 endif
             endif 
         endif
         if (j eq N-1) then begin
             if (i mod 2 eq 0) then begin
                 if NOT (KEYWORD_SET(pcorr) and i eq N-1) then begin
                     AXIS, XAXIS=0, CHARSIZE=axis_size, /XSTYLE
                 endif
             endif 
         endif

     endfor

  endfor

; close PS device?

  if KEYWORD_SET(psout) then PS_CLOSE

; ---------------------------------------------------------
; Return to user

END
