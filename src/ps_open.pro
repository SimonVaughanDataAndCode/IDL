PRO PS_OPEN, filename, portrait=portrait, colour=colour, $
             position=position, charsize=charsize, thick=thick, $
             cthick=cthick, ctab=ctab

; ----------------------------------------------------------
;+
; NAME:
;       PS_OPEN
;
; PURPOSE:
;       Open and configure a PS device for graphical output
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       PS_OPEN,output.ps
;
; INPUTS:
;       filename  - (string) name of output file (default 'idl.ps')
;
; OPTIONAL INPUTS:
;       portrait  - (logical) portrait orientation if set
;       position  - (vector) boundaries (x0,y0,x1,y1)
;       colour    - (logical) switch on colour output
;       charsize  - (float) character scaling (default = 1.4)
;       thick     - (float) line thickness (default = 3.0)
;       cthick    - (float) character line thickness (default = THICK)       
;       ctab      - (integer) colour table
;
; OUTPUTS:
;       NONE
;
; DETAILS:
;       This simple procedure opens a PostScript graphics
;       device and configures it ready to produce standardised
;       publication quality plots. The PostScript graphics
;       device should then be closed with the PS_CLOSE
;       procedure.
;
; HISTORY:
;       24/01/2007  -- v1.0 -- first working version
;       14/07/2007  -- v1.1 -- added PORTRAIT keyword
;                              added CHARSIZE keyword
;                              added THICK keyword
;                              added CTHICK keyword
;       04/08/2008  -- v1.2 -- added CTAB keyword
;       20/07/2012  -- v1.3 -- added explicit DECOMPOSED=0 option
;                               when /COLOUR is used.
;       04/07/2013  -- v1.4 -- changed to DECOMPOSED=1
;       03/10/2013  -- v1.5 -- minor bug fix (removed
;                               KEYWORD_SET from check of lanscape)
;                               Added more info to common block
; EXAMPLE USAGE:
;       IDL> x = findgen(10)
;       IDL> y = x^2.0 - 4.0*x - 4.0
;       IDL> PS_OPEN,'myplot.ps'
;       IDL> PLOT,x,y,xtitle="position",ytitle="height"
;       IDL> PS_CLOSE
; or
;       IDL> z = SHIFT(DIST(64),32,32)
;       IDL> z = EXP(-(0.1*z)^2)
;       IDL> LOADCT,1
;       IDL> PS_OPEN,/colour,charsize=3
;       IDL> SHADE_SURF,z,xstyle=1,ystyle=1,pixels=1000
;       IDL> PS_CLOSE
;
; NOTES: 
;       Loosely based on the PSON.PRO procedure by Liam E Gumley.
;

;-
; ----------------------------------------------------------

; watch out for errors
  on_error, 2

; ----------------------------------------------------------
; Check the arguments

  if (n_elements(filename) eq 0) then filename='idl.ps'

  if (n_elements(position) eq 0) then position=[0.0,0.0,1.0,1.0]

  if keyword_set(portrait) then begin
      landscape = 0
  endif else begin
      landscape = 1
  endelse

  if (n_elements(charsize) eq 0) then charsize = 1.8

  if (n_elements(thick) eq 0) then thick = 3.0

  if (n_elements(cthick) eq 0) then cthick = thick

  if (n_elements(ctab) eq 0) then ctab=1

; ----------------------------------------------------------
; Check that PS device is not already configured

  if (!d.name eq 'PS') then begin
      print, '** PS device already open in PS_OPEN.'
      return
  endif

; ----------------------------------------------------------
; Remember the current settings (place in common block)

 common ps_block, cpd, plot_thick, x_thick, y_thick, $
   char_thick, char_size, font_type, lands, filen

 cpd = !d.name
 plot_thick = !p.thick
 x_thick = !x.thick
 y_thick = !y.thick
 char_thick = !p.charthick
 char_size = !p.charsize
 font_type = !p.font
 lands = landscape
 filen = filename

; ----------------------------------------------------------
; Open and configure the PS device

; change the system variables for plots
  !p.thick = thick
  !x.thick = thick
  !y.thick = thick
  !z.thick = thick
  !p.charthick = cthick
  !p.charsize  = charsize
  !p.font      = -1

; Size, shape of A4 paper
  width=21.0          
  height=29.7

  case landscape of 
      0 : begin
          xsize = (position[2] - position[0]) * width
          ysize = (position[3] - position[1]) * height
          xstart = position[0] * width
          ystart = position[1] * height
      end
      1 : begin
          ysize = (position[3] - position[1]) * width
          xsize = (position[2] - position[0]) * height
          xstart = position[1] * width
          ystart = (1.0 - position[0]) * height 
      end
  endcase

; open the plot device
  SET_PLOT, 'PS'
   DEVICE, FILE=filename, LANDSCAPE=landscape

; set options
  DEVICE, XSIZE=xsize, YSIZE=ysize, XOFFSET=xstart, YOFFSET=ystart
  if KEYWORD_SET(colour) THEN BEGIN
    DEVICE, DECOMPOSED=1, /COLOR, BITS_PER_PIXEL=8
    LOADCT, ctab
  ENDIF

; 'trick' the system into using "Complex Roman" font from now on...
  XYOUTS, '!6'

END
