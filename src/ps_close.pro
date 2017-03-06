PRO PS_CLOSE

; ----------------------------------------------------------
;+
; NAME:
;       PS_CLOSE
;
; PURPOSE:
;       Close PS device and return to previous settings
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       PS_CLOSE
;
; INPUTS:
;       NONE
;
; OPTIONAL INPUTS:
;       NONE
;
; OUTPUTS:
;       NONE
;
; DETAILS:
;       This simple procedure closes a PostScript graphics
;       device opened by PS_OPEN.PRO and returns the
;       previous graphics settings.
;
; HISTORY:
;       24/01/2007  -- v1.0 -- first working version
;       03/10/2013  -- v1.1 -- added call to PSLANDFIX
;
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

  if (!d.name ne 'PS') then begin
      print, '** PS device not open in PS_CLOSE.'
      return
  endif

; finish with the PS device

  device, /close_file

; access common block used by PS_OPEN

  common ps_block, cpd, plot_thick, x_thick, y_thick, $
    char_thick, char_size, font_type, lands, filen

  if (n_elements(cpd) eq 0) then begin
      print,'** Missing information in PS_CLOSE.'
      print,'** Ensure that PS_OPEN has been run properly.'
      return
  endif

; fix the upside-down landscape .ps file 

  if (lands eq 1) then begin
;    pslandfix, filen
  endif

; return to original plot device and settings

  SET_PLOT, cpd

  !p.thick = plot_thick
  !x.thick = x_thick
  !y.thick = y_thick
  !p.charthick = char_thick
  !p.charsize = char_size
  !p.font = font_type

END
