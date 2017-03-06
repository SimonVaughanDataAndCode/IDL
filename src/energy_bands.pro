FUNCTION energy_bands, N_CHAN, E_min=E_MIN, E_MAX=E_max, $
                       FILE=file, CHANNELS=channels, DE=de, LIN=lin, $
                       CHATTER=chatter

; ----------------------------------------------------------
;+
; NAME:
;       ENERGY_BANDS
;
; PURPOSE:
;       Calculate logarithmically spaced PI bands
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       result = ENERGY_BANDS(50)
;
; INPUTS:
;       N_chan       - Number of bands to use
;
; OPTIONAL INPUTS:
;       E_min         - Minimum channel to use (default = 200)
;       E_max         - Maximum channel to use (default = 10000)
;       file          - Write list to an ASCII file (bands.dat)?
;       channels      - Write list of channels to ASCII file
;                       (channels.dat)?
;       dE            - coarseness of the PI channels
;       LIN           - (logical) use linear (not log) spacing?
;       CHATTER       - (logical) control amount of output to screen
;
; OUTPUTS:
;       result        - array of lower/upper PI bounds
;
; OPTIONAL OUTPUTS:
;       none
;
; DETAILS:
;       Divide the PI range from E_MIN to E_MAX into N_CHAN
;       non-overlapping bands. Output the list of lower/upper band
;       boundaires as an array (and optionally output to as ASCII file
;       called bands.dat).
;
;       The PI boundaries are discretized into multiples of dE, since
;       the pn detector bins are 5 eV (PI) wide and the MOS are 15 eV
;       (PI). The default is dE = 5 which is massively below the
;       resolution anyway. But with large N_chan the lower energy
;       bands may be smaller than dE, which can cause problems. If
;       this happens set dE=1.
;
; EXAMPLE USAGE:
;      e_bounds = ENERGY_BANDS(50, E_min=200, E_max=10000)
;
; HISTORY:
;      15/06/2010 - v1.0 - first working version
;      16/07/2010 - v1.1 - made PI channels multiples of dE.
;      09/08/2010 - v1.2 - added LIN keyword
;      19/07/2012 - v1.3 - added a couple of catches and error
;                           warnings in case of N_chan = 1 and
;                           E_max <= E_min. Use FLOOR rather than 
;                           ROUND to ensure lowest/highest bins
;                           are not smaller than expected.
;      19/07/2013 - v2.0 - changed the way the main loop works.
;                           Now ensures each (output) band contains 
;                           at least one detector (input) channel, 
;                           while trying to preserve lin/log spacing.
;                           This means the number of channels can be
;                           smaller than requested N_CHAN (typically
;                           when N_CHAN > 36). N_CHAN is adjusted 
;                           accordingly. Added CHATTER keyword
;       29/01/2014 - v2.1 - Minor fix. Added ROUND statement to if...then
;                            inside the main loop. Ensures N_CHAN=1 
;                            returns a single channel.
;       10/03/2014 - v2.2 - Ensures INTEGER format output files
;
; NOTES:
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)

  COMPILE_OPT idl2, HIDDEN

; watch out for errors

  ON_ERROR, 0

; ---------------------------------------------------------
; Check arguments

  if (N_PARAMS() eq 0) then N_chan = 50

  if NOT KEYWORD_SET(E_min) then E_min = 200

  if NOT KEYWORD_SET(E_max) then E_max = 10000

  if ((E_max - E_min) le 0) THEN BEGIN
      PRINT, '** Invalid channel ranges in ENERGY_BANDS.'
      RETURN, !NULL
  endif

  if NOT KEYWORD_SET(dE) then dE = 5

  if NOT KEYWORD_SET(N_chan) then N_chan = 1

; set level of output to screen
   
  IF NOT KEYWORD_SET(chatter) THEN chatter = 0

  IF (CHATTER gt 0) THEN PRINT, '-- ENERGY_BANDS v2.0'

; ---------------------------------------------------------
; Main routine

; enumerate all the valid PI channels

  N_E = ROUND((E_max - E_min)/dE)

  x = INDGEN(N_E + 1, /LONG)*dE + E_min

; calculate the step size (in log(PI)?)

  if KEYWORD_SET(lin) then begin
    binwidth = (MAX(x) - MIN(x)) / N_chan
  endif else begin
    binwidth = (ALOG10(MAX(x)/MIN(x))) / N_chan
  endelse

; assign each PI channel to a bin
; First establish arrays to contain the boundaries.
; The BIN array lists the partition points for the binning, i.e. the
; index of the highest X point in each new (output) bin.

  bin = MAKE_ARRAY(N_chan+3, /LONG)
  u_x = MAKE_ARRAY(N_chan+3, /FLOAT)
  l_x = MAKE_ARRAY(N_chan+3, /FLOAT)

; Define the lower edge of the lowest (output) bin

  j = 0
  bin[j] = 0
  l_x[j] = x[0]
  minbin = 1
  no_in_bin = 0
  lastbin = x[0]
  for i = 0, N_E do BEGIN
    if KEYWORD_SET(lin) then begin
      nextbin = l_x[j] + binwidth
    endif else begin
      nextbin = l_x[j] * (10.0^binwidth) 
    endelse
    no_in_bin = no_in_bin + 1
    if (no_in_bin lt minbin) then CONTINUE               
    if (x[i] le ROUND(nextbin)) then CONTINUE      
    no_in_bin = 0
    u_x[j] = x[i] - 1
    l_x[j+1] = u_x[j] + 1
    j = j + 1
    bin[j] = i
   endfor

; Define the upper edge of the highest (output) bin

  u_x[j] = x[N_E]

; Remove unused elements of the arrays

  bin = bin[0:j]
  u_x = u_x[0:j]
  l_x = l_x[0:j]

; check for missing channels

  if (N_chan gt N_ELEMENTS(bin)) then begin
     N_chan = N_ELEMENTS(bin)         
     if (chatter ge 0) then begin
       PRINT, '-- N_CHAN reduced due to empty bands in ENERGY_BANDS:', N_chan
     ENDIF
  ENDIF

; re-arrange the output as a matrix

  data = TRANSPOSE(([[l_x],[u_x]]))

; save to a file

  if KEYWORD_SET(file) then begin
      WRITE_TABLE, ROUND(data), 'bands.dat', /INTEGER
  endif

; convert from PI channels to response 
; channels (PN only)
;
; for pn: PI = [0, 20479]
; and detector response = [0, 4095] = pi/5
;
; for MOS: PI = [0, 11999]
; and detector response = [0, 799] = pi/15

  if KEYWORD_SET(channels) then begin
      pi = data/5
      channels = MAKE_ARRAY(N_chan+2, 3, /INTEGER)
      channels[*,0] = [0, REFORM(pi[0,*]), REFORM(pi[1,N_chan-1])+1]
      channels[*,1] = [REFORM(pi[0,*])-1, REFORM(pi[1,N_chan-1]), 4095]
      channels[*,2] = channels[*,1] - channels[*,0] + 1
      WRITE_TABLE, TRANSPOSE(channels), 'channels.dat', /INTEGER
  endif

; ---------------------------------------------------------
; Return to user

  RETURN, data

END
