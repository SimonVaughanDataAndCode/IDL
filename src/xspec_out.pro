PRO xspec_out, X, DX, CHAN_LIST=chan_list, OTYPE=otype, $
               FRAC_ERR=frac_err, FILENAME=filename

; ----------------------------------------------------
;+
; NAME:
;       XSPEC_OUT
;
; PURPOSE:
;       Write data to file ready for XSPEC
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       XSPEC_OUT, X
;
; INPUTS:
;       X           - (vector) data to write
;
; OPTIONAL INPUTS:
;       DX          - (vector) errors on data to write
;       CHAN_LIST   - (array) the list of PI channel ranges
;       OTYPE       - (integer) output type (see details)
;       FRAC_ERR    - (float) fractional error to assign to data
;       FILENAME    - (string) name of output file
;
; OUTPUTS:
;       none
;
; OPTIONAL OUTPUTS:
;       none
;
; DETAILS:
;       Takes as input a vector of data (X), optionally with errors
;       (DX) and a listing of channel boundaries (CHAN_LIST) and
;       writes as output an ASCII file listing channels, data (and
;       errors) in a format that can easily be converted for use in
;       XSPEC. The final ASCII to XSPEC conversion can be made using
;       the FTOOLS ASCII2PHA or FLX2XSP depending on the output type
;       (see below).
;
;       If used, CHAN_LIST contains the upper and lower bounardies of
;       each channel (e.g. output from ENERGY_BANDS.PRO) as a [2,M] array,
;       where M is the number of channels used. If this is not
;       supplied we simply assign values using ENERGY_BANDS(M).
;
;       If errors (DX) are not supplied as input we assign a small
;       fractional error (0.01*x) simply so that XSPEC does not hit a
;       "divide by zero" problem. The size of the fractional error can
;       be changed using the FRAC_ERR keyword (default = 0.01).
;
;         OTYPE=0 (default):
;
;       The output file can be constructed in two different ways. The
;       default (OTYPE=0) constructs a file in terms of detector
;       channels, i.e. we list just the data:
;
;         x, dx
;
;       With channels running 0,1,2,...,M,M+1 where M is the number of
;       channels input in the data, X. We include two extra channels,
;       one at each end (0 and M+1), with data values of zero, to
;       allow for energies in  
;       the detector response that were not included in the
;       data. These can be ignored prior to fitting in XSPEC. Such a
;       file can be converted to XSPEC format using ASCII2PHA. One
;       then needs to generate a response matrix (RMF * ARF - e.g. 
;       using MARFRMF) that contains M+2 detector channels: channel 0
;       includes all unused low energies, channels 1-M as the 'good'
;       energies, and channel M+1 includes all the unused high
;       energies. This can be created using the FTOOL RBNRMF to rebin
;       a response matrix using the channel-PI listings output by
;       ENERGY_BANDS.PRO. 
;
;         OTYPE=1:
;
;       The second type of output assumes you have a diagonal and flat
;       response, i.e. the data are in 'physical' rather than
;       'detector' units. In this case the output is a list:
;
;         E_low, E_high, x*dE, dx*dE
;
;       where E_low and E_high are the lower and upper energies
;       included in each channel, dE = E_high - E_low is the channel
;       width and the data (and errors) are multiplied by dE to
;       account for the fact that XSPEC automatically divides by the
;       bin width dE. (It assumes the data are in ct or ct/s units,
;       but performs fitting in ct/s/keV and so divides out the bin
;       widths.) 
;
; EXAMPLE USAGE:
;
; HISTORY:
;       06/08/2010 - v1.0 - first working version
;       06/07/2013 - v1.1 - fix orientation of CHAN_LIST array
;
; PROCEDURES CALLED:
;       WRITE_TABLE, ENERGY_BANDS
;
; NOTES:
;
;
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)

  COMPILE_OPT idl2, HIDDEN

; watch out for errors

;  on_error,2

; ---------------------------------------------------------
; Check arguments

  N = N_ELEMENTS(x)

; add dummy errors if not supplied

  if NOT KEYWORD_SET(frac_err) then frac_err = 0.01
  if (N_ELEMENTS(dx) eq 1) then err = MAKE_ARRAY(N, VALUE=dx)
  if (N_ELEMENTS(dx) eq 0) then err = frac_err * ABS(x)
  if (N_ELEMENTS(dx) eq N) then err = dx

; set default output file type

  if NOT KEYWORD_SET(otype) then otype=0

; if CHAN_LIST is not set then provide a default

  if NOT KEYWORD_SET(chan_list) then begin
      chan_list = ENERGY_BANDS(N)
  endif

; set defailt filename

  if NOT KEYWORD_SET(filename) then filename = "xspec_out.dat"

; ---------------------------------------------------------
; Begin main routine

  if (otype eq 0) then begin

      y = TRANSPOSE(x)
      dy = TRANSPOSE(err)
      data = ([[0, 0], [y, dy], [0, 0]])
      WRITE_TABLE, data, filename

  endif

  if (otype eq 1) then begin

      e_lo = REFORM(chan_list[0, *])
      e_hi = [REFORM(chan_list[0, 1:(N-1)]), REFORM(chan_list[1, (N-1)])]
      dE = e_hi - e_lo
      y = x*dE
      dy = err*dE
      data = [[e_lo], [e_hi], [y], [dy]]
      WRITE_TABLE, TRANSPOSE(data), filename

  endif

; ---------------------------------------------------------
; Return to user

  RETURN

END
