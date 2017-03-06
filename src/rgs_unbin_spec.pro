FUNCTION RGS_UNBIN_SPEC, FILE=file, RESP=resp, CHATTER=chatter, TEXP=texp, $
                     QUAL=qual, BACKSCAL=backscal, EFF_AREA=spcar, $
                     W_U=w_u, W_L=w_l, NORAND=norand

; ----------------------------------------------------------
;+
; NAME:
;       RGS_UNBIN_SPEC
;
; PURPOSE:
;      Load and RGS spectrum from a FITS file and 'unbin' 
;
; AUTHOR:
;      Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;      spec = RGS_UNBIN_SPEC(src, resp)
;
; INPUTS:
;      file     - names of FITS file of spectra (full path)
;      resp     - name of FITS file of response matrix (full path)
;      norand   - (logical) switch off wavelength randomisation?
;
; OPTIONAL INPUTS:
;      chatter  - (integer) level of on-screen output (default=0)
;
; OUTPUTS:
;      spec     - wavelengths of each photon in the spectrum (Angstrom)
;
; OPTIONAL OUTPUTS:
;      qual     - quality flag (0=good, 1=bad) for each channel
;      backscal - BACKSCAL value for each channel from the spectral file
;      w_l      - wavelength (Angstrom) for lower boundary of each channel
;      w_u      - wavelength (Angstrom) for upper boundary of each
;                 channel
;      texp     - Expsure time of the spectrum
;      eff_area - (vector) Effective area (cm^2) in each energy channel
;
; DETAILS:
;      Called by RGS_PLOT.
;      Performs an 'unbinning' of finely binned spectral data in a
;      FITS file ('src') given a response matrix ('resp'). Currently
;      the columns of the FITS file that are used are those suitable
;      for processing RGS data, but could be adjusted if needed. The
;      'unbinning' works by redistributing each count in the spectrum
;      within the wavelength range of the channel it was collected
;      in. The redistribution is random and uniform between w[i] and
;      w[i]+dw[i] where w[i] and dw[i] are the lower boundary and
;      width of the wavelength bin of the ith count in Angstroms. This
;      of course assumes the counts are distributed unformly within
;      each spectral channel, an assumption that should be true if the
;      input data are given with small bin widths. In the case of the
;      RGS we recommend using 3400 channels. 
;
;      Here's exactly what happens:
;       - load the data from the FITS file as channel, counts, quality
;       - load the corresponding response file FITS
;       - use channel-energy table to convert spectral channel to energy
;       - convert data wavelength: wavelength, bin width, counts, quality
;       - create 'unbinned' spectrum by assigning each count to a
;          wavelength uniformly distributed in [w, w+dw]
;       - keep track of 'quality' flag
;       
;
; EXAMPLE USAGE:
;          src = 'spec.FIT'
;          resp = 'resp.FIT'
;          spec = RGS_UNBIN_SPEC(file=src, resp=resp)
;          plot, HISTOGRAM(spec, BINSIZE=0.05), PSYM=10
;
; PROCEDURES CALLED:
;          READFITS, TBGET, SXPAR, RD_OGIP_RMF [MRDFITS]
;
;
; HISTORY:
;      26/11/2009 - v1.0 - first working version
;      15/12/2009 - v1.1 - added TEXP output
;      23/09/2013 - v2.0 - removed SPEC_READFITS call, 
;                           replaced with explicit use of FITS I/O.
;                           Needed to cope with real/fake spectral files.
;      04/10/2013 - v2.1 - added call to RD_OGIP_RMF to extract the
;                           effective area vs. channel curve.
;                           Use this to flag bad QUAL channels also.
;      12/08/2014 - v2.2 - added CHATTER keyword; added NORAND keyword
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

; define the level out output to screen

  IF NOT KEYWORD_SET(chatter) THEN chatter = 0

; check number of spectra to process

  if (N_ELEMENTS(file) ne 1) then begin
      PRINT, '** RGS_UNBIN_SPEC requies one file for FILE input.'
      RETURN, !NULL
  endif

; check number of response files to process

  if (N_ELEMENTS(resp) ne 1) then begin
      PRINT, '** RGS_UNBIN_SPEC requies one file for RESP input.'
      RETURN, !NULL
  endif

; ----------------------------------------------------------
; Main routine

; load the spectrum from FITS file

; is the file name defined?

  if (N_ELEMENTS(file) eq 0) then begin
      file = ''
      read,'-- Enter file name (ENTER to list current directory): ', file
      if (file eq '') then begin
          list = FINDFILE()
          PRINT, list
          READ, '-- Enter file name: ', file
      endif
  endif

; call READFITS to read the data from FITS file

  if (chatter lt 2) then begin
    filedata = READFITS(file, htab, EXTEN_NO=1, /SILENT)
  endif else begin
    filedata = READFITS(file, htab, EXTEN_NO=1)
  endelse
 
; check for errors from READFITS

  if (N_ELEMENTS(filedata) le 1) then begin
      if (filedata eq -1) then begin
          PRINT, '** Error calling READFITS in RGS_UNBIN_SPEC'
          PRINT, '** '+!ERROR_STATE.MSG
          RETURN, !NULL
      endif
  endif
  
; extract the exposure time

  if (chatter lt 2) then begin
    texp = SXPAR(htab, "EXPOSURE", /SILENT)
  endif else begin
   texp = SXPAR(htab, "EXPOSURE")
  endelse

; For a normal RGS spectrum the columns are:
;             1 = channel
;             2 = counts/bin
;             3 = quality (0=good, 1=bad)
;             4 = AREASCAL (0.0 - 1.0)
;             5 = BACKSCAL (arb. units)
; In this case we want the counts, quality and BACKSCAL only.
; But simulated spectrum (e.g. from XSPEC) have only
;             1 = channel
;             2 = counts/bin
;             3 = quality
;             4 = grouping
; And so we want columns 2 and 3 only and assume BACKSCAL=1

; extract the vectors of counts/bin, quality and BACKSCAL

  counts = TBGET(htab, filedata, 'COUNTS')
  qual = TBGET(htab, filedata, 'QUALITY')

  ncols = SXPAR(htab, 'TFIELDS')
  backscal = MAKE_ARRAY(N_ELEMENTS(counts), VALUE=1.0D)
  if (ncols ge 5) then backscal = TBGET(htab, filedata, 'BACKSCAL')
 
; Load the channel-energy table from the response matrix

; is the file name defined?

  if (N_ELEMENTS(resp) eq 0) then begin
      resp = ''
      read,'-- Enter RESP file name (ENTER to list current directory): ', resp
      if (resp eq '') then begin
          list = FINDFILE()
          PRINT, list
          READ, '-- Enter file name: ', resp
      endif
  endif

; call RD_OGIP_RMF to read the response matrix

  verbose = 0
  if (chatter gt 1) then verbose = 5
  rmf = RD_OGIP_RMF(resp, EFFAR=effar, EQVAR=eqvar, SPCAR=spcar, VERBOSE=verbose)
  qual = (qual gt 0) or (eqvar eq 0)

; channel boundaries (initially in keV units)

  e_l = rmf.EMN
  e_u = rmf.EMX

; convert channel boundary energies to wavelengths (Angstrom)

  apkev = 12.3984428D

  w_l = apkev/e_u
  w_u = apkev/e_l

; channel widths

  dw = w_u - w_l

; total number of counts and channels in the spectrum

  n = TOTAL(counts)
  n_chan = N_ELEMENTS(counts)
  n_chan2 = N_ELEMENTS(e_l)

; cumulative number of counts in channels <i

  cumsum = TOTAL(counts, /CUMULATIVE)

; display data parameters

  if (chatter gt 0) then begin
      PRINT,'-- Number of channels in data:     ', n_chan
      PRINT,'-- Number of channels in response: ', n_chan2
      PRINT,'-- Number of counts:               ', LONG(n)
      PRINT,'-- '
  endif

; prepare array for output

  spec = MAKE_ARRAY(n, /DOUBLE)

; form the 'unbinned' spectrum by assigning each count a wavelength
; uniformly distributed in [w_l[i], w_u[i]] where i is the channel
; that the count was collected in.

  for i=0, n_chan-1 do begin

      if (counts[i] eq 0) then continue

      rand = RANDOMU(seed, counts[i])

      if KEYWORD_SET(norand) then begin
        rand = (INDGEN(counts[i])+1) / FLOAT(counts[i]+1)
      endif

      j = cumsum[i] - counts[i] 
      k = j + counts[i] - 1
;      PRINT,i,j,k,counts[i], cumsum[i], n
      spec[j:k] = rand * dw[i] + w_l[i]

  endfor

; ----------------------------------------------------------
; return to calling routine

  RETURN, spec

END

