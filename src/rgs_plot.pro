PRO rgs_PLOT, OBS=obs, FILENAME=filename, CENT=cent, VLIM=vlim, $
              Z=z, PS=ps, DV=dv, AREA=area, CHATTER=chatter, $
              BINV=binv, SPEC=net_biny, ERR_SPEC=err, NORAND=norand

; ----------------------------------------------------------
;+
; NAME:
;       RGS_PLOT
;
; Purposex:
;      Plot RGS data in velocity space
;
; AUTHOR:
;      Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;      rgs_PLOT, obs=[1, 2, 3]
;
; INPUTS:
;      obs      - (array) list of which observations to process, as listed
;                         in text files in the current directory.
;
; OPTIONAL INPUTS:
;      file     - (string) name of file listing RGS files (default 'rgsfiles.txt')
;      cent     - rest-frame wavelengths of lines to plot (Angstrom)
;      vlim     - velocity limits to plot (in km/s), e.g. [-5000, 2000]
;      z        - redshift of source (default = 0.0)
;      ps       - TRUE/FALSE: produce a PostScript file output (velocity.ps)
;      dv       - (scalar) velocity bin width (km/s)
;      area     - (logical) if set, normalise by effective area
;      CHATTER  - (integer) more or less feedback to screen?
;      norand   - (logical) switch off wavelength randomisation?
;
; OUTPUTS:
;      velocity.qdp - QDP file showing final results
;
; OPTIONAL OUTPUTS:
;      none
;
; DETAILS:
;      Produce a 'composite' line profile in velocity space from raw
;      RGS spectral files (FITS). The input is a set of files for the
;      source spectra, the background spectra, and response
;      matricies. Each source spectral file should have a
;      corresponding response (and background if used), and should be
;      the non-background subtracted RGS spectrum in high resolution
;      (e.g. 3400 channels). If there are multiple spectra, e.g. one
;      for each of RGS1 and RGS2, or multiple observations, these may
;      be combined. For example:
; 
;        rgs_PLOT, obs=[1, 2], cent=[16.006], z=0.002336, vlim=[-7000, 7000]	
;
;      This is combine the first and second observations as listed in 
;      the file 'rgsfiles.txt' in the current working directory.
;      
;      We cannot do the more simple trick of summing the counts in
;      each channel over all spectra because the channel-energy
;      conversion for each spectrum can be (and usually is) different
;      for each RGS spectrum.  A similar problem applies to the
;      BACKSCAL and QUALITY information, which is channel
;      dependent. The quality flag shows whether the channel is good
;      or not (e.g. bad due to bad column or missing chip) and this
;      multiplied by the overall exposure gives an effective exposure
;      per channel. The BACKSCAL values are also per channel. As the
;      channel boundaries (in wavelengths) are different from spectrum
;      to spectrum, these cannot simply be averaged by averaging the
;      data in each channel. What we do is re-sample the exposure and
;      BACKSCAL vectors onto a much finer grid of wavelengths
;      (specified if necesaary using the XXX keywords) that is the
;      same for all spectral files.
;
;      Each spectrum (and background) is 'unbinned' using the
;      RGS_UNBIN_SPEC function. Briefly, this randomised the observed counts within
;      the wavelength range of each spectral channel, to produce a
;      list of unbinned 'observed' wavelengths for each count. The
;      counts for all the source spectra are then combined (likewise
;      the background if used). The channel-wavelength conversion
;      comes from the response files. 
;
;      The photon wavelengths are then de-redshifted. For each
;      rest wavelength specified with the CENT input the photon
;      wavelengths are converted to velocity shifts (conversion
;      assuming v << c), and those counts in the velocity range
;      specified by the VLIM input are collected. This is repeated for
;      each line, and for the source and background spectra
;      individually. The velocities of each photon from each line are
;      then collected together (likewise for background if used) and
;      are binned in even velocity units (specified by the DV
;      input). For example, one may wish to examine the velocity range
;      -3000,+3000 km/s around the lines at 16.006. This is done by
;      setting cent=[16.006, 18.967] and vlim=[-3000,3000] and then
;      specifiying e.g. dv=300 (km/s) for the final binning of the
;      velocity spectrum. If a background is specified the background
;      is subtracted from the source spectrum. The resulting spectrum
;      is in units of counts per velocity bin, with velocity in units
;      of km/s. Errors come from propagating the usual Poisson
;      variances (e.g. error = sqrt(N1+N2) for N1 source counts and N2
;      background counts). The result is a composite of all the counts
;      for each spectrum and for each line (shifted into velocity
;      shifts relative to that line). 
;
;      The source and background region sizes (from the BACKSCAL
;      vector) and the exposure time per channel are then
;      combined. The area is calculated for each channel by an
;      exposure time weighted average of the areas of the contributing
;      spectral files. The ratio of source and background regions
;      provides the wavelength dependent background scaling factor
;      required for the background subtraction.
;
;      An example of the file 'rgsfiles.txt' is:
;      
;        rev0278_rgs12_o1_src.pi rev0278_rgs12_o1_bkg.pi rev0278_rgs12_o1.rsp
;        rev0830_rgs12_o1_src.pi rev0830_rgs12_o1_bkg.pi rev0830_rgs12_o1.rsp
;        rev1471_rgs12_o1_src.pi rev1471_rgs12_o1_bkg.pi rev1471_rgs12_o1.rsp
;
;      This lists source, background and response files in columns. If we 
;      wanted to examine only the first and second of these we set obs=[1,2].
;
; EXAMPLE USAGE:
;
;
; PROCEDURES CALLED:
;          SPEC_READFITS, RGS_UNBIN_SPEC, RGS_UNBIN, PS_OPEN, PS_CLOSE,
;          WRITE_TABLE, PLOT_ERR, READFITS, TBGET, SXPAR, SEQ
;
; HISTORY:
;      12/11/2009 - v0.1 - first working version
;      26/11/2009 - v0.2 - modified to work on raw (total, background) 
;                          fits files, multiple spectra, and combine
;                          the counts from each spectrum and each line
;                          then bin to produce a single "composite"
;                          velocity profile in counts/bin against
;                          velocity. 
;      22/12/2009 - v0.3 - extended to cope with differing
;                          wavelength-channel conversions in each
;                          spectrum (using the RGS_UNBIN function).
;      27/02/2012 - v0.4 - Fixed minor bug in de-redshift conversion
;      14/03/2012 - v0.5 - extended wavelength range to 5.5 - 38 A
;                          replaced v(z) with relativistic formula
;      07/04/2014 - v0.6 - Added printout of wavelength range(s)
;      04/08/2014 - v1.0 - replace input system with one based on a text file,
;                           by default 'rgsfiles.txt'
;      14/08/2014 - v1.1 - moved RGS_UNBIN and RGS_UNBIN_SPEC to separate
;                           files; Added AREA and CHATTER keywords; Added legend to
;                           diagnostic plot; use VELOCITY_SHIFT function to 
;                           compute wavelength-velocity transformations
;      20/08/2014 - v1.2 - replace RGS_UNBIN with more general REGRID function;
;                           added NORAND keyword 
;
; NOTES:
;
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)

  COMPILE_OPT idl2

; watch out for errors

  ON_ERROR, 3

; ----------------------------------------------------------

; define the level out output to screen

  IF NOT KEYWORD_SET(chatter) THEN chatter = 0

; set NORAND?

  IF NOT KEYWORD_SET(norand) THEN norand = !NULL

; Use a look-up table to map the numbers in the OBS input to to the
; files RGS files in the current directory

  CD, CURRENT=pwd
  PRINT, '-- RGS_PLOT'
  if (not KEYWORD_SET(filename)) then filename = 'rgsfiles.txt'
  PRINT, '-- Collecting files from:', pwd+PATH_SEP()+filename
  
; load the list of RGS data files

  filelist = READ_TABLE(filename, /TEXT)

; extract the relevant files from this list

  size_dat = SIZE(filelist, /DIMENSION)
  m = size_dat[0]
  if (m ne 3) then begin
    PRINT, '** There must be three columns in file: ', filename
    RETURN
  endif

  if (N_ELEMENTS(obs) eq 0) then obs=INDGEN(N_obs)+1
  N_obs = N_ELEMENTS(obs)
  src = filelist[0, obs-1]
  bkg = filelist[1, obs-1]
  resp = filelist[2, obs-1]
      
; test to see if the file exist
  
  FOR i = 0, N_obs-1 do begin    
    FOR j = 0, 2 DO BEGIN
      file = (FILE_SEARCH(filelist[j, obs[i]-1], /TEST_READ, /TEST_REGULAR))
      IF (file[0] eq '') THEN BEGIN
        PRINT,'** File not found:', filelist[j, obs[i]-1]
        RETURN
      ENDIF
    ENDFOR
  ENDFOR    
      
  PRINT, '-- Loading files:'
  FOR i = 0, N_obs-1 DO PRINT, '-- ', src[i]
 
; ----------------------------------------------------------
; Check arguments

; set the centroid if not specified

  if (N_ELEMENTS(cent) eq 0) then cent = 21.602

; number of lines to combine

  n_line = N_ELEMENTS(cent)

; set the velocity limit if not specified

  if (N_ELEMENTS(vlim) eq 0) then begin
      PRINT, '-- VLIM parameter not given in RGS_PLOT. Using default setting.'
      vlim = [-5000.0, 5000.0]
  endif
  if (N_ELEMENTS(vlim) eq 1) then begin
      PRINT, '-- Only one value of VLIM parameter given in RGS_PLOT.'
      PRINT, '-- Will use [-VLIM, +VLIM] as suitable range.'
      vlim = [-ABS(vlim), ABS(vlim)]
  endif
  if (vlim[0] gt vlim[1]) then vlim = vlim[1,0]

; set the redshift if not specified

  if (N_ELEMENTS(z) ne 1) then z = 0.0

; set the velocity bin width (km/s) if not specified

  if (N_ELEMENTS(dv) ne 1) then dv = 300.0

; ----------------------------------------------------------
; construct the output velocity grid 

  binv = SEQ(vlim[0], vlim[1], dv)
  n_bin = N_ELEMENTS(binv)

; ----------------------------------------------------------
; correct for redshift

  cent = cent * (1.0D + z)

; display velocity range(s) to be examined

  wave = MAKE_ARRAY(2, N_line)
  FOR k = 0, N_line-1 do begin
    wave[*, k] = VELOCITY_SHIFT(vlim, cent[k], /VELOCITY)   
    PRINT, '-- Wavelength range: ', STRTRIM(MIN(wave[*,k]), 2), ' - ', $
                                    STRTRIM(MAX(wave[*,k]), 2), ' A'
  ENDFOR

; ----------------------------------------------------------
; Loop over each spectrum file (FITS). For each one 'unbin' the
; spectrum using the RGS_UNBIN_SPEC function. The result is a
; list of wavelengths for each count in the original (binned)
; spectrum. 
;
; Also calculte the total exposure time, summed over all spectra,
; along with the exposure time per channel (exposure multiplied by
; quality flag), exposure time-weighted area. These are mapped onto a
; fine grid of wavelengths.

; contruct output wavelength grid, u

  u_min = 5.5
  u_max = 38.0
  u_range = (u_max - u_min)
  du = 0.002
  n_u = ROUND(u_range/du)
  u = INDGEN(n_u)*du + u_min

; loop over the source files

  totime = MAKE_ARRAY(n_u)
  src_backscal = MAKE_ARRAY(n_u)
  src_effarea = MAKE_ARRAY(n_u)
  src_spec = []

  src_t = 0.0
  FOR i = 0, N_obs-1 do begin
      PRINT,'-- Processing ', src[i]
      spec_i = RGS_UNBIN_SPEC(FILE=src[i], RESP=resp[i], CHATTER=chatter, $
                              W_L=w_l, W_U=w_u, TEXP=texp, QUAL=qual, $
                              BACKSCAL=backscal, EFF_AREA=effarea, NORAND=norand)
      src_spec = [src_spec, spec_i]
      qual = ROUND(1 - qual)
      exptime = qual * texp
      src_t = src_t + texp

;      u_time = RGS_UNBIN(w_l, w_u, exptime, u, du)
;      u_backscal = RGS_UNBIN(w_l, w_u, backscal, u, du)
;      u_effarea = RGS_UNBIN(w_l, w_u, effarea, u, du)

      u_time = REGRID([w_l, MAX(w_u)], exptime, [u, MAX(u)+du])
      u_backscal = REGRID([w_l, MAX(w_u)], backscal, [u, MAX(u)+du])
      u_effarea = REGRID([w_l, MAX(w_u)], effarea, [u, MAX(u)+du])

      src_backscal = src_backscal + u_backscal * u_time
      src_effarea = src_effarea + u_effarea * u_time
      totime = totime + u_time

  endfor

; Expsure time and (time weighted) area per channel for source

  src_exp = totime
  mask = WHERE(src_exp gt 0.0, COMPLEMENT=cmask, count)
  if (count GT 0) then begin
    src_backscal[mask] = src_backscal[mask] / src_exp[mask]
    src_effarea[mask] = src_effarea[mask] / src_exp[mask]
  endif
  
  if (count LT n_u) then begin
    src_backscal[cmask] = 0.0
    src_effarea[cmask] = 0.0
  endif

; ----------------------------------------------------------
; Do the same for the background 

  totime = MAKE_ARRAY(n_u)
  bkg_backscal = MAKE_ARRAY(n_u)
  bkg_spec = []

  bkg_t = 0.0
  for i = 0, N_obs-1 do begin
    PRINT,'-- Processing ',bkg[i]
    spec_i = RGS_UNBIN_SPEC(FILE=bkg[i], RESP=resp[i], CHATTER=chatter, $
                            W_L=w_l, W_U=w_u, TEXP=texp, QUAL=qual, $
                            BACKSCAL=backscal, NORAND=norand)
    bkg_spec = [bkg_spec,spec_i]
    qual = round(1 - qual)
    exptime = qual * texp
    bkg_t = bkg_t + texp
          
;    u_time = RGS_UNBIN(w_l, w_u, exptime, u, du)
;    u_backscal = RGS_UNBIN(w_l, w_u, backscal, u, du)
 
   u_time = REGRID([w_l, MAX(w_u)], exptime, [u, MAX(u)+du])
   u_backscal = REGRID([w_l, MAX(w_u)], backscal, [u, MAX(u)+du])

   bkg_backscal = bkg_backscal + u_backscal * u_time
   totime = totime + u_time
          
  endfor

; Expsure time and (time weighted) area per channel for background

  bkg_exp = totime
  mask = WHERE(bkg_exp gt 0.0, COMPLEMENT=cmask, count)
  if (count GT 0) then bkg_backscal[mask] = bkg_backscal[mask] / bkg_exp[mask]
  if (count LT n_u) then bkg_backscal[cmask] = 0.0

; The source/background region scaling factor, per channel

  scale = MAKE_ARRAY(N_ELEMENTS(src_backscal))
  mask = WHERE(bkg_backscal gt 0.0, COMPLEMENT=cmask)
  scale[mask] = src_backscal[mask] / bkg_backscal[mask]
  scale[cmask] = !VALUES.F_NAN

; ----------------------------------------------------------
; make a few quick diagnostic plots
; open the PS device if requested

  if KEYWORD_SET(ps) then PS_OPEN, "velocity.ps", /COLOUR, $
                          POSITION=[0.1,0.1,0.9,0.9]

; plot the count spectra

  dw = ABS(cent[0] - VELOCITY_SHIFT(dv, cent[0], /VELOCITY))
  hist = HISTOGRAM(src_spec, BINSIZE=dw, locations=binw)
  PLOT, binw, hist, PSYM=10, XTITLE="Wavelength (!6!sA!r!u!9 %!6!n)", $
    YTITLE="counts/bin", TITLE="Raw counts data"
  hist_bkg = HISTOGRAM(bkg_spec, BINSIZE=dw, LOCATIONS=binw_bkg)
  OPLOT, binw_bkg, hist_bkg, PSYM=10, COLOR='FF0000'x
  FOR k = 0, N_line-1 do begin
    xrange_k = [MIN(wave[*,k]), MAX(wave[*,k])]
    PLOT, binw, hist, PSYM=10, XTITLE="Wavelength (!6!sA!r!u!9 %!6!n)", $
      YTITLE="counts/bin", TITLE="Raw counts data", XRANGE=xrange_k
    OPLOT, binw_bkg, hist_bkg, PSYM=10, COLOR='FF0000'x
  ENDFOR

; plot the BACKSCALs

  yrange = [0, MAX([src_backscal, bkg_backscal])]
  PLOT, u, src_backscal, PSYM=10, XTITLE="Wavelength (!6!sA!r!u!9 %!6!n)", $
        YTITLE="Area", TITLE="src (bkg) BACKSCAL", YRANGE=yrange
  OPLOT, u, bkg_backscal, PSYM=10, COLOR='F09090'x
  FOR k = 0, N_line-1 do begin
    xrange_k = [MIN(wave[*,k]), MAX(wave[*,k])]
    PLOT, u, src_backscal, PSYM=10, XTITLE="Wavelength (!6!sA!r!u!9 %!6!n)", $
          YTITLE="Area", TITLE="src (bkg) BACKSCAL", XRANGE=xrange_k, YRANGE=yrange
    OPLOT, u, bkg_backscal, PSYM=10, COLOR='F09090'x
  ENDFOR

; plot the src/bkg BACKSCAL ratio

  PLOT, u, scale, PSYM=10, XTITLE="Wavelength (!6!sA!r!u!9 %!6!n)", $
        YTITLE="src/bkg area", TITLE="src/bkg BACKSCAL"
  FOR k = 0, N_line-1 do begin
    xrange_k = [MIN(wave[*,k]), MAX(wave[*,k])]
    PLOT, u, scale, PSYM=10, XTITLE="Wavelength (!6!sA!r!u!9 %!6!n)", $
          YTITLE="src/bkg area", TITLE="src/bkg BACKSCAL", XRANGE=xrange_k
  ENDFOR

; plot the effective area

  PLOT, u, src_effarea, PSYM=10, XTITLE="Wavelength (!6!sA!r!u!9 %!6!n)", $
        YTITLE="Area (cm^2)", TITLE="src effective area"
  FOR k = 0, N_line-1 do begin
    xrange_k = [MIN(wave[*,k]), MAX(wave[*,k])]
    PLOT, u, src_effarea, PSYM=10, XTITLE="Wavelength (!6!sA!r!u!9 %!6!n)", $
          YTITLE="Area (cm^2)", TITLE="src effective area", XRANGE=xrange_k
  ENDFOR

; plot the exposure time (*qual)

  PLOT, u, src_exp, PSYM=10, XTITLE="Wavelength (!6!sA!r!u!9 %!6!n)", $
        YTITLE="Exposure", TITLE="Exposure time"
  FOR k = 0, N_line-1 do begin
    xrange_k = [MIN(wave[*,k]), MAX(wave[*,k])]
    PLOT, u, src_exp, PSYM=10, XTITLE="Wavelength (!6!sA!r!u!9 %!6!n)", $
          YTITLE="Exposure", TITLE="Exposure time", XRANGE=xrange_k
  ENDFOR

; ----------------------------------------------------------
; loop over each line to consider
; For each one collect the counts within vlim[0] <= v <= vlim[1]
; for both source and background. Once counts are assigned to 
; a line they are removed from the list of counts to avoid
; double counting.

  src_n = []
  bkg_n = []

  for i = 0, n_line-1 do begin

; calcuate Doppler shift velocities

     src_v = VELOCITY_SHIFT(src_spec, cent[i])
     bkg_v = VELOCITY_SHIFT(bkg_spec, cent[i])

; extract src and bkg data data from velocities in range 
; vlim[0] - vlim[1] only

     mask = WHERE(src_v LE vlim[1] and src_v GE vlim[0], count, $
                  COMPLEMENT=unmask, /NULL)
     IF (N_ELEMENTS(mask) gt 0) then src_i = src_v[mask]
     IF (N_ELEMENTS(unmask) gt 0) then src_v = src_v[unmask]

     mask = WHERE(bkg_v LE vlim[1] and bkg_v GE vlim[0], count, $
                  COMPLEMENT=unmask, /NULL)
     IF (N_ELEMENTS(mask) gt 0) then bkg_i = bkg_v[mask]
     IF (N_ELEMENTS(unmask) gt 0) then bkg_v = bkg_v[unmask]

; include the current line in the composite array

     src_n = [src_n,src_i]
     bkg_n = [bkg_n,bkg_i]

  endfor

; sort counts from all lines into velocity order

  indx = SORT(src_n)
  src_n = src_n[indx]
  indx = SORT(bkg_n)
  bkg_n = bkg_n[indx]
 
; ----------------------------------------------------------
; Bin the counts by velocity to get a binned spectrum
; in ct/bin units

  src_biny = HISTOGRAM(src_n, BINSIZE=dv, MIN=binv[0], NBINS=n_bin)
  bkg_biny = HISTOGRAM(bkg_n, BINSIZE=dv, MIN=binv[0], NBINS=n_bin)
  yrange = [0.0, 1.2*MAX(src_biny)]

; plot the histograms to check quality

  title = 'spectra before exposure corrections'
  PLOT, binv, src_biny, PSYM=10, XRANGE=vlim, /XSTYLE, YRANGE=yrange, $ 
    /YSTYLE, XTITLE="Velocity (km s!U-1!N)", YTITLE="counts/bin", TITLE=title, $
    /NODATA
  PLOT_HIST, src_biny, X0=binv, X1=binv+dv, /NOPLOT, /NOLINES
  PLOT_HIST, bkg_biny, X0=binv, X1=binv+dv, COLOR='FF0000'x, /NOPLOT, /NOLINES

; ----------------------------------------------------------
; loop over each line to consider. Find channels that 
; correspond to vlim[0] <= v <= vlim[1] for each line.
; Calculate the mean exposure time and src/bkg area within
; each velocity bin. (The mean is the average over all channels, from
; all lines, contributing to a given velocity bin. Where only a
; faction of a channel falls with a given velocity bin, the
; contribution to the mean is weighted accordingly.)

  src_exp_vbin = MAKE_ARRAY(n_bin)
  src_backscal_vbin = MAKE_ARRAY(n_bin)
  src_effarea_vbin = MAKE_ARRAY(n_bin)
  bkg_exp_vbin = MAKE_ARRAY(n_bin)
  bkg_backscal_vbin = MAKE_ARRAY(n_bin)

  u_l = u 
  u_u = (u+du) 

  for i = 0, n_line-1 do begin

; convert channel boundary wavelengths to velocity in units (km/s)
; NB: because the RGS channels are ordered in decreasing energy,
; i.e. increasing wavelength, the channel-velocity conversion
; produces velocities in acending order.

      v_l = VELOCITY_SHIFT(u_l, cent[i])
      v_u = VELOCITY_SHIFT(u_u, cent[i])

; inner loop: over each velocity bin within each line

      for j = 0, n_bin-1 do begin

; extract channels from velocities in range of Jth velocity bin

          lo_v = binv[j]
          hi_v = binv[j] + dv
          mask = WHERE(v_u GT lo_v and v_l LT hi_v, count)

; how much of the first and last of these channels are actually
; included in the velocity range of Jth velocity bin?

          if (count EQ 0) then CONTINUE

          ch_1 = mask[0]
          frac_1 = (v_u[ch_1] - lo_v) / (v_u[ch_1] - v_l[ch_1])

          ch_N = mask[count-1]
          frac_N = (hi_v - v_l[ch_N]) / (v_u[ch_N] - v_l[ch_N])

          n_channels = count - 2 + frac_1 + frac_N
 
          src_exp_vbin_i = TOTAL(src_exp[mask]) $
            - (1-frac_1)*src_exp[ch_1] $
            - (1-frac_N)*src_exp[ch_N] 

          src_backscal_vbin_i = TOTAL(src_backscal[mask]) $
            - (1-frac_1)*src_backscal[ch_1] $
            - (1-frac_N)*src_backscal[ch_N] 

          src_effarea_vbin_i = TOTAL(src_effarea[mask]) $
            - (1-frac_1)*src_effarea[ch_1] $ 
            - (1-frac_N)*src_effarea[ch_N] 

; XXX need to average properly over lines (i...) XXX

          src_exp_vbin[j]  = src_exp_vbin[j] + (src_exp_vbin_i / n_channels)
          src_backscal_vbin[j] = src_backscal_vbin[j] + src_backscal_vbin_i
          src_effarea_vbin[j] = src_effarea_vbin[j] + src_effarea_vbin_i

          bkg_exp_vbin_i = TOTAL(bkg_exp[mask]) $
            - (1-frac_1)*bkg_exp[ch_1] $
            - (1-frac_N)*bkg_exp[ch_N] 

          bkg_backscal_vbin_i = TOTAL(bkg_backscal[mask]) $
            - (1-frac_1)*bkg_backscal[ch_1] $
            - (1-frac_N)*bkg_backscal[ch_N] 

          bkg_exp_vbin[j]  = bkg_exp_vbin[j] + (bkg_exp_vbin_i / n_channels)
          bkg_backscal_vbin[j] = bkg_backscal_vbin[j] + bkg_backscal_vbin_i
       endfor
  endfor

  mask = WHERE(src_exp_vbin eq 0, count)
  if (count gt 0) then src_exp_vbin[mask] = !VALUES.F_NAN
  mask = WHERE(src_effarea_vbin eq 0, count)
  if (count gt 0) then src_effarea_vbin[mask] = !VALUES.F_NAN
  mask = WHERE(src_backscal_vbin eq 0, count)
  if (count gt 0) then src_backscal_vbin[mask] = !VALUES.F_NAN
  
; ----------------------------------------------------------
; convert to count rate and then subtract the (area scale corrected)
; background spectrum 
; Calculate the error by propogating the Poisson [sqrt(N)] error on the
; scaled source and background counts

  net_biny = MAKE_ARRAY(n_bin)
  bkg_rate = MAKE_ARRAY(n_bin)
  err      = MAKE_ARRAY(n_bin)
  backscal = MAKE_ARRAY(n_bin)

; If we are using the AREA keyword to divide out the effective area
; we do not need to also correct for channel variations in exposure
; (this would result in an over-correction of the data)

;  if KEYWORD_SET(area) then begin
;    src_exp_vbin[*] = MAX(src_exp_vbin)
;    bkg_exp_vbin[*] = MAX(bkg_exp_vbin)
;  endif

  mask = WHERE(src_exp_vbin GT 0 and bkg_backscal_vbin GT 0 and bkg_exp_vbin GT 0)

  backscal[mask] = src_backscal_vbin[mask] / bkg_backscal_vbin[mask]
  bkg_rate[mask] = bkg_biny[mask]/bkg_exp_vbin[mask] * backscal[mask]
  net_biny[mask] = src_biny[mask]/src_exp_vbin[mask] - bkg_rate[mask]
  err[mask] = SQRT(src_biny[mask]/src_exp_vbin[mask]^2 + $
                   bkg_biny[mask]/bkg_exp_vbin[mask]^2 * backscal[mask]^2)

; ----------------------------------------------------------
; apply effective area scaling if requested
; output will then be ct/s/cm^2/bin

  IF NOT KEYWORD_SET(area) then src_effarea_vbin[*] = 1.0
 
  mask = WHERE(src_backscal_vbin gt 0, count, COMPLEMENT=unmask)

  src_biny = FLOAT(src_biny)

  net_biny[mask] = net_biny[mask] / src_backscal_vbin[mask] / src_effarea_vbin[mask]
  net_biny[unmask] = 0.0
  src_biny[mask] = src_biny[mask] / src_backscal_vbin[mask] / src_effarea_vbin[mask]
  src_biny[unmask] = 0.0
  bkg_rate[mask] = bkg_rate[mask] / src_backscal_vbin[mask] / src_effarea_vbin[mask]
  bkg_rate[unmask] = 0.0
  err[mask] = err[mask] / src_backscal_vbin[mask] / src_effarea_vbin[mask]
  err[unmask] = 0.0
 
; ----------------------------------------------------------
; plot the resulting histogram(s)

  title = "Wavelengths"
  for i = 0, n_line-1 do begin
      title = title + STRING(cent[i]/(1.0+z), FORMAT='(f7.3)')
      title = title +"!6!sA!r!u!9 %!6!n "
  endfor

  yrange = [0.0, 1.2*MAX(src_biny)]

  PLOT, binv, src_biny, PSYM=10, XRANGE=vlim, /XSTYLE, YRANGE=yrange, $ 
    /YSTYLE, XTITLE="Velocity (km s!U-1!N)", YTITLE="counts/bin", TITLE=title
  OPLOT, binv, bkg_rate, PSYM=10, LINESTYLE=3
    
  yrange = [0.0, 1.2*MAX(net_biny)]

  ytitle="counts/s/bin"
  if KEYWORD_SET(area) then ytitle="counts/s/bin/cm!U2!N"
  PLOT, binv, net_biny, PSYM=10, XRANGE=vlim, /XSTYLE, YRANGE=yrange, $ 
    /YSTYLE, XTITLE="Velocity (km s!U-1!N)", YTITLE=ytitle, $
    TITLE=title
    
  PLOT_ERR, binv, net_biny, err

; alternative x-axis

;  AXIS, XAXIS=1, XTICKC=xv, XTICKNAME=xlab

; plot the exposure time per channel from zero to max

  y = src_exp_vbin/MAX(src_exp_vbin)*yrange[1]*0.975
  OPLOT, binv, y, PSYM=10, THICK=10, COLOR='FF0000'x

; plot the area scaling per channel from zero to max

  y = src_backscal_vbin/MAX(src_backscal_vbin)*yrange[1]*0.95
  OPLOT, binv, y, PSYM=10, THICK=10, COLOR='00FF00'x

; plot the effective area  per channel from zero to max

  y = src_effarea_vbin/MAX(src_effarea_vbin)*yrange[1]*0.95
  OPLOT, binv, y, PSYM=10, THICK=10, COLOR='0000FF'x

; add labels

  xspan = vlim[1] - vlim[0]
  yspan = yrange[1] - yrange[0]
  dxx = [vlim[1] - xspan*0.2, vlim[1] - xspan*0.15]
  dyy = [yrange[1] - yspan*0.2, yrange[1] - yspan*0.2]
  OPLOT, dxx, dyy, THICK=10, COLOR='00FF00'x
  XYOUTS, dxx[1]*1.01, dyy[1], 'backscal'

  dyy = [yrange[1] - yspan*0.3, yrange[1] - yspan*0.3]
  OPLOT, dxx, dyy, THICK=10, COLOR='FF0000'x
  XYOUTS, dxx[1]*1.01, dyy[1], 'exposure'

  dyy = [yrange[1] - yspan*0.1, yrange[1] - yspan*0.1]
  OPLOT, dxx, dyy, THICK=10, COLOR='0000FF'x
  XYOUTS, dxx[1]*1.01, dyy[1], 'area'

  OPLOT, binv, bkg_rate, PSYM=10, LINESTYLE=3
  dyy = [yrange[1] - yspan*0.4, yrange[1] - yspan*0.4]
  OPLOT, dxx, dyy, LINESTYLE=3
  XYOUTS, dxx[1]*1.01, dyy[1], 'bkg'

; close the PS device if requested

  if KEYWORD_SET(ps) then PS_CLOSE

; ----------------------------------------------------------
; Make a 'pretty' plot

; open the PS device if requested

  if KEYWORD_SET(ps) then PS_OPEN, "velocity-pretty.ps", $
                          CHARSIZE=1.4, POSITION=[0.1,0.1,0.9,0.9]

  title = "Wavelengths"
  for i = 0, n_line-1 do begin
      title = title + STRING(cent[i]/(1.0+z), FORMAT='(f7.3)')
      title = title + "!6!sA!r!u!9 %!6!n "
  endfor

  yrange = [0.0, 1.2*MAX(net_biny)]

  PLOT, binv, net_biny, PSYM=10, XRANGE=vlim, /XSTYLE, YRANGE=yrange, $ 
    /YSTYLE, XTITLE="Velocity (km s!U-1!N)", YTITLE="counts s!U-1!N bin!U-1!N", TITLE=title
    
  PLOT_ERR, binv, net_biny, err

  if KEYWORD_SET(bkg) then begin
      OPLOT, binv, bkg_rate, PSYM=10, LINESTYLE=3
  endif

; close the PS device if requested

  if KEYWORD_SET(ps) then PS_CLOSE

; ----------------------------------------------------------
; open a file for the output
; and copy the QDP header information

  OPENW, lun, "velocity.qdp", /GET_LUN
    PRINTF, lun, "READ SERR 2"
    PRINTF, lun, "!    velocity            flux          error"
    PRINTF, lun, "!       (km/s)       (counts)       (counts)"
  FREE_LUN, lun

; store the output in a file
; after clipping to required velocity range
; note, the velocity is the CENTRE of the bin

  v = binv + 0.5*dv
  y = net_biny
  dy = err

  data = TRANSPOSE([[v],[y],[dy]])

  WRITE_TABLE, data, "velocity.dat"

  IF (!VERSION.OS eq 'Win32' or !VERSION.OS eq 'Win64') then begin
    PRINT, '-- output in velocity.dat'
  ENDIF ELSE BEGIN
    SPAWN, 'cat velocity.dat >> velocity.qdp'
    PRINT, '-- output in velocity.qdp'
  ENDELSE

; ----------------------------------------------------------

END
