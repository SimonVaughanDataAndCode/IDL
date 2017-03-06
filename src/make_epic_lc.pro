FUNCTION make_epic_lc, ROOTNAME, SCALE, PI_LIST=pi_list, $
                       ROOTPATH=rootpath, DT=dt, CHATTER=chatter, $
                       ERR=err, T0=t0, TIME=time, BKG=bkg, FRAC_I=frac_i, $
                       POIS=pois, PATT_LIM=patt_lim, SEED=seed, $
                       MINF=minf, B_ERR=b_err, T_CLIP_START=t_clip_start, $ 
                       T_CLIP_end=t_clip_end, NOINTERP=nointerp, $
                       NOBACK=noback, CCDNR=ccdnr

; ----------------------------------------------------------
;+
; NAME:
;       MAKE_EPIC_LC
;
; PURPOSE:
;       Extract time series from XMM/EPIC event data
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       net = MAKE_EPIC_LC("mydata", pi_list=pi_list)
;
; INPUTS:
;       ROOTNAME     - (string) root of generic filenames 
;       SCALE        - (array) list of source/background scalings
;
; OPTIONAL INPUTS:
;       PI_LIST      - (array) list of low-high PI channels
;       ROOTPATH     - (string) path to event files
;       DT           - (float) time bin size
;       CHATTER      - (integer) amount of output to screen
;       POIS         - (logical) randomise counts for interpolated bins?
;       PATT_LIM     - (array) set lower and upper PATTERN ranges
;       MINF         - (float) minimum fraction exposure to allow
;       T_CLIP_START - (float) length of time (sec) to 'clip' from the
;                        start of (the pn) GTI
;       T_CLIP_END   - (float) ...end of the GTI
;       NOINTERP     - (logical) switch off interpolation
;       NOBACK       - (logical) switch off backound subtraction?
;       CCDNR        - (two integers) CCD numbers (pn, MOS)
;
; OUTPUTS:
;       NET          - (array) combined, net count rate
;
; OPTIONAL OUTPUTS:
;       TIME         - (array) start time of each bin
;       ERR          - (array) error on each bin
;       T0           - (float) zero time for the time series
;       BKG          - (array) combined, scaled background count rate
;       B_ERR        - (array) error on background output
;       FRAC_I       - (array) fractional exposure per bin for each instrument
;       SCALE        - (array) src/bkg scale coefficients for pn, MOS1, MOS2
;
; DETAILS:
;       This function extracts a set of time series from EPIC event
;       files. The event files for each camera (PN, M1, M2) must have
;       the same root name. E.g.: ROOTNAME="m101" then we expect to
;       find files with names "m101_pn_src.ds" for the pn source data
;       (and likewise for M1, M2 and background files, "bkg"). These
;       should be found in the directory given by the optional
;       ROOTPATH keyword.
; 
;       For each camera the event file is loaded and the time and PI
;       channel of each event stored. Then, for each set of PI values
;       we select the events within this range and use HISTOGRAM to
;       produce a binned time series. We repeat this for the
;       background event file, and subtract the scaled background
;       (in counts/bin) from the source to get the net (background
;       subtracted) counts/bin. The time series are generated over
;       exactly identical time bins (so that all time series - from
;       source and background and different cameras are
;       synchronised properly). 

;       The scaling factor is supplied by the
;       keyword SCALE which should list three values: the
;       source/background scale factors for the pn, M1, M2. If any of
;       these are set to zero then we ignore that camera. The scaling
;       factors are the ratios of source/background good area
;       (e.g. taken from the BACKSCAL keywords of suitable FITS
;       files). 
;
;       For each energy channel range the net counts/bin is
;       calculated as follows:
;
;         net = sum_{i=0}^{2} SRC[i] - SCALE[i] * BKG[i]
;
;       where i=0,1,2 represent the three EPIC cameras (pn, M1,
;       M2). The error on the count/bin comes from:
;
;         err^2 = sum_{i=0}^{2} SRC[i] + SCALE[i]^2 * BKG[i]
;
;       Since the variance on the sum of Poisson distributed variables
;       is just the sum of their expectations. (The SCALE^2 weighting
;       comes from the scaling of the background regions.)
;       And what's RETURNed is the count rate net/dt and error
;       sqrt(err2)/dt. 
;
;       The PATT_LIM keyword is passed along to the routine LOAD_EVENTS.
;       The PATT_LIM input can be used to specify the lower and upper
;       event pattern ranges that are to be collected. Single pixel
;       events are PATTERN=0, double pixel events are 1-4, triples are
;       5-8 and quadruples are 9-12. By default we use all PATTERNS in
;       the input files (typically 0-12). But to use only singles and
;       doubles set PATT_LIM=[0,4], which uses only events for which 
;       PATTERN is in the range 0-4. TO use only single events set
;       PATT_LIM=[0,0]. 
;
;       (Updated 12/07/12) The T_CLIP_START/END keywords are used to 
;       'prune' the pn exposure. By default T_CLIP_START = 10 sec 
;       and T_CLIP_END = 100s, data are ignored from the first 10s 
;       after the start of the irst good time interval (GTI) and the 
;       final 100s before the end of the final GTI. For this purpose we 
;       use the GTI of the first event file to be loaded (usually the pn).
;       This removes problems caused by
;       the exposures occasionally starting/ending inside the boundary 
;       set by the GTI as tabulated. This seems to happen (by a few 10's 
;       of sec) for the pn (but not MOS as far as I know).
;         
;       The CHATTER keyword switches of the screen output showing the
;       data information. The settings are:
;        -1 - minimal output
;         0 - obsid level summary only
;         1 - info for each event file
;         2 - info for each FITS header 
;         
;       Gaps in sampling are filled by interpolation (using GTI_FIX).
;       The interpolated data is randomised using the appropriate
;       Poisson distribution if the POIS keyword it set.
;       If the NOINTERP keyoword is set the interpolated data
;       is set to NaN and treated as missing data.
;
;       ADDED 29/01/2014: You may specify the CCDNR for the pn and MOS
;       targets. Passed to LOAD_EVENTS. See documentation there.
;
; EXAMPLE USAGE:
;       pi_list = [[200,500,1000], [500,1000,10000]]
;       scale = [0.136071, 0.0, 1.57629]
;       x = MAKE_EPIC_LC("rev1721", scale, $
;                        rootpath="/data/49/sav2/xmm/ngc4051/events/", $
;                        dt=100.0, pi_list=pi_list, time=time, err=err)
;       plot, time, x[0,*], psym=10, /xstyle
;       plot_err, time, x[0,*], err[0,*]
;
; HISTORY:
;       04/06/2010 - v1.0 - first working version
;       16/06/2010 - v1.1 - added BKG output
;                            fixed bug with selecting global good time
;                            using bin_0:bin_N
;       09/07/2010 - v1.2 - added POIS keyword, which is passed to
;                            GTI_FIX, and PATT_LIM which is passed to
;                            LOAD_EVENTS. 
;       20/01/2011 - v1.3 - removed need for GTI filenames, following
;                            improvement in LOAD_EVENTS
;       21/12/2011 - v1.4 - added missing line err = err[*, mask] from
;                            last step
;                            added SCALE output
;       26/04/2012 - v1.5 - added SEED input/output keyword, for
;                            controlling randomisation step in
;                            GTI_FIX, useful for debugging.
;                            Minor fix in MASK computation, now used 
;                            "and pi_s le e_hi" where previous this was
;                            "lt" (exclusive) not "le" (inclusive).
;                            Added MINF input keyword. Use double precision 
;                            arrays for time series output and exposure 
;                            calculations.
;                            Include 1/FRAC rescaling in the calculation of 
;                            source, background and net count rate errors 
;                            (by adding src2_ij array to keep track of variance). 
;       02/07/2012 - v1.6 - added B_ERR (background error) output
;       11/07/2012 - v1.7 - bug fix - removed duplicated lines
;       13/07/2012 - v1.8 - Added T_CLIP_START and T_CLIP_END keywords 
;                            to allow finer control of acceptable time ranges.
;                            This was previously handled 'internally' by 
;                            LOAD_EVENTS.PRO but now can be controlled at a 
;                            'higher' level. Removed clipping for first and last
;                            data points (no longer needed). Added missing line 
;                            to trim FRAC_I in line with src_ij etc. using 
;                            bin_0:bin_N. Changed SILENT keyword for CHATTER
;                            allowing more control over output.
;       16/07/2012 - v1.9 - Added check for DT/FRAME_TIME incompatibility
;       19/07/2012 - v2.0 - Added check/warning for N_BIN < 2.
;                            Changed BIN_N to N_BIN-2 to 
;                            remove the final histogram bin, which is always empty
;                            (due to the way histogram produces NBINS too large)
;       23/07/2012 - v2.1 - Few changes to improve memory usage, i.e.
;                            no need to duplicate HIST array, reset the t_s 
;                            and t_b arrays after each j-loop.
;       04/07/2013 - v2.2 - Added NOINTERP and NOBACK keywords
;       29/01/2014 - v2.3 - Added CCDNR input keywords (passed to LOAD_EVENTS)
;
; PROCEDURES CALLED:
;          LOAD_EVENTS [FITS I/O], GTI_FIX [POIDEV]
;
; NOTES:
;
;-
; ----------------------------------------------------------
; options for compilation (recommended by RSI)

  COMPILE_OPT idl2  
  
; watch out for errors

  ON_ERROR, 0

; ---------------------------------------------------------
; Check arguments

  IF NOT KEYWORD_SET(rootname) THEN BEGIN
      PRINT,"** No filename root specified in MAKE_EPIC_LC"
      RETURN, !NULL
  ENDIF

; set the number of cameras to 3

  N_inst = 3

; check we have some SCALEs to use

  IF (N_ELEMENTS(scale) ne N_inst) THEN BEGIN
      PRINT,"** SCALE has wrong number of elements in MAKE_EPIC_LC"
      PRINT,"** Found", N_ELEMENTS(scale), "instead of", N_inst
      RETURN, !NULL
  ENDIF

; set the dT to default

  IF (N_ELEMENTS(dt) ne 1) THEN dt = 100.0D

; set the PI_LIST to default (200,10000) if not given on input

  IF (N_ELEMENTS(pi_list) eq 0) THEN BEGIN
      pi_list = MAKE_ARRAY(1, 2, /LONG)
      pi_list[0,0] = 200L
      pi_list[0,1] = 10000L
  ENDIF

  N_size = SIZE(pi_list, /DIMENSIONS)
  N_chan = N_size[0]

  IF NOT KEYWORD_SET(rootpath) THEN rootpath = ""

; set PATT_LIM (PATTERN selection) to default (0-4) if not given on
; input 

  IF (N_ELEMENTS(patt_lim) eq 0) THEN patt_lim = [0, 4]

; have we been given minf?

  IF (N_ELEMENTS(minf) ne 1) THEN minf = 0.3D

; set the GTI 'pruning' parameters if not set already

  IF NOT KEYWORD_SET(t_clip_start) THEN t_clip_start = 10.0D
  IF NOT KEYWORD_SET(t_clip_end) THEN t_clip_end = 100.0D
  
; set level of output to screen
   
  IF not KEYWORD_SET(chatter) THEN chatter = 0

; set other processing options

  if NOT KEYWORD_SET(nointerp) THEN nointerp = 0
  if NOT KEYWORD_SET(noback) THEN noback = 0

; set the CCDs to extract 04 for pn and 01 for MOS1/2 by default

  IF NOT KEYWORD_SET(ccdnr) THEN ccdnr = [4,1]

; ---------------------------------------------------------
; Main routine

; construct source filenames for each camera

  file_pn = rootname + "_pn_src.ds"
  file_m1 = rootname + "_M1_src.ds"
  file_m2 = rootname + "_M2_src.ds"
  file_src = [file_pn, file_m1, file_m2]

; construct background filenames for each camera

  bfile_pn = rootname + "_pn_bkg.ds"
  bfile_m1 = rootname + "_M1_bkg.ds"
  bfile_m2 = rootname + "_M2_bkg.ds"
  file_bkg = [bfile_pn, bfile_m1, bfile_m2]

; set absolute lower/upper channel bounds

  chan_lim = [MIN(pi_list), MAX(pi_list)]

; -------------------------------------
; load data from event files, for each camera
;  i = 0 --> pn
;  i = 1 --> M1
;  i = 2 --> M2

; COUNTER is set to zero to signify the first event file, which is used
; to define arrays needed during processing of subsequent event
; files. Once this is done we set COUNTER=1.

  counter = 0

  FOR i = 0, N_inst-1 DO BEGIN

; IF no data for camera i (indicated by scale=0) then skip to next one

      IF (scale[i] eq 0.0) THEN CONTINUE

; load the source event file for instrument i

      t_s = LOAD_EVENTS(file_src[i], fileroot=rootpath, t0=t0_i, $
                        chan_lim=chan_lim, gti=gti, pi=pi_s, $
                        CHATTER=chatter, PATT_LIM=patt_lim, $
                        FRAME_TIME=frame_time, CCDNR=ccdnr)

; check the frame time and time bin size are not incompatible

      binup = dt / frame_time
      IF (binup lt 1.0) THEN BEGIN
          PRINT, '** DT is smaller than FRAME_TIME in MAKE_EPIC_LC'
          RETURN, !NULL
      ENDIF
      IF (binup lt 10) THEN BEGIN
          leftover = dt/frame_time - FLOOR(dt/frame_time)
          IF (ABS(leftover) gt 0.01) THEN BEGIN
              PRINT, '** DT =',dt,' and FRAME_TIME =',frame_time,' are not well matched'
              PRINT, '** Better to pick DT an integer multiple of FRAME_TIME.'
              RETURN, !NULL
          ENDIF
      ENDIF
  
; If this is the first event file, use the start time of the first GTI
; as the zero point. Define the overall start/end times.

      IF (counter eq 0) THEN BEGIN
      
          t0 = t0_i 
          tstart = gti[0]
          tend = gti[-1]
          
          ; apply 'pruning' of exposure
         
          tstart = tstart + t_clip_start
          tend = tend - t_clip_end
          n_timebins = FLOOR((tend-tstart) / dt)
;          tend = tstart + (n_timebins * dt) 

          IF (n_timebins le 1) THEN BEGIN
              PRINT, '** Too few time bins in MAKE_EPIC_LC.'
              PRINT, '** Try different DT.'
              RETURN, !NULL
          ENDIF

      ENDIF ELSE BEGIN

; If this is not the first event file, correct the zero time
; of all time series to the 'global' start time defined previously

          toffset = (t0_i - t0)
          t_s = t_s + toffset
          gti = gti + toffset

      ENDELSE

; load the background event file for instrument i

      t_b = LOAD_EVENTS(file_bkg[i], fileroot=rootpath, t0=bt0_i, $
                        chan_lim=chan_lim, pi=pi_b, $
                        CHATTER=chatter, PATT_LIM=patt_lim)

; correct for different zero times of source and background

      toffset = (bt0_i - t0)
      t_b = t_b + toffset

; loop over each PI (energy) range

      FOR j = 0, N_chan-1 do BEGIN

          ; select low/high PI channels

          e_lo = pi_list[j, 0]
          e_hi = pi_list[j, 1]

          ; select only source counts within channel bounds
 
          mask = WHERE(pi_s ge e_lo and pi_s le e_hi, count)
          t_j = t_s[mask]
          
          ; use HISTOGRAM to make a time series, i.e. counts per time bin

          hist = HISTOGRAM(t_j, BINSIZE=dt, MIN=tstart, MAX=tend, LOCATIONS=t_ij)

          ; store the results in a suitable array

          IF (counter eq 0) THEN BEGIN
              N_bin = N_ELEMENTS(t_ij)
              src_ij = MAKE_ARRAY(N_inst, N_chan, N_bin, /DOUBLE)   ; src rate
              src2_ij = MAKE_ARRAY(N_inst, N_chan, N_bin, /DOUBLE)  ; src var
              frac_i = MAKE_ARRAY(N_inst, N_bin, /DOUBLE)           ; frac expo
              bkg_ij = src_ij                                       ; bkg rate
              bkg2_ij = src_ij                                      ; bkg var

              ; select the range of 'good' data for this event list
              ; N_bin - 2 is use because the last bin produduced by the
              ; HISTOGRAM command is empty. 

              bin_0 = 0                                             ; first good bin
              bin_N = N_bin - 2                                     ; last good bin
          ENDIF

          ; correct for dropouts in the GTI
          ; make sure to use double prec.

 ;         x = DOUBLE(hist)  -- removed 23/07/2012
          src_ij[i, j, *] = GTI_FIX(t_ij, hist, gti, DT=dt, FRAC=frac, $
                                    POIS=KEYWORD_SET(pois), SEED=seed, $
                                    INTP_NUM=intp_num, MINF=minf, $
                                    WHICH_RESCALE=which_rescale, $
                                    WHICH_INTP=which_intp)

          ; compute the variance on the counts. Following Poisson stats.
          ; where data have been rescaled (using 1/FRAC) then apply this 
          ; rescaling to variance. Where data have been interpolated there 
          ; is no need if POIS keyword was used (since the result is Poisson)
          
          src2_ij[i, j, *] = src_ij[i, j, *]
          if (which_rescale ne !NULL) then begin
              src2_ij[i, j, which_rescale] = src_ij[i, j, which_rescale] / (frac[which_rescale])^2
          ENDIF  

          if KEYWORD_SET(nointerp) THEN BEGIN
            if (nointerp ne 0) THEN BEGIN
              if (which_intp ne !NULL) then begin
                 src_ij[i, j, which_intp] = !VALUES.F_NAN
                 src2_ij[i, j, which_intp] = !VALUES.F_NAN
              ENDIF  
            ENDIF
          ENDIF
          
          ; selection and binning of background events

          mask = WHERE(pi_b ge e_lo and pi_b le e_hi, count)
          t_j = t_b[mask]   

          ; use HISTOGRAM to make a time series, i.e. counts per time bin
      
          hist = HISTOGRAM(t_j, BINSIZE=dt, MIN=tstart, MAX=tend)

          ; correct for dropouts in the GTI
          ; make sure to use double prec.

;          x = DOUBLE(hist)  -- removed 23/07/2012
          bkg_fix = GTI_FIX(t_ij, hist, gti, DT=dt, MINF=minf, $
                            FRAC=frac, WHICH_RESCALE=which_rescale, $
                            POIS=KEYWORD_SET(pois), SEED=seed)

          ; apply background region scaling
          ; apply twice to bkg2 for error calculation
          ; (see error formula given above)

          bkg_ij[i, j, *] = bkg_fix * scale[i]

          ; compute the variance on the counts. Following Poisson stats.
 
          bkg2_ij[i, j, *] = bkg_fix * (scale[i])^2
          if (which_rescale ne !NULL) then begin
              bkg2_ij[i, j, which_rescale] = bkg_ij[i, j, which_rescale] / (frac[which_rescale])^2
          ENDIF  

          if KEYWORD_SET(nointerp) THEN BEGIN
            if (nointerp ne 0) THEN BEGIN
              if (which_intp ne !NULL) then begin
                bkg_ij[i, j, which_intp] = !VALUES.F_NAN
                bkg2_ij[i, j, which_intp] = !VALUES.F_NAN
              ENDIF  
            ENDIF
          ENDIF
          
          ; now we have processed the first dataset, 
          ; change the COUNTER to 1

          counter = 1

      ENDFOR

; free memory

      t_s = 0                 
      t_b = 0

; store the FRAC vector

      frac_i[i, *] = frac

; find first/last really 'good' bins
;      
      bin_0_i = MIN(WHERE(frac ge minf, count))
      bin_N_i = MAX(WHERE(frac ge minf, count)) 
;
; choose the first/last really 'good' bins over all cameras
; i.e. remove time bins at the start/end for which one
; of the (selected) cameras is underexposed.

      bin_0 = MAX([bin_0, bin_0_i])
      bin_N = MIN([bin_N, bin_N_i])

; end of main loop; next camera...

  ENDFOR

; -------------------------------------
; merge the time series where all good data overlap

  src_ij = src_ij[*, *, bin_0:bin_N]
  src2_ij = src2_ij[*, *, bin_0:bin_N]
  bkg_ij = bkg_ij[*, *, bin_0:bin_N]
  bkg2_ij = bkg2_ij[*, *, bin_0:bin_N]
  t_ij = t_ij[bin_0:bin_N]
  frac_i = frac_i[*, bin_0:bin_N]
  N_bin = N_ELEMENTS(t_ij)

; -------------------------------------
; combine all EPIC cameras within each energy band

  net = MAKE_ARRAY(N_chan, N_bin, /DOUBLE)
  err = MAKE_ARRAY(N_chan, N_bin, /DOUBLE)
  bkg = MAKE_ARRAY(N_chan, N_bin, /DOUBLE)
  b_err = MAKE_ARRAY(N_chan, N_bin, /DOUBLE)

  FOR j=0, N_chan-1 DO BEGIN

; sum source counts over all cameras
; and sum the variances

      src = REFORM( TOTAL(src_ij[*, j, *], 1, /NAN) )
      src2 = REFORM( TOTAL(src2_ij[*, j, *], 1, /NAN) )

; sum the scaled background counts over all cameras
; and sum the variances

      bg = REFORM( TOTAL(bkg_ij[*, j, *], 1, /NAN) )
      bg2 = REFORM( TOTAL(bkg2_ij[*, j, *], 1, /NAN) )

; find the net (= src - scale * bkg) counts per bin
; calculate the 'error'  on net count/bin (sqrt sum of variance)

      IF (NOBACK eq 0) THEN BEGIN
        net[j,*] = src - bg
        err[j,*] = SQRT( src2 + bg2 ) 
      ENDIF ELSE BEGIN
        net[j,*] = src    
        err[j,*] = SQRT( src2 ) 
      ENDELSE

; (scaled) background data

      bkg[j,*] = bg
      b_err[j,*] = SQRT( bg2 )

  ENDFOR

; convert from count/bin to count/sec 

  net = TEMPORARY(net) / DOUBLE(dt)
  err = TEMPORARY(err) / DOUBLE(dt)
  bkg = TEMPORARY(bkg) / DOUBLE(dt)
  b_err = TEMPORARY(b_err) / DOUBLE(dt)
  time = t_ij

; ---------------------------------------------------------
; Return to user

  RETURN, net

END
