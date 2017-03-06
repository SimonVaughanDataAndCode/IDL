FUNCTION load_events, filename, fileroot=fileroot, t0=t0, texp=texp, $
    chan_lim=chan_lim, gti=gti, CHATTER=chatter, $
    pi=pi, PATT_LIM=patt_lim, FRAME_TIME=frame_time, CCDNR=ccdnr
    
; ----------------------------------------------------------
;+
; NAME:
;       LOAD_EVENTS
;
; PURPOSE:
;       Read XMM EPIC data from filtered events file (FITS)
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       t = LOAD_EVENTS("mydata.fits")
;
; INPUTS:
;       FILENAME  - (string) file name
;
; OPTIONAL INPUTS:
;       CHATTER      - (integer) amount of output to screen
;       CHAN_LIM     - (array) set lower and upper PI channel boundaries
;       PATT_LIM     - (array) set lower and upper PATTERN ranges
;       CCDNR        - (two integers) CCD numbers (pn, MOS)
;
; OUTPUTS:
;       T         - vector of event arrival times
;
; OPTIONAL OUTPUTS:
;       T0        - (double) time (MET sec) of first event
;       TEXP      - (float) Expsure time (from LIVETIME keyword)
;       GTI       - (array) two column table of good start/stop times
;       PI        - (vector) list of event PI (energy in eV)
;       FRAME_TIME - (float) frame time of event file (in sec)
;
; PROCEDURES CALLED:
;       READFITS, TBGET, SXPAR, FITS_INFO, GTIMERGE
;
; ADDITION DEPENENCIES:
;       GETTOK, VALID_NUM, TBINFO, TAG_EXIST, FITS_INFO, STRN,
;       MRD_SKIP, GTITRIM
;
; DETAILS:
;       Read columns of data from a FITS file. The input is a filtered
;       XMM EPIC event list. This output is a list of event arrival
;       times that can be used to generate a light curve. Optional
;       outputs include a list (same length as main output) of the
;       event PIs (energy in eV), which can be used to filter on
;       event energy.
;       A GTI (Good Time Interval) table is also read. This shows
;       the start and stop times of each GTI (to the nearest sec). The
;       MET (Mission Elapsed Time) in sec of the first GTI start is
;       assigned to T0 - the main output uses this as zero time, i.e. t[0]
;       = 0. The "live" exposure time is also passed via TEXP.
;
;       The CHATTER keyword switches of the screen output showing the
;       data information. The settings are:
;        -1 - minimal output
;         0 - obsid level summary only
;         1 - info for each event file
;         2 - info for each FITS header
;
;       The CHAN_LIM input can be used to specify the lower and upper
;       PI channel ranges that are to be collected. Events with PI
;       valuess outside this range will be ignored. (The limit is
;       inclusive of the end points.)
;
;       The PATT_LIM input can be used to specify the lower and upper
;       event pattern ranges that are to be collected. Single pixel
;       events are PATTERN=0, double pixel events are 1-4, triples are
;       5-8 and quadruples are 9-12. By default we use all PATTERNS in
;       the input files (typically 0-12). But to use only singles and
;       doubles set PATT_LIM=[0,4], which uses only events for which
;       PATTERN is in the range 0-4. TO use only single events set
;       PATT_LIM=[0,0].
;
;       GTI handling:
;       The routine automatically collects the GTI tables
;       from each event file. The GTIs are contained in extensions to
;       the FITS event list file. There should be one standard GTI
;       (STDGTInn) for each operating CCD. These are produce by running
;       the standard EPIC chain/proc SAS tasks to produce the event
;       list. We want the GTI only for the CCD containing our target.
;       For the pn this is assumed to be CCDNR=4 (counting 1-12) and
;       for both MOS this is assumed to be CCDNR=1 (counting 1-7).
;       Hence, we look for STDGTI04 (pn) and STDGTI01 (MOS) by default.
;       The FXPOSIT routine is used to find the extension number in each case.
;
;       If necessary you may specify the CCDNR for the pn and MOS
;       targets. E.g. if the source is at the aimpoint this should be
;       CCDNR=4 (running 1-12) for pn and CCDNR=1 (running 1-7) for
;       MOS1/2. Hence the default is CCDNR = [4,1]. If your source lies on
;       a different chip then you need to set CCDNR appropriately so that
;       the correct GTIs are applied.
;
;       There may be additional GTI tables in other extensions
;       if further (user-specified) time filtering has been applied
;       to the event file. These are called GTI0nnxx, where nn
;       is the CCDNR now starting from 00, and xx is a number starting
;       from 02. For some reason the nn used in GTI0nnxx is 1 less
;       than the nn used in STDGTInn. So GTI003yy and GTI000yy are the
;       extensions for the aimpoint chips of the pn and MOS,
;       respectively. Except in pn SW or TIMING mode, when CCDNR=0!
;       The number yy indicates how many levels of GTI filtering
;       have been applied.
;
;       Where there are GTI0nnyy tables available we load all those with
;       the right nn. Then all GTI tables are merged. This is done using
;       the GTIMERGE function to form a new GTI table from the intersection
;       of all GTIs. I.e. to count as a GTI a given time needs to be "good"
;       in all GTI tables.
;
;       ADDED 30/06/2013: Automatically detecting TIMING data and use
;       appropriate columns from event data table. Now both IMAGING and
;       TIMING (timing and burst [pn only]) modes are supported.
;
; Example calls:
;
; IDL> t = load_events("/data/49/sav2/xmm/ngc4051/events/rev0263_pn_src.ds", $
;                       CHAN_LIM = [300,1000], T0=t0, CHATTER=2)
; IDL> plot, histogram(t, BINSIZE=100.0), PSYM=10
;
; HISTORY:
;       23/04/10 - v1.0 - First basic version
;       09/07/10 - v1.1 - Added PATT_LIM keyword for PATTERN
;                          selection
;       01/10/10 - v1.2 - Added T_CLIP keyword
;       19/01/11 - v1.3 - Automatically find GTI extension
;                          event files (using FXPOSIT routine) which
;                          resulted in a large reduction in the amount
;                          of code needed for the MOS GTI
;                          collection. Removed redundant GT_FILE
;                          keyword. Use double precision
;                          arrays for time series output.
;       12/07/12 - v1.4 - Removed T_CLIP as time interval 'pruning'
;                          is now handled explicitly by MAKE_EPIC_LC.
;                          T0 now defined as first entry in GTI table.
;                          Changed SILENT to CHATTER keyword, with
;                          more levels of control over output.
;       16/07/12 - v1.5 - Added FRAME_TIME output
;       25/07/12 - v1.6 - Removed call to FXPOSIT (which was slow) and
;                          replaced with simpler call to FITS_INFO
;       29/06/13 - v1.7 - Added check for DATAMODE. Use correct columns
;                          for event data in TIMING mode.
;       22/01/14 - v1.8 - Used STRCMP to match names of FITS file extensions
;                          containing GTIs to the standard names. Look for
;                          GTI000* first. If not present use STDGTI??
;       29/01/14 - v1.9 - Further refinement of GTI extension selection rules.
;       30/01/14 - v2.0 - Improved GTI handling. Now combined all appropriate GTIs
;-
; ----------------------------------------------------------
    
  ; options for compilation (recommended by RSI)
    
  COMPILE_OPT idl2
  
  ; watch out for errors
  
  ON_ERROR, 0
  
  ; ----------------------------------------------------------
  ; Check the arguments
  
  ; is the file name defined?
  
  IF (n_elements(filename) eq 0) THEN BEGIN
    filename = ''
    READ,'-- Enter file name (ENTER to list current directory): ', filename
    IF (filename eq '') THEN BEGIN
      list = findfile()
      PRINT, list
      READ,'-- Enter file name: ', filename
    ENDIF
  ENDIF
  
  ; set level of output to screen
  
  IF not KEYWORD_SET(chatter) THEN chatter = 0
  IF (chatter le 1) THEN SILENT = 1
  
  ; set the CCDs to extract 04 for pn and 01 for MOS1/2 by default
  
  IF NOT KEYWORD_SET(ccdnr) THEN ccdnr = [4,1]
  
  ; ----------------------------------------------------------
  ; call READFITS to read the data from FITS file
  
  IF KEYWORD_SET(fileroot) THEN BEGIN
    file = fileroot+filename
  ENDIF ELSE BEGIN
    file = filename
  ENDELSE
  
  filedata = READFITS(file, htab, EXTEN_NO=1, SILENT=silent)
  
  ; check for errors from READFITS
  
  IF (N_ELEMENTS(filedata) le 1) THEN BEGIN
    IF (filedata eq -1) THEN BEGIN
      PRINT,'** Error calling READFITS in LOAD_EVENTS'
      PRINT,'** '+!ERROR_STATE.MSG
      RETURN, !NULL
    ENDIF
  ENDIF
  
  ; extract header information
  
  target_name = SXPAR(htab, "OBJECT")
  t_duration = SXPAR(htab, "TELAPSE")
  t_ontime = SXPAR(htab, "ONTIME")
  t_exp = SXPAR(htab, "LIVETIME")
  name_tscope = SXPAR(htab, "TELESCOP")
  name_inst = STRTRIM(SXPAR(htab, "INSTRUME"))
  n_rev = SXPAR(htab, "REVOLUT")
  datamode = STRTRIM(SXPAR(htab, "DATAMODE"), 2)
  submode = STRTRIM(SXPAR(htab, "SUBMODE"), 2)
  pi_name = STRTRIM(SXPAR(htab, "OBSERVER"), 2)
  
  ; now convert three columns -
  ;   TIME (col 1)
  ;   PI (~energy; col 9 or 5)
  ;   PATTERN (col 11 or 7)
  ; - from binary to floating point.
  ; columns are different for IMAGING and TIMING mode data
  ; See the "XMM-Newton Data Files Handbook"
  
  if (datamode eq "IMAGING") then col_list = [1, 9, 11]
  if (datamode eq "TIMING") then col_list = [1, 5, 7]
  if (datamode ne "IMAGING" and datamode ne "TIMING") then begin
    PRINT, '** DATAMODE neither TIMING not IMAGING in LOAD_EVENTS.'
    RETURN, !NULL
  endif
  
  time = TBGET(htab, filedata, col_list[0], /NOSCALE)
  pi = TBGET(htab, filedata, col_list[1])
  pattern = TBGET(htab, filedata, col_list[2])
  
  ; ----------------------------------------------------------
  ; Find a GTI table in the event file header
  
  ; get a list of the extension names from the event list file
  
  IF (chatter ge 2) THEN BEGIN
    FITS_INFO, file, EXTNAME=extname
  ENDIF ELSE BEGIN
    FITS_INFO, file, /SILENT, EXTNAME=extname
  ENDELSE
  
  ; how many GTI extensions in the file?
  
  num_gti = TOTAL( STREGEX(extname, 'GTI', /BOOLEAN) )
  IF (num_gti eq 0) THEN BEGIN
    PRINT, '** Unable to find any GTI extensions in LOAD_EVENTS:', file
    RETURN, !NULL
  ENDIF
  
  ; define GTI extension names for pn, MOS
  
  ccdnr2 = ccdnr - 1
  IF (name_inst eq "EPN") THEN BEGIN
    IF (datamode eq "TIMING") THEN ccdnr2[0] = 0
    IF (submode eq "PrimeSmallWindow") THEN ccdnr2[0] = 0
    ccd = ccdnr[0]
    ccd2 = ccdnr2[0]
  ENDIF ELSE BEGIN
    ccd = ccdnr[1]
    ccd2 = ccdnr2[1]
  ENDELSE
  
  ; build-up the likely names of the GTIs
  
  yy = 2
  gti0_root = 'GTI0'
  stdgti_root = 'STDGTI'
  IF (ccd2 lt 10) THEN gti0_root = gti0_root + '0'
  IF (ccd lt 10) THEN stdgti_root = stdgti_root + '0'
  gti0_name = gti0_root + STRTRIM(ccd2, 2)
  stdgti_name = stdgti_root + STRTRIM(ccd, 2)
  
  ; find extensions matching gti0_name or stdgti_name
  
  gti0_mask = WHERE(STRCMP(extname, gti0_name, 6), count, /NULL)
  stdgti_mask = WHERE(STRCMP(extname, stdgti_name, 8), count, /NULL)
  gti_mask = [gti0_mask, stdgti_mask]
  gti_name = extname[gti_mask]
  
  ; warning if no matching STDGTInn or GTI0nnyy extension names found
  
  N_gti = N_ELEMENTS(gti_mask)
  IF (N_gti eq 0) THEN BEGIN
    PRINT, '** Unable to find any matching GTI extensions in LOAD_EVENTS:', file
    RETURN, !NULL
  ENDIF
  
  ; --------------------
  ; loop over each matching GTI extensions
  
  FOR i = 0, N_gti-1 DO BEGIN
  
    ; call READFITS to read the data from FITS file
  
    filedata = READFITS(file, htab, EXTEN_NO=gti_mask[i], SILENT=silent)
    
    ; check for errors from READFITS
    
    IF (N_ELEMENTS(filedata) le 1) THEN BEGIN
      if (filedata eq -1) THEN BEGIN
        PRINT,'** Error calling READFITS in LOAD_EVENTS'
        PRINT,'** '+!ERROR_STATE.MSG
        RETURN, !NULL
      ENDIF
    ENDIF
    
    ; how many rows are there in the file?
    
    n_row = SXPAR(htab, "NAXIS2")
    
    ; prepare data array for output
    
    gti_i = MAKE_ARRAY(n_row, 2, /DOUBLE)
    
    ; now convert two columns - start, stop times -
    ; from binary to floating point
    
    x = TBGET(htab, filedata, 1)
    gti_i[0,0] = x[*]
    x = TBGET(htab, filedata, 2)
    gti_i[0,1] = x[*]
    x = 0
    
    ; merge with previous GTIs using GTIMERGE (and GTITRIM)
    
    IF (i eq 0) THEN BEGIN
      gti = gti_i
    ENDIF ELSE BEGIN
      gti = TRANSPOSE(GTIMERGE(TRANSPOSE(gti), TRANSPOSE(gti_i), /INTERSECT))
    ENDELSE
    
  ENDFOR
  
  ; --------------------
  
  ; define t0 as start of first listed GTI
  ;  and subtract this from the time arrays
  
  t0 = gti[0,0]
  gti = gti - t0
  time = time - t0
  
  ; extract header information
  
  frame_time = SXPAR(htab, "FRMTIME") / 1000.0D
  
  ; ----------------------------------------------------------
  ; apply the channel limits if required
  
  IF KEYWORD_SET(chan_lim) THEN BEGIN
  
    mask = WHERE(pi ge chan_lim[0] and pi le chan_lim[1], count)
    time = time[mask]
    pi = pi[mask]
    pattern = pattern[mask]
    
  ENDIF
  
  ; ----------------------------------------------------------
  ; apply the pattern limits if required
  
  IF KEYWORD_SET(patt_lim) THEN BEGIN
  
    mask = WHERE(pattern ge patt_lim[0] and pattern le patt_lim[1], count)
    time = time[mask]
    pi = pi[mask]
    
  ENDIF
  
  ; ----------------------------------------------------------
  
  ; PRINT header information
  
  IF (chatter ge 1) THEN BEGIN
  
    PRINT,"--------------------------------"
    PRINT,"-- FITS file:", filename
    PRINT,"-- Observation details:"
    PRINT,"--   Target:     ", target_name
    PRINT,"--   Mission:    ", name_tscope
    PRINT,"--   Instrument: ", name_inst
    PRINT,"--   Sub-mode:   ", submode
    PRINT,"--   Data mode:  ", datamode
    PRINT,"--   Observer:   ", pi_name
    PRINT,"--   GTI used:   ", gti_name
    PRINT, format = '("--   Duration:   ", F20.8)', t_duration
    PRINT, format = '("--   Exposure:   ", F20.8)', t_exp
    PRINT, format = '("--   Revolution: ", I20)', n_rev
    PRINT, format = '("--   Frame time: ", F20.8)', frame_time
    PRINT, format = '("--   PI range    ", I10, I10)', MIN(pi), MAX(pi)
    PRINT,"--   GTI range   ", MIN(gti), MAX(gti)
    PRINT, format = '("--   Start time: ", F20.8)', t0
    PRINT, format = '("--   First event:", F20.8)', MIN(time)
    PRINT,"--------------------------------"
    
  ENDIF
  
  ; ----------------------------------------------------------
  
  RETURN, time
  
END

