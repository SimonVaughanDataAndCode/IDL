FUNCTION EPIC_LOAD_LC, ROOTPATH=rootpath, DT=dt, N_CHAN=N_chan, $
                       OBS_LIST=obs_list, DATA_T=data_t, ERROR=error, $
                       SEG_LIST=seg_list, T_START=t_start, $
                       DATA_B=data_b, ERR_B=err_b, $
                       CHAN_LIST=chan_list, OBS_NAME=obs_name, $
                       INST_MASK=inst_mask, PATT_LIM=patt_lim, PI_LIST=pi_list, $
                       SEED=seed, MINF=minf, FRAC=frac, CHATTER=chatter, $
                       T_CLIP_START=t_clip_start, T_CLIP_end=t_clip_end, $
                       MAX_GAP=max_gap, GAP_LIST=gap_list, NOINTERP=nointerp, $
                       NOBACK=noback, CCDNR=ccdnr
                       
; ----------------------------------------------------
;+
; NAME:
;       EPIC_LOAD_LC
;
; PURPOSE:
;       Load into an array multiple XMM/EPIC time series 
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       result = EPIC_LOAD_LC()
;
; INPUTS:
;       ROOTPATH    - (string) path to the event and index files
;
; OPTIONAL INPUTS:
;       DT           - (float) time bin size
;       N_CHAN       - (integer) Number of (log spaced) frequency bands
;       OBS_LIST     - (vector) which observations to process
;       INST_MASK    - (vector) use to ignore PN, M1, M2 cameras
;       PATT_LIM     - (array) set lower and upper PATTERN ranges
;       PI_LIST      - (array) list of PI energy intervals to use
;       MINF         - (float) minimum fraction exposure to allow
;       SEED         - (integer) optional seed for random number generation
;       CHATTER      - (integer) more or less feedback to screen?
;       T_CLIP_START - (float) length of time (sec) to 'clip' from the
;                       start of (the pn) GTI
;       T_CLIP_END   - (float) ...end of the GTI
;       MAX_GAP      - (integer/float) maximum acceptable period of bad data
;       NOINTERP     - (logical) switch off interpolation?
;       NOBACK       - (logical) switch off backound subtraction?
;       CCDNR        - (two integers) CCD numbers (pn, MOS)
;
; OUTPUTS:
;       result       - (array) combined, net EPIC count rate in
;                     [N_chan, N_time] 
;
; OPTIONAL OUTPUTS:
;       DATA_T       - (vector) start time of each output bin (sec)
;       ERROR        - (array) error on each output count rate
;       SEG_LIST     - (array) [n_seg, 2] array listing first and last
;                       element numbers of the combined time series
;       T_START      - (vector) list of the start times of each obs
;       DATA_B       - (array) background count rate in each
;                       time/energy bin
;       ERR_B        - (array) error on bkg data
;       CHAN_LIST    - (array) the list of PI channel ranges
;       OBS_NAME     - (vector) list of the observation names
;       FRAC         - (array) combined fractional exposure per bin 
;       GAP_LIST     - (array) [N_gap, 2] array listing gap indices
;
; DETAILS:
;       This is essentially a wrapper routine for MAKE_EPIC_LC, which
;       loads EPIC data from event files (and auxillary files),
;       combined data over all instruments (if specified) correcting
;       for GTI losses and background (if specified) with appropriate
;       scaling (contained in auxillary files).
;
;       ROOTPATH must specify path (directories) where the necessary 
;       files can be found: event lists and index file.
;
;       Input are the time and energy resolution for the output time
;       series. Time bins size DT is in sec. The energy resolution N_CHAN
;       is number of logarithmically space channels between PI
;       200-10000 (i.e. 0.2-10 keV). 
;       The CHAN_LIM input can be used to specify the lower and upper
;       PI channel ranges that are to be collected. Events with PI
;       valuess outside this range will be ignored. (The limit is
;       inclusive of the end points.). E.g. 
;
;         pi_list = [[200,500], [1000,3000], [5000,10000]]
;
;       Using PI_LIST takes precident of N_CHAN, in which case
;       PI_LIST has dimensions (2, N_CHAN).
;
;       Output is an array of dimensions [N_chan, N_time] where N_time
;       is the total number of time bins and contains the count rates
;       in each energy band for each time bin. These are combined over
;       all available EPIC cameras (PN, M1, M2) and background
;       subtracted and corrected for GTI losses. The actual time of each bin
;       is contained in DATA_T, and SEG_LIST records the start:stop
;       elements which define each time series. ERROR has the same
;       dimensions as RESULT and contains the uncertainties on the
;       count rates. See MAKE_EPIC_LC, LOAD_EVENTS and GTI_FIX for
;       more details of the extraction and correction procedures. The
;       output array T_START list the start times of each observation
;       (time for first bin) in absolute (spacecraft seconds) units.
;
;       If you wish to ignore a camera completely use the INST_MASK
;       vector. By default this is set to [1,0,0] meaning pn only, M1 and
;       M2 are ignored. Set a value to 0 to ignore than camera
;       completely. For example [1,1,1] will use pn+M1+M2, but
;       [0,1,1] will use only M1+M2 data.
;
;       In addition to the event files we need one ASCII (text) file 
;       (also in the directory ROOTPATH): 'scaling.txt'
;       This has one header line at the top
;       and then 7 columns listing for each of the N_OBS
;       observations an identifer name, and the source and background
;       region areas (e.g. from the BACKSCAL keywords) for each of the
;       three cameras. For example:
;
;          obs     pn_src  pn_nkg  M1_src  M1_bkg M2_src  M2_bkg
;          rev1721 1975222 7452800 1532800 832000 1542400 816000
;          rev1722 2001628 7436800 1532800 820800 1503976 832000
;              ...     ...     ...     ...    ...     ...    ...
;
;       And for each observation identifier (e.g. "rev1721") we need
;       some event files to process. These should be of the form:
;
;               <obsid>_<inst>_<src/bkg>.ds
;
;       where <obsid> is the observation identifier. <inst> is the
;       instrument ("pn", "M1" or "M2"). <src/bkg> is either "src" or
;       "bkg" to indicate source or background events. Asuming all
;       three cameras are available for "rev1721" we expect to find
;       the following files:
;
;               rev1721_M1_bkg.ds
;               rev1721_M1_src.ds
;               rev1721_M2_bkg.ds
;               rev1721_M2_src.ds
;               rev1721_pn_bkg.ds
;               rev1721_pn_src.ds
;
;       If there are no data for a given camera we set the area
;       scalings (in the 'scaling.txt' file) to zero. This
;       will then be skipped over.
;
;       For each obsid listed in the 'scaling.txt' files we
;       therefore expect 6 files (3 src events, 3 bkg events) assuming
;       all three cameras were used. 
;       
;       The vector OBS_LIST can be used to select only a few observations. 
;       For example, if 5 observations are listed in the index file
;       (<obsid>_<inst>_<src/bkg>.ds) the default value is 
;       OBS_LIST = [1,2,3,4,5] (note: start counting from 1).
;       But if you only wish to take data from observation 4 and 5 set
;       OBS_LIST = [4,5].
;       
;       The FRAC output keyword shows the fractional exposure per bin 
;       (0.0-1.0) averaged over the cameras used.
;       
;       The CHATTER keyword switches of the screen output showing the
;       data information. The settings are:
;        -1 - minimal output
;         0 - obsid level summary only (default)
;         1 - info for each event file
;         2 - info for each FITS header 
;
;       Gaps in sampling are filled by interpolation (using GTI_FIX).
;       The interpolated data is randomised using the appropriate
;       Poisson distribution if the POIS keyword it set.
;       If the NOINTERP keyoword is set the interpolated data
;       is set to NaN and treated as missing data.
;
;       NOTE: There must be event files for source and background.
;       However, in some cases (such as very bright sources) you
;       may prefer not to subtract the background data. In such
;       cases set the NOBACK keyword. This causes the area-scaled
;       background time series to be created, but not subtracted 
;       from the source time series. The background can therefore
;       still be used to check for background flares, etc.
;
;       ADDED 29/01/2014: You may specify the CCDNR for the pn and MOS
;       targets. Passed to LOAD_EVENTS. See documentation there.
; 
; EXAMPLE USAGE:
;       rootpath = '/data/49/sav2/xmm/ngc4051/events/'
;       x = EPIC_LOAD_LC(rootpath=ROOTPATH, N_chan=5, data_t=time)
;       plot, time, x[0,*], psym=10
;
; HISTORY:
;       18/06/2010 - v1.0 - first working version
;       23/06/2010 - v1.1 - added OBS_NAME output
;       24/06/2010 - v1.2 - replaced FLAG with SEG_LIST output
;       09/07/2010 - v1.3 - added INST_MASK keyword
;       16/07/2010 - V1.4 - Added PATT_LIM keyword
;       10/11/2010 - v1.5 - Added PI_LIST keyword
;       25/04/2012 - v1.6 - added SEED input/output keyword, for
;                            controlling randomisation step in
;                            GTI_FIX, useful for debugging
;                            Added MINF input keyword and FRAC 
;                            output keyword. Use double precision 
;                            arrays for time series output and 
;                            exposure calculations
;       26/04/2012 - v1.7 - Changed default values for ROOTPATH;
;                            N_OBS now defines itself based on 
;                            the length of the event index list.
;                            New OBS input can use used to select 
;                            specific observations. Added QUIET
;                            keyword
;       02/07/2012 - v1.8 - Added ERR_B output
;       15/07/2012 - v1.9 - No longer uses three scaling files, 
;                            just one: 'scaling.txt'
;                            Neater output to screen.
;                            Added T_CLIP_START/END keywords.
;                            Added CHATTER keyword (replaces QUIET).
;                            Neater format of on-screen output
;                            (for CHATTER >= 0). BACKSCAL values stored
;                            in double prec arrays throughout.
;      16/07/2012 - v1.91 - Added check for DT/FRAME_TIME incompatibility
;                            Added check for low ct/bin levels
;      17/07/2012 - v1.92 - Changed definition of N_chan to catch bug.
;      20/07/2012 - v2.0  - Set ROOTPATH default to be the current working
;                            directory. Added a better warning if no index
;                            file is found. 
;      23/07/2012 - v2.1  - Added GTI_CHECK step before finish.
;                            Upon error return !NULL
;      02/10/2012 - v2.2  - Added clause for single-obs case, when 'scaling.txt'
;                            lists only one line of observations
;      04/07/2013 - v2.3  - Added NOINTERP and NOBACK keywords
;      29/01/2014 - v2.4 - Added CCDNR input keywords (passed to
;                           LOAD_EVENTS)
;      03/02/2014 - v2.5 - Changing loading of 'scaling.txt' file to
;                           avoid problems due to trailing blank lines.
;                            
; PROCEDURES CALLED:
;          ENERGY_BANDS, READ_TABLE, GTI_CHECK
;          MAKE_EPIC_LC [LOAD_EVENTS, GTI_FIX]
;          [these depend on FITS I/O routines]
;
; NOTES:
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)
 
  COMPILE_OPT idl2

; watch out for errors

  ON_ERROR, 3

; ---------------------------------------------------------
; Check arguments

; set path to data files. 
; Defaults to the current directory. 
; NB: Use PATH_SEP to determine the correct (OS-dependent)
; file path segment separator character 

  CD, CURRENT = cwd
  IF NOT KEYWORD_SET(rootpath) THEN rootpath = cwd + PATH_SEP()

; set channel number

  IF (N_ELEMENTS(N_chan) ne 1) THEN N_chan = 1

; make list of N_CHAN log-spaced energy bands

  IF (N_ELEMENTS(pi_list) eq 0) THEN BEGIN
      pi_list = (ENERGY_BANDS(N_chan))
  ENDIF
  xpi_list = TRANSPOSE(pi_list)
  N_chan = N_ELEMENTS(pi_list)/2

; set time bin size to default 

  IF (N_ELEMENTS(dt) ne 1) THEN dt = 100.0D

; set INST_MASK to default [1,0,0], i.e. only pn

  IF NOT KEYWORD_SET(inst_mask) THEN inst_mask = [1,0,0]

; set PATT_LIM (PATTERN selection) to default (0-4) if not given on
; input

  IF (N_ELEMENTS(patt_lim) eq 0) THEN patt_lim = [0, 4]

; check the SEED keyword, if not set assign to !NULL and it will be
; automatically defined when the ranndom number generator is first
; called 

  IF NOT KEYWORD_SET(seed) THEN seed = !NULL

; have we been given minf?

  IF (N_ELEMENTS(minf) ne 1) THEN minf = 0.3D

; define the level out output to screen

  IF NOT KEYWORD_SET(chatter) THEN chatter = 0

; set the GTI 'pruning' parameters if not set already

  if NOT KEYWORD_SET(t_clip_start) then t_clip_start = 10.0D
  if NOT KEYWORD_SET(t_clip_end) then t_clip_end = 100.0D

; set the 'bad time interval' length
; (used by GTI_CHECK)

  if NOT KEYWORD_SET(max_gap) then max_gap = 100.0D
  
; set other processing options

  if NOT KEYWORD_SET(nointerp) THEN nointerp = 0
  if NOT KEYWORD_SET(noback) THEN noback = 0

; set the CCDs to extract 04 for pn and 01 for MOS1/2 by default

  IF NOT KEYWORD_SET(ccdnr) THEN ccdnr = [4,1]

; ---------------------------------------------------------
; Begin main routine

; ----------------------------------------------------
; Load the file 'scaling.txt' which contains the following
; vital information about the event files
;  1. the rootname of the dataset OBS_NAME
;  2. the src and bkg BACKSCAL (extraction region areas) for 
;      each camera (pn_src, pn_bkg, M1_src, M1_bkg, M2_src, M2_bkg)
; From these calculate the src/bkg scaling factor, SCALE. 

  scale_name = 'scaling.txt'
  scale_data = READ_TABLE(rootpath+scale_name, head=1, /TEXT)

; check the file was picked up

  IF (N_ELEMENTS(scale_data) le 1) THEN BEGIN
    IF (chatter ge 0) THEN BEGIN
      PRINT, '** Problem with index file:', rootpath + scale_name
      PRINT, '** Check ROOTPATH keyword and try again.'
    ENDIF
    RETURN, !NULL
  ENDIF

; define N_OBS, and OBS_LIST if not alredy defined

  N_obs = TOTAL(scale_data[0, *] ne "", /INTEGER)
  scale_data = scale_data[0:6, 0:(N_obs-1)]

  IF NOT KEYWORD_SET(obs_list) THEN obs_list = INDGEN(N_obs)+1

  IF (MAX(obs_list) gt N_obs) THEN BEGIN
    IF (chatter ge 0) THEN PRINT, '** OBS_LIST exceeds N_obs in EPIC_LOAD_LC.'
    RETURN, !NULL
  ENDIF
  N_obs = N_ELEMENTS(obs_list)

; extract the src and bkg areas

  src = DOUBLE(scale_data[[1,3,5], *])
  bkg = DOUBLE(scale_data[[2,4,6], *])

; remove unwanted observations

  src = src[*, obs_list-1]
  bkg = bkg[*, obs_list-1]

; compute the src/bkg area scaling factors

  mask = WHERE(bkg gt 0.0, count)
  scale = MAKE_ARRAY(3, N_obs, /DOUBLE)
  scale[mask] = src[mask] / bkg[mask]
  scale = TRANSPOSE(scale)

; Extract the list of OBS_NAMES, and trim the array to one dimension

  obs_name = scale_data[0, obs_list-1]
  obs_name = REFORM(obs_name)

; ignore some cameras according to INST_MASK

  mask = WHERE(inst_mask eq 0, count)
  IF (count gt 0) THEN scale[*, mask] = 0.0D
  
; ----------------------------------------------------
; now load the event files for each dataset
; for each dataset call MAKE_EPIC_LC to produce a merged EPIC time
; series from the relevant event files. Record the start time and
; store the count rate and error in an array. Then concatinate these
; arrays over all datasets. Note that the arrays contain
; [energy, time], i.e. they contain time series in N_chan energy bins
; and N_time time intervals in an array of dimensions [N_chan,
; N_time]. Use the FLAG array to keep track of which elements in the
; output array are from which dataset, by listing the elements of the
; output are the start:stop times of each observation

  seg_list = MAKE_ARRAY(N_obs, 2, /LONG)
  t_start = MAKE_ARRAY(N_obs, /DOUBLE)

  str = '-- EPIC_LOAD_LC processing' + STRING(N_obs) + $
        ' observations using PATTERNs:' + STRING(PATT_LIM[0]) + $
        '-' + STRING(PATT_LIM[1])
  PRINT, STRCOMPRESS(str)

  FOR j=0, N_obs-1 DO BEGIN

      scale_j = REFORM(scale[j, *])
      file_j = obs_name[j]

; extract a binned time series, summed over all active cameras

      x = MAKE_EPIC_LC(file_j, scale_j, ROOTPATH=rootpath,  T0=t0, DT=dt, /POIS, $
                       PI_LIST=xpi_list, TIME=time, ERR=err, CHATTER=chatter, BKG=bkg, $
                       PATT_LIM=PATT_LIM, FRAC_I=frac_i, SEED=seed, MINF=minf, $
                       B_ERR=b_err, T_CLIP_START=t_clip_start, T_CLIP_end=t_clip_end, $
                       NOINTERP=nointerp, NOBACK=noback)

; check for errors 

      if (N_ELEMENTS(x) eq 0) THEN RETURN, !NULL

; sum over all energy ranges, and compute mean of this total rate

      mean_x = MEAN(TOTAL(x, 1))
      mean_bkg = MEAN(TOTAL(bkg, 1))
      mask = WHERE(inst_mask eq 1)
      loss = WHERE(frac_i[mask,*] lt minf, count)
      loss = FLOAT(count)/N_ELEMENTS(frac_i[mask,*])

; average the FRAC (fractional exposure per bin) over all active cameras

      frac_i = TOTAL(frac_i[mask,*], 1) / N_ELEMENTS(mask)

; print (to screen) a summary of the data

      IF (chatter ge 0) THEN BEGIN

        IF (j eq 0) THEN BEGIN
          PRINT,'--          Name Inst.       Duration  Start time           ', $
                'Src rate   Bkg rate  Intp frac Intp no  ', $
                'Scale (pn, MOS1, MOS2)'
          ENDIF

         qinst = ''
         IF (scale_j[0] ne 0) THEN qinst = qinst + 'pn' ELSE qinst = qinst + '  '
         IF (scale_j[1] ne 0) THEN qinst = qinst + ' M1' ELSE qinst = qinst + '   '
         IF (scale_j[2] ne 0) THEN qinst = qinst + ' M2' ELSE qinst = qinst + '   '

         PRINT, format = '("-- ", A13, TR1, A8, F12.2, F20.8, 2(F11.5), F9.5, I11, TR3, 3(F10.5))', $
                 STRTRIM(file_j,2), qinst, max(time), $
                 t0, mean_x, mean_bkg, loss, count, scale_j[*]

      ENDIF

; merge the arrays over multiple observations

      t_start[j] = t0
      N = N_ELEMENTS(time)

      IF (j eq 0) THEN BEGIN

          data_x = x
          data_t = time
          error  = err
          err_b = b_err
          seg_list[0,*] = [0, N-1]
          data_b = bkg
          frac = frac_i

      ENDIF ELSE BEGIN

          t_offset = t_start[j] - t_start[0]
          seg_list[j,*] = [0, N-1] + N_ELEMENTS(data_t)
          data_x   = [[data_x], [x]]
          data_t   = [data_t, time+t_offset]
          error    = [[error], [err]]
          data_b   = [[data_b], [bkg]]
          err_b    = [[err_b], [b_err]]
          frac     = [[frac], frac_i]

      ENDELSE

  ENDFOR

  IF (chatter ge 0) THEN BEGIN
    IF (MIN(data_x + data_b)*dt le 5) THEN BEGIN
      PRINT, '-- WARNING: some data have low count/bin.'
    ENDIF
  ENDIF

; check for large "bad time intervals"

  gap_list = GTI_CHECK(frac, SEG_LIST=seg_list, MINX=MINF, MAX_GAP=max_gap, DT=dt)

  IF (chatter ge 0) THEN PRINT, '-- EPIC_LOAD_LC finished.'

; ---------------------------------------------------------
; Return to user

  chan_list = pi_list
  RETURN, data_x

END
