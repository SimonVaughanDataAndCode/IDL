FUNCTION gti_check, X, SEG_LIST=seg_list, MAX_GAP=max_gap, $
                   MINX=minx, CHATTER=chatter, DT=dt

; ----------------------------------------------------------
;+
; NAME:
;       GTI_CHECK
;
; PURPOSE:
;       Check for patches of 'bad' data in multi-observation time series 
;       array
;
; AUTHOR:
;       Simon Vaughan (U.Leicester) 
;
; CALLING SEQUENCE:
;       gap_list = GTI_CHECK(x, SEG_LIST=seg_list)
;
; INPUTS:
;       X           - (vector) the time series to check
;           
; OPTIONAL INPUTS:  
;       SEG_LIST    - (array) [n_seg, 2] array listing first and last
;                      element numbers of the combined time series
;       MINX        - (float) minimum value of X to allow
;       MAX_GAP     - (integer/float) maximum acceptable period of bad data
;       CHATTER     - (integer) more or less feedback to screen?    
;       DT          - (float) time bin size (if MAX_GAP given in time units)
;
; OUTPUTS:
;       GAP_LIST    - (array) [N_gap, 2] array listing gap indices
;
; OPTIONAL OUTPUTS:  
;       none
;
; DETAILS:
;       Takes a time series and reports continuous patches where
;       X <= MINX. Patches are only considered if they are longer 
;       than MAX_GAP points in length. 
;       
;       MAX_GAP is assumed to be in units of the step size, i.e. 
;       number of consecutive data points. If the DT keyword is set
;       the MAX_GAP is assumed to be in units of sec and the
;       segments are flagged if longer than ROUND(MAX_SEG/DT) point.
;
;       The SEG_LIST array is an [N_seg, 2] array given as input 
;       used when X is built from N_seg distinct time 
;       series merged together, e.g. output from EPIC_LOAD_LC.PRO. 
;       In this case the process of searching for periods of 'bad' 
;       data is carried out separately within each segment. 
;
;       The CHATTER keyword switches of the screen output showing the
;       data information. The settings are:
;        -1 - minimal output
;         0 - standard output (default)
;
; EXAMPLE USAGE:
;       gap_list = GTI_CHECK(x, SEG_LIST=seg_list, MAX_GAP=100)
;
; USES:
;       INDEX_CLUMPS
;
; HISTORY:
;       23/07/2012 - v1.0 - first working version
;       04/07/2013 - v1.1 - bug fix: use /LONG in INDGEN()
;       19/07/2013 - v1.2 - bug fix: use 'gt' not 'ge' in WHERE
;                            search for long gaps. Before this it
;                            would crash if gap_list was empty.
;       11/07/2014 - v1.3 - removed /NOSINGLE keyword from call to 
;                            INDEX_CLUMPS; changed defaul MAX_GAP
;                            to -1 (no filtering on length)
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

; check size of the x array

  s = SIZE(x)
  IF (s[0] ne 1) THEN BEGIN
    PRINT, '** X is not 1-dimensional array in GTI_CHECK.'
    RETURN, -1
  ENDIF
  N = s[1]
  
; check SEG_LIST or set default value

  IF NOT KEYWORD_SET(seg_list) THEN seg_list = [[0], [N-1]]
  IF (MAX(seg_list) gt N) THEN BEGIN
    PRINT, '** SEG_LIST extends beyond X array in GTI_CHECK.'
    RETURN, -1
  ENDIF

  s = SIZE(seg_list)
  IF (s[0] ne 2) THEN BEGIN
    PRINT, '** SEG_LIST is not [N_seg, 2] array in GTI_CHECK.'
    RETURN, -1
  ENDIF
  N_seg = s[1]
  seg_len = seg_list[*, 1] - seg_list[*, 0]
  
; set MINX to defaul if not set already

  if NOT KEYWORD_SET(minx) THEN MINX = 1.0D

; define the level out output to screen

  IF NOT KEYWORD_SET(chatter) THEN chatter = 0

; define the largest acceptable gap size
; (i.e. larger gaps are the ones flagged as 'bad')

  IF NOT KEYWORD_SET(max_gap) THEN max_gap = -1

; if the time bin size is set, then reset MAX_GAP
; into dimensionless units (step sizes)

  IF NOT KEYWORD_SET(dt) THEN dt = 1
  test_gap = ROUND(max_gap / dt)

; ---------------------------------------------------------
; Main routine

; define output array, empty to begin with.

  gap_list = []
  bad = 0

  FOR i = 0, N_seg-1 DO BEGIN

; consider each 'segment' of the x time series

    indx = INDGEN(seg_len[i]+1, /LONG) + seg_list[i, 0]  
    mask = WHERE(x[indx] lt MINX, count)

; escape quick if no further work needed

    IF (count eq 0) THEN CONTINUE

    bad = bad + count 
    
    INDEX_CLUMPS, mask, lo, hi, COUNT=N_gap
    IF (N_gap eq 0) THEN CONTINUE 
    lo = mask[lo]
    hi = mask[hi]
    gap_list_i = [[lo],[hi]] + seg_list[i, 0]
    gap_list = [gap_list, gap_list_i]
    
  ENDFOR

  gap_len = 0
  IF (N_ELEMENTS(gap_list) gt 0) THEN BEGIN
    gap_len = gap_list[*, 1] - gap_list[*, 0]
  ENDIF

  ; remove the gaps shorter than TEST_GAP 

  N_gap = 0
  mask = WHERE(gap_len gt test_gap, count, /NULL)   ;  changed 'ge' to 'gt' [19/07/2013]
  IF (count gt 0) THEN BEGIN
    gap_list = gap_list[mask, *]
    gap_len = gap_len[mask]
    N_gap = N_ELEMENTS(gap_len)
  ENDIF
 
  ; display results (depending on CHATTER setting)
  
  IF (chatter ge 0) THEN BEGIN
    PRINT, FORMAT='("-- fraction of bad bins (%): ",F7.3)', FLOAT(bad)/N*100
    IF (N_gap gt 0) THEN $
    PRINT, FORMAT='("-- Largest gap:              ", I7)', MAX(gap_len)
    PRINT, FORMAT='("-- Number of gaps:           ", I7)', N_gap
  ENDIF
  
; ---------------------------------------------------------
; Return to user

  RETURN, gap_list

END
