FUNCTION epic_segment, X, SEG_LIST=seg_list, N_SEG=n_seg, $
                       MEANS=means, CHATTER=chatter, $
                       INDX=indx, GAP_LIST=gap_list

; ----------------------------------------------------------
;+
; NAME:
;       EPIC_SEGMENT
;
; PURPOSE:
;       Take output from EPIC_LOAD_LC and divide into equal length
;       segments 
;
; AUTHOR:
;       Simon Vaughan (U.Leicester) 
;
; CALLING SEQUENCE:
;        X_SEG = EPIC_SEGMENT(x, error, seg_list, N_SEG=100, DY_OUT=err_seg)
;
; INPUTS:
;        X        - (array) merged time series with
;                     dimensions [N_chan, N_time] 
;
; OPTIONAL INPUTS:  
;        SEG_LIST - (array) [N_obs, 2] array listing first and last
;                      indices of the combined time series
;        N_SEG    - (integer) segment size to use for splitting
;        CHATTER  - (integer) more or less feedback to screen?
;        GAP_LIST - (array) [N_gap, 2] array listing first and last
;                      indices of any large patches of "bad" data
;
; OUTPUTS:
;        X_OUT    - (array) trimmed time series array with each
;                   segment length seg_list[i,0]:seg_list[i,1] cut
;                   down to a multiple of N_SEG
;
; OPTIONAL OUTPUTS:  
;        MEANS    - (array) mean values within each segment
;        INDX     - (vector) indices of input array X that are
;                             used to form the output array
;
; DETAILS:
;        Used for 'trimming' a (multiple) time series array X into 
;        blocks of length N_SEG.
;        
;        The input array comprises [N_chan, N_time] values, with N_time
;        observations in N_chan channels, as produced by EPIC_LOAD_LC.
;        This may contain multiple observations (N_obs > 1) concatinated
;        into a single array. An accompanying array SEG_LIST [N_obs, 2]
;        lists the rows of X corresponding to the first and last elements
;        of the each observation.    
;
;        Each observation can have arbitrary length. But for analyses
;        such as Fourier analysis, it is often useful to treat time series
;        in equal length segments (non-overlapping). To do this we
;        extract a subset of consecutive data points from each
;        observation, of length M*N_seg where M is an integer.
;        The output array X_OUT is therefore a subset of X of length
;        K*N_seg, with each block of N_seg consecutive points taken
;        from within (not crossing) each observations. I.e. we respect
;        the observation boundaries given by SEG_LIST.
;        For example, SEG_LIST may look like:
;        
;                  0        9248       18313       24700 
;               9247       18312       24699       30571    
;        
;        The entire time series is 30572 points in length, but
;        this is the merger for four complete observations.
;        Observation 1 is found in X[*, 0:9247], and so on 
;        for other observations. 
;        
;        If SEG_LIST is not supplied we assume the input array X is 
;        entirely from one observation. 
;        
;        GAP_LIST provides additional restristions on how the data
;        should be segmented. GAP_LIST lists the indices corresponding
;        to the start and stop times of "bad" data. This is an [N_gap, 2]
;        array similar in format to SEG_LIST. If a gap falls with a
;        block of length N_seg then the block starting point is shifted
;        forwards to avoid the "bad" interval.
;
; EXAMPLE USAGE:
;        x_seg = EPIC_SEGMENT(x, SEG_LIST=seg_list, N_SEG=100)
;
; HISTORY:
;        15/07/2010 - v1.0 - first working version
;        11/08/2010 - v1.1 - added IF...CONTINUE get-out clause for
;                             when there are not enough data for one
;                             segment. Also added MEANS output
;        16/02/2012 - v1.2 - use LONG type for SEG_LIST array
;        17/07/2012 - v1.3 - added CHATTER keyword.
;                             added capability for 1-dimensional input
;        01/08/2012 - v2.0 - Completely re-structured the main algorithm
;                             added INDX output, removed ERROR input.
;                             Added GAP_LIST input.
;        06/05/2013 - v2.1 - Fixed error in MEANS output. Thanks to WNA.
;        03/07/2013 - v2.2 - Added /LONG option to all INDGEN calls. This
;                             fixes a bug that occurred with long time 
;                             series. Fixed a typo. Added warning when there
;                             are no good segments.
;
; PROCEDURES CALLED:
;        INDEX_CLUMPS
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

; define the level out output to screen

  IF NOT KEYWORD_SET(chatter) THEN chatter = 0

; check the shape of the input array
; if it is 1-dimensional re-shape to 2-dimensional

  s = SIZE(x)
  shape = s[0]
  IF (shape eq 0) THEN BEGIN
    PRINT, '** Data array X not defined in EPIC_SEGMENT.'
    RETURN, !NULL
  ENDIF
  IF (shape eq 1) THEN BEGIN
    N = s[1]
    N_chan = 1
    z = x
    x = MAKE_ARRAY(1, N, TYPE=s[2])
    x[0,*] = z
  ENDIF
  IF (shape eq 2) THEN BEGIN
    N = s[2]
    N_chan = s[1]
  ENDIF

; if N_SEG not defined, set default

  IF (N_ELEMENTS(n_seg) eq 0) then n_seg=256

; check the size of the input array

  IF (N le n_seg) THEN BEGIN
      PRINT, '** Array X is too small in EPIC_SEGMENT'
      RETURN, !NULL
  ENDIF

; if SEG_LIST not defined, set default

  IF (N_ELEMENTS(seg_list) lt 2) THEN BEGIN
      seg_list = MAKE_ARRAY(1, 2, /LONG)
      seg_list[0, 0] = 0
      seg_list[0, 1] = N-1
  ENDIF
  size_seg = SIZE(seg_list)
  IF (size_seg[0] ne 2) THEN BEGIN
    PRINT, '** Array SEG_LIST is wrong shape in EPIC_SEGMENT'
    RETURN, !NULL
  ENDIF

; if GAP_LIST not defined, set default

  IF (N_ELEMENTS(gap_list) lt 2) THEN gap_list = !NULL
  size_gap = SIZE(gap_list)
  IF (size_gap[0] ne 2 and size_gap[0] ne 0) THEN BEGIN
    PRINT, '** Array GAP_LIST is wrong shape in EPIC_SEGMENT'
    RETURN, !NULL
  ENDIF  
  N_gap = size_gap[1]

; ----------------------------------------------------------
; Main routine

; Determine how many time series are to be processed

  N_obs = size_seg[1]

; expand GAP_LIST into a vector GAP_INDX listing all 
; bad data indices

  IF (N_gap eq 0) THEN BEGIN

    gap_indx = -1

  ENDIF ELSE BEGIN

    gap_indx = []
    gap_len = gap_list[*, 1] - gap_list[*, 0]

    FOR i = 0, N_gap-1 DO BEGIN
      gap_i = INDGEN(gap_len[i]+1, /LONG) + gap_list[i, 0]
      gap_indx = [gap_indx, gap_i]
    ENDFOR

    ; this is a clever trick to simultaneously sort GAP_INDX
    ; into accending order and remove duplicated entries

    hist = HISTOGRAM(gap_indx, OMIN=mindx)
    gap_indx = WHERE(hist gt 0) + mindx

  ENDELSE

; prepare for main loop

  indx = []
  
  FOR i = 0, N_obs-1 DO BEGIN
  
      seg0 = seg_list[i, 0]
      seg1 = seg_list[i, 1]
      obs_len = seg1 - seg0 + 1

      ; skip this observation if it is too short (obs_len < n_seg)

      IF (obs_len lt n_seg) THEN CONTINUE

      ; look for any bad data points within this observation

      mask = WHERE(gap_indx ge seg0 and gap_indx le seg1, count)

      ; if no bad data here, add to the list all data points up to M*N_seg

      IF (count eq 0) THEN BEGIN      

        num = obs_len / n_seg
        indx_i = INDGEN(num * n_seg, /LONG) + seg0
        indx = [indx, indx_i]
        CONTINUE

      ENDIF ELSE BEGIN

        ;  which are the bad points?
        
        gap_mask = gap_indx[mask]

        ; use HISTOGRAM to find the good points of this observation
        
        hist = HISTOGRAM( gap_mask, MIN=seg0, MAX=seg1 )
        keep = WHERE( hist eq 0, count) + seg0
        IF (count lt N_seg) THEN CONTINUE

        ; find 'clumps' of good data = GTIs
        
        INDEX_CLUMPS, keep, lo, hi, /NOSINGLE, COUNT=count
        IF (count eq 0) THEN CONTINUE
        lo = keep[lo]
        hi = keep[hi]

        ; treat each GTI 
        
        len_gti = hi - lo + 1
        good_blocks = WHERE(len_gti ge N_seg, N_blocks)
        IF (N_blocks eq 0) THEN CONTINUE
        FOR j = 0, N_blocks-1 DO BEGIN

          num = len_gti[good_blocks[j]] / N_seg
          indx_i = INDGEN(num * N_seg, /LONG) + lo[good_blocks[j]]
          indx = [indx, indx_i] 
        
        ENDFOR

      ENDELSE
      
  ENDFOR

; display progress

  N_blocks = N_ELEMENTS(indx) / n_seg
  IF (chatter ge 0) THEN BEGIN
    PRINT, '-- EPIC_SEGMENT'
    PRINT, '-- Observations to process:  ', N_obs
    PRINT, '-- Channels to process:      ', N_chan
    PRINT, '-- Gaps to deal with:        ', N_gap
    PRINT, '-- No. output blocks:        ', N_blocks
  ENDIF

  IF (N_blocks lt 1) THEN BEGIN
    PRINT, '** No good segments found in EPIC_SEGMENT
    RETURN, !NULL
  ENDIF

; find the mean of each block

  IF ARG_PRESENT(means) THEN BEGIN
  
    means = MAKE_ARRAY(N_chan, N_blocks, VALUE=0.D)
    block = INDGEN(N_seg, /LONG)
    FOR i = 0, N_blocks-1 DO BEGIN
      mask = block + i*N_seg
      mask = indx[mask]
      means[*, i] = MEAN( x[*, mask], DIMENSION=2)
    ENDFOR
  
  ENDIF

; ----------------------------------------------------------
; Return data to user

  IF (shape eq 1) THEN RETURN, z[indx]

  RETURN, x[*, indx]

END
