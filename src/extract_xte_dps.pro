PRO extract_xte_dps, XDF_FILE, I=i, OUTFILE=outfile, QUIET=quiet, $
                     T_SEG=t_seg, F_NYQ=f_nyq, CHAN_LIM=chan_lim


; ----------------------------------------------------------
;+
; NAME:
;       EXTRACT_XTE_DPS
;
; PURPOSE:
;       Extract and export dynamical power spectra from RXTE
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       EXTRACT_XTE_DPS, fits_files.xdf, i=0
;
; INPUTS:
;       XDF_FILE    - (string) name of ASCII file that lists the
;                     (full) paths to the individual RXTE/PCA event
;                     files 
;
; OPTIONAL INPUTS:
;       I           - (integer) which file from the list to use
;       T_SEG       - (float) length of time intervals to use (sec)
;       F_NYQ       - (float) upper freq for periodograms
;       CHAN_LIM    - (2-array) lower/upper detector channels to
;                     include
;       QUIET       - (logical) run RADPS in 'quiet' mode?
;
; OUTPUTS:
;       <file>      - GZIP'd ASCII file
;
; OPTIONAL OUTPUTS:
;       none
;
; DETAILS:
;       This is essentially a wrapper function from Craig Markwardt's
;       XTE/PCA data extract routine RADPS. 
;
;       Input is the name of a file (XDF_FILE) that lists RXTE/PCA event
;       files. From this list file I is selected. RADPS is then run to
;       extract dynamic power spectra (DPS) from the file. The DPS is
;       an array of size N_f * N_t, where N_f is the number of
;       frequencies and N_t is the number of time intervals. For each
;       time interval (t=0,1,2,...,N_t-1) it gives the periodogram as
;       a function of frequency f_j (j=0,1,2,...,N_f-1). 
;
;       Output is a single file, named using the OUTFILE keyword, or a
;       default name. The file is a GZIP'd ASCII file. It contains one
;       'header' line followed by the DPS array. The header line lists
;       the following information:
;         N_f, N_t, file_name, CHAN_LIM, T_SEG
;       for example:
;         1024 7396 /data/.../pca/FS37_41549c0-4155836.gz [19,109] dt=0.50
;
;       The data array that follows has dimensions N_f+2 * N_t. The
;       first two 'colums' (i.e. DATA[0,*] and DATA[1,*]) give the
;       time and mean count rate for each time interval. The remaining
;       entries of DATA contain the N_f*N_t values of the DPS.
;
;       Output filename: if the optional keyword OUTFILE is not
;       specified then a filename for the output will be
;       generated. This is built-up from the basic parameters:
;        '<target>_<obsid>_seg<i>_<f_Nyq>Hz_<t_seg>s_psd.dat'
;        where <target> is the target name as specified in the event
;        file, <obsid> is the ObsID from the event file, <f_Nyq> and
;        <t_seg> are the DPS parameters as described above.
;
;       The output file is initially plain ASCII but then compressed
;       using GZIP via a SPAWN command. (This relies on GZIP being
;       installed in the background environment.)
;
;       Time/frequency resolution: The time resolution and Nyquist
;       (upper) frequency of the DPS can be specified using the
;       optional keywords: T_SEG and f_Nyq. T_SEG gives the length
;       (in seconds) of segments from which to compute periodograms,
;       it therefore defines both the time resolution and frequency
;       resolution of the DPS, by Delta_f = 1/T_SEG. The number of
;       time bins withing each segment is therefore N_t =
;       T_SEG/Delta_T; the number of frequencies is N_f =
;       N_t/2. delta_t specifies the time bins size for the raw time
;       series, and is set by the Nyquist frequency f_Nyq =
;       1/(2*Delta_T). Units are seconds/Hz unless stated otherwise.
;
;       Energy channel selection: the optional keyword CHAN_LIM
;       defines the (raw) energy channels from which to extract
;       data. E.g. CHAN_LIM = [19,109] means extract from detector
;       channels 19-109 (inclusive). The range is 0-255 abd is raw
;       channels, not bins per data mode.
;
;       Normalisation: By default the individual periodograms are normalised
;       using the Leahy convention, such that Poisson noise level is
;       2.0 for counting data.
;
; EXAMPLE USAGE:
;       extract_xte_dps, "fits_files_1.xdf", I=1, T_SEG=0.5D, F_NYQ=2048D, CHAN_LIM=[18,109]
;
; PROCEDURES CALLED:
;       RADPS, WRITE_TABLE, STRNSIGNIF
;
; HISTORY:
;       13/04/2012 - v1.0 - first working version
;       18/05/2012 - v1.1 - added channel ranges to output filename
;
; NOTES:
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)
  COMPILE_OPT idl2

; watch out for errors
  on_error, 2

; -------------------------------------------
; check the arguments

  if NOT KEYWORD_SET(xdf_file) then begin
      filename = ''
      READ,'-- Enter file name (ENTER to list current directory): ',filename
      if (filename eq '') then begin
          list = FINDFILE()
          PRINT, list
          READ,'-- Enter file name: ',filename
      endif
  endif else begin
      filename = xdf_file
  endelse

; select first file, if I not set
  if (N_ELEMENTS(i) eq 0) then i=0

; set default values for time intervals, Nyquist freq and energy channels.
  if (N_ELEMENTS(t_seg) ne 1) then t_seg=1D 
  if (N_ELEMENTS(f_Nyq) ne 1) then f_Nyq = 1024D
  if (N_ELEMENTS(chan_lim) ne 2) then chan_lim=[19, 109]

; do a check that delta_t is multiple of 1/2^N
; xxx

; -------------------------------------------
; main routine

; load the list of event files
  file_list = READ_TABLE(filename, /TEXT)
  file_list = file_list + ".gz"

; display data parameters
  if NOT KEYWORD_SET(quiet) then begin
      PRINT, '-- Event file:             ', file_list[i]
      PRINT, '-- Interval length (sec):  ', t_seg
      PRINT, '-- Nyquist frequency (Hz): ', f_Nyq
      PRINT, '-- Energy channels:        ', chan_lim
   endif

; create the DPS
  RADPS, file_list[i], T_SEG, DPS, NYQUISTFREQ=f_nyq, $
         CLIMITS=chan_lim, P0=means, TIME=tt, $
         FREQAVG=f, EXPOSURE=exposure, QUIET=KEYWORD_SET(quiet)

; extract further information from events file
  filedata = READFITS(file_list[MIN(i)], htab, EXTEN_NO=1, /SILENT)
  target_name = STRTRIM(SXPAR(htab, "OBJECT"), 2)
  if NOT KEYWORD_SET(quiet) then PRINT, '-- Target:                 ', target_name
  obsid = STRTRIM(SXPAR(htab, "OBS_ID"), 2)
  if NOT KEYWORD_SET(quiet) then PRINT, '-- ObsID:                  ', obsid

; remove under-exposed data
  frac_exp = exposure / t_seg
  mask = WHERE(frac_exp gt 0.9, count)
  tt = tt[mask]
  means = means[mask]
  dps = dps[*,mask]

; get the size of the DPS array
  data_size = SIZE(dps, /DIMENSION)
  nf = data_size[0]
  nt = data_size[1]

; reset the time intervals
  t_start = tt[0]
  tt = tt - t_start

; format the time and count rate arrays
  ti = reform(tt, 1, n_elements(tt))
  mi = reform(means, 1, n_elements(means))

; combine time, count rate, dynamic PSD arrays
  dps = TEMPORARY([ti, mi, dps])

; -------------------------------------------
; construct header info for file
  chan_string = STRTRIM(chan_lim[0],2) + "," + STRTRIM(chan_lim[1],2)
  dim_string = STRTRIM(nf,2) + " " + STRTRIM(nt,2)
  header = dim_string + " " + file_list[MIN(i)] + " [" + chan_string + "] dt=" + STRTRIM(t_seg,2) 

; construct output file name
  if NOT KEYWORD_SET(outfile) then begin
    dt_string = STRNSIGNIF(t_seg, 2)
    fn = STRSPLIT(STRING(f_nyq), ".", /EXTRACT)
    fn = STRTRIM(fn[0], 2)
    chan_string = "c" + STRTRIM(chan_lim[0],2) + "-" + STRTRIM(chan_lim[1],2)
    outfile = target_name + "_" + $
              obsid + "_seg" + STRTRIM(i+1,2) + "_" + chan_string + "_" + $
              fn + "Hz_" + STRTRIM(dt_string, 2) + "s_psd.dat"

  endif

  if NOT KEYWORD_SET(quiet) then begin
    PRINT, '-- Header info: ', header
    PRINT, '-- Output file: ', outfile
  endif

; save the output to an ASCII file
  WRITE_TABLE, dps, outfile, WIDTH=100000L, HEADER=header

; use GZIP to compress output file
  SPAWN, 'ls -sh '+outfile, result
  if NOT KEYWORD_SET(quiet) then PRINT, '-- File before compression: ', result

  SPAWN, 'gzip -f '+outfile

  SPAWN, 'ls -sh '+outfile+'.gz', result
  if NOT KEYWORD_SET(quiet) then PRINT, '-- File after  compression: ', result

; -------------------------------------------
; return to user

  if NOT KEYWORD_SET(quiet) then PRINT, '-- All done.'

END
