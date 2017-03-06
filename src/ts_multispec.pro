FUNCTION TS_MULTISPEC, filename, n_seg=n_seg, rebin=rebin, $
          dp=dp, freq=freq, f0=f0, f1=f1, pois=pois, $ 
          deadtime=deadtime, rms=rms, text=text, cols=cols, $
          meandata=meandata, head=head, bin_n=bin_n, binmap=binmap, $
          minbin=minbin

; ----------------------------------------------------------
;+
; NAME:
;       TS_MULTISPEC
;
; PURPOSE:
;       Calculate average periodogram from list of files
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       pow = TS_MULTISPEC('filename.txt',n_seg=1024,/pois)
;
; INPUTS:
;       filename  - (string) name of file listing FITS time series
;       n_seg     - (integer) number of points in each segment
;
; OPTIONAL INPUTS:
;       rebin     - (scalar) frequency rebinning factor
;       pois      - (logical) subtract Poisson noise level?
;       deadtime  - (logical) correct for detector deadtime?
;       rms       - (logical) use [rms/mean]^2 norm for spectrum?
;       text      - (logical) data in ASCII files? (default FITS)
;       cols      - (vector) two columns to load 
;       head      - (integer) number of lines to skip in header
;       minbin    - (integer) minimum no. data points per bin
;
; OUTPUTS:
;       pow       - averaged power density spectrum
;
; OPTIONAL OUTPUTS:
;       dp        - (vector) error bar on averaged spectrum
;       freq      - (vector) frequencies (centre of bin)
;       f0        - (vector) lower bound of frequency bin
;       f1        - (vector) upper bound of frequency bin
;       meandata  - (scalar) mean level of input data 
;       bin_n     - (vector) number of frequencies/bin
;       binmap    - (vector) bin number for each frequency
;
; DETAILS:
;       Estimate the power spectrum of a time series contained
;       in several files. The files can be either FITS or
;       ASCII (specified using the /text keyword).
;       The names of the files to process are listed in a single 
;       ASCII file specified by the FILENAME input parameter.
;
;       If the time series files as ASCII format, the COLS
;       parameter specifies which two columns to use for time
;       and value. Defaul is COLS=[0,1] meaning column 0=time
;       and column 1=value.
;
;       Each series is broken into segments that are contiguous
;       and from each is computed a periodogram.
;       These are ensamble averaged, and then rebinned 
;       in frequency to give a consistent estimate
;       of the power spectrum (with error bars).
; 
;       The amount of rebinning is controlled with the REBIN
;       parameter. If a positive integer, we rebin with REBIN
;       points per bin. If a negative floating point number, we
;       rebin logarimically every factor of (1-REBIN). For example
;       REBIN = -0.1 will bin over ranges f -> 1.1f.
;       If REBIN=1 then no frequency binning is applied (only 
;       ensamble averaging). The default value is REBIN = 50.
;
;       The error bars on the power densities are computed from
;       the chi-square distribution of powers, i.e. error is
;       1/sqrt{M} where M is the number of periodogram ordinates
;       averaged in a given frequency bin (either ensamble or
;       frequency averaging). In the limit of large M (>50)
;       the averaged power densities are normally distributed.
;       
;       The periodograms may be optionally corrected for 
;       Poisson noise, in which case it is assumed the time series are 
;       in units of ct/s and dt is in units of s.
;
; PROCEDURES CALLED:
;       TS_READFITS, TS_SEGMENT, DYNAMIC_PSD, LIN_REBIN, LOG_REBIN
;       READ_TABLE, STATUSLINE
;
; HISTORY:
;       02/05/07  - v1.0 - first working version
;       03/05/07  - v1.1 - added option for no-rebinning
;                          bug-fix to pow normalisation
;       04/05/07  - v1.2 - added MEANDATA output
;       17/05/07  - v1.3 - used STATUSLINE for progress report
;       05/06/09  - v1.4 - added HEAD keyword
;       08/06/09  - v1.5 - bug fix. REBIN keyword was being changed 
;                          if LIN_REBIN was called. 
;                          Added MINBIN keyword
;
; NOTES:
;       + Add deadtime correction (to DYNAMIC_PSD)
;       + WARNING: Errors are wrong when POIS keyword used!
;         (spotted 05/06/09)
;
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments

; if N_SEG not defined, set default
  if (n_elements(n_seg) eq 0) then n_seg=128

; if REBIN not defined, set default
  if (n_elements(rebin) eq 0) then rebin=50

; check FILENAME has been set
  if (n_elements(filename) eq 0) then begin
      print,'** FILENAME not set in TS_MULTISPEC'
      return,0
  endif

; is the minimum number of points per bin supplied?
  if (n_elements(minbin) eq 0) then minbin=1

; make sure it's an integer
  minbin=round(minbin)

; ----------------------------------------------------------
; Main part of procedure

; Load list of FITS file names from the file
  filelist = read_table(filename,/text)

; Determine how many FITS files are to be processed
  N = n_elements(filelist)
  print,'-- Files to process',N

; keep track of number of segments (COUNTER) and 
; total number of periodograms calculated (NDATA)
; and mean level of data array (MEANDATA)
  counter = 0
  ndata = 0
  meandata = 0.0

; Loop over each FITS file to process
  for i=0,N-1 do begin

; Read the file
      if keyword_set(text) then begin
          if (n_elements(cols) ne 2) then begin
              tcol = 0
              xcol = 1
          endif else begin
              tcol = cols[0]
              xcol = cols[1]
          endelse
          if keyword_set(head) then begin
              data = read_table(filelist[0,i], /double, nmax=200000L, head=head)
          endif else begin
              data = read_table(filelist[0,i], /double, nmax=200000L)
          endelse
          if (n_elements(data) le 1) then continue
          x = reform(data[xcol,*])
          time = reform(data[tcol,*])
          dt = time[1]-time[0]
      endif else begin
          x = ts_readfits(filelist[0,i],t=time,dt=dt)
      endelse

; Break time series into contiguous segments 
      seglist = ts_segment(time,dx=dt,Nseg=Nseg,minseg=n_seg)

      if (Nseg eq 0) then continue

; Loop over each segment
      for j=0,Nseg-1 do begin

          string = '-- Processing file '+strtrim(strupcase(i+1),2)
          string = string+'/'+strtrim(strupcase(N),2)
          string = string+' segment '+strtrim(strupcase(j+1),2)
          string = string+'/'+strtrim(strupcase(Nseg),2)
          statusline,string

; number of data points in segment j

          seglength=seglist[j,1]-seglist[j,0]

; if too few data for variance calculation (1/N > 1/n_seg) then skip

          if (seglength lt n_seg) then continue

; extract only the j-th segment of data to process
          data = x[seglist[j,0]:seglist[j,1]]

; note the start time for segment j of file i
          t_start = time[seglist[j,0]]

; Compute periodograms
          dpsd = DYNAMIC_PSD(data,n_seg,dt=dt,t=t,f=f,df=df,mean=ave, $
               pois=keyword_set(pois),rms=keyword_set(rms))

; Add all periodograms together
          if (counter eq 0) then begin
              pow = total(dpsd,1)
          endif else begin
              pow = temporary(pow) + total(dpsd,1)
          endelse

; calculate number of files (COUNTER) and number of periodograms (NDATA)
; and mean data level
          counter = counter + 1
          ndata = ndata + (size(dpsd))[1]
          meandata = meandata + total(ave)

; End loop over segments
      endfor

; End loop over FITS files
  endfor

 ; black line
  print

; normalise by number of data segments (to get ensamble average power)
  pow = temporary(pow) / float(ndata)
  meandata = meandata / float(ndata)

; ----------------------------------------------------------
; if REBIN is positive then use linear frequency bins
  if (rebin gt 1) then begin
      binwidth = floor(rebin)*df
      bin_p = lin_rebin(f, pow, binwidth=binwidth, bin_x=freq, $ 
                        bin_l=f0, bin_u=f1, bin_n=bin_n, binmap=binmap, $
                        minbin=minbin)
  endif

; if REBIN is negative use logarithmic binning
  if (rebin lt 0) then begin
      bin_p = log_rebin(f, pow, binfactor=1.0-rebin, bin_x=freq, $ 
                        bin_l=f0, bin_u=f1, bin_n=bin_n, binmap=binmap, $
                        minbin=minbin)
  endif

; otherwise do not rebin
  if ((rebin ge 0) and (rebin le 1)) then begin
      bin_n = intarr(n_elements(f)) + 1
      bin_p = pow
      freq = f
      f0 = f - 0.5*df
      f1 = f + 0.5*df
  endif

; calculate total number of periodogram ordinates in each frequency
; bin
  bin_n = bin_n * ndata

; calculate error on average power in each frequency bin
  dp = bin_p / sqrt(bin_n)

; ----------------------------------------------------------
; Return the data array to the user

  return,bin_p

END
