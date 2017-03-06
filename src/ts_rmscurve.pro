PRO TS_RMSCURVE, filename, n_seg=n_seg, freq=freq, $
                 p_n=p_n, pois=pois, deadtime=deadtime, var=var

; ----------------------------------------------------------
;+
; NAME:
;       TS_RMSCURVE
;
; PURPOSE:
;       Calculate time series for mean flux and rms
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       TS_RMSCURVE,'filename.txt'
;
; INPUTS:
;       filename  - (string) name of file listing FITS time series
;       n_seg     - (integer) number of points in each rms segment
;
; OPTIONAL INPUTS:
;       pois      - (logical) subtract Poisson noise level?
;       deadtime  - (logical) correct for detector deadtime?
;       freq      - (vector) list of two frequencies defining
;                   the range over which to calculate the rms
;       var       - (logical) output variance rather than RMS?
;
; OUTPUTS:
;       <files>   - series of ASCII files 
;       rmscurve.list - ASCII file listing the output filenames
;
; OPTIONAL OUTPUT:
;       p_n       - noise level of the variance periodogram
;
; DETAILS:
;       Calculate mean flux and rms in segments of length n_seg
;       for several time series in FITS files. The names of the
;       FITS files to process are listed in a single ASCII file
;       specified by the filename input parameter.
;
;       Each series is broken into segments that are contiguous
;       and from each is computed a time series of the mean
;       flux and the rms. The mean is calculated in the
;       standard fashion, and the rms is calculated by integrating
;       the periodogram over the frequency range freq=[f1,f2].
;       If the time (dt) is in units of sec, then f1,f2 should
;       be in units of sec^-1 = Hz.
;
;       The periodograms may be optionally corrected for 
;       Poisson noise, in which case it is assumed the time series are 
;       in units of ct/s and dt is in units of s.
;
;       It is assumed that the input values are in time order
;       in intervals of dt, with no gaps.
;
; Example call:
; IDL> TS_RMSCURVE,'output.list',n_seg=65536,freq=[2.0,20.0],/pois
;
; PROCEDURES CALLED:
;       TS_READFITS, TS_SEGMENT, TS_RMSFLUX, WRITE_TABLE
;
; HISTORY:
;       27/04/07  - v1.0 - first working version
;       01/05/07  - v1.1 - added check for seglength > n_seg
;                          added keyword_set(pois) check in call to TS_RMSFLUX
;                          added output of filename list (rmscurve.list)
;       03/05/07  - v1.2 - added calculation of variance-of-variance
;                          term (P_N)
;
; NOTES:
;       + Add deadtime correction (to TS_RMSFLUX)
;       + Strictly, the average periodogram computed as P is the
;         unweighted mean of the time-averaged periodograms from
;         each segment. These should be weighted by the amount of
;         time (data) in each segment
;
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments

; if n_seg not defined, set default
  if (n_elements(n_seg) eq 0) then n_seg=128

; check filename has been set
  if (n_elements(filename) eq 0) then begin
      print,'** FILENAME not set in TS_RMSCURVE'
      return
  endif

; if freq=[f1,f2] not defined then set to zero
; (signals to subsequent routines that it is not set)
  if (n_elements(freq) ne 2) then begin
      freq=[0,0]
  endif 

; ----------------------------------------------------------
; Main part of procedure

; Load list of FITS file names from the file
  filelist = read_table(filename,/text)

; Determine how many FITS files are to be processed
  N = n_elements(filelist)
  print,'-- Files to process',N

; Make array to store output filenames
  outlist = strarr(1000)

; initialise noise level
  p_n = 0.0

; Loop over each FITS file to process
  counter = 0
  for i=0,N-1 do begin

; Read the FITS file
      x = ts_readfits(filelist[0,i],t=time,dt=dt)

; Break time series into contiguous segments 
      seglist = ts_segment(time,Nseg=Nseg,minseg=n_seg)

      if (Nseg eq 0) then continue

; Loop over each segment
      for j=0,Nseg-1 do begin

          print,'-- Processing file/segment',i+1,j+1

; number of data points in segment j
          seglength=seglist[j,1]-seglist[j,0]

; if too few data for variance calculation (1/N > 1/n_seg) then skip
          if (seglength lt n_seg) then continue

; extract only the j-th segment of data to process
          data = x[seglist[j,0]:seglist[j,1]]

; note the start time for segment j of file i
          t_start = time[seglist[j,0]]

; Compute time series for mean and rms
          rms = ts_rmsflux(data,dt=dt,n_seg=n_seg,freq=freq,mean=mean,t=t, $
                           pois=keyword_set(pois), var=keyword_set(var),pow=pow,df=df)

          if (counter eq 0) then begin
              p = pow
          endif else begin
              p = temporary(p) + pow
          endelse

; Generate a filename for the output
          file_string = strsplit(filelist[0,i],'/.',count=nstring,/extract)
          obsid = file_string[nstring-3]
          flag_file = file_string[nstring-2]
          flag_segm = strtrim(strupcase(j+1),2)
          filename = obsid+'/'+flag_file+'-'+flag_segm+'.rms'

; add the start time back onto the bin times
          t=double(t)+t_start

; put time (t), mean flux (mean) and rms (rms) in data array
          data = transpose([[t],[mean],[rms]])
          if (n_elements(t) eq 1) then begin
              data = [t[0],mean[0],rms[0]]
              data = reform(data,3,1)
          endif

; Save the output as an ASCII file
          write_table,data,filename,/double
          outlist[counter] = filename
          counter = counter + 1
                  
; End loop over segments
      endfor

; End loop over FITS files
  endfor

; normalise to get average power spectrum
  p = p / float(counter)

; variance-of-variance term : V[sigma] = \sum p(f)^2 df^2
  f = (indgen(n_elements(p))+1)*df
  mask = where((f ge freq[0]) and (f le freq[1]))
  E_sigma = total(p[mask]) * df
  V_sigma = total(p[mask]^2) * df^2 

; convert from variance to expected power density
; (formula is different depending on whether VAR or RMS is measured)
  if keyword_set(var) then begin
      p_n = V_sigma * (2.0*n_seg*dt)
  endif else begin
      p_n = 0.25 * V_sigma / E_sigma * (2.0*n_seg*dt)
  endelse


; output the list of filenames 
  write_table,outlist[0:counter-1],'rmscurve.list'

; ----------------------------------------------------------
; Return the data array to the user

  return

END
