FUNCTION TS_MULTIBISPEC, filename, n_seg=n_seg, err=err, f=f, $
                      pois=pois, img=img, Pj=Pj, text=text, cols=cols, $
                      tbin=tbin, rebin=rebin

; ----------------------------------------------------------
;+
; NAME:
;       TS_MULTIBISPEC
;
; PURPOSE:
;       Calculate the bi-coherence from list of files
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       bicoh = TS_MULTIBISPEC('filename.txt',n_seg=1024,/pois)
;
; INPUTS:
;       filename  - (string) name of file listing FITS time series
;       n_seg     - (integer) number of points in each segment
;
; OPTIONAL INPUTS:
;       dt        - (float) the sampling rate (default = 1.0)
;       n_seg     - (int) number of points in each segment
;       pois      - (logical) correct for Poisson noise level?
;       img       - (logical) if true the return full 2d bicoherence
;       text      - (logical) data in ASCII files? (default FITS)
;       cols      - (vector) two columns to load 
;       tbin      - (integer) bin in time into dt*tbin bin size
;       rebin     - (scalar) frequency rebinning factor
;
; OUTPUTS:
;       b         - array listing bicoherence b(j) or b(j,k)
;
; OPTIONAL OUTPUTS:  
;       err   - error on bicoherence b(j) or b(j,k)
;       f     - frequency (valid for f_j and f_k)
;       Pj    - power spectrum P(f)
;
; DETAILS:
;       Computes the bicoherence of a set of univariate time series.
;       The series are containted in either FITS or ASCII
;       (specified using the /text keyword).
;       The names of the files to process are listed in a single 
;       ASCII file specified by the FILENAME input parameter.
;       
;       The bicoherence is the third-order spectrum. Whereas the
;       power spectrum is a second order statistics, formed from
;       X'(f)*X(f), where X(f) is the Fourier transform of x(t), 
;       the bispectrum is a third order statistic formed from
;       X(f_j)*X(f_k)*X'(f_j+f_k). The bispectrum is therfore
;       a function of a pair of frequencies (f_j, f_k). It is also
;       a complex-valued function. The (normalised) square amplitude
;       is called the bicoherence (by analogy with the coherence
;       from the cross-spectrum).       
;
;       The bispectrum is calculated by dividing the time series
;       into M segments of length N_seg, calculating their
;       Fourier transforms and biperiodogram, then averaging
;       over the ensamble.
;
;       Although the bicoherence is a function of two frequencies
;       the default output of this function is a one dimensional
;       output, the bicoherence rebinned as a function of only
;       the sum of the two frequencies. This loses some information
;       but provides a less noisy result in many cases. If the full
;       two dimensional output is required, it can be obtained by
;       setting the IMG keyword.
;
;           k ^               The quarter of the j,k plane
;             |               contained within the limits 
;     nf/2    |   /\          k <= j and j+k <= nf is called the
;             |  /  \         "Inner Triangle" of the "principal
;             | /    \        domain" of the bifrequencies j,k.
;             |/      \
;             +-----------> j
;                     nf 
;
;       The rest of the bifrequency plane is filled with various
;       reflections or aliased versions of the data contained in the
;       Inner Triangle. Because of this redundancy, the bicoherence
;       is measured only within the IT.
;
;       The amount of rebinning is controlled with the REBIN
;       parameter. If a positive integer, we rebin with REBIN
;       points per bin, in each axis. If a negative floating point 
;       number, we rebin logarimically every factor of (1-REBIN) 
;       along each axis. For example:
;       REBIN = -0.1 will bin over ranges f -> 1.1f.
;       If REBIN=1 then no frequency binning is applied (only
;       ensamble averaging). The default value is REBIN = -0.1.
;
;       For more information on bispectra see: 
;       * Brillinger D. R., 1994, Signal Processing, 36, 239
;       * Elgar S., Guza R. T., 1988, IEEE Trans. Acoustics
;                    Speech, and Signal Processing, 36, 1667
;       * Fackrell J., 1996, PhD Thesis, Univ. Edinburgh
;       * Kim Y. C., Powers, E. J., 1979, IEEE Trans. 
;                                     Plasma Science, 7, 120
;       * Maccarone T. J., Coppi, P. S., 2002, MNRAS, 336, 817
;       * Mendel J. M., 1991, Proc. IEEE, 79, 278
;
; PROCEDURES CALLED:
;       TS_READFITS, TS_SEGMENT, LIN_REBIN, LOG_REBIN
;       READ_TABLE, STATUSLINE
;
; HISTORY:
;       18/05/2007  - v1.0 -  First working version
;       26/06/2007  - v1.1 -  added TBIN and REBIN input
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments and input

; if N_SEG not defined, set default
  if (n_elements(n_seg) eq 0) then n_seg=128

; check FILENAME has been set
  if (n_elements(filename) eq 0) then begin
      print,'** FILENAME not set in TS_MULTIBISPEC'
      return,0
  endif

; check TBIN keyword
  if (n_elements(tbin) eq 0) then tbin=1
  if ((n_seg MOD tbin) ne 0) then begin
      print,'** TBIN not factor of N_SEG in TS_MULTIBISPEC'
      return,0
  endif

; check REBIN keyword
  if (n_elements(rebin) eq 0) then rebin = -0.1

; ----------------------------------------------------------
; Load the file of files

; Load list of FITS file names from the file
  filelist = read_table(filename,/text)

; Determine how many FITS files are to be processed
  N = n_elements(filelist)
  print,'-- Files to process',N

; ----------------------------------------------------------
; Main part of procedure

; no. +ve frequencies
  nf = (n_seg/tbin)/2 

; Prepare arrays for storing components of bicoherence
; bicoh = |A|^2 / (B*C)
; where A = X_j * X_k * X'_{j+k}   - complex
;       B = |X_j * X_k|^2          - real
;       C = |X_{j+k}|^2            - real
;      Pj = |X_j|^2                - real

  A = make_array(nf+1,nf+1,/complex,value=0.0)
  B = make_array(nf+1,nf+1,/float,value=0.0)
  C = make_array(nf+1,nf+1,/float,value=0.0)
  Pj = make_array(nf+1,/float,value=0.0)

; keep track of number of segments (COUNTER) and 
; total number of periodograms calculated (NDATA)
  counter = 0
  ndata = 0
  Pn = 0.0

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
          data = read_table(filelist[0,i],/double,nmax=200000L)
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
      for seg=0,Nseg-1 do begin

; number of data points in segment j
          seglength=seglist[seg,1]-seglist[seg,0]

; if too few data for FFT (1/N > 1/n_seg) then skip
          if (seglength lt n_seg) then continue

; frequency resolution
          df = 1.0/(dt*n_seg)

; frequency array
          f = (findgen(nf+1)) * df 

; ----------------------------------------------------------
; Calculate the periodogram for each segment using FFT command

; extract only the j-th segment of data to process
          data = x[seglist[seg,0]:seglist[seg,1]]

; Break the time series into M segments of N_seg length each
          m = seglength/n_seg
          n_cut = m*n_seg
          data = reform(data[0:n_cut-1],n_seg,m,/overwrite)

; calculate mean value (DC component) for each segment
          meanx = total(data,1)/float(n_seg)

; subtract mean (DC component)
          data = temporary(data) - transpose(rebin(meanx,m,n_seg,/sample))

; bin along time dimension by factor TBIN
          if (tbin gt 1) then begin
              rebindata = rebin(data,n_seg/tbin,m)
              data = rebindata
          endif

; Calculate the Fourier transform of each segment (row)
; this gives DFT(f_j,t_i)
          data=fft(data,1,dimension=1,/overwrite)

; extract only positive frequencies f_j : j=0,1,2,...,nf
          dft = data[0:nf,*] 

; calculate the index matrix {j,k} = j+k
          indx = indgen(nf+1)
          j = rebin(indx,nf+1,nf+1,/sample)
          k = transpose(j)
          jk = j + k

; remove the triangle j+k > nf
          mask = (jk le nf)
          jk = jk * mask

; Loop over each sub-segment II = 0,1,...,M-1

          for ii=0,m-1 do begin

              string = '-- Processing file '+strtrim(strupcase(i+1),2)
              string = string+'/'+strtrim(strupcase(N),2)
              string = string+' segment '+strtrim(strupcase(seg+1),2)
              string = string+'/'+strtrim(strupcase(Nseg),2)
              string = string+' subsegment '+strtrim(strupcase(ii+1),2)
              string = string+'/'+strtrim(strupcase(m),2)+'    '
              statusline,string

; define Fourier transform of the ith segment
              Xj = dft[*,ii] 

; Normalise so that |DFT|^2 is in [rms/mean]^2 Hz^-1 units
; XXX ???    norm = sqrt(2.0 * dt / float(n_seg) / meanx[i]^2)
              norm = sqrt(2.0 * (dt*tbin) / (n_seg/tbin) / mean(x)^2 )
              Xj = temporary(Xj) * norm

; calculate the matrix X_j*X_k
              XjXk = Xj # Xj

; calculate the matrix X_{j+k}
              Xjk = Xj[jk] 
      
; calculate the triple term X_j *X_k * X_{j+k}
              XjXkXjk =  XjXk * conj(Xjk)

; sum over the ensamble (of segments)
              A = temporary(A) + XjXkXjk
              B = temporary(B) + abs(XjXk)^2
              C = temporary(C) + abs(Xjk)^2

; define Power Spectrum 
              Pj = temporary(Pj) + abs(Xj)^2

; and its Poisson noise level
              Pn = Pn + 2.0/meanx[ii]

; end loop over K sub-segments
          endfor

          counter = counter + 1
          ndata = ndata + m

; end loop over segments SEG
      endfor

; end loop over files I
  endfor

; blank line
  print

; normalise the sums to get means
  A = temporary(A)/float(ndata)
  B = temporary(B)/float(ndata)
  C = temporary(C)/float(ndata)
  Pj = temporary(Pj)/float(ndata)
  Pn = Pn/float(ndata)

; Calculate the bicoherence from |A|^2/(B*C) 
  A = abs(A)^2
  BC = B * C
 
; use mask to aviod divide-by-zero errors
  mask = where(BC ne 0.0,count)

; prepare array for 2-dimensional bicoherence
  bicoh = make_array(nf+1,nf+1,/float,value=0.0)

; normalise |bispectrum|^2 to get bicoherence
  if (count gt 0) then begin
      bicoh[mask] = A[mask] / BC[mask]
  endif

; ----------------------------------------------------------
; subtract 1/M bias (see references)
  bicoh[mask] = bicoh[mask] - 1.0/float(ndata)

; now apply de-noising normalisation
  if keyword_set(pois) then begin
      bicoh[mask] = bicoh[mask]*BC[mask]
      B = temporary(B) - (Pj[j]*Pn + Pj[k]*Pn - Pn*Pn)
      C = temporary(C) - Pn
      BC = B*C
      mask = where(BC ne 0.0,count)
      if (count gt 0) then bicoh[mask] = bicoh[mask] / BC[mask]
  endif

; inner triangle
  mask = (k le j) * (j+k le nf)

; apply masking of principal domain
  bicoh = temporary(bicoh) * mask

; ----------------------------------------------------------
; return this if user requested image (2-d bicoherence)
  if keyword_set(img) then begin

; remove zero frequencies
      f=f[1:nf]
      bicoh=bicoh[1:nf,1:nf]

; rebin the image 
      if (rebin gt 1) then binf=lin_rebin(f,binwidth=rebin*df,binmap=binmap)
      if (rebin lt 0) then binf=log_rebin(f,binfactor=1.0-rebin,binmap=binmap)
      if ((rebin gt 1) or (rebin lt 0)) then begin
          nbin = max(binmap)+1
          binb = make_array(nbin,nbin,/float,value=0.0)
          binn = make_array(nbin,nbin,/L64,value=0)
          for j=0,nf-1 do begin
              for k=0,nf-1 do begin
                  binj = binmap[j]
                  bink = binmap[k]
                  binn[binj,bink] = binn[binj,bink] + 1
                  binb[binj,bink] = binb[binj,bink] + bicoh[j,k]
              endfor
          endfor
          bicoh = binb/(binn > 1)
      endif
      f = binf
      return, bicoh
  endif

; sum over Inner Triangle into bin(j+k)
; by looping over each k

  b = fltarr(nf+1)
  for k=0,nf/2 do begin
      b[2*k:nf] = b[2*k:nf] + bicoh[k:nf-k,k]
  endfor

; normalise to get mean
  b = b / float(indgen(nf+1)/2+1)

; ----------------------------------------------------------
; Return the result to the user

  return, b
  
END
