FUNCTION BI_SPECTRUM, x, n_seg=n_seg, dt=dt, err=err, f=f, $
                      pois=pois, img=img, Pj=Pj

; ----------------------------------------------------------
;+
; NAME:
;       BI_SPECTRUM
;
; PURPOSE:
;       Calculate the bi-coherence of a univariate time series.
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       bicoh = BI_SPECTRUM(x)
;
; INPUTS:
;       x - (float vector) univariate time series
;
; OPTIONAL INPUTS:
;       dt    - (float) the sampling rate (default = 1.0)
;       n_seg - (int) number of points in each segment
;       pois  - (logical) correct for Poisson noise level?
;       img   - (logical) if true the return full 2d bicoherence
;
; OUTPUTS:
;       b     - array listing bicoherence b(j) or b(j,k)
;
; OPTIONAL OUTPUTS:  
;       err   - error on bicoherence b(j) or b(j,k)
;       f     - frequency (valid for f_j and f_k)
;       Pj    - power spectrum P(f)
;
; DETAILS:
;       Computes the bicoherence of a univariate time series x.
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
;
; HISTORY:
;       09/02/2007  - v1.0 -  First working version
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments and input

; is the time series data array well-defined?
  n = n_elements(x)
  if (n lt 100) then begin
      print,'** Not enough data in BI_SPECTRUM.'
      return,0
  endif

; is the segment length supplied?
  if (n_elements(n_seg) eq 0) then n_seg=n/20

; is the sampling period dt supplied?
  if (n_elements(dt) eq 0) then dt=1.0

; ----------------------------------------------------------
; Break the time series into M segments of N_seg length each

  m = n/n_seg
  if (m lt 1) then begin
      print,'** SIZE(DATA) < N_SEG in BI_SPECTRUM'
      return,0
  endif
  n_cut = m*n_seg
  data = reform(x[0:n_cut-1],n_seg,m)

; no. +ve frequencies
  nf = n_seg/2                      

; frequency resolution
  df = 1.0/(dt*n_seg)

; frequency array
  f = (findgen(nf+1)) * df 

; ----------------------------------------------------------
; Calculate the periodogram for each segment using FFT command

; calculate mean value (DC component) for each segment
  meanx = total(data,1)/float(n_seg)

; subtract mean (DC component)
  data = temporary(data) - transpose(rebin(meanx,m,n_seg,/sample))

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
  Pn = 0.0

; Loop over each segment 0,1,...,M-1

  for i=0,m-1 do begin

; define Fourier transform of the ith segment
      Xj = dft[*,i] 

; Normalise so that |DFT|^2 is in [rms/mean]^2 Hz^-1 units
;      norm = sqrt(2.0 * dt / float(n_seg) / meanx[i]^2)
      norm = sqrt(2.0 * dt / float(n_seg) / mean(x)^2)
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
      Pn = Pn + 2.0/meanx[i]

      print,'--',i+1,'/',m

  endfor

; normalise the sums to get means
  A = temporary(A)/float(m)
  B = temporary(B)/float(m)
  C = temporary(C)/float(m)
  Pj = temporary(Pj)/float(m)
  Pn = Pn/float(m)

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

; subtract 1/M bias (see references)
  bicoh[mask] = bicoh[mask] - 1.0/float(m)

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

; return this if user requested image (2-d bicoherence)
  if keyword_set(img) then return, bicoh

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
