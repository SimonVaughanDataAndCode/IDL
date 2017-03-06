FUNCTION normal, SEED, COVAR=COVAR, $
                 N=N, mean=mean, samp_cov=samp_cov

; ----------------------------------------------------------
;+
; NAME:
;       NORMAL
;
; PURPOSE:
;       Draw array of random deviates from multidimensional Gaussian
;
; AUTHOR:
;       Simon Vaughan (U.Leicester) 
;
; CALLING SEQUENCE:
;       x = normal(seed,covar,100)
;
; INPUTS:
;       seed          - (integer) seed for random number generator
;       covar         - (array) M*M covariance matrix for Gaussian
;       n             - (integer) number of draws to make
;
; OPTIONAL INPUTS:  
;       mean          - (vector) vector of mean values
;
; OUTPUTS:
;       x             - (array) M*N array of random deviates
;
; OPTIONAL OUTPUTS:  
;       samp_cov      - (array) M*M sample covariance matrix
;
; DETAILS:
;       Draw N vectors (of length M) from an M-dimensional
;       normal (Gaussian) distribution with a covariance matrix 
;       given by COVAR (which much be an M*M symmetric array].
;
;       The output is an M*N dimensional array. Each of
;       the N 'rows' is an M-dimensional 'column' vector
;       drawn from the M-dimentional normal.
;
;       Uses the Cholesky method. We first generate an
;       M*N array of independent Gaussian deviates with the
;       RANDOMN procedure. Then apply a matrix transformation to 
;       this such that the transformed array has the required
;       covariance matrix. 
;
;       If MEAN is given on input, then the output deviates
;       have a mean (vector) of MEAN and a covariance matrix
;       of COVAR. Otherwise MEAN = 0.
;
;       If requested, also compute the sample covariance matrix
;       from the random deviates.
;
; EXAMPLE USAGE:
;
;       IDL> covar = [[2.0,1.2],[1.2,5.0]]
;       IDL> x = normal(seed,covar=covar,n=500)
;       IDL> plot,x[0,*],x[1,*],psym=1,xrange=[-10,10], $
;                         /xstyle,yrange=[-10,10],/ystyle
;
; METHOD:
;       The Cholesky methods works as follows:
;
;       1. Cholesky decomposition (matrix square root) of the
;          covariance matrix S is C (C is a lower left triangle)
;
;              S = C C^T  
;
;        2. Generate an array of standard normal deviates Z
;
;              Z ~ N(0,I)
;
;            where 0 is the zero mean and I is the identity matrix
;
;        3. Apply the matrix multiplication 
;
;              X = C Z
;
;        4. Add the mean [X' = X + mu] if necessary
;
;        The result should be distributed as X ~ N(mu,S)
;
;        To see that this works consider that the covariance
;        matrix of X (zero mean) is the expectation value:
;
;        cov[X] = E[X X^T] = E[C Z (C Z)^T] = E[C Z Z^T C^T]
;
;        using X=CZ. Now by construction the covariance matrix
;        of Z = I since we generated it from Z ~ N(0,I). Therefore
;
;        cov[Z] = E[Z Z^T] = I
;
;        which means cov[X] = E[C I C^T] = E[C C^T] = C C^T = S
;
;        And so we have our intended covariance matrix. 
;
; HISTORY:
;       03/12/2007  --  v1.0  -- first working version
;
; NOTES:
;       The peculiar [column,row] notation of IDL means
;       that some of the matrix algebra above looks like
;       it has been implimented backwards in the code below. 
;       
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)
  COMPILE_OPT idl2  
  
; watch out for errors
  on_error,2

; ---------------------------------------------------------
; Check arguments

; if N is not set, use default
  if n_elements(n) eq 0 then n=1

; make sure enough arguments were input
  if n_params() lt 1 then begin
      print, '** Missing argument in NORMAL(seed,[covar=],[n=])'
      return, 0
  endif

; if COVAR is not set, use default.
  if n_elements(covar) eq 0 then covar=[[1,0],[0,1]]

; check dimensions of COVAR match a matrix
  if (size(covar,/N_DIMENSIONS)) ne 2 then begin
      print, '** COVAR is not a matrix in NORMAL.'
      return, 0
  endif

; check shape of COVAR matrix
  m = (size(covar))[1]
  if (size(covar))[2] ne m then begin
      print, '** COVAR is not a square matrix in NORMAL.'
      return, 0
  endif

; if MEAN is not set then use default
  if n_elements(mean) eq 0 then mean = make_array(m,value=0.0d)

; check if right size
  if n_elements(mean) ne m then begin
      print, '** MEAN vector length does not match COVAR size in NORMAL.'
      return, 0
  endif

; ---------------------------------------------------------
; Main routine

; generate N*M array of normal deviates, Z
  Z = randomn(seed, m, n)

; use dummy array for storage of covariance matrix
  C = covar

; perform Cholesky decomposition of covariance matrix
  la_choldc, C, /double           

; remove the upper-right triangle part of matrix
; (which is left over from COVAR) leaving only 
; lower-left and diagonal elements non-zero.
  ii = rebin(indgen(m),m,m)
  jj = transpose(ii)
  mask = where(ii gt jj)
  C[mask] = 0.0d

; apply transformation to get indended covariance matrix
  x = matrix_multiply(C,Z,/ATRANSPOSE)

; add the mean values to the deviates
  x = x + rebin(mean,m,n)

; compute sample covariance matrix from X
  if arg_present(samp_cov) then $
    samp_cov = correlate(x,/covariance,/double)

; ---------------------------------------------------------
; Return to user

  return,x

END
