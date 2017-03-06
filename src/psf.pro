FUNCTION PSF, npsf, r0=r0, r1=r1, pscale=pscale, $
              diff=diff, eff=eff, tri=tri

; ----------------------------------------------------------
;+
; NAME:
;       PSF
;
; PURPOSE:
;       Generate a simple Point Spread Function (PSF) image
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       img = PSF(100,r0=4,r1=100,pscale=0.5)
;
; INPUTS:
;       npsf      - (integer) size of output array (npsf^2)
;
; OPTIONAL INPUTS:
;       r0        - (scalar) short axis of PSF (arcmin)
;       r1        - (scalar) long axis of PSF (arcmin)
;       pscale    - (scalar) pixel scale (arcmin/pixel)
;       eff       - (vector) 3 elements listings the
;                            relative effective areas of
;                             core, both cross arms, and diffuse
;       tri       - (logical) if set then use triangular 
;                             cross-section for unfocussed light
;
; OUTPUTS:
;       img       - (array) npsf*npsf image array
;
; OPTIONAL OUTPUTS:
;       diff      - (scalar) average diffuse light brightness
;                            in units of (pixel^-1)
;
; DETAILS:
;       Generate an image of a simple Point Spread Function.
;       The PSF is a 2-dimensional shape comprising the
;       sum of four components:
;        1. central focus
;        2. vertical cross-arm
;        3. horizontal cross-arm
;        4. diffuse (unfocussed) light
;       We have approximated (1) as a cone with a 
;       half width of R0 (in units of arcminutes).
;
;       We have approximated the cross-arms as 
;       triangular in shape with a half width R1 in the
;       long axis and half width R0 in the short axis.
;
;       The diffuse light is spread over the entire
;       image, and is returned as a per pixel value
;       in the output keyword DIFF. If TRI = TRUE
;       then the diffuse light has a pyramid shape,
;       otherwise is constant within the box 2*R1 on
;       a side.
;
;       The pixel scale is given by the PSCALE keyword,
;       in units of arcminute/pixel. The functions
;       used to generate the PSF are integrated over
;       each pixel using the QSIMP routine from the
;       IDL Astronomy Library. The triangular function
;       is specified by the TRIANGLE function that has
;       has one dependent variable and one parameter
;       (R - which defines with full width).
;
;       The optional input array EFF defines the relative
;       effective areas of the core, the combined cross-arms,
;       and the diffuse light. Default is 1.0,2.0,4.0, meaning
;       the core contains 1/7 of the light, each cross-arm
;       also 1/7 and the diffuse light 4/7.
;
; PROCEDURES CALLED:
;       QSIMP.PRO, TRAPZD.PRO, ZPARCHECK.PRO, TRIANGLE.PRO
;
; HISTORY:
;       19/07/07  - v1.0 - first working version
;       01/08/07  - v1.1 - default value of R0=2 (after NPB)
;                          added diffuse light evenly over PSF
;       02/08/07  - v1.2 - Added diffuse light only in square
;                          that is (2*r1) on a side.
;       03/08/07  - v1.3 - Added TRI keyword for diffuse pyramid
;       08/08/07  - v1.4 - Corrected centre position (thanks AMR)       
;
; NOTES:
;       none
;
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments

; if NPSF not defined, set default
  if (n_elements(npsf) eq 0) then npsf=1024

; if R0 not defined, set default
  if (n_elements(r0) eq 0) then r0=2

; if R1 not defined, set default
  if (n_elements(r1) eq 0) then r1=150

; if PSCALE not defined, set default
  if (n_elements(pscale) eq 0) then pscale=0.5

; if EFF not defined, set default
  if (n_elements(eff) ne 3) then eff=[1.0,2.0,4.0]

; ----------------------------------------------------------
; Main part of procedure

; convert R0, R1 from arcminutes to pixels

  pr0 = r0/pscale
  pr1 = r1/pscale

; create an empty image for the output PSF

  psf = fltarr(npsf,npsf)

; define the centre of the image in X, Y

  cenx = floor(npsf/2.0) + 0.5
  ceny = floor(npsf/2.0) + 0.5

; define index array [0,1,...,npsf-1]

  indx = indgen(npsf)

; define x as distance from centre in pixels

  x = indx - cenx

; Define two triangular functions of x
;   f_short - triangle of width R0
;   f_long  - triangle of width R1
; These are computed by integrating triangular
; function TRIANGLE between pixels i and i+1
; using the QSIMP integration routine.

  f_short = fltarr(npsf)
  f_long = fltarr(npsf)
  for i = 0L,npsf-2 do begin

      if (x[i] gt -pr0) and (x[i+1] lt pr0) then begin
          qsimp, 'triangle', x[i], x[i+1], s, r=pr0
      endif else begin
          s = 0.0
      endelse
      f_short[i] = s

      if (x[i] gt -pr1) and (x[i+1] lt pr1) then begin
          qsimp, 'triangle', x[i], x[i+1], s, r=pr1
      endif else begin
          s = 0.0
      endelse
      f_long[i] = s

  endfor

; generate two images each of size NPSF*NPSF
; one contains all the X values and one contains
; all the Y values. 
;
;  IX = 0 1 2 3 ...   IY = 0 0 0 0 ...
;       0 1 2 3 ...        1 1 1 1 ...
;       0 1 2 3 ...        2 2 2 2 ...
;       0 1 2 3 ...        3 3 3 3 ...

  ix = rebin(indgen(npsf),npsf,npsf,/sample)
  iy = transpose(ix)

; generate image of PSF core using
; PSF(x,y) = f_short(x)*f_short(y)
; normalise to integral = eff[0]

  psf_core = f_short(ix) * f_short(iy)
  psf_core = psf_core/total(psf_core)*eff[0]

; generate image of PSF vertical cross-arm
; PSF(x,y) = f_long(x)*f_short(y)
; normalise to integral = 1

  psf_vert = f_long(ix) * f_short(iy)
  psf_vert = psf_vert/total(psf_vert)*eff[1]*0.5

; generate image of PSF horizontal cross-arm
; PSF(x,y) = f_short(x)*f_long(y)
; normalise to integral = 1

  psf_horz = f_short(ix) * f_long(iy)
  psf_horz = psf_horz/total(psf_horz)*eff[1]*0.5

; compute final PSF of focussed light as
; sum of these three components

  psf = psf_core + psf_vert + psf_horz

; add diffuse light evenly in square of size (2*R1)x(2*R1)

  if keyword_set(tri) then begin

;     If TRI=TRUE then forms square-base pyramid.

      cenj = npsf/2 - 0.5
      j = triangle(indgen(npsf) - cenj,r=pr1)
      x = rebin(j,npsf,npsf,/sample)
      z = (x < transpose(x))
      z = z/total(z) 
      psf = psf + z * eff[2]

  endif else begin

;     If TRI=FALSE then this is spread evenly over square

      nn = 4*pr1*pr1
      px0 = ( (cenx - pr1) > 0)
      px1 = ( (cenx + pr1) < npsf-1)
      py0 = ( (ceny - pr1) > 0)
      py1 = ( (ceny + pr1) < npsf-1)
      psf[px0:px1,py0:py1] = psf[px0:px1,py0:py1] + eff[2]/nn

  endelse

; normalise so that sum[PSF] = 1

  psf = psf/total(psf)

; ----------------------------------------------------------
; Return the data array to the user

  return, psf

END
