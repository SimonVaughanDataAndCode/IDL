FUNCTION DENSITY,x,y, plot=plot, $
                 minx=minx, maxx=maxx, nx=nx, dx=dx, $
                 miny=miny, maxy=maxy, ny=ny, dy=dy, $
                 binx=binx, biny=biny, $
                 _EXTRA=EXTRA_KEYWORDS
; ----------------------------------------------------------
;+
; NAME:
;       DENSITY
;
; PURPOSE:
;       Empirical one or two dimensional densities
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       dens = density(x,y)
;
; INPUTS:
;       x        - one dimensional data array
;
; OPTIONAL INPUTS:
;       y        - (array) one dimensional data array
;       minx     - (float) starting value for binning along x 
;       maxx     - (float) finishing value for binning along x 
;       miny     - (float) starting value for binning along y
;       maxy     - (float) finishing value for binning along y 
;       nx       - (integer) number of bins along x axis
;       ny       - (integer) number of bins along y axis
;       dx       - (float) bin size along x axis
;       dy       - (float) bin size along y axis
;       plot     - (logical) produce a simple plot?
;
; OUTPUTS:
;       dens     - one or two dimensional array of densities
;
; OPTIONAL OUTPUTS:
;       binx     - locations of bins along x axis
;       biny     = location of bins along y axis
;
; DETAILS:
;       Calculates empirical density distribution of
;       one or two dimensional data. The output is
;       one or two dimensional depending on the input.
;       If only array X is input then work in 1d
;       and produce the histogram (density).
;       If X and Y arrays are input then work in 2d
;       and calculate the joint histogram (density).
;
;       The MINX,MAXX [MINY,MAXY] parameters are
;       used for setting the upper,lower boundaries
;       of the histogram and will default to 
;       min(X),max(X), and similarly for Y.
;
;       The number of bins is determined from the 
;       bin width DX [DY] if given on input, otherwise
;       from NX [NY]. If no values are given default 
;       for NX [NY] is 11.
;
; EXAMPLE USAGE:
;       IDL> x = randomn(seed,1000)
;       IDL> y = randomn(seed,1000)
;       IDL> dens = density(x,y,nx=20,ny=20,/plot)
;       IDL> oplot,x,y,psym=1
;
; PROCEDURES CALLED:
;       SEQ
;
; HISTORY:
;
;
; NOTES:
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)
  COMPILE_OPT idl2, HIDDEN

; watch out for errors
  on_error,2

; ---------------------------------------------------------
; Check arguments

; is there enough data in X array

  n_x = n_elements(x)
  if (n_x lt 4) then begin
      print,'** Not enough data in DENSITY.'
      return,0
  endif

; compare size of X and Y arrays
; only use arrays of matching size

  n_y = n_elements(y)
  if (n_y eq 0) then begin
      dim = 1
  endif else begin
      dim = 2
      if (n_y lt n_x) then begin
          xx = x[0:n_y-1]
          yy = y
      endif
      if (n_y gt n_x) then begin
          xx = x
          yy = y[0:n_x-1]
      endif
      if (n_y lt n_x) then begin
          xx = x
          yy = y
      endif
  endelse

; set upper/lower ranges 

  if (n_elements(minx) eq 0) then minx = min(x,max=newmax)
  if (n_elements(maxx) eq 0) then begin
      if (n_elements(newmax) ne 0) then maxx = newmax
      if (n_elements(newmax) eq 0) then maxx = max(x)
  endif
  if (dim eq 2) then begin
      if (n_elements(miny) eq 0) then miny = min(y,max=newmax)
      if (n_elements(maxy) eq 0) then begin
          if (n_elements(newmax) ne 0) then maxy = newmax
          if (n_elements(newmax) eq 0) then maxy = max(y)
      endif
  endif

; check ranges make sense

  if (minx ge maxx) then begin
      print,'** MINX > MAXX in DENSITY'
      return,0
  endif
  if (dim eq 2) then begin
      if (miny ge maxy) then begin
          print,'** MINY > MAXY in DENSITY'
          return,0
      endif
  endif

; calculate ranges

  rangex = maxx - minx
  if (dim eq 2) then rangey = maxy - miny

; store old values in case they are updated

  if keyword_set(nx) then old_nx = nx
  if keyword_set(ny) then old_ny = ny

; determine NX [NY] from DX [DY] value if given

  if (n_elements(dx) eq 1) then nx = floor(rangex/dx) + 1
  if (dim eq 2) then begin
      if (n_elements(dy) eq 1) then ny = floor(rangey/dy) + 1
  endif

; make sure NX [NY] exists (either from input NX [NY] or DX [DY])
; if not, use default

  if (n_elements(nx) eq 0) then nx=11
  if (dim eq 2) then begin
      if (n_elements(ny) eq 0) then ny=11
  endif

; make sure DX [DY] exist

  if (n_elements(dx) eq 0) then dx = rangex/(nx - 0.5)
  if (dim eq 2) then begin
      if (n_elements(dy) eq 0) then dy = rangey/(ny - 0.5)
  endif

; ---------------------------------------------------------
; Main routine

  if (dim eq 1) then begin
      dens = histogram(x,binsize=dx,min=minx,max=maxx,locations=binx)
  endif

; produce the 2D density

  if (dim eq 2) then begin
      dens = hist_2d(x,y,bin1=dx,bin2=dy,min1=minx,min2=miny,max1=maxx,max2=maxy)
  endif

; normalise from frequency to relative frequency (density estimate)

  dens = dens/total(dens)

; calculate the bin locations

  binx = seq(minx,maxx,n=nx)+dx/2
  if (dim eq 2) then biny = seq(miny,maxy,n=ny)+dy/2

; make a plot

  if keyword_set(plot) then begin
 
      if (dim eq 1) then begin
          plot,binx,dens,psym=10,xrange=[minx,maxx],/xstyle,_EXTRA=EXTRA_KEYWORDS$
      endif

      if (dim eq 2) then begin

; determine dynamic range: large = logscale, small = linscale
          lo = min(dens[where(dens gt 0)])
          hi = max(dens)
          nlevels = floor(alog(total(dens))) + 1
          nlevels = max([nlevels,6])

          if (hi/lo gt 20) then begin
              levels = exp(seq(alog(lo),alog(hi),n=nlevels))
          endif else begin
              levels = seq(lo,hi,n=nlevels)
          endelse

          contour,dens,binx,biny,levels=levels,/xstyle,/ystyle,xrange=[minx,maxx],yrange=[miny,maxy],_EXTRA=EXTRA_KEYWORDS$
      endif

  endif

; ---------------------------------------------------------
; Return to user

  if (n_elements(old_nx) ne 0) then nx = old_nx
  if (n_elements(old_ny) ne 0) then ny = old_ny
  return,dens

END
