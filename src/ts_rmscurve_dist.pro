PRO TS_RMSCURVE_DIST, filename, rmsplot=rmsplot, log=log, ps=ps

; ----------------------------------------------------------
;+
; NAME:
;       TS_RMSCURVE_DIST
;
; PURPOSE:
;       Plot the flux/rms distribution from TS_RMSCURVE
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       TS_RMSCURVE_DIST,'rmscurve.list'
;
; INPUTS:
;       filename  - (string) name of file listing FITS time series
;
; OPTIONAL INPUTS
;       rmsplot   - (logical) set to plot distribution
;       log       - (logical) use logarithmic x-axis
;       ps        - (logical) send output to PostScript file
;
; DETAILS:
;       The filename given on input should point to an
;       ASCII file listing the names of a series of files
;       produced by TS_RMSCURVE. These will be read into
;       memory, concatonated into a single dataset and the
;       histogram of fluxes (or rms) will be plotted
;
;       If the RMSPLOT keyword is set this will plot the
;       binned rms histogram instead of the fluxes.
;
; Example call:
; IDL> TS_RMSCURVE_DIST,'rmscurve.list',/rms
;
; PROCEDURES CALLED:
;       READ_TABLE, WRITE_TABLE, PLOT_ERR, PS_OPEN, PS_CLOSE
;
; HISTORY:
;       08/05/07  - v1.0 - first working version
;
; NOTES:
;
;-
; ----------------------------------------------------------

; watch out for errors
  on_error,2

; ----------------------------------------------------------
; Check the arguments

; check filename has been set
  if (n_elements(filename) eq 0) then begin
      print,'** FILENAME not set in TS_RMSCURVE_DIST'
      return
  endif

; ----------------------------------------------------------
; Load the files 

; Load list of FITS file names from the file
  filelist = read_table(filename,/text)

; Determine how many FITS files are to be processed
  N = n_elements(filelist)
  print,'-- Files to process',N

; Loop over each file to process
  for i=0,N-1 do begin

      print,'-- Loading file',i+1

; Read the file
      x = read_table(filelist[0,i],/double)

; add to the data arrays

      if (i eq 0) then begin
          flux = reform(x[1,*])
          rms  = reform(x[2,*])
          flux = flux/mean(flux)
          rms  = rms/mean(rms)
      end else begin
          iflux = reform(x[1,*])
          irms  = reform(x[2,*])
          iflux = iflux/mean(iflux)
          irms  = irms/mean(irms)
          flux = [flux,iflux]
          rms  = [rms,irms]
      endelse

; End loop over files

  endfor

; total number of data points
  n = n_elements(flux)
  print,'-- No data:    ',n

; ----------------------------------------------------------
; Calculate empirical density function

  if keyword_set(rmsplot) then begin
      x = rms
  endif else begin
      x = flux
  endelse

; convert to log units if needed
  if keyword_set(log) then x = alog10(x)

; make historgram
  binwidth=0.01
  hist = lin_rebin(x,binwidth=binwidth,bin_dy=err,bin_x=bin,/hist,minbin=25)

; normalise histogram from frequency to density
  norm = total(hist)*binwidth
  hist = hist/norm
  err = err/norm

; ----------------------------------------------------------
; Fit the empirical density function with model

; remove zero point
  mask = where(bin ne 0)
  hist = hist[mask]
  bin = bin[mask]
  n = n_elements(bin)

; fit model to data
  weights = 1.0/hist
  a = [-0.037,0.27,0.0]
  fita = [1,1,1]
  itmax = 20
  function_name = 'plognormal'

; fit model to distribution data
  model = curvefit(bin, hist, weights, a, chisq=chisq, /double, $
                   /noderivative, iter=iter, itmax=itmax, status=status, fita=fita, $
                   function_name = function_name)

; display outcome of fitting
  case status of
      1: print,'** Fit not converging in TS_RMSCURVE_DIST'
      2: print,'** ITMAX reached prior to convergence in TS_RMSCURVE_DIST'
      else: print,'-- Iterations used:',iter
  endcase

; calculate residuals
  res = (hist-model)/err

; display results of fit
  chisq = total(res*res)

; calculate degrees-of-freedom
  dof = n - total(fita)

; and p-value for chi-square
  p = 1.0 - chisqr_pdf(chisq,dof)

  print,'-- Chi-square: ', chisq
  print,'-- DOF:        ',dof
  print,'-- p-value:    ',p
  print,'--       mu               sigma          t'
  print,'--',a

; ----------------------------------------------------------
; Make 2-panel plot

; open the PS device if requested
  if keyword_set(ps) then ps_open

; ------ panel 1 ------

; define X-axis range 
  xoffset=1.1
  xrange=[0,max(bin)*xoffset]
  if keyword_set(log) then begin
      xrange[0] = min(bin)/xoffset
  endif

; define Y axis range 
  yrange=max(hist)-min(hist)
  yoffset=1.1
  yrange=[-0.05*yrange,max(hist)*yoffset]

; define axis label
  if keyword_set(rmsplot) then begin
      xtitle = "normalised rms"
  endif else begin
      xtitle = "normalised flux"
  endelse

; define plotting symbol (filled circle)

  A = FINDGEN(17) * (!PI*2/16.)  
  USERSYM, COS(A), SIN(A), /FILL  

; Plot histogram
  plot,bin,hist,psym=10,symsize=0.5,position=[0.15,0.35,0.95,0.95], $
    xrange=xrange, xstyle=1, yrange=yrange, ystyle=1, $
    ytitle="density", xtickname=REPLICATE(' ', 30)

; with errors
  plot_err,bin,hist,err

; plot model
  plots,bin,model,color=150,thick=2.0

; ------ panel 2 ------

; plot residuals
  yrange=[-4.5,4.5]
  plot,bin,res,psym=10,symsize=0.5,position=[0.15,0.15,0.95,0.35], $
    xrange=xrange, xstyle=1, yrange=yrange, ystyle=1, $
    ytitle="(x-!7l!X)/!7r!X", xtitle=xtitle, /noerase, xticklen=0.1

; with errors (+/-1 sigma)
;  plot_err,bin,res

; draw zero line
  plots,xrange,[0,0],linestyle=2

; close the PS device if requested
  if keyword_set(ps) then ps_close

; ----------------------------------------------------------
; Return control to the user

  return

END
