PRO TS_RMSCURVE_PLOT, filename, rmsflux=rmsflux, ps=ps

; ----------------------------------------------------------
;+
; NAME:
;       TS_RMSCURVE_PLOT
;
; PURPOSE:
;       Plot the output from TS_RMSCURVE
;
; AUTHOR:
;       Simon Vaughan (U.Leicester)
;
; CALLING SEQUENCE:
;       TS_RMSCURVE_PLOT,'rmscurve.list'
;
; INPUTS:
;       filename  - (string) name of file listing FITS time series
;
; OPTIONAL INPUTS
;       rmsflux   - (logical) set to plot rms-flux correlation
;       ps        - (logical) send output to PostScript file
;
; OUTPUTS:
;       <files>   - series of ASCII files 
;
; DETAILS:
;       The filename given on input should point to an
;       ASCII file listing the names of a series of files
;       produced by TS_RMSCURVE. These will be read into
;       memory, concatonated into a single dataset and the
;       two time series (one for mean flux, one for rms)
;       will be plotted.
;
;       If the RMSFLUX keyword is set this will plot the
;       binned rms vs. flux relation instead of the two
;       data vectors as time series.
;
; Example call:
; IDL> TS_RMSCURVE_PLOT,'rmscurve.list',/rmsflux
;
; PROCEDURES CALLED:
;       READ_TABLE, WRITE_TABLE, PLOT_ERR, PS_OPEN, PS_CLOSE
;
; HISTORY:
;       27/04/07  - v1.0 - first working version
;       02/05/07  - v1.1 - added RMSFLUX and PS options
;       15/02/09  - v1.2 - added sorting of data into time 
;                          order, and saving of the flux and
;                          rms data as time series in ASCII file.
;
; NOTES:
;       + Add deadtime correction (to TS_RMSFLUX)
;
;-
; ----------------------------------------------------------

; watch out for errors
;  on_error,2

; ----------------------------------------------------------
; Check the arguments

; check filename has been set
  if (n_elements(filename) eq 0) then begin
      print,'** FILENAME not set in TS_RMSCURVE_PLOT'
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
          time = reform(x[0,*])
          flux = reform(x[1,*])
          rms  = reform(x[2,*])
      end else begin
          itime = reform(x[0,*])
          iflux = reform(x[1,*])
          irms  = reform(x[2,*])

          time = [time,itime]
          flux = [flux,iflux]
          rms  = [rms,irms]
      endelse

; End loop over files

  endfor

; subtract the start time

  time = time - min(time)

; sort into time order ???

  indx = sort(time)
  time = time[indx]
  flux = flux[indx]
  rms  = rms[indx]

  sdata = transpose([[time],[flux],[rms]])
  write_table,sdata,"temp.lc"

; ----------------------------------------------------------
; Plot the output

; define plotting symbol (filled circle)

  A = FINDGEN(17) * (!PI*2/16.)  
  USERSYM, COS(A), SIN(A), /FILL  

; open the PS device if requested
  if keyword_set(ps) then ps_open

  if (keyword_set(rmsflux) eq 0) then begin

; make plot of two time series (flux and rms)

; define X-axis range (TIME)
      xoffset=0.1
      xspan=max(time)-min(time)
      xrange=[min(time)-xoffset*xspan,max(time)+xoffset*xspan]

; define Y axis range (FLUX)
      yoffset=0.2
      yspan=max(flux)-min(flux)
      yrange=[min(flux)-yoffset*yspan,max(flux)+yoffset*yspan]

; Plot FLUX time series
      plot,time,flux,psym=8,symsize=0.5,position=[0.15,0.55,0.95,0.95], $
        xtickname=REPLICATE(' ', 30), noclip=0, xrange=xrange, xstyle=1, $
        yrange=yrange, ystyle=1, ytitle="flux (ct s!U-1!N)"

; Add a legend
      xpos=min(time)-xspan*0.0
      ypos=max(flux)-yspan*0.0
      xyouts,xpos,ypos,'flux'

; define Y axis range (RMS)
      yspan=max(rms)-min(rms)
      yrange=[min(rms)-yoffset*yspan,max(rms)+yoffset*yspan]

; Plot RMS time series
      plot,time,rms,psym=8,symsize=0.5,position=[0.15,0.15,0.95,0.55],/noerase, $
        noclip=0, xrange=xrange, xstyle=1, yrange=yrange, ystyle=1, $
        xtitle="time (s)", ytitle="rms (ct s!U-1!N)"

; Add a legend
      xpos=min(time)-xspan*0.0
      ypos=max(rms)-yspan*0.0
      xyouts,xpos,ypos,'2-20 Hz rms'

  endif else begin

; bin according to flux
      bin_r = lin_rebin(flux,rms,binwidth=100,minbin=50,/varerr, $
                        bin_x=bin_f, bin_dy=bin_dr)

; define X-axis range (FLUX)
      xoffset=0.1
      xspan=max(bin_f)-min(bin_f)
      xrange=[min(bin_f)-xoffset*xspan,max(bin_f)+xoffset*xspan]

; define Y axis range (RMS)
      yoffset=0.1
      yspan=max(bin_r)-min(bin_r)
      yrange=[min(bin_r)-yoffset*yspan,max(bin_r)+yoffset*yspan]
      yrange[0]=0.0

; Plot the binned rms-flux correlation
      plot,bin_f,bin_r,xtitle="flux (ct s!U-1!N)",ytitle="rms (ct s!U-1!N)", $
        psym=8,symsize=0.5,position=[0.15,0.15,0.95,0.95], $
        xrange=xrange, xstyle=1, yrange=yrange, ystyle=1

; add error bars (from variance within each flux bin)
      plot_err,bin_f,bin_r,bin_dr

; Add a legend
      xpos=min(bin_f)-xspan*0.0
      ypos=max(bin_r)-yspan*0.0
      xyouts,xpos,ypos,'2-20 Hz rms'

  endelse

; close the PS device if requested
  if keyword_set(ps) then ps_close

; ----------------------------------------------------------
; Return control to the user

  return

END
