
data = READ_TABLE("for_simon.txt")

f = REFORM(data[1,*])
r = REFORM(data[2,*])
flux = REFORM(data[3,*])

r_abs = r/100.0*flux

plot, flux, r_abs, /xstyle, XRANGE=[0,2250], /ystyle, YRANGE=[0,1e5]

bin_y = LIN_REBIN(flux, r_abs, /SORT, MINBIN=50, bin_x = bin_x, $
                 /VARERR, bin_dy=bin_dy)

ps_open

plot, bin_x, bin_y, /xstyle, XRANGE=[1500,2500], /ystyle, YRANGE=[0,5e2], $
  xtitle="Flux (ct/s/PCU)", ytitle="QPO rms (ct/s/PCU)"

plot_err, bin_x, bin_y, bin_dy

ps_close

