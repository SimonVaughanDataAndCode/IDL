PRO LOAD_EPIC

  pi_list = ENERGY_BANDS(15)
  rootpath = "/data/49/sav2/xmm/ngc4051/events/"

  N_obs = 17
  N_inst = 3
  scale = MAKE_ARRAY(N_obs, N_inst)
  scale_name = ['scaling.txt', 'scaling-m1.txt', 'scaling-m2.txt']

; loop over all instruments and load the src/bkg scaling factors

  for i = 0, N_inst - 1 do begin

      scale_data = READ_TABLE(rootpath+scale_name[i], head=1, /TEXT)
      src = FLOAT(scale_data[1, *])
      bkg = FLOAT(scale_data[2, *])
      mask = WHERE(src gt 0.0 AND bkg gt 0.0, count)
      scale[mask, i] = src[mask] / bkg[mask]
      obs_names = scale_data[0, *]
 
  endfor

  j = 1
  x = MAKE_EPIC_LC(obs_name[j], scale[j, *], rootpath=rootpath,  dt=100.0, pi_list=pi_list, time=time, err=err, /SILENT)

  plot, time, x[0,*], psym=10, /xstyle
  plot_err, time, x[0,*], err[0,*]

END
