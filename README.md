# IDL
IDL codes (mostly legacy)

These routines are made available for you to use as you wish. 
They have all been tried and tested on IDL 8.1, to a greater or lesser extent, but undoubtedly some bugs will remain. 
All routines include some (standard format) header documentation. 

## Data analysis, statistics

lin_rebin.pro	Rebin one dimensional data into equal width bins

log_rebin.pro	Rebin one dimensional data into equal logarithmic bins

regroup.pro	Rebin (regroup) binned/spaced one dimensional data

chi2_pdf.pro	Compute density for chi-sq distribution

clip_range.pro	Find the lower and upper per cent bounds on input data

excess_variance.pro	Compute the 'excess' variance of a sample

lincorr.pro	Linear correlation coefficient with p-value

lognormal.pro	Compute density of lognormal function

normal.pro	Draw array of random deviates from multidimensional Gaussian

normal_pdf.pro	Compute density for normal distribution

parabound.pro	Calculate confidence bounds on quadratic model

randomp.pro	Generate array of random deviates with power law distribution

summary.pro	Standard five number summary of an array.

## Utilities

int_lorentz.pro	Calculate the definite integral of a Lorentzian function

int_powerlaw.pro	Find definite integral of power law function

seq.pro	Generate sequence of evenly spaced real numbers

arrow_plot.pro	Plot scatter diagram using arrows of given gradients

plot_err.pro	Plot error bars on x-y plot (in y and x).

plot_hist.pro	Make a nice histogram plot

ps_close.pro	Close PS device and return to previous settings

ps_open.pro	Open and configure a PS device for graphical output

read_table.pro	Read an ASCII table into a data array.

splot.pro	Produce an NxN matrix of scatter plots

write_table.pro	Write a data array to an ASCII table file.

xspec_out.pro	Write data to file ready for XSPEC

## Time series analysis

bi_spectrum.pro	Calculate the bi-coherence of a univariate time series.

cross_spectrum.pro	Calculate coherence and phase lags of two time series

dynamic_psd.pro	Calculate the dynamic (time-resolved) power spectrum.

endmatch.pro	Perform 'end matching' of a data series

extract_xte_dps.pro	Extract and export dynamical power spectra from RXTE

periodogram.pro	Calculate periodogram (mod-square Fourier transform) of a univariate time series

rms_mean_fit.pro	Calculate, plot and fit the rms-mean relation given input time series

ts_gen.pro	Generate a random time series from a power spectrum model

ts_multispec.pro	Calculate average periodogram from list of FITS files

ts_readfits.pro	Read two columns data (e.g. time, flux) from a FITS file

ts_rmscurve.pro	Calculate time series for mean flux and rms

ts_rmscurve_dist.pro	Plot the flux/rms distribution from TS_RMSCURVE

ts_rmscurve_plot.pro	Plot combined output from TS_RMSCURVE

ts_rmsflux.pro	Calculate mean and rms at even intervals in a series

ts_segment.pro	Define contiguous segments of an input series

ts_spec_plot.pro	Plot a binned periodogram, from e.g. TS_MULTISPEC

## XMM-Newton tasks [XMM_Extract.pdf](XMM_Extract.pdf)

rgs_plot.pro	Plot RGS data in velocity space, for multiple lines and/or datasets

make_epic_lc.pro	Extract time series from XMM/EPIC event data

load_events.pro	Read XMM EPIC data from filtered events file (FITS)

gti_check.pro	Check for patches of 'bad' data in multi-observation time series

gti_fix.pro	Correct for gaps (GTI) in EPIC time series by interpolation

energy_bands.pro	Calculate logarithmically spaced PI bands

epic_load_lc.pro	Load into an array multiple XMM/EPIC time series

epic_segment.pro	Take output from EPIC_LOAD_LC and divide into equal length segments
