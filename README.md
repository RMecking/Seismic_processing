# Seismic_processing
misc. scripts (python/ matlab) for seismic/ GPR data manipulation

-----------------------------------------------------------------
Traceinterp2D:
Doubles the traces in a dataset and divides the trace distance by 2. Usefull in cases where spatial aliasing prohibits the effective use of FK-filtering or similar. At the moment the itnerpolated tarces are inserted exactly in the centre between two adjacent traces.
Needs segy-input via the SU-matlab toolbox but can be easiy modified to accept similar data. Output SU-files.

How it works: A window with a window length dependent on the maximum wavelength in the dataset is moving along the traces and calculates the time shift of maximum correlation between two traces for this given window. The two traces are shifted by half this time towards each other and the mean value of the two traces at zero lag is calculated. This is done for every time window along the traces. Due to the overlapping of the time window, manipulated by the shift of the window per iteration ((@param: gapsamp), multiple interpolated values are generated for each sample. The mean value of this sample creates the interpolated sample. Increasing the gap between time windows (gapsamp) decreases the number of values per sample giving way to errors due to false maximum correlation time shifts. It does decrease however the computation time significantly. I found ~ 10 values/ per sample to be a sufficient amount. 
NOTE for future updates: The usage of a standard deviation criterium may be used to elimit outliers and decrease the number of time windows.

------------------------------------------------------------------
