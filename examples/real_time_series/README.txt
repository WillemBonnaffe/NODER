@ real_time_series
Contains scripts for the NODE analysis of the hare and lynx real time series with time dependence
(3 oscillations).

@ 1-functions.R
- contains all function definitions for the fit, analysis, and
visualisation, of a 2-dimensional NODE system approximating p.c. growth rates of hare and lynx 

@ 2-initiate.R
- Loads functions
- Loads real time series from (hare_lynx_data.csv)

@ 3-optimise.R
- runs 2-initiate
- trains a pilot ensemble of 200 elements by running a BFGS optimisation
  routine starting from random parameters draw 
- stores pilot ensemble in out/ folder

@ 4-analysis.R
- filter 25% best elements of the pilot ensemble
- compute NODE predictions 
- compute direct effects and contributions to hare and lynx p.c. growth rates
- visualise predictions, effects, and contributions (figure 4)
- visualise p.c. growth rates with different hare and lynx density levels (figure 5)
- stores figures (4 and 5) in out/plots/ folder

@ main.R 
- runs and analyses the ensemble (scripts 3 and 4)

@ out/ensemble.csv
