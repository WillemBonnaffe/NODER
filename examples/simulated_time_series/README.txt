@ simulated_time_series
Contains scripts for the NODE analysis of the hare and lynx simulated time series
(3 oscillations).

@ 1-functions.R
- contains functions to simulate time series from Lotka-Volterra equations
- contains all function definitions for the fit, analysis, and
visualisation, of a 2-dimensional NODE system approximating p.c. growth rates of hare and lynx 

@ 2-initiate.R
- Loads functions
- simulates time series from Lotka-Volterra equations

@ 3-optimise.R
- runs 2-initiate
- trains a pilot ensemble of 200 elements by running a BFGS optimisation
  routine starting from random parameters draw 
- stores pilot ensemble in out/ folder

@ 4-analysis.R
- filter 25% best elements of the pilot ensemble
- compute NODE predictions 
- compute direct effects and contributions to hare and lynx p.c. growth rates
- visualise predictions, effects, and contributions
- visualise p.c. growth rates with different hare and lynx density levels
- stores figures (2 and 3) in out/plots/ folder

@ main.R 
- runs and analyses the ensemble (scripts 3 and 4)