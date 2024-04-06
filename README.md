Data, model, and figure generation code to accompany the manuscript:

R.M. Strickland and M.D.Covington (2024). The formation of glacier dirt cones, submitted to Journal of Glaciology. 

Running the model scripts requires Landlab. Please go to landlab.github.io for download instructions. 

# Description of Contents:

## Scripts:

* time-dependent-experiments.py - Model script used for time-dependent cone development experiments

* time-dependent-experiments-3d.py - Same as above, but saves data for 3d plots; user must specify timesteps to save data

* volume-experiments.py - Model script used for debris volume vs. cone height experiments

* ell-experiments_Sc09.py - Model script used to examine cone geometry vs. characteristic length experiments with critical slope S_c = 0.9

* ell-experiments_Sc115.py - Model script used to examine cone geometry vs. characteristic length experiments with critical slope S_c = 1.15

* ell-experiments_Sc14.py - Model script used to examine cone geometry vs. characteristic length experiments with critical slope S_c = 1.4

* infinite-pits-transects.py - Model script used to examine the growth of cones from debris-filled pits of "infinite" depth

* time-dependent-experiments/time-dependent-plots.ipynb - iPython notebook for creating the time dependent cone development figures

* time-dependent-experiments/transect-plots.ipynb - iPython notebook for creating cone transect plots for the time dependent experiments

* 3d-plot/3dplots.ipynb - iPython notebook to create the 3d plot of cone development
  
* infinite-pits-transects/infinite-pit-analysis.ipynb - iPython notebook for creating infinite pit plots 

* ell-experiments/ell-experiment-plots.ipynb - iPython notebook for creating plots that examined cone geometry vs. the characteristic length

* volume-experiments/volume-plots.ipynb - iPython notebook for creating plot examining cone heights vs. initial debris volume

## Directories:

* time-dependent-experiments/ - Data from time-dependent experiments

* infinite-pits-transects/ - Data from the infinite pit experiments

* ell-experiments/ - Data from the cone geometry vs. characteristic length experiments

* volume-experiments/ - Data from debris volume vs. cone height experiments

* 3d-plot/ - Data from simulation used to create the 3d plot
