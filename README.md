HypernovaeGit
Description of the files contained in HypernovaeGit
=============

Bremsstrahlung.ipynb

IJulia notebook which serves as a test ground for the different bremstrahlung functions and features. Test all new bremsstrahlung contributions here. 
-----------------------------------

Source\ Physics.ipynb 

Ijulia notebook which serves as a test ground for the source physics of a shock produced by a hypernova event. Test all new source physics contributions here. 
-----------------------------------

TestHypernovae.c

??? 
-----------------------------------

brem_funcs.jl 

Contains all of the relevant functions for calculating the emissivity of a plasma given its temperature and the frequency of light that will be emitted by the plasma. Both relativistic and non-relativistic plasmas (with respect to the electron rest mass) are contained here. The former is calculated using the integration of the functions found in Rybicki Lightman , while the later uses the analytical approximation found in the same text. 

Results are output with all normalizations except for  the factor (Z^2 n_e n_i) which is the ion charge, electron density, and ion density respectively). 
-----------------------------------

run_brem.jl 

Contains code that will produce a .dat file containing the bremsstrahlung emissivity of a plasma for a range of frequencies (Hz) and temperatures (K). The emissivities can be calculated using a parallel machine by calling 

julia -p n run_brem.jl 

where n is the number of ADDITIONAL processors to be created. Results are stored in a file "parallel_bremsstrahlung_emissivity.dat" with each line corresponding to the emissivity at a given temperature. Timing results are also given. 
