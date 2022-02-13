[![DOI](https://zenodo.org/badge/317295428.svg)](https://zenodo.org/badge/latestdoi/317295428)

# Hydro-osmotic chain
**Hydro-osmotic** chain is a Python-based code for simulations of coarsening in a 1-dimensional network of hydroosmotic micro-lumens. It was developped to study the contributions of osmotic pressure, permeation and active pumping to the dynamics of micro-lumens, and to investigate the cynamical collective effects of a chain of micro-lumens.
The micro-lumens are modelled as pressurized cavities with concentration, connected by bridges of variable length. Water and solutes can be exchanged with the external medium via the membrane, through permeation or active pumping.
The code solves a system of coupled non-linear equations for the area and concentration dynamics of each micro-lumen, by computing the flux between them. Topological events such as collapse and coalescence are handled during the simulation. The simulation runs until one micro-lumen is left, either by absorbing all neighboring micro-lumens via coalescence or coarsening, or by global collapse of the micro-lumens of the chain.

Please cite the original publication if you use this code :
Le Verge--Serandour, M and Turlier, H, biorXiv, 2020


## Requirements
This code runs in Python3, but can be used with Python2.7 (not updated since Jan. 2020).
* Python3 or Python2.7
    * getpass
    * os
    * time
    * socket
    * subprocess
    * sys
    * warnings
* [Scipy](https://www.scipy.org/) :
* [Integrate](https://docs.scipy.org/doc/scipy/reference/integrate.html#module-scipy.integrate) : sub-package for numerical integration.
    * [Numpy](https://numpy.org) : arrays ad algebraic operations on arrays
    * [Matplotlib](https://matplotlib.org) : graphical representation
* os, sys, time, subprocess, warnings, getpass, socket

## Content
* _ressources/ : contains the ressources used to simulate a chain of lumens
	* analysis_tools/
		* analysis.py : file used to extract exploitable data (concentrations, pressures, ...) of the chain ;
		* distribution.py : library of functions to compute the size distribution of lumens.
	* chain.py : main file, script used to run a simulation of a chain, given a .conf file ;
        * compress.py : script compressing directories in .tar.xvf ;
	* configreader.py : script that read configuration files (.conf) ;
	* flux.py : library used to perform flux calculations ;
	* gen_config.py : script that generates configurations files (.conf) from templates (.conf.tpl) ;
	* lumen_class.py : library that contains the definition of the objects chain, lumen and bridge ;
	* network.py : library of functions to work on graphs (nearest neighbors, etc) :
	* submit.py : script sending simulations to a cluster using SLURM language ;
	* tools.py : library of diverse useful functions
	* topology.py : library of functions handling the topology of the graph

* _configfiles/ : contains some configuration files for example or frequently used.

* _notebook/ : contains notebooks used to plot figures

## Instructions
### Quick launch
1. Go to the _data folder `cd ~/git/hydroosmotic_chain/_data`
2. Change the parameters that you what in the configuration file (config.conf)
3. Run the command `../_ressources/chain.py config.conf`
4. Outputs will be stored in the out/ folder.

### Options of simulation
The list of options is accessible using chain.py -h. This print the help of the script.
1. chain.py config.conf -v : print the state of the chain every frames.
2. chain.py config.conf -r : records the outputs only at the end of the simulation
3. chain.py config.conf -e : write the last event of the simulation for classification.

### Configuration file
The configuration is stored in a configuration file, ending by `.conf` extension. Parameters are sorted in categories starting with brackets []. Each entry is detailled below.
1. [sim] Generic category that applies for the simulation itself.
	* path : path where lumens.dat and bridges.dat are located. Useful for a pre-existing chain.
	* nlumens : int, number of lumens if path is not specified.
	* recording : boolean, if True, records all steps of the integration in .dat output files.
	* rec_distrib : boolean, if True, records the distribution of lengths and number of ions every frame.
	* chain_type : hydroosmotic or hydraulic. Type of the chain.
	* pumping : boolean, if True, allows pumping activity.
	* nb_frames : int, number of frames.
	* savefig : boolean, if True, save pictures of the chain every frame.
	* outdir : path, folder where output files are stored.
	* seed : int, if specified, it is the seed used for RNG.
2. [topology] Category for the topology of the chain : sizes, topological lengths, etc.
	* l_merge : float, length below which two lumens fuse. Usually 0.001
	* l_dis : float, length below which a lumen disappears from the chain. Usually 0.1 
	* avg_size : float, average area of the lumens.
	* std_size : float, standard deviation of the area of the lumens.
	* avg_dist : float, average lengt of the bridges.
	* std_dist : float, standard deviation of the length of the bridges.
	* dist_to_left : float, distance of the left-most lumen from the left-border
	* dist_to_right : float, distance of the right-most lumen from right-border
	* eps : float, Laplace/Osmotic pressure ratio. Usually, eps = 1e-3
	* merge : boolean, if True, allows for merging lumens.
	* pattern : uniform or normal. If normal, the chain is generated with normal distributions for the areas of lumens and lengths of the bridges, otherwise, they all have the same size (avg_size) and bridge lengths (avg_dist).
3. [integration] Category that accounts for the integration parameters.
	* solver : rkf45 (adaptative time-step). Another option (less stable) is rk45.
	* max_step : int, maximal number of time-steps allowed.
	* alpha : float. Arbitrary constant used to give the initial time-step. Usually, alpha = 1e-5
	* tolerance : float. Tolerance allowed to calculate the adaptative time-step. Usually, tolerance = 1e-10 to 1e-20.
4. [hydroosmotic] Category containing parameters for a hydro-osmotic chain.
	* equilibrium : boolean, if True, the chin starts at osmotic equilibrium, if False, the chain is generated with gaussian distribution
	* nions_avg : float, average number of ions (gaussian distribution)
	* nions_std : float, standard deviation of the number of ions (gaussian distribution)
	* chis : float, initial rescaled concentration screening length
	* chiv : float, initial rescaled pression screening length
	* taus : float, solute relaxation time
	* tauv : float, solvent relaxation time
5 [pumping] Category for the pumping parameters
	* pattern : None, constant, gaussian, sigmoide, rectangular. Profile of the active pumping along the chain.
	* param : float, values of the parameters of the pumping profile.

### Outputs
Output files are stored in the outdir folder specified in the configuration file. 
* Basic output files :
	* init_chain.dat : Initial chain informations.
	* end_chain.dat :  State of the chain at the end of the simulation.
	* sim_L_avg.dat : Average length of the lumens of the chain through time.
	* sim_ell_avg.dat : Average length of the bridges of the chain through time.
	* sim_nlum.dat : Number of lumens in time.
	* log.txt : simulation informations
	* events.txt : list of the integration anf topological events that happen during the simulation
* If rec_distrib = True :
	* distrib_length.dat : distribution of the lumen lengths at every frame
	* distrib_nion.dat : distribution of the number of ions in the lumens at every frame
* If recording = True :
	* sim_all.dat : lengths of all lumen at every time step
	* sim_bridges.dat : lengths of all bridges at every time step

## About
### Authors
Developped by Mathieu Le Verge--Serandour and Hervé Turlier, 2019
