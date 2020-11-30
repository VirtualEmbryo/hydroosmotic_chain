
# Requirements

# Content
* _ressources/ : contains the ressources used to simulate a chain of lumens
	* analysis_tools/
		* analysis.py : file used to extract exploitable data (concentrations, pressures, ...) of the chain ;
	* chain.py : main file, used to launch a simulation of a chain, given a .conf file ;
	* configreader.py : reads the configuration files (.conf) ;
	* flux.py : library used to perform flux calculations ;
	* gen_config.py : generates configurations files (.conf) from templates (.conf.tpl) ;
	* lumen_class.py : library that contains the definition of the objects chain, lumen and bridge ;
	* network.py : library of useful functions to work on graphs (nearest neighbors, etc) :
	* submit.py : sends simulations to a cluster using SLURM ;
	* tools.py : library of diverse useful functions
	* topology.py : library of useful functions to work on the topology of the graph (merge and disappear functions are defined there)

* _configfiles/ : contains some configfiles for example or frequently used

# Instructions


# About
## Authors
Developed by Mathieu Le Verge--Serandour and Hervé Turlier, 2020
