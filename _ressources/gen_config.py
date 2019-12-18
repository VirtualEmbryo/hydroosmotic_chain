#!/usr/bin/env python
# gen_config.py

"""
gen_config.py config.conf.tpl [options]

    Options
    -------
    nsim : int, optional, default : 1
        Number of simulations.

    outdir : float, optional default : out
        
"""

import os, sys

module_path = os.path.abspath(os.path.join('..', 'chain_lumen/'))

if module_path not in sys.path :
    sys.path.append(module_path)

try :
    import _ressources.configreader as configreader
except :
    import configreader as configreader

import numpy.random as rand

outdir = 'run'
datadir = '~/git/chain_lumen/_data'

def create_directory(dirname, zfill=4, nmin=0) :
    
    foldername = dirname + str(nmin).zfill(zfill)
    while os.path.isdir( dirname + str(nmin).zfill(zfill) ) :
        nmin += 1
        foldername = dirname + str(nmin).zfill(zfill) + '/'
    print(foldername)
    os.mkdir(foldername)
    return os.path.abspath(foldername)
    
def cp_config(config, foldername, seed=True) :
    
    config.add('sim', 'outdir', foldername)
    if seed :
        config.add('sim', 'seed', str(int(rand.random()*1e9)))
        #config['sim']['seed'] = str(int(rand.random()*1e9))
    
    #with open(os.path.join(foldername, 'config.conf'), 'w') as configfile:
    #    config.write(configfile)
    config.write(os.path.join(foldername, 'config.conf'))
    return ;
    

def main(conf_name, args) :
    global datadir
    # default args
    dirname = outdir
    nsim = 1
    seed = True
    
    # Read arglist
    for arg in args :
        if arg.startswith('nsim=') :
            nsim = int(arg[len('nsim='):])
        elif arg.startswith('outdir=') :
            dirname = arg[len('outdir='):]
        elif arg.startswith('seed=') :
            seed = eval(arg[len('seed='):])

    # Read config file        
    config = configreader.Config()
    config.read(conf_name)
    
    # Create directories
    for n in range(nsim) :
        foldername = create_directory(dirname)
        print(foldername)
        cp_config(config, foldername, seed=seed)
        
    return ;

if __name__ == '__main__' :
    if len(sys.argv) < 2:
        print('[network_simulation.py] Error: missing config file, type help for instructions.')
    elif sys.argv[1]=='help' or sys.argv[1]=='-h':
        print(__doc__)
    # first argument should be a readable config file:
    elif os.access(sys.argv[1], os.R_OK) and sys.argv[1].endswith('.conf.tpl'):
        conf_name = sys.argv[1]
        args = []
        try :
            args = sys.argv[2:]
        except :
            pass
            
        main(conf_name, args)
    else :
        print('Error : no readable config file found, must be a .conf.tpl file.')
        print(sys.argv[1])
    sys.exit()