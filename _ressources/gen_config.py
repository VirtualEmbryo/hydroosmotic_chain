#!/usr/bin/env python3
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
import numpy as np

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

keys_order = []
keys_order += [['sim', ['path', 'nlumens', 'e0', 'recording', 'chain_type', 'pumping', 'nb_frames', 'savefig', 'outdir', 'seed']]]
keys_order += [['topology', ['l_merge', 'l_dis', 'avg_size', 'std_size', 'avg_dist', 'std_dist', 'dist_toleft', 'dist_toright', 'eps', 'merge', 'pattern']]]
keys_order += [['integration', ['solver', 'max_step', 'alpha', 'tolerance']]]
keys_order += [['hydroosmotic', ['equilibrium', 'nions_avg', 'nions_std', 'chis', 'chiv', 'taus', 'tauv', 'leaks']]]
#keys_order += [['hydraulic', ['kappa', 'tau', 'gamma']],]
keys_order += [['hydraulic', ['kappa']],]
keys_order += [['pumping', ['pattern', 'param_1', 'param_2', 'param_3', 'param_4']]]


def create_directory(dirname, zfill=4, nmin=0) :
    
    foldername = dirname + str(nmin).zfill(zfill)
    while os.path.isdir( dirname + str(nmin).zfill(zfill) ) :
        nmin += 1
        foldername = dirname + str(nmin).zfill(zfill) + '/'
    os.mkdir(foldername)
    return os.path.abspath(foldername)
    
def cp_config(config, foldername, seed=True, config_path='') :
    global keys_order
    config.add('sim', 'outdir', foldername)
    if seed :
        config.add('sim', 'seed', str(int(rand.random()*1e9)))
        #config['sim']['seed'] = str(int(rand.random()*1e9))
    
    if len(config_path) > 0 :
        config.set_item('sim', 'path', config_path)
    config.write(os.path.join(foldername, 'config.conf'), keys_order=keys_order)
    return ;
    
def write_config_files(config_path, foldername, params=None) :
    os.mkdir(foldername)
    
    lumens = np.loadtxt(os.path.join(config_path, 'lumens.dat'))
    if len(lumens[0]) == 4 :
        header_l = 'index\tpos\tlength\tca'
        header_b = 'index\tlum1\tlum2\tlength\tca'
    elif len(lumens[0]) == 5 :
        header_l = 'index\tpos\tlength\tnb_ions\tca'
        header_b = 'index\tlum1\tlum2\tlength\tca'
        
    bridges = np.loadtxt(os.path.join(config_path, 'bridges.dat'))
    
    if params != None :
        key, value = params
    
        for k in range(len(key)) :
            if key[k].startswith('l') :
                index, attribute = key[k][1:].split(':')        # Split the key at the ':' excluding the object type (l or b)
                v = value[k]
                i_lum = np.argwhere(lumens[:, 0] == int(index))[0, 0]

                if attribute == 'pos' :
                    lumens[i_lum, 1] = v
                elif attribute == 'length' :
                    lumens[i_lum, 3] = v
                elif attribute == 'nb_ions' :
                    lumens[i_lum, 3] = v
                elif attribute == 'ca' :
                    lumens[i_lum, -1] = v

            
            elif key[k].startswith('b') :
                index, attribute = key[k][1:].split(':')        # Split the key at the ':' excluding the object type (l or b)
                v = value[k]
                i_br = np.argwhere(bridges[:, 0] == int(index))[0, 0]
                
                if attribute == 'length' :
                    bridges[i_br, 3] = v
                elif attribute == 'ca' :
                    bridges[i_br, -1] = v
        
        np.savetxt(os.path.join(foldername, 'lumens.dat'), lumens, header=header_l, fmt = '%.5f', delimiter='\t')
        np.savetxt(os.path.join(foldername, 'bridges.dat'), bridges, header=header_b, fmt = '%.5f', delimiter='\t')
        
    else :
        fw_lumens = open(os.path.join(foldername, 'lumens.dat'), 'w')
        fw_lumens.write(lumens)
    
        fw_bridges = open(os.path.join(foldername, 'bridges.dat'), 'w')
        fw_bridges.write(bridges)
    
        fw_lumens.close()
        fw_bridges.close()
    return ;

def main(conf_name, args) :
    global datadir
    # default args
    dirname = outdir
    nsim = 1
    seed = True
    params = None
    
    # Read arglist
    for arg in args :
        if arg.startswith('nsim=') :
            nsim = int(arg[len('nsim='):])
        elif arg.startswith('outdir=') :
            dirname = arg[len('outdir='):]
        elif arg.startswith('seed=') :
            seed = eval(arg[len('seed='):])
        elif arg.startswith('params=') :
            params = arg[len('params='):]
    
    # Read config file
    config = configreader.Config()
    config.read(conf_name)
    config_path = config.get_item('sim', 'path')
    
    # Create directories
    if params == None :
        for n in range(nsim) :
            foldername = create_directory(dirname)
            if len(config_path) > 0 :
                sim_config_path = os.path.join(foldername, 'config')
                cp_config(config, foldername, seed=seed, config_path=sim_config_path)
                write_config_files(config_path, sim_config_path, params='')
            else :
                cp_config(config, foldername, seed=seed)
    else :
        param_array = np.loadtxt(params)
        header = str.split(open(params, 'r').readlines()[0])[1:]
        
        for n in range(len(param_array)) :
            foldername = create_directory(dirname)
            if len(config_path) > 0 :
                sim_config_path = os.path.join(foldername, 'config')
                cp_config(config, foldername, seed=seed, config_path=sim_config_path)
                write_config_files(config_path, sim_config_path, params=(header, param_array[n]))
            
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