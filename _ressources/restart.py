#!/usr/bin/env python
# restart.py

import numpy as np
import matplotlib.pyplot as plt

import os
import sys

import chain
import tools
import lumenclass as lc
import configreader

def search_value(filename, value_name) :
    my_file = open(filename).readlines()
    for line in my_file :
        if line.find(value_name) != -1 :
            value = line.split(' ')[-1]
    return value

def import_config(old_dirname) :
    

    config = configreader.Config()
    new_conf = config.read(os.path.join(old_dirname, 'config.conf'))
    end_file = open(os.path.join(old_dirname, 'end_chain.dat')).readlines()

    nb_lumens = int(search_value(os.path.join(old_dirname, 'end_chain.dat'), 'Nb lumens'))

    chain_type = new_conf['sim']['chain_type']
    e0 = float(new_conf['sim']['e0'])
    l_dis = float(new_conf['topology']['l_dis'])
    l_merge = float(new_conf['topology']['l_merge'])

    theta = np.pi/3.
    eps = float(new_conf['topology']['eps'])
    pbc = False


    if chain_type == 'hydraulic' :
        new_ch = lc.Chain(nb_lumens=nb_lumens, e0 = e0, theta=theta, l_dis=l_dis, l_merge=l_merge, pbc=pbc)
    elif chain_type == 'hydroosmotic' :
        new_ch = lc.Osmotic_Chain(nb_lumens=nb_lumens, e0 = e0, theta=theta, l_dis=l_dis, l_merge=l_merge, pbc=pbc)
        new_ch.xis = float(search_value(os.path.join(old_dirname, 'end_chain.dat'), 'xi_s'))
        new_ch.xiv = float(search_value(os.path.join(old_dirname, 'end_chain.dat'), 'xi_v'))
        new_ch.taus = float(search_value(os.path.join(old_dirname, 'end_chain.dat'), 'tau_s'))
        new_ch.tauv = float(search_value(os.path.join(old_dirname, 'end_chain.dat'), 'tau_v'))
        new_ch.pumping = eval(new_conf['sim']['pumping'])
        new_ch.merge = eval(new_conf['topology']['merge'])
    
    new_ch.lumens_dict, new_ch.bridges_dict = {}, {}
    new_ch.time = float(search_value(os.path.join(old_dirname, 'end_chain.dat'), 'Current Time'))
    new_ch.total_length = float(search_value(os.path.join(old_dirname, 'end_chain.dat'), 'Total length'))

    for line in end_file :
        if line.startswith('Lumen') :
            lum = line.split(' ')
            #print(lum)
            index = int(lum[1])
            init_pos = float(lum[5])
            init_length = float(lum[8])
            init_nb_ions = float(lum[10])
            ca = float(lum[14])
            new_ch.lumens_dict[index] = lc.Osmotic_Lumen(index, init_pos, init_length, init_nb_ions, theta, eps, ca)
        
    for line in end_file :
        if line.startswith('Bridge') :
            bridge = line.split(' ')
            index = int(bridge[1])
            lumen1 = float(bridge[3][1:-1])
            lumen2 = float(bridge[4][:-1])
            length = float(bridge[7])
            ca = float(bridge[10])
            new_ch.bridges_dict[index] = lc.Osmotic_Bridge(index, lumen1, lumen2, length, ca)
    
    return new_ch
    

def main(configname, args) :
    # Initialize arguments
    state_simul = False
    frame_reg = True
    new_dirname = 'bis'
    pics_dirname = os.path.join(new_dirname, 'pics')
    
    # Import arguments from the console (if any...)
    if len(args) > 0 :
        for arg in args :
            if arg.startswith('-v') :
                state_simul = True
            if arg.startswith('-r') :
                frame_reg = False
            if arg.startswith('dirname=') :
                new_dirname = arg[len('dirname='):]
                
    # Make dir bis
    try : 
        os.mkdir(new_dirname)
        os.mkdir(os.path.join(new_dirname, 'pics'))
    except : pass
    
    config, old_chain = chain.load_config(configname)
    
    # Parameters of the simulation
    max_step = int(config['integration']['max_step'])
    alpha = float(config['integration']['alpha'])
    recording = eval(config['sim']['recording'])
    tolerance = float(config['integration']['tolerance'])
    nb_frames = int(config['sim']['nb_frames'])
    solver = config['integration']['solver']
    chain_type = config['sim']['chain_type']
    # Try to see if savefig is valid for saving figures
    if chain.MATPLOTLIB_BOOL == True :
        savefig = eval(config['sim']['savefig'])
    else :
        savefig = False
    
    # Import the previous chain
    new_chain = import_config(old_dirname=os.getcwd())
    
    # Run Simulation
    end = chain.run(new_chain, max_step = max_step, alpha = alpha, recording=recording, tolerance=tolerance, nb_frames=nb_frames, solver=solver, savefig=savefig, state_simul=state_simul, dir_name=new_dirname,  pics_dirname=pics_dirname, frame_reg=frame_reg)
    # Add the ending event to events.log
    tools.write_ending_event(end, chain=new_chain, eventfilename='events.log')
    
    # Save the final chain in the end_chain.dat file
    new_chain.__save__(os.path.join(new_dirname, 'end_chain.dat'))

if __name__ == '__main__' :
    if len(sys.argv) < 2:
        print('[network_simulation.py] Error: missing config file, type help for instructions.')
    elif sys.argv[1]=='help' or sys.argv[1]=='-h':
        print(__doc__)
    # first argument should be a readable config file:
    elif os.access(sys.argv[1], os.R_OK):
        conf_name = sys.argv[1]
        args = []
        try :
            args = sys.argv[2:]
        except :
            pass
            
        main(conf_name, args)
    else :
        print('Error : no readable config file found')
        print(sys.argv[1])
    sys.exit()