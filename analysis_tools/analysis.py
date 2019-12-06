#!/usr/bin/env python
# analysis.py

"""
python3 

    Options
    -------
    main_dir
    tpl
    filename
    outfile

"""

import numpy as np
try : import matplotlib.pyplot as plt
except : pass

import os, sys

def get_listdir(main_dir, tpl) :
    L = os.listdir(main_dir)
    
    listdir = []
    for elem in L :
        if elem.startswith(tpl) :
            listdir += [elem]
    return listdir

def get_initialconditions(filename, mu=0.6105653703843762, eps = 1e-3) :
    f = open(filename, 'r')
    s = f.readlines()
    f.close()
    
    xis = float(s[5].split(' ')[-1])
    xiv = float(s[6].split(' ')[-1])
    
    L1 = float(s[13].split(' ')[8])
    N1 = float(s[13].split(' ')[10])
    L2 = float(s[14].split(' ')[8])
    N2 = float(s[14].split(' ')[10])
    
    ell_12 = float(s[18].split(' ')[7])
    W = int(s[21].split(' ')[-1])
    
    L0 = L1+L2+ell_12
    P1 = L0*eps/L1
    P2 = L0*eps/L2
    
    C1 = N1 * mu / L1**2
    C2 = N2 * mu / L2**2
    
    return xis, xiv, P1, C1, P2, C2, W

def get_outputs(filename, list_dir, main_dir = '_data/') :
    L = []
    for elem in list_dir :
        L += [get_initialconditions(os.path.join(main_dir, elem, filename))]
    return np.array(L)

def make_file_outputs(outputs, filename, header='') :
    np.savetxt(filename, outputs, delimiter='\t', header=header)
        
def main(args) :
    for arg in args :
        if arg.startswith('main_dir=') :
            main_dir = arg[len('main_dir='):]
        elif arg.startswith('tpl=') :
            tpl = arg[len('tpl='):]
        elif arg.startswith('filename=') :
            filename = arg[len('filename='):]
        elif arg.startswith('outfile=') :
            outfile = arg[len('outfile='):]
            
    list_dir = get_listdir(main_dir, tpl)
    outputs = get_outputs(filename, list_dir, main_dir)
    
    make_file_outputs(outputs, os.path.join(main_dir, outfile))
    
    return ;

if __name__ == '__main__' :
    if len(sys.argv) < 2:
        print('[network_simulation.py] Error: missing args.')
    elif sys.argv[1]=='help' or sys.argv[1]=='-h':
        print(__doc__)
    # first argument should be a readable config file:
    else :
        args = sys.argv[1:]
        main(args)
    sys.exit()
    
    