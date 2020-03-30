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
    if len(s) == 22 :
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
    else : 
        xis = float(s[5].split(' ')[-1])
        xiv = float(s[6].split(' ')[-1])
    
        L1 = float(s[13].split(' ')[8])
        N1 = float(s[13].split(' ')[10])
        L2 = float(s[14].split(' ')[8])
        N2 = float(s[14].split(' ')[10])
    
        ell_12 = float(s[18].split(' ')[7])
        W = 0
    
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
    
    print(main_dir, tpl, filename, outfile)
    
    list_dir = get_listdir(main_dir, tpl)
    outputs = get_outputs(filename, list_dir, main_dir)
    
    make_file_outputs(outputs, os.path.join(main_dir, outfile))
    
    return ;

def distrib(line, nbins=10, Lmax=None) :
    """Calculate the distribution of a configuration given a time step in the shape of a line
    
    line = [time, L1, L2, L3, ...]
    
    
    """
    time = float(line.split('\t')[0])
    s = line.split('\t')[1:]
    values = []
    for elem in s :
        if elem != '' and elem != '\n' :
            if float(elem) > 0.1 :
                if Lmax != None :
                    if float(elem) <= Lmax :
                        values += [float(elem)]
                else :
                    values += [float(elem)]
    x, y = np.histogram(values, bins=nbins)
    y = 0.5*(y[1:]+y[:-1])
    return time, [x, y]

def batch_window(data, wmin, wmax, nwindow) :
    window = np.logspace(wmin, wmax, nwindow)
    time = np.cumsum(window)
    batch = []
    for i in range(len(time)) :
        indices = np.argwhere(np.abs(data[:, 0] - time[i]) <= window[i])[:, 0]
        batch += [data[indices]]
    return batch    

def batch_average(batchlist) :
    B_avg = []
    B_std = []

    for i in range(len(batchlist[0])) :
        Lavg = []
        Lstd = []
        for b in batchlist :
            Lavg += [np.average(b[i], axis=0)]
            Lstd += [np.std(b[i], axis=0)]
    
        tavg = np.nanmean([Lavg[j][0] for j in range(len(Lavg))])
        navg = np.nanmean([Lavg[j][1] for j in range(len(Lavg))])
        
        tstd = np.nanstd([Lstd[j][0] for j in range(len(Lstd))])
        nstd = np.nanstd([Lstd[j][1] for j in range(len(Lstd))])
        
        B_avg += [[tavg, navg]]
        B_std += [[tstd, nstd]]
        
    B_avg = np.array(B_avg)
    B_std = np.array(B_std)
    
    return B_avg, B_std

def batch(data_dict, wmin, wmax, nwindow) :
    window = np.logspace(wmin, wmax, nwindow)
    time = np.cumsum(window)
    dat_batch_list = []
    for k in data_dict.keys() :
        dat_batch_list += [batch_window(data_dict[k], wmin=wmin, wmax=wmax, nwindow=nwindow)]
        print(k, end='\r')
    print('End of import !')
    B_avg, B_std = batch_average(dat_batch_list)
    return B_avg, B_std

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
    
    