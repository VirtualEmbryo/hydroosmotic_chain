#!/usr/bin/env python
# chain.py

"""
python3 chain.py config.conf [options]


    Options
    -------
    -h : help
    -v : print the steps
    -r : do not register at frames


    Contains
    --------

    Required
    --------
    

"""
import os, sys
import time

module_path = os.path.abspath(os.path.join('..', 'chain_lumen/'))

if module_path not in sys.path :
    sys.path.append(module_path)

import numpy as np

# homemade libraries
try :
    import _ressources.flux as flux
    import _ressources.tools as tools
    import _ressources.network as net
    import _ressources.topology as tplg
    import _ressources.lumenclass as lc
    import _ressources.configreader as configreader
except :
    import flux as flux
    import tools as tools
    import network as net
    import topology as tplg
    import lumenclass as lc
    import configreader

try : 
    import matplotlib.pyplot as plt
    MATPLOTLIB_BOOL = True
except : 
    MATPLOTLIB_BOOL = False
    pass        

pow_law = 3.

# Ending conditions
one_lumen_end = True
#one_lumen_end = False

# ========================================================================
# ========================= LOAD CONFIG ==================================
# ========================================================================

def load_config(filename) :
    global my_chain
    
    # Import config
    conf = configreader.Config()
    config = conf.read(filename)    
    path = config['sim']['path']
    
    lumen_type = config['sim']['chain_type']
        
    # Import pre-written configuration if specified
    if len(path) > 0 :
        print('Import config from ' + path)      
        lumens_array = np.loadtxt(os.path.join(path,'lumens.dat'), delimiter = '\t', comments = '#')
        bridges_array = np.loadtxt(os.path.join(path,'bridges.dat'), delimiter = '\t', comments = '#')
        
        if lumen_type == 'hydroosmotic' :
            my_chain = lc.Osmotic_Chain(nb_lumens = len(lumens_array)-2, e0=float(config['sim']['e0']), l_merge=float(config['topology']['l_merge']), l_dis=float(config['topology']['l_dis']))
            my_chain.__import_config__(lumens_array, bridges_array, eps = float(config['topology']['eps']))
            
            chis = float(config['hydroosmotic']['chis'])
            chiv = float(config['hydroosmotic']['chiv'])
            my_chain.taus = float(config['hydroosmotic']['taus'])
            my_chain.tauv = float(config['hydroosmotic']['tauv'])
        
        elif lumen_type == 'hydraulic' :
            my_chain = lc.Chain(nb_lumens = len(lumens_array)-2, e0=float(config['sim']['e0']), l_merge=float(config['topology']['l_merge']), l_dis=float(config['topology']['l_dis']))
            my_chain.__import_config__(lumens_array, bridges_array, eps = float(config['topology']['eps']), kappa=float(config['topology']['kappa']))
            
            my_chain.tau = float(config['hydraulic']['tau'])
            my_chain.kappa = float(config['hydraulic']['kappa'])
            my_chain.gamma = float(config['hydraulic']['gamma'])
            
        my_chain.pumping = config['pumping']['pattern']
        
    else :
        # Seed
        if conf.has_option('sim', 'seed') and len(config['sim']['seed']) > 0 :
            np.random.seed(int(config['sim']['seed']))
        
        # Chain type generation 
        if lumen_type == 'hydroosmotic' :
            my_chain = lc.Osmotic_Chain(
                            nb_lumens = int(config['sim']['nlumens']), 
                            e0=float(config['sim']['e0']), 
                            l_merge=float(config['topology']['l_merge']), 
                            l_dis=float(config['topology']['l_dis']))
            
            equilibrium = config['hydroosmotic']['equilibrium']
            pattern = config['topology']['pattern']
            nions_avg = float(config['hydroosmotic']['nions_avg'])
            nions_std = float(config['hydroosmotic']['nions_std'])
            
            chis = float(config['hydroosmotic']['chis'])
            chiv = float(config['hydroosmotic']['chiv'])
            my_chain.taus = float(config['hydroosmotic']['taus'])
            my_chain.tauv = float(config['hydroosmotic']['tauv'])
            
        elif lumen_type == 'hydraulic' :
            my_chain = lc.Chain(
                nb_lumens = int(config['sim']['nlumens']), 
                e0=float(config['sim']['e0']), 
                l_merge=float(config['topology']['l_merge']), 
                l_dis=float(config['topology']['l_dis']))
            
            my_chain.tau = float(config['hydraulic']['tau'])
            #my_chain.kappa = float(config['hydraulic']['kappa'])
            #my_chain.gamma = float(config['hydraulic']['gamma'])
            
        # Pumping
        if conf.has_option('pumping', 'pattern') : 
            my_chain.pumping = config['pumping']['pattern']
            if config['pumping']['pattern'] == 'normal' :
                
                br_avg, br_std = float(config['pumping']['bridge_1']), float(config['pumping']['bridge_2'])
                lum_avg, lum_std = float(config['pumping']['lumen_1']), float(config['pumping']['lumen_2'])
                
                ca_bridge_list = [np.random.normal(br_avg, br_std) for b in range(my_chain.nb_lumens+1)]
                ca_lumen_list = [np.random.normal(lum_avg, lum_std) for m in range(my_chain.nb_lumens+2)]
            
            elif config['pumping']['pattern'] == 'normal_abs' :
                
                br_avg, br_std = float(config['pumping']['bridge_1']), float(config['pumping']['bridge_2'])
                lum_avg, lum_std = float(config['pumping']['lumen_1']), float(config['pumping']['lumen_2'])
                
                ca_bridge_list = [np.abs(np.random.normal(br_avg, br_std)) for b in range(my_chain.nb_lumens+1)]
                ca_lumen_list = [np.abs(np.random.normal(lum_avg, lum_std)) for m in range(my_chain.nb_lumens+2)]
            
            elif config['pumping']['pattern'] == 'uniform' :
                br_low, br_high = float(config['pumping']['bridge_1']), float(config['pumping']['bridge_2'])
                lum_low, lum_high = float(config['pumping']['lumen_1']), float(config['pumping']['lumen_2'])
                
                ca_bridge_list = [np.random.uniform(br_low, br_high) for b in range(my_chain.nb_lumens+1)]
                ca_lumen_list = [np.random.uniform(lum_low, lum_high) for m in range(my_chain.nb_lumens+2)]
                
            elif config['pumping']['pattern'] == 'gradient' :
                br_1, br_2 = float(config['pumping']['bridge_1']), float(config['pumping']['bridge_2'])
                lum_1, lum_2 = float(config['pumping']['lumen_1']), float(config['pumping']['lumen_2'])
                
                ca_bridge_list = [b*(br_2 - br_1)/(my_chain.nb_lumens+1) + br_1 for b in range(my_chain.nb_lumens+1)]
                ca_lumen_list = [m*(lum_2 - lum_1)/(my_chain.nb_lumens+2) + lum_1 for m in range(my_chain.nb_lumens+2)]
            
            elif config['pumping']['pattern'] == 'constant' :
                br_1 = float(config['pumping']['bridge_1'])
                lum_1 = float(config['pumping']['lumen_1'])
                
                ca_bridge_list = [br_1 for b in range(my_chain.nb_lumens+1)]
                ca_lumen_list = [lum_1 for m in range(my_chain.nb_lumens+2)]
            
            else :
                my_chain.pumping = 'None'
                ca_lumen_list = [0. for i in range(my_chain.nb_lumens+2)]
                ca_bridge_list = [0. for i in range(my_chain.nb_lumens+1)]
                
        else :
            my_chain.pumping = 'None'
            ca_lumen_list = [0. for i in range(my_chain.nb_lumens+2)]
            ca_bridge_list = [0. for i in range(my_chain.nb_lumens+1)]
        
        # Generate the graph
        if lumen_type == 'hydroosmotic' :
            my_chain.__gen_network_lumen_object__(
                avg_size=float(config['topology']['avg_size']), std_size=float(config['topology']['std_size']), 
                avg_dist=float(config['topology']['avg_dist']), std_dist=float(config['topology']['std_dist']), 
                dist_toleft=float(config['topology']['dist_toleft']), dist_toright=float(config['topology']['dist_toright']), 
                eps = float(config['topology']['eps']), 
                equilibrium=equilibrium, 
                pattern=pattern, 
                nions_avg=nions_avg, nions_std=nions_std, 
                ca_lumen_list=ca_lumen_list, ca_bridge_list=ca_bridge_list)
        elif lumen_type == 'hydraulic' :
            my_chain.__gen_network_lumen_object__(
                avg_size=float(config['topology']['avg_size']), std_size=float(config['topology']['std_size']), 
                avg_dist=float(config['topology']['avg_dist']), std_dist=float(config['topology']['std_dist']), 
                dist_toleft=float(config['topology']['dist_toleft']), dist_toright=float(config['topology']['dist_toright']), 
                gamma = float(config['hydraulic']['gamma']), 
                ca_lumen_list=ca_lumen_list, ca_bridge_list=ca_bridge_list, 
                kappa = float(config['hydraulic']['kappa']))
    
    if lumen_type == 'hydroosmotic' :
        #my_chain.xis = chis*my_chain.bridges_dict[1].length
        #my_chain.xiv = chiv*my_chain.bridges_dict[1].length
        if 1 :
            average_length_bridges = np.average([my_chain.bridges_dict[k].length for k in my_chain.bridges_dict.keys()])
            my_chain.xis = chis*average_length_bridges
            my_chain.xiv = chiv*average_length_bridges
            
        my_chain.leaks = eval(config['hydroosmotic']['leaks'])
            
        #else :
        #    my_chain.xis = chis*float(config['topology']['avg_dist'])
        #    my_chain.xiv = chiv*float(config['topology']['avg_dist'])
        
        #my_chain.xis = chis*my_chain.total_length
        #my_chain.xiv = chiv*my_chain.total_length
    
    # Merging
    if conf.has_option('topology', 'merge') :
        my_chain.merge = eval(config['topology']['merge'])
        if not my_chain.merge : print('Merging not allowed')
    else :
        my_chain.merge = True
        
    #print(my_chain)
    return config, my_chain

# ========================================================================
# ======================== Runge-Kutta ===================================
# ========================================================================

def calc_new_timestep(h, error, tolerance, secure=0.8, cst_tstep=0) :

    if cst_tstep :
        return h
    else :
        if error == 0. :
            #print('Error !', error)
            return h*secure
        else :
            if error > tolerance :
                s = secure*(tolerance / error)**(0.2)
                #ratio = (tolerance / (2.*error))**0.5
                #s = secure*min(2., max(0.3, ratio))
            else :
                s = secure*(tolerance / error)**(0.25)
                #ratio = (tolerance / (2.*error))**0.5
                #s = secure*min(2., max(0.3, ratio))
        #print(s)
        return h*s

def calc_K_i(T, h, L, ell, chain, K_list=[], coeff_list=[]) :
    """
    T = t_k + a*h

    K_list = (K1, K2, ..., K_{X-1})
    """
    K = {}
    # corresponds to k1
    if len(K_list) == 0 and len(coeff_list) == 0 :
        for j in chain.lumens_dict.keys() :
            if j != 0 and j != -1 :
                K[j] = h*flux.func_Lj_hydraulic(j, T, L, ell, chain)
                #K[j] = flux.func_Lj_hydraulic(j, T, L, ell, chain)     ## test
            else :
                K[j] = 0.
    
    # corresponds to k2, k3, ...
    else :
        if len(coeff_list) != len(K_list) :
            print('ERROR ! Note the same lengths of coefficient or K vectors !')
            
        new_L, new_ell = {}, {}
        
        for j in chain.lumens_dict.keys() :
            if j != 0 and j != -1 :
                new_L[j] = L[j]
                for n in range(len(coeff_list)) :
                    new_L[j] += coeff_list[n]*K_list[n][j]
            else :
                new_L[j] = 0. 
                
        for b in chain.bridges_dict.keys() :
            x, y = chain.bridges_dict[b].lumen1, chain.bridges_dict[b].lumen2
            new_ell[b] = ell[b] + (L[x] - new_L[x]) + (L[y] - new_L[y])
        
        for j in chain.lumens_dict.keys() :
            if j != 0 and j != -1 :
                K[j] = h*flux.func_Lj_hydraulic(j, T, new_L, new_ell, chain)
                #K[j] = flux.func_Lj_hydraulic(j, T, new_L, new_ell, chain)     ##test
            else :
                K[j] = 0.
    return K

def calc_KQ_i(T, h, L, N, ell, chain, K_list=[], Q_list=[], coeff_list=[]) :
    """
    T = t_k + a*h

    K_list = (K1, K2, ..., K_{X-1})
    """
    K, Q = {}, {}
    if len(K_list) == 0 and len(Q_list) == 0 and len(coeff_list) == 0 :
        # corresponds to k1
        for j in chain.lumens_dict.keys() :
            if j != 0 and j != -1 :
                K[j] = h*flux.func_Lj(j, T, L, N, ell, chain)
                Q[j] = h*flux.func_Nj(j, T, L, N, ell, chain)
            else :
                K[j] = 0.
                Q[j] = 0.
    else :
        if len(coeff_list) != len(K_list) != len(Q_list) :
            print('ERROR ! Note the same lengths of coefficient or K, Q vectors !')
            
        # corresponds to k2, k3, ...
        new_L, new_N, new_ell = {}, {}, {}
        
        for j in chain.lumens_dict.keys() :
            if j != 0 and j != -1 :
                new_L[j] = L[j]
                new_N[j] = N[j]
                for n in range(len(coeff_list)) :
                    new_L[j] += coeff_list[n]*K_list[n][j]
                    new_N[j] += coeff_list[n]*Q_list[n][j]
            else :
                new_L[j] = 0. 
                new_N[j] = 0.
                
        for b in chain.bridges_dict.keys() :
            x, y = chain.bridges_dict[b].lumen1, chain.bridges_dict[b].lumen2
            new_ell[b] = ell[b] + (L[x] - new_L[x]) + (L[y] - new_L[y])
        
        for j in chain.lumens_dict.keys() :
            if j != 0 and j != -1 :    
                K[j] = h*flux.func_Lj(j, T, new_L, new_N, new_ell, chain)
                Q[j] = h*flux.func_Nj(j, T, new_L, new_N, new_ell, chain)
            else :
                K[j] = 0.
                Q[j] = 0.
    return K, Q

def test_DeltaK(L_vec, ell_vec, Ky, chain) :
    repeat = False
    index_L_list = []
    index_ell_list = []
    for j in L_vec.keys() :
        if L_vec[j]+Ky[j] < 0 :
            #print(j)
            repeat = True
            index_L_list += [j]
    for b in ell_vec.keys() :
        
        b_1, b_2 = chain.bridges_dict[b].lumen1, chain.bridges_dict[b].lumen2
        if b_1 != 0 and b_1 != -1 and b_2 != 0 and b_2 != -1 :
            if ell_vec[b] - Ky[b_1] - Ky[b_2] < 0 :
                repeat = True
                index_ell_list += [b]
    
    if len(index_L_list) == 0 and len(index_ell_list) == 0:
        index_L_list = None
        index_ell_list = None
        
    return repeat, index_L_list, index_ell_list
        
def rk45_step(t0, chain, h, alpha) :
    global pow_law
    cp_chain = chain.__copy__()
    
    L_vec, ell_vec = chain.__L_vec__(), chain.__ell_vec__()
    if chain.lumen_type == 'hydroosmotic' :
        N_vec = chain.__N_vec__()
        
    repeat = True
    while repeat :
        if chain.lumen_type == 'hydroosmotic' :
            ### RK4 - Step 1
            K1, Q1 = calc_KQ_i(t0, h, L_vec, N_vec, ell_vec, chain)

            ### RK4 - STEP 2
            K2, Q2 = calc_KQ_i(t0+0.5*h, h, L_vec, N_vec, ell_vec, chain, K_list=[K1], Q_list=[Q1], coeff_list=[0.5])
    
            ### RK4 - STEP 3
            K3, Q3 = calc_KQ_i(t0+0.5*h, h, L_vec, N_vec, ell_vec, chain, K_list=[K2], Q_list=[Q2], coeff_list=[0.5])
    
            ### RK4 - STEP 4
            K4, Q4 = calc_KQ_i(t0+h, h, L_vec, N_vec, ell_vec, chain, K_list=[K3], Q_list=[Q3], coeff_list=[1.])
    
    
            ### FINAL
            K = {j: (K1[j]+2.*K2[j]+2.*K3[j]+K4[j])/6. for j in chain.lumens_dict.keys()}
            Q = {j: (Q1[j]+2.*Q2[j]+2.*Q3[j]+Q4[j])/6. for j in chain.lumens_dict.keys()}
            
            ### Test if configuration is allowed
            repeat, index_L_list, index_ell_list = test_DeltaK(L_vec, ell_vec, K, chain)
        
        elif chain.lumen_type == 'hydraulic' :
            ### RK4 - Step 1

            K1 = calc_K_i(t0, h, L_vec, ell_vec, chain)

            ### RK4 - STEP 2
            K2 = calc_K_i(t0+0.5*h, h, L_vec, ell_vec, chain, K_list=[K1], coeff_list=[0.5])
    
            ### RK4 - STEP 3
            K3 = calc_K_i(t0+0.5*h, h, L_vec, ell_vec, chain, K_list=[K2], coeff_list=[0.5])
    
            ### RK4 - STEP 4
            K4 = calc_K_i(t0+h, h, L_vec, ell_vec, chain, K_list=[K3], coeff_list=[1.])
    
    
            ### FINAL
            K = {j: (K1[j]+2.*K2[j]+2.*K3[j]+K4[j])/6. for j in chain.lumens_dict.keys()}
            if chain.lumen_type == 'hydroosmotic' :
                Q = {j: (Q1[j]+2.*Q2[j]+2.*Q3[j]+Q4[j])/6. for j in chain.lumens_dict.keys()}
            
            ### Test if configuration is allowed
            repeat, index_L_list, index_ell_list = test_DeltaK(L_vec, ell_vec, K, chain)
            
        count = 0
        
        if repeat :
            count += 1
            h = 0.5*h
            chain.events += 'ERROR Time : ' + str(chain.time) + ' : length(s) of lumen(s) ' + str(index_L_list) + ' is negative. Time step is divided.\n'
            
            if count >= 10 : print('More than '+str(count)+' trials without convergence. There might be a problem.')
    
    
    # UPDATE CHAIN
    for j in chain.lumens_dict.keys() :
        chain.lumens_dict[j].length  += K[j]
        if chain.lumen_type == 'hydroosmotic' :
            chain.lumens_dict[j].nb_ions += Q[j]
    
    net.calc_ell_list(chain)
    new_time_step = alpha/(len(chain.lumens_dict.keys()) -2.)**(pow_law)
    
    return new_time_step

def rkf45_step(t0, chain, h, tolerance = 1e-6) :
    cp_chain = chain.__copy__()
    
    L_vec, ell_vec = chain.__L_vec__(), chain.__ell_vec__()
    
    if chain.lumen_type == 'hydroosmotic' :
        N_vec = chain.__N_vec__()
    
    repeat = True
    
    count = 0
    
    while repeat :
        if chain.lumen_type == 'hydroosmotic' :
            
            ### RK4 - Step 1
            K1, Q1 = calc_KQ_i(t0, h, L_vec, N_vec, ell_vec, chain)

            ### RK4 - STEP 2
            coeff_list = [0.25]
            K2, Q2 = calc_KQ_i(t0+0.25*h, h, L_vec, N_vec, ell_vec, chain, K_list=[K1], Q_list=[Q1], coeff_list=coeff_list)
    
            ### RK4 - STEP 3
            coeff_list = [3./32, 9./32]
            K3, Q3 = calc_KQ_i(t0+(3./8)*h, h, L_vec, N_vec, ell_vec, chain, K_list=[K1, K2], Q_list=[Q1, Q2], coeff_list=coeff_list)
    
            ### RK4 - STEP 4
            coeff_list = [1932./2197, -7200./2197, 7296./2197]
            K4, Q4 = calc_KQ_i(t0+(12./13)*h, h, L_vec, N_vec, ell_vec, chain, K_list=[K1, K2, K3], Q_list=[Q1, Q2, Q3], coeff_list=coeff_list)
    
            ### RK4 - STEP 5
            coeff_list = [439./216, -8., 3680./513, -845./4104]
            K5, Q5 = calc_KQ_i(t0+h, h, L_vec, N_vec, ell_vec, chain, K_list=[K1, K2, K3, K4], Q_list=[Q1, Q2, Q3, Q4], coeff_list=coeff_list)
    
            ### RK4 - STEP 6
            coeff_list = [-8./27, 2., -3544./2565, 1859./4104, -11./40]
            K6, Q6 = calc_KQ_i(t0+0.5*h, h, L_vec, N_vec, ell_vec, chain, K_list=[K1, K2, K3, K4, K5], Q_list=[Q1, Q2, Q3, Q4, Q5], coeff_list=coeff_list)
    
            ### FINAL
            Ky = {j: (25./216)*K1[j] + (1408./2565)*K3[j] + (2197./4104)*K4[j] -(1./5)*K5[j] for j in chain.lumens_dict.keys()}
            Qy = {j: (25./216)*Q1[j] + (1408./2565)*Q3[j] + (2197./4104)*Q4[j] -(1./5)*Q5[j] for j in chain.lumens_dict.keys()}
    
            Kz = {j: (16./135)*K1[j] + (6656./12825)*K3[j] + (28561./56430)*K4[j] - (9./50)*K5[j] + (2./55)*K6[j] for j in chain.lumens_dict.keys()}
            Qz = {j: (16./135)*Q1[j] + (6656./12825)*Q3[j] + (28561./56430)*Q4[j] - (9./50)*Q5[j] + (2./55)*Q6[j] for j in chain.lumens_dict.keys()}
            
            ### Test if configuration is allowed
            repeat, index_L_list, index_ell_list = test_DeltaK(L_vec, ell_vec, Kz, chain)
            
        elif chain.lumen_type == 'hydraulic' :
            ### RK4 - Step 1
            K1 = calc_K_i(t0, h, L_vec, ell_vec, chain)

            ### RK4 - STEP 2
            coeff_list = [0.25]
            K2 = calc_K_i(t0+0.25*h, h, L_vec, ell_vec, chain, K_list=[K1], coeff_list=coeff_list)
    
            ### RK4 - STEP 3
            coeff_list = [3./32, 9./32]
            K3 = calc_K_i(t0+(3./8)*h, h, L_vec, ell_vec, chain, K_list=[K1, K2], coeff_list=coeff_list)
    
            ### RK4 - STEP 4
            coeff_list = [1932./2197, -7200./2197, 7296./2197]
            K4 = calc_K_i(t0+(12./13)*h, h, L_vec, ell_vec, chain, K_list=[K1, K2, K3], coeff_list=coeff_list)
    
            ### RK4 - STEP 5
            coeff_list = [439./216, -8., 3680./513, -845./4104]
            K5 = calc_K_i(t0+h, h, L_vec, ell_vec, chain, K_list=[K1, K2, K3, K4], coeff_list=coeff_list)
    
            ### RK4 - STEP 6
            coeff_list = [-8./27, 2., -3544./2565, 1859./4104, -11./40]
            K6 = calc_K_i(t0+0.5*h, h, L_vec, ell_vec, chain, K_list=[K1, K2, K3, K4, K5], coeff_list=coeff_list)
    
            ### FINAL
            Ky = {j: (25./216)*K1[j] + (1408./2565)*K3[j] + (2197./4104)*K4[j] -(1./5)*K5[j] for j in chain.lumens_dict.keys()}
    
            Kz = {j: (16./135)*K1[j] + (6656./12825)*K3[j] + (28561./56430)*K4[j] - (9./50)*K5[j] + (2./55)*K6[j] for j in chain.lumens_dict.keys()}

            ### Test if configuration is allowed
            repeat, index_L_list, index_ell_list = test_DeltaK(L_vec, ell_vec, Kz, chain)
        
        if repeat :
            count += 1
            h = 0.5*h
            chain.events += 'ERROR Time : ' + str(chain.time) + ' : length(s) of lumen(s) ' + str(index_L_list) + ' or bridge(s) ' + str(index_ell_list) + ' is negative. Time step is divided.\n'
        
            if count >= 10 : print('More than '+str(count)+' trials without convergence. There might be a problem.')
    
    # UPDATE CHAIN
    for j in chain.lumens_dict.keys() :
        chain.lumens_dict[j].length  += Kz[j]
        if chain.lumen_type == 'hydroosmotic' :
            chain.lumens_dict[j].nb_ions  += Qz[j]
        
    net.calc_ell_list(chain)
    
    # Update Time Step
    error_list = []
    for j in Ky.keys() :
        error_list += [abs(Ky[j]-Kz[j])]
        
    error = max(error_list)
        
    new_time_step = calc_new_timestep(h, error, tolerance, secure=0.9, cst_tstep=False)
    
    return new_time_step
    
# ========================================================================
# ========================= SIMULATION ===================================
# ========================================================================

def system(chain, h=1e-2, recording = False, method = 'rk45', tolerance=1e-10, alpha=1e-2) :
    stop = False
    t0 = chain.time
    if method == 'rk45' :
        new_tstep = rk45_step(t0, chain, h, alpha=alpha)
    elif method == 'rkf45' :
        new_tstep = rkf45_step(t0, chain, h, tolerance=tolerance)
    else :
        print('Method not recognized')
    
    stop_cause = ''
    
    #if len(tplg.check_emptylumens(chain)) > 0 or len(tplg.check_merginglumens(chain)) > 0 or len(tplg.empty_ions(chain)) > 0:
    #    1
        
    if len(tplg.check_emptylumens(chain)) > 0 :
        stop = True
        stop_cause = 'Empty'
        
    elif len(tplg.check_merginglumens(chain)) > 0 :
        stop = True
        stop_cause = 'Empty'
        
    elif chain.lumen_type == 'hydroosmotic' and len(tplg.empty_ions(chain)) > 0 :
        stop = True
        stop_cause = 'Empty_ions'
        
    elif chain.nb_lumens == 1 and np.abs(2*chain.lumens_dict[max(chain.lumens_dict.keys())].length - chain.total_length) < chain.l_merge :
        stop = True
        stop_cause = 'full_size'
        print('FULL SIZE')
        
    elif chain.nb_lumens <= 1 :
        stop = True
        stop_cause = 'end_simul'
        
    if recording :
        chain.__record_state__()

    if method == 'rk45' :
        return stop, stop_cause, new_tstep 
    elif method == 'rkf45' :
        return stop, stop_cause, new_tstep

def make_Nfile(N0, filename, folder) :
    file_N = open(os.path.join(folder, filename), 'w')
    file_N.write('#t\tN(t)\n')
    file_N.write(str(0.)+'\t'+str(N0)+'\n')
    file_N.close()

def make_ellfile(ell_avg, filename, folder) :
    file_ell = open(os.path.join(folder, filename), 'w')
    file_ell.write('#t\tell(t)\n')
    
    file_ell.write(str(0.)+'\t'+str(ell_avg)+'\n')
    file_ell.close()
    
def make_Lfile(L_avg, filename, folder) :
    file_L = open(os.path.join(folder, filename), 'w')
    file_L.write('#t\tL(t)\n')
    
    file_L.write(str(0.)+'\t'+str(L_avg)+'\n')
    file_L.close()

def save_distribfile(chain, filename_l, filename_nion, folder) :
    file_length = open(os.path.join(folder, filename_l), 'a')
    s_length = str(chain.time) + '\t'
    
    if chain.lumen_type == 'hydroosmotic' :
        file_nion = open(os.path.join(folder, filename_nion), 'a')
        s_nion = str(chain.time)+ '\t'
    
    for k in chain.lumens_dict.keys() :
        if k != 0 and k != -1 :
            s_length += str(chain.lumens_dict[k].length) + '\t'
            if chain.lumen_type == 'hydroosmotic' : 
                s_nion += str(chain.lumens_dict[k].nb_ions) + '\t'
    s_length += '\n'
    if chain.lumen_type == 'hydroosmotic' : 
        s_nion += '\n'
    
    file_length.write(s_length)
    file_length.close()
    
    if chain.lumen_type == 'hydroosmotic' : 
        file_nion.write(s_nion)
        file_nion.close()

def save_N(t, Nt, filename) :
    file_N = open(filename, 'a')
    file_N.write(str(t)+'\t'+str(Nt)+'\n')
    file_N.close()
    
def save_ell(t, ellt, filename) :
    file_ell = open(filename, 'a')
    file_ell.write(str(t)+'\t'+str(ellt)+'\n')
    file_ell.close()
    
def save_L(t, Lt, filename) :
    file_L = open(filename, 'a')
    file_L.write(str(t)+'\t'+str(Lt)+'\n')
    file_L.close()
    
# ========================================================================
# ===========================  RUN  ======================================
# ========================================================================

def run(chain, max_step=1000, alpha=1e-4, savefig=False, nb_frames=1000, dir_name = 'out', recording=True, tolerance=1e-9, solver='rkf45', state_simul=False, pics_dirname='pics', frame_reg=True) :
    """
    
    Returns
    -------
    end : int
        0 : a problem occured during simulation
    
        10 : max step reached
        11 : full size reached
    """
    # Default parameters
    stop = False            # to stop simulation
    step = 0                # simulation step number
    N0 = chain.nmax         # number of lumens
    h = alpha / N0          # integration constant
    x = np.linspace(0, chain.total_length, 1001)     # x-axis positions

    make_Nfile(N0=N0, filename='sim_nlum.dat', folder = dir_name)
    make_ellfile(ell_avg=chain.__calc_ell_avg__(), filename='sim_ell_avg.dat', folder = dir_name)
    make_Lfile(L_avg=chain.__calc_L_avg__(), filename='sim_L_avg.dat', folder = dir_name)
    
    if recording :
        save_distribfile(chain, filename_l = 'distrib_length.dat', filename_nion = 'distrib_nion.dat', folder = dir_name)
    
    #if 1 :
    try :
        for i in range(max_step) :
            step += 1
            # make a step

            stop, stop_cause, h = system(chain, h=h, recording = recording, method = solver, tolerance=1e-10, alpha=alpha)
            chain.time += h
            
            # check topology
            if stop_cause != 'full_size' :
                tplg.topology(chain)
    
            # Save the number of lumens
            save_N(chain.time, len(chain.lumens_dict)-2, os.path.join(dir_name, 'sim_nlum.dat'))
            
            # Save bridges average length
            ellt_avg = chain.__calc_ell_avg__()
            save_ell(chain.time, ellt_avg, os.path.join(dir_name, 'sim_ell_avg.dat'))
            
            # Save lumens average length
            Lt_avg = chain.__calc_L_avg__()
            save_L(chain.time, Lt_avg, os.path.join(dir_name, 'sim_L_avg.dat'))
            
            # Save distributions
            if recording :
                save_distribfile(chain, filename_l = 'distrib_length.dat', filename_nion = 'distrib_nion.dat', folder = dir_name)
            
            if stop == 1 :
                if stop_cause == 'Empty' or stop_cause == 'Fusion' :
                    tplg.topology(chain)
                    stop = 0
                    
                elif stop_cause == 'end_simul':
                    
                    # One lumen left : return 2 ; otherwise return 1
                    if len(chain.lumens_dict) - 2 == 1 and one_lumen_end :
                        end = 2
                        print('End simulation : 1 Lumen left')
                        break ;
                        
                        #chain.events += 'Time ' + "{:4.6f}".format(chain.time) + ' : winner (' + end_event + ') is lumen ' + str(int(chain.winner))
                    
                    elif len(chain.lumens_dict) - 2 == 0 :
                        end = 1
                        #chain.events += 'Time ' + "{:4.6f}".format(chain.time) + ' : winner (' + 'D' + ') is lumen ' + str(0)
                        print('End simulation : 0 Lumen left')
                        break ;
                
                elif stop_cause == 'full_size' :
                    print('End simulation : system is fully occupied')
                    end = 11
                    break ;

                
            if savefig == True and step % nb_frames == 0 :
                tools.plot_profile(x, chain, centers=False, lw=1.5, show=False, savefig=True, savename=os.path.join(pics_dirname, 'pic'+str(step).zfill(8)+'.png'))
    
            if step % nb_frames == 0 :
                if state_simul :
                    print('Step : ', step, ' ; Time : ', "{:4.6f}".format(chain.time), ' ; Nb Lumens : ', len(chain.lumens_dict)-2)
                
                # Records the details of the chain
                save_data = recording
                if save_data and frame_reg :
                    #print('save', chain.lumen_type)
                    #tools.save_recording(chain, filename='sim_all.dat', filename_events='events.log', folder=dir_name, chain_type = chain.lumen_type, erase=True)
                    
            if step == max_step :
                print('\n\nEnd simulation : max step reached')
                tools.save_recording(chain, filename='sim_all.dat', filename_events='events.log', folder=dir_name)
        
        
        # Last step
        if i == max_step-1 :    # reached maximum step
            end = 10
            
        if savefig :
            tools.plot_profile(x, chain, centers=False, lw=1.5, show=0, savefig=True, savename=os.path.join(pics_dirname, 'pic'+str(step).zfill(8)+'.png'))        

        save_data = recording
        if save_data :
            tools.save_recording(chain, filename='sim_all.dat', filename_events='events.log', folder=dir_name, chain_type = chain.lumen_type, erase=False)
    #elif 0 :
    except RuntimeWarning :
        # RuntimeWarning is raised only in the flux.py library
        # If a RuntimeWarning exception occurs, it is raised as an exception.
        tools.save_recording(chain, filename='sim_all.dat', filename_events='events.log', folder=dir_name, erase=False)
        print('\nSimulation stopped before the end, flow error.')
        end = 0
    #else :
    except :
        tools.save_recording(chain, filename='sim_all.dat', filename_events='events.log', folder=dir_name, erase=False)
        print('\nSimulation stopped before the end...')
        end = 0
    
    return end ;

# ========================================================================
# ===========================  MAIN  =====================================
# ========================================================================
        
def main(configname, args) :
    # Initialize arguments
    state_simul = False
    frame_reg = True
    
    # Import arguments from the console (if any...)
    if len(args) > 0 :
        for arg in args :
            if arg.startswith('-v') :
                state_simul = True
            if arg.startswith('-r') :
                frame_reg = False
    
    # Load Configuration file
    config, my_chain = load_config(configname)
    
    # Parameters of simulation
    max_step = int(config['integration']['max_step'])
    alpha = float(config['integration']['alpha'])
    recording = eval(config['sim']['recording'])
    tolerance = float(config['integration']['tolerance'])
    nb_frames = int(config['sim']['nb_frames'])
    solver = config['integration']['solver']
    chain_type = config['sim']['chain_type']
    
    
    dir_name = config['sim']['outdir']
    
    # Clean dir_name if not empty
    if len(os.listdir(dir_name)) > 0 :
        for elem in os.listdir(dir_name) :
            if not elem.endswith('.conf') and not os.path.isdir(elem) :
                os.remove(os.path.join(dir_name, elem))
    
    # Try to see if savefig is valid for saving figures
    if MATPLOTLIB_BOOL == True :
        savefig = eval(config['sim']['savefig'])
    else :
        savefig = False
    
    # Save initial state of the chain (for reconstruction)
    my_chain.__save__(os.path.join(dir_name, 'init_chain.dat'))
    
    # If needed, create a folder for pictures filees
    pics_dirname=''
    if savefig :
        pics_dirname = os.path.join(dir_name,'pics')
        
        if os.path.isdir(pics_dirname) :
            if os.listdir(pics_dirname) > 0 :
                print('Remove files from ' + pics_dirname)
                for elem in os.listdir(pics_dirname) :
                    os.remove(os.path.join(pics_dirname, elem))

        else :
            os.mkdir(pics_dirname)
        print(pics_dirname)
        x = np.linspace(0, my_chain.total_length, 1001)
        tools.plot_profile(x, my_chain, centers=False, lw=1.5, show=False, savefig=True, savename=os.path.join(pics_dirname, 'pic'+str(0).zfill(8)+'.png'))
    
    # Run Simulation
    end = run(my_chain, max_step = max_step, alpha = alpha, recording=recording, tolerance=tolerance, nb_frames=nb_frames, solver=solver, savefig=savefig, state_simul=state_simul, dir_name=dir_name,  pics_dirname=pics_dirname, frame_reg=frame_reg)
    # Add the ending event to events.log
    tools.write_ending_event(end, chain=my_chain, eventfilename='events.log')
    
    # Save the final chain in the end_chain.dat file
    my_chain.__save__(os.path.join(dir_name, 'end_chain.dat'))
    return ;
    
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
    
    
#