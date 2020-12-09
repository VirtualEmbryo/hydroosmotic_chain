#!/usr/bin/env python3
# chain.py

"""
python3 chain.py config.conf [options]

    Main script. Import a configuration file, integrate it and save output.
    

    Options
    -------
    -h : help
    -v : print the steps (verbose)
    -r : do not register at frames
    -e : register the last event (useful for fate diagram)

    Example
    -------
    1. Run a simulation of a chain with configuration file config.conf
        chain.py config.conf
    2. Run a simulation and print some information on the console
         chain.py config.conf -v

    Contains
    --------
        Load config
    load_config         : Load a configuration from a file and creates the chain.
    set_pumping         : Sets the active pumping coefficients for the lumens and bridges.
    calc_pumping        : Calculate the active pumping for the positions of the chain given a profile.
    gen_pumping_list    : Generates the active pumping list of coefficients.
        Runge-Kutta
    calc_new_timestep   : Calculate a new time-step given the rule of Press & Teukolsky, Comput. Phys., 1992.
    calc_K_i            : Calculate the coefficient K_i for RK scheme.
    calc_KQ_i           : Calculate the coefficient K_i and Q_i for RK scheme.
    test_DeltaK         : Test whether the new configuration is allowed.
    rk45_step           : Calculate an integration step using RK45 numerical scheme.
    rkf45_step          : Calculate an integration step using RKF45 numerical scheme.
        Simulation
    system              : Integrates the chain for one step. Checks whether the simulation must stop or not.
    make_Nfile          : Create a file for the number of lumens versus time.
    make_ellfile        : Create a file for the average length of bridges versus time.
    make_Lfile          : Create a file for the average length of lumens versus time.
    save_distribfile    : Save the distribution of lengths and number of ions of the chain.
    save_N              : Save the number of lumens of the chain.
    save_ell            : Save the average length of bridges of the chain.
    save_L              : Save the average and mean-field length of lumens of the chain.
        Run
    run                 : Run the full numerical integration of the chain.
        Main
    main                : Load the configuration, create a chain and integrate it.
        
        
    Required
    --------
        Python libraries
    os
    sys
    time
    numpy (np)
        Homemade libraries
    flux
    tools
    network (net)
    topology (tplg)
    lumenclass (lc)
    configreader
    functions

Mathieu Le Verge--Serandour, 2020
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
    import _ressources.functions as functions
except :
    import flux as flux
    import tools as tools
    import network as net
    import topology as tplg
    import lumenclass as lc
    import configreader
    import functions

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
    """
    load_config(filename)
    
        Load the configuration file and make a chain from it.
    
        Parameters
        ----------
        filename : string
            Name of the configuration file.
        Returns
        -------
        config : dictionnary, configuration
            Configuration stored in a dictionnary. For informations about configuration, please read configreader.py
            
        my_chain : chain_object, global
            The chain
    """
    global my_chain
    
    # Import config
    conf = configreader.Config()
    config = conf.read(filename)    
    path = config['sim']['path']
    
    lumen_type = config['sim']['chain_type']
    
    # ======================== IMPORT CONFIG FROM PATH ===============
    # Import pre-written configuration if specified
    if len(path) > 0 :
        print('Import config from ' + path)      
        lumens_array = np.loadtxt(os.path.join(path,'lumens.dat'), delimiter = '\t', comments = '#')
        bridges_array = np.loadtxt(os.path.join(path,'bridges.dat'), delimiter = '\t', comments = '#')
        
        if lumen_type == 'hydroosmotic' :
            my_chain = lc.Osmotic_Chain(nb_lumens = len(lumens_array)-2, l_merge=float(config['topology']['l_merge']), l_dis=float(config['topology']['l_dis']))
            my_chain.__import_config__(lumens_array, bridges_array, eps = float(config['topology']['eps']))
            
            chis = float(config['hydroosmotic']['chis'])
            chiv = float(config['hydroosmotic']['chiv'])
            my_chain.taus = float(config['hydroosmotic']['taus'])
            my_chain.tauv = float(config['hydroosmotic']['tauv'])
        
        elif lumen_type == 'hydraulic' :
            my_chain = lc.Chain(nb_lumens = len(lumens_array)-2,  l_merge=float(config['topology']['l_merge']), l_dis=float(config['topology']['l_dis']))
            my_chain.__import_config__(lumens_array, bridges_array, eps = float(config['topology']['eps']))
            
            
        my_chain.pumping = config['pumping']['pattern']
        
    # ======================== IMPORT PARAMETERS =====================   
    else :
        # Seed
        if conf.has_option('sim', 'seed') and len(config['sim']['seed']) > 0 :
            np.random.seed(int(config['sim']['seed']))
        
        # Chain type generation 
        if lumen_type == 'hydroosmotic' :
            my_chain = lc.Osmotic_Chain(
                            nb_lumens = int(config['sim']['nlumens']), 
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
                l_merge=float(config['topology']['l_merge']), 
                l_dis=float(config['topology']['l_dis']))
        
        # GENERATE THE GRAPH
        if lumen_type == 'hydroosmotic' :
            my_chain.__gen_network_lumen_object__(
                avg_size=float(config['topology']['avg_size']), std_size=float(config['topology']['std_size']), 
                avg_dist=float(config['topology']['avg_dist']), std_dist=float(config['topology']['std_dist']), 
                dist_toleft=float(config['topology']['dist_toleft']), dist_toright=float(config['topology']['dist_toright']), 
                eps = float(config['topology']['eps']), 
                equilibrium=equilibrium, 
                pattern=pattern, 
                nions_avg=nions_avg, nions_std=nions_std
                )
                
        elif lumen_type == 'hydraulic' :
            my_chain.__gen_network_lumen_object__(
                avg_size=float(config['topology']['avg_size']), std_size=float(config['topology']['std_size']), 
                avg_dist=float(config['topology']['avg_dist']), std_dist=float(config['topology']['std_dist']), 
                dist_toleft=float(config['topology']['dist_toleft']), dist_toright=float(config['topology']['dist_toright'])
                )
                
        # PUMPING
        if conf.has_option('pumping', 'pattern') : 
            my_chain.pumping = config['pumping']['pattern']
            
            gen_pumping_list(my_chain, config, my_chain.nb_lumens)            
            
        else :
            my_chain.pumping = 'None'

    if lumen_type == 'hydroosmotic' :
        average_length_bridges = my_chain.__calc_ell_avg__()
        average_L = my_chain.__calc_L_avg__()
        
        my_chain.xis = chis*average_length_bridges
        my_chain.xiv = chiv*average_length_bridges
            
        my_chain.leaks = False      # Allowds for leaks from the borders. Not implemented properly.
    
    # Merging
    if conf.has_option('topology', 'merge') :
        my_chain.merge = eval(config['topology']['merge'])
        if not my_chain.merge : print('Merging not allowed')
    else :
        my_chain.merge = True
        
    return config, my_chain
    
def set_pumping(chain, ca_lumen_list, ca_bridge_list) :
    """
    set_pumping(chain, ca_lumen_list, ca_bridge_list)
    
        Set the active pumping values for the lumens and the bridges
    
        Parameters
        ----------
        chain : chain object
        ca_lumen_list : list
            List of the active pumping for the lumens. Must be of length N, where N is the total number of lumens (borders excluded).
        ca_bridge_list : list
            List of the active pumping for the bridges. Must be of length N+1, where N is the total number of lumens (borders excluded).
    
        NB : Modify directly the chain.
    """
    for b in chain.bridges_dict.keys() :
            chain.bridges_dict[b].ca = ca_bridge_list[b]
        
    for k in chain.lumens_dict.keys() :
        if k != 0 and k != -1 :
            chain.lumens_dict[k].ca = ca_lumen_list[k]

def calc_pumping(chain, func, args):
    """
    calc_pumping(chain, func, args)
    
        Calculate the pumping profile given the specified positions and input profile.
    
        Parameters
        ----------
        chain : chain-object
        func : string
            Function that specifies the profile. Can be choosen between
            'None', 'constant', 'linear', 'gaussian', 'rectangular', 'sigmoide' 
            If None, no profile is calculated.
        args : dict
            Dictionnary containing the arguments of the profile function. The keys and values depend on the input function.
        
    """
    for k in chain.lumens_dict.keys() :
        if k != 0 and k != -1 :
            pos, length = chain.lumens_dict[k].pos, chain.lumens_dict[k].length
            x1, x2 = pos - length, pos + length
            chain.lumens_dict[k].ca = functions.integrate(func, x1, x2, args)
            
    for b in chain.bridges_dict.keys() :
        i, j = chain.bridges_dict[b].lumen1, chain.bridges_dict[b].lumen2
        x1, x2 = chain.lumens_dict[i].pos, chain.lumens_dict[j].pos
        chain.bridges_dict[b].ca = functions.integrate(func, x1, x2, args)

def gen_pumping_list(chain, config, nb_lumens) :
    '''
    gen_pumping_list(chain, config, nb_lumens)
        
        Generates a pumping list for lumens and bridges given the configuration
        
        Parameters
        ----------
        chain : chain-object
        config : dictionnary (configuration)
            Configuration of the chain
        nb_lumens : int
            Number of lumens of the chain.
    '''
    Ltot = chain.total_length
    func = config['pumping']['pattern']
    chain.pumping_func = func
    
    if func == 'constant' :
        chain.pumping_args = {'value' : float(config['pumping']['param_1'])}
        calc_pumping(chain, func, chain.pumping_args)
    # ========== RANDOM ============
    elif func == 'normal' :
        lum_avg, lum_std = float(config['pumping']['param_1']), float(config['pumping']['param_2'])
        br_avg, br_std = float(config['pumping']['param_3']), float(config['pumping']['param_4'])
        
        ca_bridge_list = [np.random.normal(br_avg, br_std) for b in range(nb_lumens+1)]
        ca_lumen_list = [np.random.normal(lum_avg, lum_std) for m in range(nb_lumens+2)]
        
        set_pumping(chain, ca_lumen_list, ca_bridge_list)
    
    elif func == 'normal_abs' :
        lum_avg, lum_std = float(config['pumping']['param_1']), float(config['pumping']['param_2'])
        br_avg, br_std = float(config['pumping']['param_3']), float(config['pumping']['param_4'])
        
        ca_bridge_list = [np.abs(np.random.normal(br_avg, br_std)) for b in range(nb_lumens+1)]
        ca_lumen_list = [np.abs(np.random.normal(lum_avg, lum_std)) for m in range(nb_lumens+2)]
        
        set_pumping(chain, ca_lumen_list, ca_bridge_list)
    
    elif func == 'uniform' :
        lum_low, lum_high = float(config['pumping']['param_1']), float(config['pumping']['param_2'])
        br_low, br_high = float(config['pumping']['param_3']), float(config['pumping']['param_4'])
        
        
        ca_bridge_list = [np.random.uniform(br_low, br_high) for b in range(nb_lumens+1)]
        ca_lumen_list = [np.random.uniform(lum_low, lum_high) for m in range(nb_lumens+2)]

        set_pumping(chain, ca_lumen_list, ca_bridge_list)
        
    # ========== FUNCTIONS ============
    elif func == 'linear' :
        chain.pumping_args = {'slope' : float(config['pumping']['param_1']), 'offset' : float(config['pumping']['param_2'])}
        calc_pumping(chain, func, chain.pumping_args)
            
    elif func == 'gaussian' :
        chain.pumping_args = {'amp' : float(config['pumping']['param_1']), 'mu' : float(config['pumping']['param_2'])*Ltot, 'sigma' : float(config['pumping']['param_3'])*Ltot, 'fmin' : float(config['pumping']['param_4'])}
        calc_pumping(chain, func, chain.pumping_args)

    elif func == 'rectangular' :
        chain.pumping_args = {'fmin' : float(config['pumping']['param_1']), 'fmax' : float(config['pumping']['param_2']), 'start' : float(config['pumping']['param_3'])*Ltot, 'stop' : float(config['pumping']['param_4'])*Ltot}
        calc_pumping(chain, func, chain.pumping_args)
        
    elif func == 'sigmoide' :
        chain.pumping_args = {'fmin' : float(config['pumping']['param_1']), 'fmax' : float(config['pumping']['param_2']), 'slope' : float(config['pumping']['param_3'])/Ltot, 'inflexion_point' : float(config['pumping']['param_4'])*Ltot}
        calc_pumping(chain, func, chain.pumping_args)
    
    else :
        my_chain.pumping = 'None'
        ca_lumen_list = [0. for i in range(nb_lumens+2)]
        ca_bridge_list = [0. for i in range(nb_lumens+1)]
        set_pumping(chain, ca_lumen_list, ca_bridge_list)
        
    return;
    
# ========================================================================
# ======================== Runge-Kutta ===================================
# ========================================================================

def calc_new_timestep(h, error, tolerance, secure=0.9, cst_tstep=0) :
    """
    calc_new_timestep(h, error, tolerance, secure=0.9, cst_tstep=0)
    
        Calculate the new time step from the previous one (h), given an input error and tolerance.
        Criterion is based on Press & Teukolsky, Computers in Physics, 1992.
    
        Parameters
        ----------
        h : float
            Initial time step.
        error : float
            Calculated error from the two estimations of RKF45 integration.
        tolerance : float
            Imposed tolerance. Used as a parameter of the configuration.
        secure : float
            Security factor, arbitrary value.
        cst_tstep : boolean, optional, default : False
            True if the time step must be constant.
        
        Returns
        -------
        time_step : float
            New time step
        
    """
    if cst_tstep :
        return h
    else :
        if error == 0. :
            return h*secure
        else :
            if error > tolerance :
                s = secure*(tolerance / error)**(0.2)
            else :
                s = secure*(tolerance / error)**(0.25)
                
        time_step = max(h*s, 1e-10)   # if h*s = 0, then min t_step = 1e-10
        return time_step

def calc_K_i(T, h, L, ell, chain, K_list=[], coeff_list=[]) :
    """
    calc_K_i(T, h, L, ell, chain, K_list=[], coeff_list=[])
    
        Calculate the list K_list of runge-kutta coefficients K_i at time T_new for each variable.
        The coefficients are calculated one after the others.
    
        Example : considering RK45 for a single variable y(t), the coefficient k_i are given by
            i = 1 : k_1 = h*f(t, y)
            i = 2 : k_2 = h*f(t + 0.5*h, y + 0.5*k_1)
            i = 3 : k_3 = h*f(t + 0.5*h, y + 0.5*k_2)
            i = 4 : k_4 = h*f(t + h, y + k_3)
        Only one coefficient is calculated at a time, but for a vector Y = [y_1(t), ..., y_N(t)] of N variables.
            
        For informations about the coefficients and their values, refer to Supplementary Informations.
        
        Parameters
        ----------
        T : float
             Old time
        h : float
            Time step
        L : dictionnary
            Dictionnary of the lengths of the lumens, with indices as keys, lengths as values.
        ell : dictionnary
            Dictionnary of the lengths of the bridges, with indices as keys, lengths as values.
        chain : chain-object
            
        K_list : list, optional, default : []
            List of previously calculated coefficient. For instance, K_2 needs K_1, K_3 needs K_1 and K_2, etc.
        coeff_list : []
            List of the coefficients (Butcher table). Depends on the use of RK45 or RKF45.
    
        Returns
        -------
        K : dictionnary
            Dictionnary that contains the coefficients k for each lumens.
    """
    K = {}
    # Empty coeff_list and K_list corresponds to K_1 coefficients,
    if len(K_list) == 0 and len(coeff_list) == 0 :
        for j in chain.lumens_dict.keys() :
            if j != 0 and j != -1 : # Exclude the borders
                K[j] = h*flux.func_Lj_hydraulic(j, T, L, ell, chain)
            else :
                K[j] = 0.
    
    # corresponds to k2, k3, ...
    else :
        if len(coeff_list) != len(K_list) :
            print('ERROR ! Note the same lengths of coefficient or K vectors !')
            
        new_L, new_ell = {}, {}
        
        for j in chain.lumens_dict.keys() :
            if j != 0 and j != -1 : # Exclude the borders
                new_L[j] = L[j]
                for n in range(len(coeff_list)) :
                    new_L[j] += coeff_list[n]*K_list[n][j]
            else :
                new_L[j] = 0. 

        # Closure relation must be imposed : the length L_ij = L_i + L_j + ell_ij is a constant
        for b in chain.bridges_dict.keys() :
            x, y = chain.bridges_dict[b].lumen1, chain.bridges_dict[b].lumen2
            new_ell[b] = ell[b] + (L[x] - new_L[x]) + (L[y] - new_L[y])
        
        for j in chain.lumens_dict.keys() :
            if j != 0 and j != -1 :
                K[j] = h*flux.func_Lj_hydraulic(j, T, new_L, new_ell, chain)
            else :
                K[j] = 0.
    return K

def calc_KQ_i(T, h, L, N, ell, chain, K_list=[], Q_list=[], coeff_list=[]) :
    """
    calc_K_i(T, h, L, N, ell, chain, K_list=[], Q_list=[], coeff_list=[])
    
        Calculate the lists K_list, Q_list of runge-kutta coefficients K_i and Q_i at time T_new for each variable.
        The coefficients are calculated one after the others.
    
        Example : considering RK45 for a single variable y(t), the coefficient k_i are given by
            i = 1 : k_1 = h*f(T, y)
            i = 2 : k_2 = h*f(T + 0.5*h, y + 0.5*k_1)
            i = 3 : k_3 = h*f(T + 0.5*h, y + 0.5*k_2)
            i = 4 : k_4 = h*f(T + h, y + k_3)
        Only one coefficient is calculated at a time, but for a vector Y = [y_1(t), ..., y_N(t)] of N variables.
        The coefficients K correspond to the lengths, the coefficients Q corresponds to the number of ions.
            
        For informations about the coefficients and their values, refer to Supplementary Informations.
        
        Parameters
        ----------
        T : float
             Old time
        h : float
            Time step
        L : dictionnary
            Dictionnary of the lengths of the lumens, with indices as keys, lengths as values.
        N : dictionnary
            Dictionnary of the number of ions of the lumens, with indices as keys, lengths as values.
        ell : dictionnary
            Dictionnary of the lengths of the bridges, with indices as keys, lengths as values.
        chain : chain-object
            
        K_list : list, optional, default : []
            List of previously calculated coefficient. For instance, K_2 needs K_1, K_3 needs K_1 and K_2, etc.
        Q_list : list, optional, default : []
            List of previously calculated coefficient. For instance, Q_2 needs Q_1, Q_3 needs Q_1 and Q_2, etc.
        coeff_list : []
            List of the coefficients (Butcher table). Depends on the use of RK45 or RKF45.
    
        Returns
        -------
        K : dictionnary
            Dictionnary that contains the coefficients k for each lumens.
        Q : dictionnary
            Dictionnary that contains the coefficients q for each lumens.
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
    """
    test_DeltaK(L_vec, ell_vec, Ky, chain)
    
        Test whether the new vectors L_vec and ell_vec are allowed. 
        If one of the numbers of L_vec or ell_vec is negative, return repeat=True.
    
        Parameters
        ----------
        L_vec : dictionnary
            Lengths of the lumens
        ell_vec : dictionnary
            Lengths of the bridges
        Ky : dictionnary
            Coefficients K from Runge-Kutta method to calculate the new step in lengths.
        chain : chain-object
    
        Returns
        -------
        repeat : boolean
            True if the integration needs to be repeated.
        index_L_list : list
            List of lumen indices where Lj < 0
        index_ell_list
            List of bridges indices where ell_j < 0
    """
    # Pre-set outputs
    repeat = False
    index_L_list = []
    index_ell_list = []
    
    # Check for each lumen wheter its new length is strictly positive. If not, repeat is set to True and its index is stored.
    for j in L_vec.keys() :
        if L_vec[j]+Ky[j] < 0 :
            repeat = True
            index_L_list += [j]
            
    # Check for each bridge wheter its new length is strictly positive. If not, repeat is set to True and its index is stored.
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
    """
    rk45_step(t0, chain, h, alpha)
    
        Calculate a step of integration given the previous configuration, using RK45 numerical scheme.
    
        Parameters
        ----------
        t0 : float
            Old time, from which integration is calculated
        chain : chain-object
    
        h : float
            Time-step
        alpha : float
            Constant specified in the configuration file (config['integration']['alpha']). Its value specifies the initial time step.
    
        Returns
        -------
        new_time_step : float
            New time step of the simulation, calculated as
    
                new_time_step = alpha / N(t) ** pow_law
    
            such that as the number of lumens N(t) decreases, the time-step increases.
    """
    global pow_law # Arbitrary coefficient used for an adaptative time-step for RK45 method.
    
    # For safety checks : copy the chain.
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
        
        # If the configuration is not allowed, a new length is negative or zerom thus time-step is divided by two
        # and the integration is recalculated. The incriminated lumen/bridge is flagged in the event-file.
        # If the number of recalculation exceeds 10, a flag is raised in the console (but does not stop the integration).
        if repeat :
            count += 1
            h = 0.5*h
            chain.events += 'ERROR Time : ' + str(chain.time) + ' : length(s) of lumen(s) ' + str(index_L_list) + ' is negative. Time step is divided.\n'
            
            if count >= 10 : print('More than '+str(count)+' trials without convergence. There might be a problem with time-step.')
    
    
    # UPDATE CHAIN
    for j in chain.lumens_dict.keys() :
        chain.lumens_dict[j].length  += K[j]
        if chain.lumen_type == 'hydroosmotic' :
            chain.lumens_dict[j].nb_ions += Q[j]
    # Recalculate the bridges lengths
    net.calc_ell_list(chain)
    
    # Set the new time step
    new_time_step = alpha/(len(chain.lumens_dict.keys()) -2.)**(pow_law)    ## NB : the borders are stored as lumens, so the number of lumens is the length of lumens_dict minus 2
    
    return new_time_step

def rkf45_step(t0, chain, h, tolerance = 1e-6) :
    """
    rkf45_step(t0, chain, h, alpha)
    
        Calculate a step of integration given the previous configuration, using RK45 numerical scheme.
    
        Parameters
        ----------
        t0 : float
            Old time, from which integration is calculated
        chain : chain-object
    
        h : float
            Time-step
        tolerance : float, optional, default : 1e-6
            Tolerance imposed to calculate the adaptative time step.
    
        Returns
        -------
        new_time_step : float
            New time step of the simulation. The error is is given by the max error between estimates in length.
    
        NB : the secure factor S=0.9 is imposed here.
    """
    
    # Makes a copy of the chain
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
        
        # If the configuration is not allowed, a new length is negative or zerom thus time-step is divided by two
        # and the integration is recalculated. The incriminated lumen/bridge is flagged in the event-file.
        # If the number of recalculation exceeds 10, a flag is raised in the console (but does not stop the integration).
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
        
    # Update Time Step by getting the error.
    error_list = []
    for j in Ky.keys() :
        error_list += [abs(Ky[j]-Kz[j])]
        
    error = max(error_list)
        
    # Calculate the new time-step
    new_time_step = calc_new_timestep(h, error, tolerance, secure=0.9, cst_tstep=False)
    return new_time_step
    
# ========================================================================
# ========================= SIMULATION ===================================
# ========================================================================

def system(chain, h=1e-2, recording = False, method = 'rkf45', tolerance=1e-10, alpha=1e-2) :
    """
    system(chain, h=1e-2, recording = False, method = 'rkf45', tolerance=1e-10, alpha=1e-2)
    
        Integrates the chain for a step. Checks whether the simulation must stop or not.
    
        Parameters
        ----------
        chain : chain-object
    
        h : float, optional, default : 1e-2
            Initial time-step.
        recording : boolean, optional, default : False
            If True, records every step of the chain
        method :  string, optional, default : 'rkf45'
            Numerical scheme. Can be 'rk45' or 'rkf45'
        tolerance : float, optional, default : 1e-10
            Tolerance for the adaptative time-step
        alpha : float, optional, default : 1e-2
            Integration constant, used for rk45 method.
        Returns
        -------
        stop : boolean
            If True, the integration stops.
        stop_cause : string
            Contains why the integration stopped.
        new_tstep : float
    """
    stop = False
    t0 = chain.time
    
    # INTEGRATION
    if method == 'rk45' :
        new_tstep = rk45_step(t0, chain, h, alpha=alpha)
    elif method == 'rkf45' :
        new_tstep = rkf45_step(t0, chain, h, tolerance=tolerance)
    else :
        print('Method not recognized')
    
    stop_cause = ''
    
    # CHECK IF SIMULATION STOPS
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
        
    # RECORDS 
    if recording :
        chain.__record_state__()

    # RETURNS
    if method == 'rk45' :
        return stop, stop_cause, new_tstep 
    elif method == 'rkf45' :
        return stop, stop_cause, new_tstep

def make_Nfile(N0, filename, folder) :
    """
    make_Nfile(N0, filename, folder)
    
        Creates the file containing the number of lumens as a function of time.
    
        Parameters
        ----------
        N0 : int
            Initial number of lumens
        filename : string
            Name of the file to write
        folder : string
            Folder of the file.
    """
    file_N = open(os.path.join(folder, filename), 'w')
    file_N.write('#t\tN(t)\n')
    file_N.write(str(0.)+'\t'+str(N0)+'\n')
    file_N.close()

def make_ellfile(ell_avg, filename, folder) :
    """
    make_ellfile(ell_avg, filename, folder)
    
        Creates the file containing the average bridge length as a function of time.
    
        Parameters
        ----------
        ell_avg : int
            Initial average bridge length
        filename : string
            Name of the file to write
        folder : string
            Folder of the file.
    """
    file_ell = open(os.path.join(folder, filename), 'w')
    file_ell.write('#t\tell(t)\n')
    
    file_ell.write(str(0.)+'\t'+str(ell_avg)+'\n')
    file_ell.close()
    
def make_Lfile(L_avg, L_mf, filename, folder) :
    """
    make_Lfile(ell_avg, filename, folder)
    
        Creates the file containing the average lumen length as a function of time.
    
        Parameters
        ----------
        L_avg : int
            Initial average lumen length
        filename : string
            Name of the file to write
        folder : string
            Folder of the file.
    """
    file_L = open(os.path.join(folder, filename), 'w')
    file_L.write('#t\tL(t)\tL_*(t)\n')
    
    file_L.write(str(0.)+'\t'+str(L_avg)+'\t'+str(L_mf)+'\n')
    file_L.close()

def save_distribfile(chain, filename_l, filename_nion, folder) :
    """
    save_distribfile(chain, filename_l, filename_nion, folder)
    
        Save the distribution in length (and number of ions) of the lumens of the chain 
        at the current time chain.time.
    
        Parameters
        ----------
        chain : chain-object
        filename_l : string
            Name of the distribution file for the lengths of lumens.
        filename_nion : string
            Name of the distribution file for the number of ions of lumens.
        folder : string
            Directory where the files are.        
    """
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
    """
    save_N(t, Nt, filename)
    
       Save the number of lumens as a function of time.
    
        Parameters
        ----------
        t : float
            Current time.
        Nt : int
            Number of lumens at time t
        filename : string
            Name of the file to write
    """
    file_N = open(filename, 'a')
    file_N.write(str(t)+'\t'+str(Nt)+'\n')
    file_N.close()
    
def save_ell(t, ellt, filename) :
    """
    save_ell(t, ellt, filename)
    
       Save the average bridge length as a function of time.
    
        Parameters
        ----------
        t : float
            Current time.
        ellt : float
            Average bridge length at time t
        filename : string
            Name of the file to write
    """
    file_ell = open(filename, 'a')
    file_ell.write(str(t)+'\t'+str(ellt)+'\n')
    file_ell.close()
    
def save_L(t, Lt, L_mf, filename) :
    """
    save_L(t, Lt, L_mf, filename)
    
       Save the average lumen length and mean-field length as a function of time.
    
        Parameters
        ----------
        t : float
            Current time.
        Lt : float
            Average lumen length at time t
        L_mf : float
            Mean-field lumen length at time t. The mean-field lumen length is the harmonic average.
        filename : string
            Name of the file to write
    """
    file_L = open(filename, 'a')
    file_L.write(str(t)+'\t'+str(Lt)+'\t'+str(L_mf)+'\n')
    file_L.close()
    
# ========================================================================
# ===========================  RUN  ======================================
# ========================================================================

def run(chain, max_step=1000, alpha=1e-4, savefig=False, nb_frames=1000, dir_name = 'out', recording=False, rec_distrib=False, tolerance=1e-9, solver='rkf45', state_simul=False, pics_dirname='pics', frame_reg=True) :
    """
    run(chain)
    
        Numerical integration of the chain, managing topological events and adaptative time-step. 
        Write the states in output files, stored in the specified directory.
    
        Parameters
        ----------
        chain : object-type
        max_step : int, optional, default : 1000
            Maximum number of integration steps.
        alpha : float, optional, default : 1e-4
            Integration constant. Used to calculate the initial time-step as h = alpha / N0 were N0 is the initial number of lumens.
        savefig : boolean, optional, default : False
            If True, save the chain picture after nb_frames steps.
        nb_frames : int, optional, default : 1000
            Number of steps required before saving data files or pictures of the chain.
        dir_name : string, optional, default : 'out'
            Name of the directory where the chain is. Outputs while be stored there.
        recording : boolean, optional, default : False
            If True, save the state of chain at every time step.
        rec_distrib : boolean, optional, default : False
            If True, save the distributions in length and number of ions of the chain.
        tolerance : float, optional, default : 1e-9
            Tolerance used to calculate the adaptative time step.
        solver : string, optional, default : 'rkf45'
            Method used for the numerical scheme.
        state_simul : boolean, optional, default : False
            If True, print the state of the simulation in the terminal.
        pics_dirname : string, optional, default : 'pics'
            Name of the directory where pictures of the chain are stored.
        frame_reg : boolean, optional, default : True
            If True, register data of the chain every nb_frames.
            Otherwise, this is stored at the end of the simulation.
    
        Returns
        -------
        end : int
            0 : a problem occured during simulation
            1 : one lumen left
            2 : no lumen left
            3 : full size reached
            4 : max step reached
            5 : RunTime error in flow calculation
            
    """
    # Default parameters
    stop = False            # to stop simulation
    step = 0                # simulation step number
    N0 = chain.nmax         # number of lumens
    h = alpha / N0          # integration constant
    x = np.linspace(0, chain.total_length, 1001)     # x-axis positions

    # Create files that contain data
    make_Nfile(N0=N0, filename='sim_nlum.dat', folder = dir_name)
    make_ellfile(ell_avg=chain.__calc_ell_avg__(), filename='sim_ell_avg.dat', folder = dir_name)
    make_Lfile(L_avg=chain.__calc_L_avg__(), L_mf=chain.__calc_L_mf__(), filename='sim_L_avg.dat', folder = dir_name)
    
    if rec_distrib :
        save_distribfile(chain, filename_l = 'distrib_length.dat', filename_nion = 'distrib_nion.dat', folder = dir_name)
    else : 1
    
    #try :
    if 1 :
        # Check the topology before integration.
        tplg.topology(chain)
        
        for i in range(max_step) :
            step += 1
            # make a step
            stop, stop_cause, h = system(chain, h=h, recording = recording, method = solver, tolerance=1e-10, alpha=alpha)
            chain.time += h
            
            # check topology
            if stop_cause != 'full_size' :
                tplg.topology(chain)
                
            # Save events
            tools.save_events(chain, folder=dir_name, filename_events='events.txt')
    
            # Save the number of lumens
            save_N(chain.time, len(chain.lumens_dict)-2, os.path.join(dir_name, 'sim_nlum.dat'))
            
            # Save bridges average length
            ellt_avg = chain.__calc_ell_avg__()
            save_ell(chain.time, ellt_avg, os.path.join(dir_name, 'sim_ell_avg.dat'))
            
            # Save lumens average length
            Lt_avg = chain.__calc_L_avg__()
            Lt_mf = chain.__calc_L_mf__()
            save_L(chain.time, Lt_avg, Lt_mf, os.path.join(dir_name, 'sim_L_avg.dat'))
            
            # Save distributions
            if rec_distrib and step % nb_frames == 0 :
                save_distribfile(chain, filename_l = 'distrib_length.dat', filename_nion = 'distrib_nion.dat', folder = dir_name)
            
            # If the simulation stops
            if stop == 1 :
                # If the cause of stop is topological event, then proceed to topology and reset stop boolean.
                if stop_cause == 'Empty' or stop_cause == 'Fusion' :
                    tplg.topology(chain)
                    stop = 0
                    
                elif stop_cause == 'end_simul':
                    # One lumen left : return 2 ; otherwise return 1
                    
                    # No lumen left
                    if len(chain.lumens_dict) - 2 == 0 :
                        end = 2
                        chain.events += 'Time ' + "{:4.6f}".format(chain.time) + ' : empty chain.\n'
                        tools.save_events(chain, folder=dir_name, filename_events='events.txt')
                        
                        print('End simulation : 0 Lumen left')
                        if recording :
                            tools.save_recording(chain, filename='sim_all.dat', folder=dir_name, chain_type = chain.lumen_type, erase=frame_reg)#
                        break ;
                    
                    # One lumen left
                    elif len(chain.lumens_dict) - 2 == 1 and one_lumen_end :
                        end = 1
                        chain.events += 'Time ' + "{:4.6f}".format(chain.time) + ' : end. One lumen left.'
                        tools.save_events(chain, folder=dir_name, filename_events='events.txt')
                        print('End simulation : 1 Lumen left')
                        
                        if recording :
                            tools.save_recording(chain, filename='sim_all.dat', folder=dir_name, chain_type = chain.lumen_type, erase=frame_reg)#
                        if savefig :
                            tools.plot_profile(x, my_chain, centers=False, lw=1.5, show=False, savefig=True, savename=os.path.join(pics_dirname, 'pic'+str(step).zfill(8)+'.png'))
                        
                        break ;
                
                elif stop_cause == 'full_size' :
                    end = 3
                    chain.events += 'Time ' + "{:4.6f}".format(chain.time) + ' : full chain.\n'
                    tools.save_events(chain, folder=dir_name, filename_events='events.txt')
                    
                    print('End simulation : system is fully occupied')
                    tools.save_recording(chain, filename='sim_all.dat', folder=dir_name, chain_type = chain.lumen_type, erase=frame_reg)#
                    break ;
    
            if step % nb_frames == 0 :
                # Print the state of the chain in terminal
                if state_simul :
                    print('Step : ', step, ' ; Time : ', "{:4.6f}".format(chain.time), ' ; Nb Lumens : ', len(chain.lumens_dict)-2)
                
                # Records the details of the chain
                if recording and frame_reg :
                    tools.save_recording(chain, filename='sim_all.dat', folder=dir_name, chain_type = chain.lumen_type, erase=frame_reg)#
                    
                # Save a picture of the chain
                if savefig == True :
                    tools.plot_profile(x, chain, centers=False, lw=1.5, show=False, savefig=True, savename=os.path.join(pics_dirname, 'pic'+str(step).zfill(8)+'.png'))
                        
            if step == max_step :
                end = 4
                print('\n\nEnd simulation : max step reached')
                chain.events += 'Time ' + "{:4.6f}".format(chain.time) + ' : max step.\n'
                tools.save_events(chain, folder=dir_name, filename_events='events.txt')
                if recording :
                    tools.save_recording(chain, filename='sim_all.dat', folder=dir_name)
            
    #except RuntimeWarning :
    elif 0 :
         # RuntimeWarning is raised only in the flux.py library
         # If a RuntimeWarning exception occurs, it is raised as an exception.
         chain.events += 'Time ' + "{:4.6f}".format(chain.time) + ' : Flow error.\n'
         tools.save_events(chain, folder=dir_name, filename_events='events.txt')

         if recording :
             tools.save_recording(chain, filename='sim_all.dat', folder=dir_name, erase=False)
         print('\nSimulation stopped before the end, flow error.')
         end = 5
        
    #except :
    else :
        chain.events += 'Time ' + "{:4.6f}".format(chain.time) + ' : error occured in integration.\n'
        tools.save_events(chain, folder=dir_name, filename_events='events.txt')
        if recording :
            tools.save_recording(chain, filename='sim_all.dat', folder=dir_name, erase=False)
    
    return end ;

# ========================================================================
# ===========================  MAIN  =====================================
# ========================================================================
        
def main(configname, args) :
    """
    main(configname, args)
    
        Main function of the script. Import a configuration, then run the simulation by integrating numerically the ODEs
        of the system.
        
        Parameters
        ----------
        configname : string
            Name of the configuration file.
        args : list
            List of arguments, that can be
                -v : verbose, print informations on the chain in the terminal.w
                -r : do not record at every frame.
                -e : save the last event.
        
    """
    # Initialize arguments
    state_simul = False                 # Output print on the terminal (time, n of lumens, etc.)
    frame_reg = True                    # Register the state of the chain when step == frame
    rec_last_event = False              # Record the last event
    
    # Import arguments from the console (if any...)
    if len(args) > 0 :
        for arg in args :
            if arg.startswith('-v') :
                state_simul = True
            if arg.startswith('-r') :
                frame_reg = False
            if arg.startswith('-e') :
                rec_last_event = True
    
    # Load Configuration file
    config, my_chain = load_config(configname)
    
    # Parameters of simulation
    max_step = int(config['integration']['max_step'])
    alpha = float(config['integration']['alpha'])
    recording = eval(config['sim']['recording'])
    rec_distrib = eval(config['sim']['rec_distrib'])
    tolerance = float(config['integration']['tolerance'])
    nb_frames = int(config['sim']['nb_frames'])
    solver = config['integration']['solver']
    chain_type = config['sim']['chain_type']
    
    
    dir_name = config['sim']['outdir']
    if not os.path.isdir(dir_name) :
        os.mkdir(dir_name)
    
    # Clean dir_name if not empty
    if len(os.listdir(dir_name)) > 0 :
        for elem in os.listdir(dir_name) :
            if not elem.endswith('.conf') and not os.path.isdir(elem) :
                os.remove(os.path.join(dir_name, elem))
                
    # Write the log.txt file
    tools.write_log(outdir=dir_name, confname=os.path.abspath(configname), args=args)
    
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
            if len(os.listdir(pics_dirname)) > 0 :
                print('Remove files from ' + pics_dirname)
                for elem in os.listdir(pics_dirname) :
                    os.remove(os.path.join(pics_dirname, elem))

        else :
            os.mkdir(pics_dirname)
        
        x = np.linspace(0, my_chain.total_length, 1001)
        tools.plot_profile(x, my_chain, centers=False, lw=1.5, show=False, savefig=True, savename=os.path.join(pics_dirname, 'pic'+str(0).zfill(8)+'.png'))
    
    # Run Simulation
    end = run(my_chain, max_step = max_step, alpha = alpha, recording=recording, rec_distrib=rec_distrib, tolerance=tolerance, nb_frames=nb_frames, solver=solver, savefig=savefig, state_simul=state_simul, dir_name=dir_name,  pics_dirname=pics_dirname, frame_reg=frame_reg)
    
    ### Classification
    # Add the ending event to events.txt
    if 'classification' not in my_chain.__dict__.keys() :
        my_chain.classification = 'collapse'
        
        
    if rec_last_event :
        tools.write_ending_event(end, chain=my_chain, eventfilename=os.path.join(dir_name, 'events.txt'))
    
    # Save the final chain in the end_chain.dat file
    my_chain.__save__(os.path.join(dir_name, 'end_chain.dat'))
    
    # Add ending time to the log.txt file
    tools.add_end_time(outdir=dir_name, end=end)
    
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