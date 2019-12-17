#!/usr/bin/env python
# chain.py

"""
python3 chain.py config.conf [options]


    Options
    -------
    -h : help
    -v : print the steps


    Contains
    --------
    

"""
import os, sys

module_path = os.path.abspath(os.path.join('..', 'chain_lumen/'))

if module_path not in sys.path :
    sys.path.append(module_path)
    #print(sys.path)

import numpy as np
try :
    import configparser
except :
    import ConfigParser

# homemade libraries
try :
    import _ressources.flux as flux
    import _ressources.tools as tools
    import _ressources.network as net
    import _ressources.topology as tplg
    import _ressources.lumenclass as lc
except :
    import flux as flux
    import tools as tools
    import network as net
    import topology as tplg
    import lumenclass as lc

try : 
    import matplotlib.pyplot as plt
    MATPLOTLIB_BOOL = True
    #print('pyplot imported')
except : 
    MATPLOTLIB_BOOL = False
    pass        
# ========================================================================
# ============================ Move ======================================
# ========================================================================
def check_switches(old_pos, new_pos) :
    switches = []
    for k in range(1, len(old_pos)) :
        old_sign = np.sign(old_pos[k, 1]-old_pos[k-1, 1])
        new_sign = np.sign(new_pos[k, 1]-new_pos[k-1, 1])
        if old_sign != new_sign :
            switches += [[new_pos[k], new_pos[k-1]]]
    return switches

def move_lumens(chain, events, motion = 'noisy', avg_displ = 0., std_displ = 1, friction=1, dt = 1., kT = 1.) :
    """
    
    """
    pos_init = chain.__give_positions__()
    
    new_pos = np.copy(pos_init)
    new_pos[:, 1] = np.zeros(len(new_pos[:, 1]))
    new_pos[-1, 1] = pos_init[-1, 1]
    
    
    for k in range(len(new_pos)) :
        index = int(new_pos[k, 0])
        if index != 0 and index != -1 :
            if motion == 'noisy' :
                new_pos[k, 1] = chain.lumens_dict[index].pos
                new_pos[k, 1] += np.random.normal(avg_displ, std_displ)
            elif motion == 'fiction' :
                new_pos[k, 1] = chain.lumens_dict[index].pos
                new_pos[k, 1] += np.random.normal(avg_displ, std_displ) / (chain.lumens_dict[index].length * friction)

            elif motion == 'langevin' :
                new_pos[k, 1] = chain.lumens_dict[index].pos
                new_pos[k, 1] += np.random.normal(avg_displ, std_displ) * np.sqrt(dt) * np.sqrt(2*kT / (friction))
                    
        chain.lumens_dict[index].pos = new_pos[k, 1]

    net.calc_ell_list(chain)
    
    switches = check_switches(pos_init, new_pos)
    for elem in switches :
        if elem[0][0] == 0 or elem[1][0] == 0 or elem[0][0] == -1 or elem[1][0] == -1 :
            if elem[0][0] == 0 or elem[0][0] == -1 :
                index = elem[1][0]
                if elem[0][0] == 0 :
                    displ = chain.lumens_dict[index].length - (chain.lumens_dict[index].pos-chain.lumens_dict[0].pos)
                elif elem[0][0] ==  -1 :
                    displ = -chain.lumens_dict[index].length + (chain.lumens_dict[-1].pos-chain.lumens_dict[index].pos)
                    
            elif elem[1][0] == 0 or elem[1][0] == -1:
                index = elem[0][0]
                if elem[1][0] == 0 :
                    displ = chain.lumens_dict[index].length - (chain.lumens_dict[index].pos-chain.lumens_dict[0].pos)
                elif elem[1][0] == -1 :
                    displ = -chain.lumens_dict[index].length + (chain.lumens_dict[-1].pos-chain.lumens_dict[index].pos)
            tplg.move_lumens_borders(index, displ, chain)
        else :
            i, j = elem[0][0], elem[1][0]
            k, e = tplg.merge_lumens(i, j, net.find_bridge(i, j, chain), chain)
            events += e

# ========================================================================
# =========================== Evolve =====================================
# ========================================================================

def evolve(chain, dt, flux_type='standard') :
    lumens_dict = chain.lumens_dict
    if chain.lumen_type == 'standard' :
        for k in lumens_dict.keys() :
            if k != 0 and k != -1 :
                lumens_dict[k].length += dt*chain.fluxes[k]
    elif chain.lumen_type == 'hydraulic' :
        for k in lumens_dict.keys() :
            if k != 0 and k != -1 :
                new_A = lumens_dict[k].area + dt*chain.fluxes[k]
                lumens_dict[k].length = np.sign(new_A)*np.sqrt(lumens_dict[k].mu*np.abs(new_A))
    elif chain.lumen_type == 'hydroosmotic' :
        for k in lumens_dict.keys() :
            if k != 0 and k!= -1 :
                lumens_dict[k].length += flux.func_Lj(lumens_dict[k], chain.fluxes[k][0], chain.tauv)*dt
                lumens_dict[k].nb_ions += flux.func_Nj(lumens_dict[k], chain.fluxes[k][1], chain.taus)*dt

def integrate(chain, dt, threshold=0.5, flux_type='standard', kappa=1, flux_val=1e-2) :
    
    chain.fluxes = flux.calc_fluxes(chain, threshold,kappa=kappa, flux_val=flux_val)
    
    evolve(chain, dt, flux_type=flux_type)
    
    net.calc_ell_list(chain)
    #chain.__update_chain__()
    
def func_ivp(t, y) :
    fl.calc_fluxes(my_chain, flux_type=my_chain.lumen_type)
    my_chain.__make_flux_vec__()
    return my_chain.flux_vec
    
def func_odeint(y, t) :
    fl.calc_fluxes(my_chain, flux_type=my_chain.lumen_type)
    my_chain.__make_flux_vec__()
    print(my_chain.flux_vec)
    return my_chain.flux_vec
    
def integrate_bis(chain, config) :
    sim = True

    nb_pts = int(config['integration']['nb_pts'])
    alpha = float(config['sim']['alpha'])
    
    t_start = 0
    t_end = float(config['integration']['time_interval'])
    t_eval = np.linspace(t_start, t_end, nb_pts)

    my_chain.__make_Y0__()

    if config['integration']['solver'] == 'ivp' :
        sol = solve_ivp(func_ivp, t_span=[t_start, t_end], y0=my_chain.y0, t_eval=t_eval)
    elif config['integration']['solver'] == 'odeint' :
        sol = odeint(func_odeint, y0=my_chain.y0, t=t_eval)
    
    return sol
    
# ========================================================================
# ============================ Simul =====================================
# ========================================================================

def run_simul(step, chain, t0, threshold=0.5, motion=False, flux_type='standard', kappa=1, friction=1, max_step=10000, flux_val=1e-2, recording=False, tolerance = 1e-4) :
    stop = False
    
    nmax = max(chain.lumens_dict.keys())
    dt = t0 / ((chain.nb_lumens)**3)
    chain.time += dt
    
    # Integrate
    chain = integrate(chain, dt, threshold, kappa=kappa, flux_val=flux_val)
    
    # Motion
    if motion :
        move_lumens(chain, '', motion = motion, avg_displ = 0., std_displ = 1, friction=friction, dt = dt, kT = 1.)
        
    # Topology
    chain.events = tplg.topology(chain, events)
    
    # Recording
    if recording :
        chain.__record_state__()
    
    total_length = chain.__calc_total_length__()
    # stop condition
    if chain.nb_lumens <= 1 or step > max_step or np.abs(total_length - chain.total_length) > tolerance :
        stop = True
        print(total_length)
        
    return stop, events

def load_config(filename) :
    global my_chain
    
    config = configparser.ConfigParser()
    config.read(filename)

    
    path = config['sim']['path']
    
    if len(path) > 0 :
        print('Import config from ' + path)      
        lumens_array = np.loadtxt(path+'lumens.dat', delimiter = '\t', comments = '#')
        bridges_array = np.loadtxt(path+'bridges.dat', delimiter = '\t', comments = '#')
        
        my_chain = lc.Osmotic_Chain(nb_lumens = len(lumens_array)-2, taus=float(config['hydroosmotic']['taus']), tauv=float(config['hydroosmotic']['tauv']), e0=float(config['sim']['e0']), l_merge=float(config['topology']['l_merge']), l_dis=float(config['topology']['l_dis']))
        
        my_chain.__import_config__(lumens_array, bridges_array, eps = float(config['topology']['eps']))
        
        my_chain.pumping = config['pumping']['pattern']
    else :
        if config.has_option('sim', 'seed') and len(config['sim']['seed']) > 0 :
            np.random.seed(int(config['sim']['seed']))
            
        my_chain = lc.Osmotic_Chain(nb_lumens = int(config['sim']['M']), taus=float(config['hydroosmotic']['taus']), tauv=float(config['hydroosmotic']['tauv']), e0=float(config['sim']['e0']), l_merge=float(config['topology']['l_merge']), l_dis=float(config['topology']['l_dis']))
        
        if config.has_option('pumping', 'pattern') :
            #print('No pumping')    
            my_chain.pumping = config['pumping']['pattern']
            if config['pumping']['pattern'] == 'normal' :
                
                br_avg, br_std = float(config['pumping']['bridge_1']), float(config['pumping']['bridge_2'])
                lum_avg, lum_std = float(config['pumping']['lumen_1']), float(config['pumping']['lumen_2'])
                
                ca_bridge_list = [np.random.normal(br_avg, br_std) for b in range(my_chain.nb_lumens+1)]
                ca_lumen_list = [np.random.normal(lum_avg, lum_std) for m in range(my_chain.nb_lumens+2)]
            
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
            
            else :
                my_chain.pumping = 'None'
                ca_lumen_list = [0. for i in range(my_chain.nb_lumens+2)]
                ca_bridge_list = [0. for i in range(my_chain.nb_lumens+1)]
                
        else :
            my_chain.pumping = 'None'
            ca_lumen_list = [0. for i in range(my_chain.nb_lumens+2)]
            ca_bridge_list = [0. for i in range(my_chain.nb_lumens+1)]
        
        equilibrium = eval(config['hydroosmotic']['equilibrium'])
        nions_avg = float(config['hydroosmotic']['nions_avg'])
        nions_std = float(config['hydroosmotic']['nions_std'])
        
        my_chain.__gen_network_lumen_object__(avg_size=float(config['topology']['avg_size']), std_size=float(config['topology']['std_size']), avg_dist=float(config['topology']['avg_dist']), std_dist=float(config['topology']['std_dist']), dist_toleft=float(config['topology']['dist_toleft']), dist_toright=float(config['topology']['dist_toright']), eps = float(config['topology']['eps']), equilibrium=equilibrium, nions_avg=nions_avg, nions_std=nions_std, ca_lumen_list=ca_lumen_list, ca_bridge_list=ca_bridge_list)
    
    chis = float(config['hydroosmotic']['chis'])
    chiv = float(config['hydroosmotic']['chiv'])
    
    my_chain.xis = chis*my_chain.bridges_dict[1].length
    my_chain.xiv = chiv*my_chain.bridges_dict[1].length
    #my_chain.xis = chis*my_chain.total_length
    #my_chain.xiv = chiv*my_chain.total_length
    
    #print(my_chain)
    return config, my_chain

# ========================================================================
# ======================== Runge-Kutta ===================================
# ========================================================================

def temporary_chain(Lj, Nj, chain, index) :
    # Create a temporary new lumen with updated values
    lumen_j = chain.lumens_dict[index].__copy__()
    #print(lumen_j.length, Lj, ';', lumen_j.nb_ions, Nj)
    lumen_j.length = Lj
    lumen_j.nb_ions = Nj
    
    # Find the neighbors
    i, k = net.find_neighbors(index, chain.bridges_dict)
    lumen_i = chain.lumens_dict[i].__copy__()
    lumen_k = chain.lumens_dict[k].__copy__()
    
    # Find the connecting bridges
    br_ij, br_jk = net.connected_bridges(index, chain.bridges_dict)
    new_br_ij = chain.bridges_dict[br_ij].__copy__()
    new_br_jk = chain.bridges_dict[br_jk].__copy__()
    
    # Calculate modified ell_ij, ell_jk 
    ell_ij = np.abs(lumen_j.pos - lumen_i.pos) - (lumen_i.length + lumen_j.length)
    ell_jk = np.abs(lumen_j.pos - lumen_k.pos) - (lumen_k.length + lumen_j.length)
    new_br_ij.length = ell_ij
    new_br_jk.length = ell_jk
    
    return lumen_j, lumen_i, lumen_k, new_br_ij, new_br_jk

def calc_new_timestep(error, tolerance, secure=0.9, cst_tstep=1) :

    if cst_tstep :
        return 1.
    else :
        if error == 0. :
            #print('Error !', error)
            return 1.
        else :
            if error > tolerance :
                #s = secure*(tolerance / error)**(0.2)
                ratio = (tolerance / (2.*error))**0.5
                s = secure*min(2., max(0.3, ratio))
            else :
                #s = secure*(tolerance / error)**(0.25)
                ratio = (tolerance / (2.*error))**0.5
                s = secure*min(2., max(0.3, ratio))
        #print(s)
        return s

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

def test_DeltaK(L_vec, Ky) :
    repeat = False
    index_list = []
    for j in L_vec.keys() :
        if L_vec[j]+Ky[j] < 0 :
            print(j)
            repeat = True
            index_list += [j]
    
    if len(index_list) == 0 :
        index_list = None
        
    return repeat, index_list
        
def rk45_step(t0, chain, h) :

    L_vec, N_vec, ell_vec = chain.__L_vec__(), chain.__N_vec__(), chain.__ell_vec__()
    repeat = True
    while repeat :
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
        repeat, index = test_DeltaK(L_vec, K)
        count = 0
        
        if repeat :
            count += 1
            h = 0.5*h
            chain.events += 'ERROR Time : ' + str(chain.time) + ' : length(s) of lumen(s) ' + str(index) + ' is negative. Time step is divided.\n'
            
        if count >= 10 : print('More than 10 trials without convergence. There might be a problem.')
    
    
    # UPDATE CHAIN
    chain.time = t0+h
    for j in chain.lumens_dict.keys() :
        chain.lumens_dict[j].length  += K[j]
        chain.lumens_dict[j].nb_ions += Q[j]
    
    net.calc_ell_list(chain)

def rkf45_step(t0, chain, h, tolerance = 1e-6) :
    cp_chain = chain.__copy__()
    
    L_vec, N_vec, ell_vec = chain.__L_vec__(), chain.__N_vec__(), chain.__ell_vec__()

    repeat = True
    
    count = 0
    while repeat :
        
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
        repeat, index = test_DeltaK(L_vec, Kz)
        
        if repeat :
            count += 1
            h = 0.5*h
            chain.events += 'ERROR Time : ' + str(chain.time) + ' : length(s) of lumen(s) ' + str(index) + ' is negative. Time step is divided.\n'
        
        if count >= 10 : print('More than 10 trials without convergence. There might be a problem.')
    
    # UPDATE CHAIN
    chain.time = t0+h
    for j in chain.lumens_dict.keys() :
        chain.lumens_dict[j].length  += Kz[j]
        chain.lumens_dict[j].nb_ions += Qz[j]
        
    net.calc_ell_list(chain)
    
    # Update Time Step
    error_list = []
    for j in Ky.keys() :
        #print(abs(Ky[j]-Kz[j]))
        error_list += [abs(Ky[j]-Kz[j]), abs(Qy[j] - Qz[j])]
    error = max(error_list)
        
    new_time_step = h*calc_new_timestep(error, tolerance, secure=0.9, cst_tstep=0)
    #print(new_time_step)
    return new_time_step
    
# ========================================================================
# ========================= SIMULATION ===================================
# ========================================================================

def system(chain, h=1e-2, recording = False, method = 'rk45', tolerance=1e-10) :
    stop = False
    t0 = chain.time
    
    if method == 'rk45' :
        rk45_step(t0, chain, h)
    elif method == 'rkf45' :
        new_tstep = rkf45_step(t0, chain, h, tolerance=tolerance)
    else :
        print('Method not recognized')
        
    stop_cause = ''
    
    if len(tplg.check_emptylumens(chain)) > 0 or len(tplg.check_merginglumens(chain)) > 0 or len(tplg.empty_ions(chain)) > 0:
        1
        
    if len(tplg.check_emptylumens(chain)) > 0 :
        stop = True
        stop_cause = 'Empty'
        
    elif len(tplg.check_merginglumens(chain)) > 0 :
        stop = True
        stop_cause = 'Empty'
        
    elif len(tplg.empty_ions(chain)) > 0 :
        stop = True
        stop_cause = 'Empty_ions'
        
    elif chain.nb_lumens == 1 and 2*chain.lumens_dict[max(chain.lumens_dict.keys())].length >= chain.total_length :
        stop = True
        stop_cause = 'end_simul'
        
    elif chain.nb_lumens <= 1 :
        stop = True
        stop_cause = 'end_simul'
        
    if recording :
        chain.__record_state__()
    
    if method == 'rk45' :
        return stop, stop_cause, h
    elif method == 'rkf45' :
        return stop, stop_cause, new_tstep
    
def save_N(t, Nt, filename) :
    file_N = open(filename, 'a')
    file_N.write(str(t)+'\t'+str(Nt)+'\n')
    file_N.close()
    
def run(chain, max_step=1000, alpha=1e-4, savefig=False, nb_frames=1000, dir_name = 'out', recording=True, tolerance=1e-9, solver='rkf45', state_simul=False, pics_dirname='pics') :
    
    stop = False
    step = 0
    N0 = chain.nmax
    h = alpha / N0
    
    x = np.linspace(0, my_chain.total_length, 1001)

    #os.rmdir(dirname)

    file_N = open(os.path.join(dir_name, 'sim_nlum.dat'), 'w')
    file_N.write('#t\tN(t)\n')
    file_N.write(str(0.)+'\t'+str(N0)+'\n')
    file_N.close()
    
    if 1 :
    #try :
        for i in range(max_step) :
            step += 1
                        
            # 
            stop, stop_cause, h = system(chain, h=h, recording = recording, method = solver, tolerance=1e-10)
            tplg.topology(chain)
    
            # Save the number of lumens
            save_N(chain.time, len(chain.lumens_dict)-2, os.path.join(dir_name, 'sim_nlum.dat'))
            
            if stop==1 :
                if stop_cause == 'end_simul':
                    # One lumen left : returns 2 ; otherwise return 1
                    if len(chain.lumens_dict) - 2 == 1 :
                        end = 2
                        print('End simulation : 1 Lumen left')
                        
                        for k in chain.lumens_dict.keys() : 
                            if k != 0 and k != -1 :
                                my_chain.winner = k
                        
                        #chain.events += 'Time : ' + "{:4.6f}".format(chain.time) + ' : winner is lumen ' + str(int(winner))
                    elif len(chain.lumens_dict) - 2 == 0 :
                        end = 1
                        print('End simulation : 0 Lumen left')
                    break ;
            
                elif stop_cause == 'Empty' or stop_cause == 'Fusion' :
                    tplg.topology(chain)
                    stop = 0

                
            if savefig == True and step % nb_frames == 0 :
                tools.plot_profile(x, chain, centers=False, lw=1.5, show=False, savefig=True, savename=os.path.join(pics_dirname, 'pic'+str(step).zfill(8)+'.png'))
    
            if step % nb_frames == 0 :
                if state_simul :
                    #print('Step : ', step, ' ; Time : ', "{:4.6f}".format(chain.time), ' ; Nb Lumens : ', len(chain.lumens_dict)-2, end='\t\r')
                    print('Step : ', step, ' ; Time : ', "{:4.6f}".format(chain.time), ' ; Nb Lumens : ', len(chain.lumens_dict)-2)
                save_data = recording
                if save_data :
                    tools.save_recording(chain, filename='sim_all.dat', filename_events='events.log', folder=dir_name)
                    
            if step == max_step :
                print('\n\nEnd simulation : max step reached')
                tools.save_recording(chain, filename='sim_all.dat', filename_events='events.log', folder=dir_name)

        
        if i == max_step-1 :
            end = 10
            
        if savefig :
            
            tools.plot_profile(x, chain, centers=False, lw=1.5, show=0, savefig=True, savename=os.path.join(pics_dirname, 'pic'+str(step).zfill(8)+'.png'))        

        save_data = recording
        if save_data :
            tools.save_recording(chain, filename='sim_all.dat', filename_events='events.log', folder=dir_name)
    else :
    #except :
        tools.save_recording(chain, filename='sim_all.dat', filename_events='events.log', folder=dir_name)
        print('\n\nSimulation stopped before the end...')
        end = 0
    
    return end ;
        
def main(configname, args) :
    # Initialize arguments
    state_simul = False
    
    # Import arguments from the console (if any...)
    if len(args) > 0 :
        for arg in args :
            if arg.startswith('-v') :
                state_simul = True
    
    # Load Configuration file
    config, my_chain = load_config(configname)
    
    # Parameters of simulation
    max_step = int(config['integration']['max_step'])
    alpha = float(config['integration']['alpha'])
    recording = eval(config['sim']['recording'])
    tolerance = float(config['integration']['tolerance'])
    nb_frames = int(config['sim']['nb_frames'])
    solver = config['integration']['solver']
    
    dir_name = config['sim']['outdir']
    
    if MATPLOTLIB_BOOL == True :
        savefig = eval(config['sim']['savefig'])
    else :
        savefig = False
        
    my_chain.__save__(os.path.join(dir_name, 'init_chain.dat'))
    
    pics_dirname=''
    if savefig :
        pics_dirname = os.path.join(dir_name,'pics')
        #print(pics_dirname)
        try :
        #if 1:
            os.mkdir(pics_dirname)
            print(pics_dirname)
        except :
            print('Remove previous pictures from ' + pics_dirname)
            for elem in os.listdir(pics_dirname) :
               os.remove(os.path.join(pics_dirname, elem))
    
    # Run Simulation
    end = run(my_chain, max_step = max_step, alpha = alpha, recording=recording, tolerance=tolerance, nb_frames=nb_frames, solver=solver, savefig=savefig, state_simul=state_simul, dir_name=dir_name,  pics_dirname=pics_dirname)
    
    if end == 2 :
        f = open(os.path.join(dir_name, 'init_chain.dat'), 'a+')
        f.write('========= END =========\n')
        f.write('Winner : ' + str(int(my_chain.winner)))
        f.close()
        
    
    
    # Move the config file into the directory
    #os.rename(configname, os.path.join(dir_name, configname))
    
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