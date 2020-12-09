"""
analysis.py is a library of functions used in analysis

    Contains
    --------
        AVERAGE
    gen_log_times           : Calculate a time-list in logspace, given tmin and tmax
    find_argmins            : Find the argument (index) of the times with minimal distances from log_times.
    gen_index_array         : Gives the indices of times corresponding to the minimal distance between choosen time-points in logscale.
        BATCH
    batch_window            : Make batches of data in a given window.
    batch_average           : Make average of a batch.
    batch                   : Calculate an average of values contained in data_dict based on "windows".
        IMPORT
    make_path_dict          : Gives the list of path of the simulation folders.
    import_osmotic          : Import sim_nlum.dat files and returns a dictionnary from it.
        PLOT
    plot_conf               : Plot the  number of lumens vs time given the choosen configuration
    plot_osmotic_pumping    : Plot the number of lumens vs time for HO-chain with active pumping.
    plot_osmotic            : Plot the number of lumens vs time for HO-chain without active pumping.

    Requirements
    ------------

Mathieu Le Verge--Serandour, 2020
"""

import numpy as np
try : import matplotlib.pyplot as plt
except : pass

import os, sys


ca_dict = {0.01 : ['ca1e-2', 'b'], 0.1 : ['ca1e-1', 'r'], 1 : ['ca1e0', 'purple'],10 : ['ca1e1', 'orange']}
lw = 1. # Linewidth
ms = 4  # Marker size
chi_dict = {(5, 5) : ['^', ms, 1., lw],
            (5, 1) : ['^', ms, 1., lw],
            (50, 5) : ['s', ms, 1., lw],
            (50, 50) : ['+', ms, 1., lw],
            (50, 1) : ['+', ms, 1., lw],
            (500, 500) : ['s', ms, 1., lw],
            (500, 1) : ['s', ms, 1., lw],
            (500, 5) : ['o', ms, 1., lw],
           }
           
# AVERAGE
def gen_log_times(tmin, tmax, npts) :
    """
    gen_log_times(tmin, tmax, npts)
    
        Calculate a time-list in logspace, given tmin and tmax
    
        Parameters
        ----------
        tmin : float
            Minimal time
        tmax : float
            Maximal time
    
        Returns
        -------
        T : list
            List of times generated in logscale.
    """
    T = np.logspace(np.log10(tmin), np.log10(tmax), npts)
    return T

def find_argmins(Nt, log_times) :
    """
    find_argmins(Nt, log_times)
        
        Find the argument (index) of the times with minimal distances from log_times.
        
        Parameters
        ----------
        Nt : list
            Array that represents the number of lumens versus time, such that Nt[:, 0] are the time-steps, Nt[:, 1] the number of lumens
        log_times : list
            List of times (in log-scale)
        Returns
        -------
        index_array : list
            list of indices
    """
    index_array = np.zeros(len(log_times), dtype=int)
    for k in range(len(log_times)) :
        index_array[k] = np.argmin(np.abs(log_times[k]-Nt[:, 0]))
    return index_array

def gen_index_array(Nt, npts) :
    """
    gen_index_array(Nt, npts)
    
        Gives the indices of times corresponding to the minimal distance between choosen time-points in logscale.
    
        Parameters
        ----------
        Nt : array
            Array that represents the number of lumens versus time, such that Nt[:, 0] are the time-steps, Nt[:, 1] the number of lumens
        npts : int
            Number of points
    
        Returns
        -------
        index_array : array
            List of indices of time-points for the simulation.
    """
    size = len(Nt)
    
    tmin = np.min(Nt[1:, 0])   # skip the first time since t=0
    tmax = np.max(Nt[1:, 0])   # skip the first time since t=0
    
    log_times = gen_log_times(tmin, tmax, npts)
    index_array = find_argmins(Nt, log_times)
    
    return index_array

# BATCH
def batch_window(data, wmin, wmax, nwindow) :
    """
    batch_window(data, wmin, wmax, nwindow)
    
        Create a list of intervals/windows in log-scale.
        Make batches of data in a given window.
    
        Parameters
        ----------
        data :  array
            List of data. Must be an array of at least 2 dimension.
        wmin : float
            Minimal time for the intervals
        wmax : float
            Maximal time for the intervals
        nwindow : int
            Number of windows
    
        Returns
        -------
        batch : list
            List of batched data
    """
    window = np.logspace(wmin, wmax, nwindow)
    time = np.cumsum(window)
    batch = []
    for i in range(len(time)) :
        indices = np.argwhere(np.abs(data[:, 0] - time[i]) <= window[i])[:, 0]
        batch += [data[indices]]
    return batch    

def batch_average(batchlist) :
    """
    batch_average(batchlist)
        
        Make average of a batch.
        
        Parameters
        ----------
        batchlist : list
            List of data to batc, such as simulations
            Example :
            batchlist = [ [[0, 0], [1, 2], [2, 4], ...],
                          [[0, 3], [1, 2], [2, 8], ...]
                        ]
        Returns
        -------
        B_avg : array
            Array of averaged data
        B_std : array
            Array of the standard-deviations of the data
    """
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
    """
    batch(data_dict, wmin, wmax, nwindow)
        
        Calculate an average of values contained in data_dict based on "windows".
        The number of windows and the min and max windows give points arround whch the windows are taken
        and all the values contained in these windows are averaged to form a single point, with standard deviation.
        
        Parameters
        ----------
        data_dict : dict
            Dictionnary of data to average.
        wmin : float
            Minimal value (in log-scale) of the averaging area
        wmax : float
            Maximal value (in log-scale) of the averaging area
        nwindow : int
            Number of windows
    
        Returns
        -------
        B_avg : list
            List of averaged values
        B_std :
            List of standard deviations associated to the averaged values
    """
    window = np.logspace(wmin, wmax, nwindow)
    time = np.cumsum(window)
    dat_batch_list = []
    for k in data_dict.keys() :
        dat_batch_list += [batch_window(data_dict[k], wmin=wmin, wmax=wmax, nwindow=nwindow)]
        print(k, end='\r')
    print('End of import !')
    B_avg, B_std = batch_average(dat_batch_list)
    return B_avg, B_std

# IMPORT
def make_path_dict(nsim, main_dir, subdir, subsubdir='') :
    """
    make_path_dict(nsim, main_dir, subdir, subsubdir=')
        Gives the list of path of the simulation folders.
    
        Parameters
        ----------
        nsim : int
            Number of simulation
        main_dir : path
        subdir : path
        subsubdir : path, optional, default : ''
            
        Returns
        -------
        pathdict :
            List of paths of each simulation.
    
    """
    pathdict = {}
    for n in range(nsim) :
        pathdict[n] = os.path.join(main_dir, subdir, subsubdir, 'run'+str(n).zfill(4))
    return pathdict

def import_osmotic(chiv, chis, path_list, ca = None) :
    """
    import_osmotic(chiv, chis, path_list, ca = None)
    
        Import sim_nlum.dat files and returns a dictionnary from it.
    
        Parameters
        ----------
        chiv : float
            Value of chis for the folders
        chis : float
            Value of chis for the folders
        pathlist : list
            List of path of the simulation folders, indexed by pathlist[(chiv, chis)][n]
        ca : string, optional, default : None
            If specified, value of active pumping.
        Returns
        -------
        Nt : dict
            Dictionnary containing all N(t) vs t arrays for the simulations
    """
    Nt = {}
    if ca == None :
        for n in path_list[(chiv, chis)].keys() :
            Nt[n] = np.loadtxt(os.path.join(path_list[(chiv, chis)][n], 'sim_nlum.dat'))
    else :
        for n in path_list[(chiv, chis)][ca].keys() :
            Nt[n] = np.loadtxt(os.path.join(path_list[(chiv, chis)][ca][n], 'sim_nlum.dat'))
    return Nt

# PLOT
def plot_conf(ca, chis, chiv, Nt, npts=20,  wmin=-6, wmax=3, rescaled=False, rescaled_pumping=False, ell0=10, L0=12000, tau=1, show_sim=False, ax=None) :
    """
    plot_conf(ca, chis, chiv, Nt, npts=20,  wmin=-6, wmax=3, rescaled=False, rescaled_pumping=False, ell0=10, L0=12000, tau=1, show_sim=False, ax=None)
        
        Plot the  number of lumens vs time given the choosen configuration
        
        chi_dict, ca_dict are colors, linestyle, etc.
    """
    global chi_dict, ca_dict
    mu, nu, eps = 0.6105653703843762, 1.2091995761561452, 1e-3
    chi = (chis, chiv)
    xiv = ell0*chiv
    T = (2*tau*L0**2)/((xiv**2)*mu*eps)
    if ca != 0 :
        T_ca = 1./(mu*nu*ca)
    else :
        print('No pumping!')
        T_ca = 1.
    if show_sim :
        for n in range(nsim) :
            index_array = gen_index_array(Nt[(chis, chiv)][ca][n], npts)
            if rescaled :
                if ax != None :
                    ax.plot(Nt[(chis, chiv)][ca][n][index_array, 0]/T, Nt[(chis, chiv)][ca][n][index_array, 1], color = ca_dict[ca][1], marker=chi_dict[chi][0], markersize=chi_dict[chi][1], alpha = 0.1)
                else :
                    plt.plot(Nt[(chis, chiv)][ca][n][index_array, 0]/T, Nt[(chis, chiv)][ca][n][index_array, 1], color = ca_dict[ca][1], marker=chi_dict[chi][0], markersize=chi_dict[chi][1], alpha = 0.1)
            elif rescaled_pumping :
                if ax != None :
                    ax.plot(Nt[(chis, chiv)][ca][n][index_array, 0]/T_ca, Nt[(chis, chiv)][ca][n][index_array, 1], color = ca_dict[ca][1], marker=chi_dict[chi][0], markersize=chi_dict[chi][1], alpha = 0.1)
                else :
                    plt.plot(Nt[(chis, chiv)][ca][n][index_array, 0]/T_ca, Nt[(chis, chiv)][ca][n][index_array, 1], color = ca_dict[ca][1], marker=chi_dict[chi][0], markersize=chi_dict[chi][1], alpha = 0.1)
            else :
                if ax != None :
                    ax.plot(Nt[(chis, chiv)][ca][n][index_array, 0], Nt[(chis, chiv)][ca][n][index_array, 1], color = ca_dict[ca][1], marker=chi_dict[chi][0], markersize=chi_dict[chi][1], alpha=0.1)
                else :
                    plt.plot(Nt[(chis, chiv)][ca][n][index_array, 0], Nt[(chis, chiv)][ca][n][index_array, 1], color = ca_dict[ca][1], marker=chi_dict[chi][0], markersize=chi_dict[chi][1], alpha=0.1)

    res_avg, res_std = an.batch(data_dict=Nt[(chis, chiv)][ca], wmin=wmin, wmax=wmax, nwindow=npts)
    if rescaled :
        if ax != None :
            ax.plot(res_avg[:, 0]/T, res_avg[:, 1], color = ca_dict[ca][1], marker=chi_dict[chi][0], markersize=chi_dict[chi][1], label = r'$c^a = ' + str(ca) + '$')
        else :
            plt.plot(res_avg[:, 0]/T, res_avg[:, 1], color = ca_dict[ca][1], marker=chi_dict[chi][0], markersize=chi_dict[chi][1], label = r'$c^a = ' + str(ca) + '$')
    elif rescaled_pumping :
        if ax != None :
            ax.plot(res_avg[:, 0]/T_ca, res_avg[:, 1], color = ca_dict[ca][1], marker=chi_dict[chi][0], markersize=chi_dict[chi][1], label = r'$c^a = ' + str(ca) + '$')
        else :
            plt.plot(res_avg[:, 0]/T_ca, res_avg[:, 1], color = ca_dict[ca][1], marker=chi_dict[chi][0], markersize=chi_dict[chi][1], label = r'$c^a = ' + str(ca) + '$')
    else :
        if ax != None :
            ax.plot(res_avg[:, 0], res_avg[:, 1], color = ca_dict[ca][1], marker=chi_dict[chi][0], markersize=chi_dict[chi][1], label = r'$c^a = ' + str(ca) + '$')
        else :
            plt.plot(res_avg[:, 0], res_avg[:, 1], color = ca_dict[ca][1], marker=chi_dict[chi][0], markersize=chi_dict[chi][1], label = r'$c^a = ' + str(ca) + '$')
    
def plot_osmotic_pumping(Nt, ca_bools, chi_bools, npts=50, rescaled = False, rescaled_pumping = False, show_sim=False, scaling_laws=True, savefig=False, savename='H0-coarsening_pumping.png', extension='png') :
    """
    plot_osmotic_pumping(chiv, chis, path_list, plot_param_list, Nt_list, rescale = False)
    
        Plot the number of lumens vs time for HO-chain with active pumping.
    """
    plt.figure(figsize=(6, 6))
    plt.xscale('log')
    plt.yscale('log')

    try :
        chis5_chiv5, chis50_chiv50, chis500_chiv500, chis5_chiv500 = chi_bools
    except :
        chis5_chiv5, chis50_chiv50, chis500_chiv500 = chi_bools
        chis_chiv500 = 0
        
    ca_0, ca_1e_2, ca_1e_1, ca_1e0, ca_1e1 = ca_bools
    
    # chis, chiv = (5, 5)
    if chis1_chiv5 :
        chis, chiv = 1, 5
        # ca = 0.
        #plot_conf(ca=0., chis=chis, chiv=chiv, Nt=Nt_N1000, npts=npts, wmin=-6, wmax=4, rescaled=rescaled, ell0=10, L0=1, tau=1, show_sim=show_sim)
        # ca = 0.01
        if ca_1e_2 :
            plot_conf(ca=0.01, chis=chis, chiv=chiv, Nt=Nt_N1000, npts=npts, wmin=-6, wmax=6, rescaled=rescaled, rescaled_pumping=rescaled_pumping, ell0=10, L0=1, tau=1, show_sim=show_sim)
        # ca = 0.1
        if ca_1e_1 :
            plot_conf(ca=0.1, chis=chis, chiv=chiv, Nt=Nt_N1000, npts=npts, wmin=-6, wmax=5, rescaled=rescaled, rescaled_pumping=rescaled_pumping, ell0=10, L0=1, tau=1, show_sim=show_sim)
        # ca = 1
        if ca_1e0 :
            plot_conf(ca=1, chis=chis, chiv=chiv, Nt=Nt_N1000, npts=npts, wmin=-6, wmax=4, rescaled=rescaled, rescaled_pumping=rescaled_pumping, ell0=10, L0=1, tau=1, show_sim=show_sim)
        # ca = 10
        if ca_1e1 :
            plot_conf(ca=10, chis=chis, chiv=chiv, Nt=Nt_N1000, npts=npts, wmin=-6, wmax=4, rescaled=rescaled, rescaled_pumping=rescaled_pumping, ell0=10, L0=1, tau=1, show_sim=show_sim)

    # chis, chiv = (50, 50)
    if chis1_chiv50 :
        chis, chiv = 1, 50
        # ca = 0.
        #plot_conf(ca=0., chis=chis, chiv=chiv, Nt=Nt_N1000, npts=npts, wmin=-6, wmax=4, rescaled=rescaled, ell0=10, L0=1, tau=1, show_sim=show_sim)
        # ca = 0.01
        if ca_1e_2 :
            plot_conf(ca=0.01, chis=chis, chiv=chiv, Nt=Nt_N1000, npts=npts, wmin=-6, wmax=6, rescaled=rescaled, rescaled_pumping=rescaled_pumping, ell0=10, L0=1, tau=1, show_sim=show_sim)
        # ca = 0.1
        if ca_1e_1 :
            plot_conf(ca=0.1, chis=chis, chiv=chiv, Nt=Nt_N1000, npts=npts, wmin=-6, wmax=5, rescaled=rescaled, rescaled_pumping=rescaled_pumping, ell0=10, L0=1, tau=1, show_sim=show_sim)
        # ca = 1
        if ca_1e0 :
            plot_conf(ca=1, chis=chis, chiv=chiv, Nt=Nt_N1000, npts=npts, wmin=-6, wmax=4, rescaled=rescaled, rescaled_pumping=rescaled_pumping, ell0=10, L0=1, tau=1, show_sim=show_sim)
        # ca = 10
        if ca_1e1 :
            plot_conf(ca=10, chis=chis, chiv=chiv, Nt=Nt_N1000, npts=npts, wmin=-6, wmax=4, rescaled=rescaled, rescaled_pumping=rescaled_pumping, ell0=10, L0=1, tau=1, show_sim=show_sim)

    # chis, chiv = (500, 500)
    if chis1_chiv500 :
        chis, chiv = 1, 500
        # ca = 0.
        if ca_0 :
            plot_conf(ca=0., chis=chis, chiv=chiv, Nt=Nt_N1000, npts=npts, wmin=-6, wmax=4, rescaled=rescaled, rescaled_pumping=rescaled_pumping, ell0=10, L0=1, tau=1, show_sim=show_sim)
        # ca = 0.01
        if ca_1e_2 :
            plot_conf(ca=0.01, chis=chis, chiv=chiv, Nt=Nt_N1000, npts=npts, wmin=-6, wmax=6, rescaled=rescaled, rescaled_pumping=rescaled_pumping, ell0=10, L0=1, tau=1, show_sim=show_sim)
        # ca = 0.1
        if ca_1e_1 :
            plot_conf(ca=0.1, chis=chis, chiv=chiv, Nt=Nt_N1000, npts=npts, wmin=-6, wmax=5, rescaled=rescaled, rescaled_pumping=rescaled_pumping, ell0=10, L0=1, tau=1, show_sim=show_sim)
        # ca = 1
        if ca_1e0 :
            plot_conf(ca=1, chis=chis, chiv=chiv, Nt=Nt_N1000, npts=npts, wmin=-6, wmax=4, rescaled=rescaled, rescaled_pumping=rescaled_pumping, ell0=10, L0=1, tau=1, show_sim=show_sim)
        # ca = 10
        if ca_1e1 :
            plot_conf(ca=10, chis=chis, chiv=chiv, Nt=Nt_N1000, npts=npts, wmin=-6, wmax=4, rescaled=rescaled, rescaled_pumping=rescaled_pumping, ell0=10, L0=1, tau=1, show_sim=show_sim)

    if scaling_laws :
        coarsening_law = True
        merging_law = True

        if coarsening_law :
            if rescaled :
                k_c, alpha_c = 3000, -2./5
                t_c = np.logspace(0.5, 4, 101)
            else :
                k_c, alpha_c = 20., -2./5
                t_c = np.logspace(-3, 0, 101)
            plt.plot(t_c, k_c*t_c**alpha_c, linestyle='--', color = 'k', linewidth=2, label = r'$t^{-2/5}$')

        if merging_law :
            if rescaled :
                k_m, alpha_m = 200., -1.
                t_m = np.logspace(0, 2, 101)
            else :
                k_m, alpha_m = 50000., -1.
                t_m = np.logspace(1.7, 3.7, 101)
            plt.plot(t_m, k_m*t_m**alpha_m, linestyle='-.', color = 'k', linewidth=2, label = r'$t^{-1}$')

    plt.legend(fontsize=12, loc=3)
    if rescaled :
        plt.xlabel(r'$t / T_h$', fontsize=12)
    elif rescaled_pumping :
        plt.xlabel(r'$t/T_p$', fontsize=12)
    else :
        plt.xlabel(r'$t/\tau_v$', fontsize=12)

    plt.ylabel('N(t)', fontsize=12)
    plt.grid()
    plt.title('HO-chain with active pumping', fontsize=18)
    if savefig :
        plt.savefig(savename, format=extension)
    plt.show()

def plot_osmotic(chiv, chis, path_list, plot_param_list, Nt_list, rescale = False) :
    """
    plot_osmotic(chiv, chis, path_list, plot_param_list, Nt_list, rescale = False)
    
        Plot the number of lumens vs time for HO-chain without active pumping.
    """
    label = False
    Nt_list[(chis, chiv)] = {}
    for n in range(nsim) :
        Nt_list[(chis, chiv)][n] = np.loadtxt(os.path.join(path_list[(chis, chiv)][n], 'sim_nlum.dat'))
        if label :
            if not rescale :
                plt.plot(Nt_list[(chis, chiv)][n][:, 0], Nt_list[(chis, chiv)][n][:, 1], color=plot_param_list[(chis, chiv)][0], marker=plot_param_list[(chis, chiv)][1], markersize=plot_param_list[(chis, chiv)][2], linewidth=plot_param_list[(chis, chiv)][3], alpha=plot_param_list[(chis, chiv)][4])
            else :
                plt.plot(Nt_list[(chis, chiv)][n][:, 0]*chiv**2, Nt_list[(chis, chiv)][n][:, 1], color=plot_param_list[(chis, chiv)][0], marker=plot_param_list[(chis, chiv)][1], markersize=plot_param_list[(chis, chiv)][2], linewidth=plot_param_list[(chis, chiv)][3], alpha=plot_param_list[(chis, chiv)][4])
        else :
            label = True
            if not rescale :
                plt.plot(Nt_list[(chis, chiv)][n][:, 0], Nt_list[(chis, chiv)][n][:, 1], color=plot_param_list[(chis, chiv)][0], marker=plot_param_list[(chis, chiv)][1], markersize=plot_param_list[(chis, chiv)][2], linewidth=plot_param_list[(chis, chiv)][3], label = r'$\chi_s$ = '+ str(chis) + r' ; $\chi_v = $' + str(chiv), alpha=1.)
            else :
                plt.plot(Nt_list[(chis, chiv)][n][:, 0]*chiv**2, Nt_list[(chis, chiv)][n][:, 1], color=plot_param_list[(chis, chiv)][0], marker=plot_param_list[(chis, chiv)][1], markersize=plot_param_list[(chis, chiv)][2], linewidth=plot_param_list[(chis, chiv)][3], label = r'$\chi_s$ = '+ str(chis) + r' ; $\chi_v = $' + str(chiv), alpha=plot_param_list[(chis, chiv)][4])
    return Nt_list[(chis, chiv)]

def plot_osmotic_pumping(chiv, chis, path_list, plot_param_list, Nt_list, rescale = False) :
    """
    plot_osmotic_pumping(chiv, chis, path_list, plot_param_list, Nt_list, rescale = False)
    
        Plot the number of lumens vs time for HO-chain with active pumping.
    """
    label = False
    Nt_list[(chis, chiv)] = {}
    for n in range(nsim) :
        Nt_list[(chis, chiv)][n] = np.loadtxt(os.path.join(path_list[(chis, chiv)][n], 'sim_nlum.dat'))
        if label :
            if not rescale :
                plt.plot(Nt_list[(chis, chiv)][n][:, 0], Nt_list[(chis, chiv)][n][:, 1], color=plot_param_list[(chis, chiv)][0], marker=plot_param_list[(chis, chiv)][1], markersize=plot_param_list[(chis, chiv)][2], linewidth=plot_param_list[(chis, chiv)][3], alpha=plot_param_list[(chis, chiv)][4])
            else :
                plt.plot(Nt_list[(chis, chiv)][n][:, 0]*chiv**2, Nt_list[(chis, chiv)][n][:, 1], color=plot_param_list[(chis, chiv)][0], marker=plot_param_list[(chis, chiv)][1], markersize=plot_param_list[(chis, chiv)][2], linewidth=plot_param_list[(chis, chiv)][3], alpha=plot_param_list[(chis, chiv)][4])
        else :
            label = True
            if not rescale :
                plt.plot(Nt_list[(chis, chiv)][n][:, 0], Nt_list[(chis, chiv)][n][:, 1], color=plot_param_list[(chis, chiv)][0], marker=plot_param_list[(chis, chiv)][1], markersize=plot_param_list[(chis, chiv)][2], linewidth=plot_param_list[(chis, chiv)][3], label = r'$\chi_s$ = '+ str(chis) + r' ; $\chi_v = $' + str(chiv), alpha=1.)
            else :
                plt.plot(Nt_list[(chis, chiv)][n][:, 0]*chiv**2, Nt_list[(chis, chiv)][n][:, 1], color=plot_param_list[(chis, chiv)][0], marker=plot_param_list[(chis, chiv)][1], markersize=plot_param_list[(chis, chiv)][2], linewidth=plot_param_list[(chis, chiv)][3], label = r'$\chi_s$ = '+ str(chis) + r' ; $\chi_v = $' + str(chiv), alpha=plot_param_list[(chis, chiv)][4])
    return Nt_list[(chis, chiv)]

##