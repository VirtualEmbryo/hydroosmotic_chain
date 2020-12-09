'''
distribution.pyis a library of functions used to calculate distributions of size given folders containing simulations.

    Contains
    --------
    find_tstep       : Find the closest time from a given time in a time-list.
    calc_avg_distrib : Calculate the averages of dictionnaries.
    make_T_dict      : Make time lists from input simulation folders.
    make_step_dict   : Make dictionnary of selected times from simulation folders.
    calc_distrib     : Calculate the distribution at a given time of simulations in a folder.

    Requirements
    ------------
        Python libraries
    numpy (np)
    os

Mathieu Le Verge--Serandour, 2020
'''


import numpy as np
import os

def find_tstep(time_list, t0) :
    """
    find_tstep(time_list, t0)
    
        Find the closest time from t0 in time_list
        
        Parameters
        ----------
        time_list : array
            List of time
        t0 : float
            Time
        Returns
        -------
        t : float
            Clostest time from t0.        
    """
    return np.argmin(np.abs(time_list - t0))

def calc_avg_distrib(X_dict, Y_dict) :
    """
    calc_avg_distrib(X_dict, Y_dict)
    
        Parameters
        ----------
        X_dict, Y_dict : dict
            Dictionnaries containing the data to average.
        Returns
        -------
        X_avg, Y_avg : arrays
            Arrays containing the averages of the dictionnaries, averaged on axis=0.
        
    """
    X_array = np.array([X_dict[k] for k in X_dict.keys()])
    X_avg = np.average(X_array, axis=0)

    Y_array = np.array([Y_dict[k] for k in Y_dict.keys()])
    Y_avg = np.average(Y_array, axis=0)
    return X_avg, Y_avg

def make_T_dict(folder, npts, nsim=20, log_tmin=-3, log_tmax=7) :
    """
    make_T_dict(folder, npts, nsim=20, log_tmin=-3, log_tmax=7)
    
        Calculate the time-list distributed linearly in log-scale between 10^log_tmin and 10^log_tmax
        From simulations contained in folder, import the time of the simulation.
        
        Parameters
        ----------
        folder : path
            Folder containing the simulations, starting by 'run' and containing the file 'distrib_length.dat'
        npts : int
            Number of points of the time-list
        nsim : int, optional, default : 20
            Number of simulations in the folder
        log_tmin : float, optional, default : -3
            Min time ot time-list, such that t_min = 10^log_tmin
            Ex : log_tmin = -3 => t_min = 0.001
        log_tmax : float, optional, default : 7
            Max time ot time-list, such that t_max = 10^log_tmax
            Ex : log_tmax = 7 => t_max = 10000000
    
        Returns
        -------
        T_dict : dict
            Dictionnary of the times of the simulations, indexed by folder number
        time_plot_list : array
            Array of time-list points        
    """
    time_plot_list = np.logspace(log_tmin, log_tmax, npts)
    T_dict = {}
    for k in range(nsim) :
        T_dict[k] = np.loadtxt(os.path.join(folder, 'run'+str(k).zfill(4)+'/distrib_length.dat'), usecols=0)
    return T_dict, time_plot_list

def make_step_dict(T_dict, time_plot_list, npts, nsim) :
    """
    make_step_dict(T_dict, time_plot_list, npts, nsim)
    
        Calculate for each time-point in time_plot_list the closest time in T_dict for each simulation.
    
        Parameters
        ----------
        T_dict : dict
            Dictionnary containing the times of simulations.
        time_plot_list : array
            List of time-points
        npts : int
            Number of time points in time_plot_list
        nsim : int
            Number of simulation in T_dict
    
        Returns
        -------
        step_dict : dict
            Dictionnary containing the closest time in T_dict from given time in time_plot_list.
        
    """
    step_dict = {}
    for n in range(npts) :
        step_dict[n] = {}
        for k in range(nsim) :
            step_dict[n][k] = np.argmin(np.abs(T_dict[k] - time_plot_list[n]))
    return step_dict

def calc_distrib(time, folder, area = False) :
    """
    calc_distrib(time, folder, area = False)
    
        Calculate the distribution of simulations stored in the folder at a given time, using numpy.histogram
        The distribution is normalized such that its integral is 1.
        The binning is calculated as the max between 10 and 10+5*log(len(SIZES)) where SIZES is a list containing concatenated sizes of all simulations.
    
        Parameters
        ----------
        time : float
            Time-point at which the distribution will be calculated. 
            Due to the adaptative time-stepping of simulations, it is not the exact same time for each simulations
            but the closest time found.
        folder : path
            Folder containing the simulations, starting by 'run'
            Ex : run0000, run0001, run0002, etc.
        area : boolean, option, default : False
            If True, the distribution will be calculated as a function of area.
            If False, the distribution will be calculated as a function of length.
    
        Returns
        -------
        distrib = [x, y] : list of distribution points at time t
        x : non-rescaled size (length or area)
        y : normalized number of lumens with given size.
    """
    mu = 0.6105653703843762
    
    dat = {}
    for elem in os.listdir(folder) :
        if elem.startswith('run') :
            tdat = np.loadtxt(os.path.join(folder, elem, 'distrib_length.dat'), usecols=0)
            step = np.argmin(np.abs(tdat-time))
            Ldat = np.genfromtxt(os.path.join(folder, elem, 'distrib_length.dat'), skip_header=step, skip_footer=len(tdat)-step-1)
            dat[int(elem[-4:])] = [step, Ldat]

    # Concatenate all simulations at a given time point into a single list new_L
    new_L = np.concatenate([dat[k][1][1:] for k in dat.keys()])
    
    #Define the number of bins
    bins = np.max([10, 10+int(np.log10(len(new_L)))*5])
    
    if area :
        new_A = new_L**2 / mu
        weights = np.ones_like(new_A)/len(new_A)
        y, x = np.histogram(new_A, bins=bins)
    else :
        # If the histogram has abscissa L / L_*
        y, x = np.histogram(new_L, bins=bins, weights=np.ones_like(new_L)/len(new_L))
        
    x = 0.5*(x[1:]+x[:-1])
    distrib = [x, y]
    return distrib



