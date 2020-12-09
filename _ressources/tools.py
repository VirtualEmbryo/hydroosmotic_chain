"""
tools.py is a library containg generic tools used for the chain : graphics, functions, etc.

    Contains
    --------
        Log file
    write_log               : Write a log.txt file for the simulation
    add_end_time            : Add ending time in log.txt
        Profile
    calc_height_up          : Calculate upper membrane profile of the chain.
    calc_height_dw          : Calculate down membrane profile of the chain.
    profile                 : Calculate the profile/height at point x.
    plot_profile            : Plot the profile.
    plot_profile2           : Another version of plot_profile.
    plot_evolution          : Plot the evolution of the chain (area, concentration, nb ions, bridge lengths)
    plot_evolution_hydraulic: Plot the evolution of the hydraulic chain (area, bridge lengths)
        Save functions
    clear_folder            : Clean the given folder.
    save_events             : Save the events in an event-file (events.txt).
    save_recording          : Save the states of the chain.
        Events
    check_slope             : Calculate the slope of length L(t) between t0, t1.
    find_winner             : Find the winner of a simulation in event-file.
    find_nmax               : Calculate n-th maxima of a list.
    check_ending_event      : Classification of last event.
    write_ending_event      : Write the last event in event-file.
    get_winner              : Get the type of last event of the simulation.
        Loading
    load_file               : Load lumens data in dictionnaries.
    load_brfile             : Load bridges data in dictionnaries.
    calc_A_tot              : Calculate the total area through time.
    
    Requirements
    ------------
        Python libraries
    os
    sys
    numpy (np)
    time
    getpass
    socket
    matplotlib.pyplot (for graphical functions only)s

        Homemade libraries
    functions

Mathieu Le Verge--Serandour, 2020
"""


import numpy as np
import os
import sys

import time
import getpass
import socket

module_path = os.path.abspath(os.path.join('..', 'chain_lumen/'))
if module_path not in sys.path :
    sys.path.append(module_path)

try :
    import _ressources.functions as functions
except :
    import functions
    
try :
    import matplotlib.pyplot as plt
    
except : 
    pass

# ========================================================================
# ============================= Log file =================================
# ========================================================================

def write_log(outdir, confname, args) :
    """
    write_log(outdir, confname, args)
    
        Write log.txt file containing informations about the simulation.
    
        Parameters
        ----------
        outdir : string, directory
            Directory of the simulation
        confname : string
            Location of the configuration file
        args : list
            Optional arguments of the simulation given in the chain.main.py function (-v, -r, etc.)
    
    """
    username = getpass.getuser()
    hostname = socket.gethostname()
    start_time = time.ctime()
    end_time = time.ctime()
    
    f = open(os.path.join(outdir, 'log.txt'), 'w')
    
    f.write('username : ' + str(username) + '\n')
    f.write('hostname : ' + str(hostname) + '\n')
    f.write('wdir : ' + str(outdir) + '\n')
    f.write('\n')
    f.write('config   : '  + str(os.path.join(outdir, confname)) + '\n')
    if len(args) > 0 :
        f.write('args     : ' + str(args) + '\n')
    else :
        f.write('args     : [-]' + '\n')
    f.write('\n')
    f.write('start    : ' + str(start_time) + '\n')
    
def add_end_time(outdir, end) :
    """
    add_end_time(outdir, end)
    
        Write the (real) end time of the simulation in log.txt
    """
    end_time = time.ctime()
    f = open(os.path.join(outdir, 'log.txt'), 'a')
    
    f.write('status   : ' + str(end) + '\n')
    f.write('stop     : ' + str(end_time))
    
# ========================================================================
# ============================ Profile ===================================
# ========================================================================

def calc_height_up(x, Li, thetai, h0, xci) :
    """
    calc_height_up(x, Li, thetai, h0, xci)
        
        Calculate the upper profile h(x) at position x of the lumen centered at xci
        
        Parameters
        ----------
        x : float
            Position where to calculate h(x)
        Li : float
            Length of lumen i
        thetai : float
            Contact angle of lumen i
        h0 : float
            Diameter of the bridge
        xci : float
            Center of lumen i
        
        Returns
        -------
        h_x : float
            Upper position of the membrane at x, given that its within the lumen i.
    """
    Ri = Li/np.sin(thetai)
    yci = h0 - np.sqrt(Ri**2 - Li**2)
    h_x = yci + np.sqrt(Ri**2 - (x-xci)**2)
    return h_x

def calc_height_dw(x, Li, thetai, h0, xci) :
    """
    calc_height_up(x, Li, thetai, h0, xci)
        
        Calculate the bottom profile h(x) at position x of the lumen centered at xci
        
        Parameters
        ----------
        x : float
            Position where to calculate h(x)
        Li : float
            Length of lumen i
        thetai : float
            Contact angle of lumen i
        h0 : float
            Diameter of the bridge
        xci : float
            Center of lumen i
        
        Returns
        -------
        h_x : float
            Bottom position of the membrane at x, given that its within the lumen i.
    """
    Ri = Li/np.sin(thetai)
    yci = - h0 + np.sqrt(Ri**2 - Li**2)
    h_x = yci - np.sqrt(Ri**2 - (x-xci)**2)
    return h_x

def profile(x, chain, theta=np.pi/3., h0=0.1) :
    """
    profile(x, chain, theta=np.pi/3., h0=0.1)
    
        Calculates the profile h(x) of the chain membrane (up and down).
    
        Parameters
        ----------
        x : float or array
            Array of positions or position
        chain : chain-object
            
        theta : float, optional, default : np.pi/3.
            Contact angle (assumed to be identical)
        h0 : float, optional, default : 0.1
            Diameter of the bridges of the chain.
    """
    L_list = [chain.lumens_dict[k].length for k in list(chain.lumens_dict.keys())]
    pos_x = [chain.lumens_dict[k].pos for k in list(chain.lumens_dict.keys())]
    h_u = np.ones(len(x))*h0
    h_d = -np.ones(len(x))*h0
    for i in range(len(x)) :
        for k in range(len(pos_x)) :
            if np.abs(x[i]-pos_x[k]) < L_list[k] :
                h_u[i] = calc_height_up(x[i], L_list[k], theta, h0, pos_x[k])
                h_d[i] = calc_height_dw(x[i], L_list[k], theta, h0, pos_x[k])
    return h_u, h_d

def plot_profile(x, chain, theta=np.pi/3., centers = True, axis = False, savefig = False, show=False, savename = 'pic.png', picformat='png', lw = 1, contour_color='k', center_color='r', xlim=[]) :
    """
    plot_profile(x, chain, ...)
        
        Plot the profile of the chain.
    
        Parameters
        ----------
        x : array
            Positions of the chain.
        chain : chain-object
        
        theta : float, optional, default : np.pi/3.
            Contact angle
        centers : boolean, optionalm default : True
            Plot the centers of the lumens
        axis : boolean, optional, default : False
            Show x and y axis
        savefig : boolean, optional, default : False
            If True, save the figure under savename
        show : boolean, optional, default : False
            If True, show the chain.
        savename : string, optional, default : 'pic.png'
            Name of the output file/picture if saved.
        picformat : string, optional, default : 'png'
            Format of the output picture
        lw : float, optional, default : 1
            Line width on the figure.
        contour_color : string, optional, default : 'k'
            Color of the membrane on the picture. 'k' is black.
        center_color : string, optional, defautl : 'r'
            Color of the centers if plotted.
        xlim : list, optional, default : []
            Limits of the x-axis.    
    """
    e0 = 0.01
    h_u, h_d = profile(x, chain, theta=theta, h0=e0)
    lw=1

    plt.plot(x, h_d, linewidth = lw, color = contour_color)
    plt.plot(x, h_u, linewidth = lw, color = contour_color)

    if len(xlim) > 0 :
        xmin, xmax = xlim[0], xlim[1]
    else :
        xmin, xmax = np.min(x), np.max(x)
    xmin, xmax = chain.total_length*0.18, chain.total_length*0.82
    
    if centers :
        for k in list(chain.lumens_dict.keys()) :
            if k == 0 or k == -1 :
                plt.scatter(chain.lumens_dict[k].pos, 0, color = 'b')
            else :
                plt.scatter(chain.lumens_dict[k].pos, 0, color = center_color)
                
    plt.axis('equal')
    plt.axis('off')
    
    plt.suptitle(r'$\bar{t}$ = '+"{:4.4e}".format(chain.time), fontsize=20)
    
    if savefig :
        plt.savefig(savename[:-4]+'.png', format='png', dpi=200)

    if show : plt.show()
    else : plt.close()

def plot_profile2(x, chain, theta=np.pi/3., centers = True, axis = True, savefig = False, show=True, savename = 'pic.png', format='png', lw = 2, contour_color='k', center_color='r') :
    fig, ax = plt.subplots(1, 1, figsize=(6, 8))
    
    try :
        args = functions.get_parameters(chain.pumping_args, chain.pumping_func)
        fpump = functions.calc_func(x, chain.pumping_func, args, chain.total_length*2)
        ax[0].plot(x, fpump)
        ax[0].plot(pos, ca, marker='o', linewidth=0)
        ax[0].set_ylim(-np.max(fpump)*5e-2, np.max(fpump)*1.1)
    except : pass
    
    #number = int(savename[-11:-4])
    number=8000
    h_u, h_d = profile(x, chain, theta=theta, h0=chain.e0)
    ###ax[1].suptitle('t = ' + "{:5.5f}".format(chain.time))
    
    plt.plot(x-2e-3*number, h_u-2e-3*number, linewidth = lw, color = contour_color)
    ax[1].plot(x, h_d, linewidth = lw, color = contour_color)
    ax[1].plot(x, h_u, linewidth = lw, color = contour_color)

    xmin, xmax = np.min(x), np.max(x)
    ax[1].set_xlim(xmin, xmax)
    ax[0].set_xlim(xmin, xmax)
    
    if centers :
        for k in list(chain.lumens_dict.keys()) :
            if k == 0 or k == -1 :
                ax[1].scatter(chain.lumens_dict[k].pos, 0, color = 'b')
            else :
                ax[1].scatter(chain.lumens_dict[k].pos, 0, color = center_color)
    #ax[1].axis('equal')
    
    if not axis :
        ax[1].axis('off')
    #format = 'eps'
    #savefig = 1
    
    #print(number)
    #savename=savename[:-4] + '.eps'
    #if savefig and number == 8000 :
    
    if savefig :
        plt.savefig(savename, format=format)

    if show : plt.show()
    else : plt.close()

def plot_evolution(L, nions, ell, show_totalarea=False, savefig=False, savename='graph.eps', figsize=(7, 7), x_logscale=False, y_logscale=False, show_meanlength = True, title='', xlim=[], nbins=0) :
    """
    plot_evolution(L, nions, ell)
    
        Plot the evolution of the hydroosmotic chain in brige length, number of ions, area, and concentration.
    
        Parameters
        ----------
        L : array
            Array of lengths, with first column as time, the others being lengths of lumen sorted by index.
        nions : array
            Array of number of ions, with first column as time, the others being lengths of lumen sorted by index.
        ell : array
            Array of bridge lengths, with first column as time, the others being lengths of bridges sorted by index.
        show_totalarea : boolean, optional, default : False
             Plot the total area.
        savefig : boolean, optional, default : False
             If True, save the figure
        savename : string, optional, default : 'graph.eps'
            Name of the figure if saved
        figsize : tuple, optional, default : (7, 7)
            Figure size
        x_logscale : boolean, optional, default : False
            If True, the x-scale is logarithmic
        y_logscale : boolean, optional, default : False
            If True, the y-scale is logarithmic
        show_meanlength : boolean, optional, default : True
            If True, plot the mean bridge length
        title : string, optional, default : ''
            Title of the figure
        xlim : list, optional, default : []
            Limit in x-axis. Must be a 1d list with 2 entries.
        nbins : int, optional, default : 0
            Number of bins of the graph. If 0, plyplot calculates the best number of bins.
    
    
        NB : for the area calculation, geometrical coefficient mu = mu(pi/3) = 0.611
    """
    fig, ax = plt.subplots(2, 2, figsize=figsize)
    if len(title) > 0 :
        plt.suptitle(title, fontsize=20)
    tmin, tmax = 0., 0.4

    if x_logscale :
        ax[0,0].set_xscale('log')
        ax[1,0].set_xscale('log')
        ax[0,1].set_xscale('log')
        ax[1,1].set_xscale('log')

    if y_logscale :
        ax[0,0].set_yscale('log')
        ax[1,0].set_yscale('log')
        ax[0,1].set_yscale('log')
        ax[1,1].set_yscale('log')

    # INTERLUMENAL LENGTHS
    ax[0, 0].set_title(r'$\ell_{ij}$', fontsize=15)
    for k in range(1, len(ell[0])) :
        ax[0, 0].plot(ell[:, 0], ell[:, k], linestyle='-', linewidth=2, label = str(k))
    mean = np.nanmean(ell[:, 1:], axis=1)
    ax[0, 0].plot(ell[1:-1, 0], mean[1:-1], color = 'k', linestyle = '--')
    ax[0, 0].grid()
    
    # Nions
    ax[0, 1].set_title(r'$N_{ions}$', fontsize=15)
    for k in range(1, len(nions[0])) :
        ax[0, 1].plot(nions[:, 0], nions[:, k], linewidth=2)
    
    ax[0, 1].grid()

    mu = 0.6105653703843762
    
    # AREAS
    ax[1, 0].set_title('Area', fontsize=15)
    for k in range(1, len(L[0])) :
        ax[1, 0].plot(L[:, 0], L[:, k]**2 / mu)
    
    if show_totalarea :
        t_a, A_tot = calc_A_tot(L)
        ax[1, 0].plot(t_a, A_tot, linestyle='--', linewidth=2, color='k')

    ax[1, 0].grid()
    ax[1, 0].set_xlabel('Time [s]')

    # Concentration
    ax[1, 1].set_title('Concentrations', fontsize=15)
    for k in range(1, len(nions[0])) :
        ax[1, 1].plot(nions[:, 0], nions[:, k]*mu / L[:, k]**2, linewidth=2)

    ax[1, 1].grid()
    ax[1, 1].set_xlabel('Time [s]')
    
    
    if len(xlim) > 0 :
        ax[0, 0].set_xlim(xlim[0], xlim[1])
        ax[0, 1].set_xlim(xlim[0], xlim[1])
        ax[1, 0].set_xlim(xlim[0], xlim[1])
        ax[1, 1].set_xlim(xlim[0], xlim[1])
        
    if nbins > 0 :
        ax[0, 0].locator_params(nbins=nbins)
        ax[0, 1].locator_params(nbins=nbins)
        ax[1, 0].locator_params(nbins=nbins)
        ax[1, 1].locator_params(nbins=nbins)
    
    if savefig :
        plt.savefig(savename, format='eps')
        
    
    plt.show()
    
def plot_evolution_hydraulic(L, ell, show_totalarea=False, savefig=False, savename='graph.eps', figsize=(7, 4), x_logscale=False, y_logscale=False, show_meanlength = True, xlim=[], nbins=0) :
    '''
    plot_evolution_hydraulic(L, ell)
    
        Plot the evolution of the hydraulic chain in brige length, number of ions, area, and concentration.
    
        Parameters
        ----------
        L : array
            Array of lengths, with first column as time, the others being lengths of lumen sorted by index.
        ell : array
            Array of bridge lengths, with first column as time, the others being lengths of bridges sorted by index.
        show_totalarea : boolean, optional, default : False
             Plot the total area.
        savefig : boolean, optional, default : False
             If True, save the figure
        savename : string, optional, default : 'graph.eps'
            Name of the figure if saved
        figsize : tuple, optional, default : (7, 7)
            Figure size
        x_logscale : boolean, optional, default : False
            If True, the x-scale is logarithmic
        y_logscale : boolean, optional, default : False
            If True, the y-scale is logarithmic
        show_meanlength : boolean, optional, default : True
            If True, plot the mean bridge length
        title : string, optional, default : ''
            Title of the figure
        xlim : list, optional, default : []
            Limit in x-axis. Must be a 1d list with 2 entries.
        nbins : int, optional, default : 0
            Number of bins of the graph. If 0, plyplot calculates the best number of bins.
    
        NB : for the area calculation, geometrical coefficient mu = mu(pi/3) = 0.611
    '''
    #fig, ax = plt.subplots(1, 2, figsize=figsize)
    fig, ax = plt.subplots(1, 2, figsize=figsize)

    tmin, tmax = 0., 0.4

    if x_logscale :
        ax[0].set_xscale('log')
        ax[1].set_xscale('log')

    if y_logscale :
        ax[0].set_yscale('log')
        ax[1].set_yscale('log')

    # INTERLUMENAL LENGTHS
    ax[0].set_title(r'$\ell_{ij}$', fontsize=15)
    for k in range(1, len(ell[0])) :
        ax[0].plot(ell[:, 0], ell[:, k], linestyle='-', linewidth=2, label = str(k))
    mean = np.nanmean(ell[:, 1:], axis=1)
    ax[0].plot(ell[1:-1, 0], mean[1:-1], color = 'k', linestyle='--')
    ax[0].grid()
    ax[0].set_xlabel('Time [s]')
    if len(xlim) > 0 :
        ax[0].set_xlim(xlim[0], xlim[1])
        
    mu = 0.6105653703843762
    # AREAS
    ax[1].set_title('Area', fontsize=15)
    for k in range(1, len(L[0])) :
        ax[1].plot(L[:, 0], L[:, k]**2 / mu)
    
    if show_totalarea :
        t_a, A_tot = calc_A_tot(L)
        ax[1].plot(t_a, A_tot, linestyle='--', linewidth=2, color='k')

    ax[1].grid()
    ax[1].set_xlabel('Time [s]')
    if len(xlim) > 0 :
        ax[1].set_xlim(xlim[0], xlim[1])        
    
    if nbins > 0 :
        ax[0].locator_params(nbins=nbins)
        ax[1].locator_params(nbins=nbins)

    if savefig :
        plt.savefig(savename, format='eps')
    plt.show()

# ========================================================================
# ======================== SAVE FUNCTIONS ================================
# ========================================================================
def clear_folder(dir_name) :
    """
    clear_folder(dir_name)
    
        Remove all existing files in dir_name
    
        Parameters
        ----------
        dir_name : string
            Directory to clean
    
        Returns
        -------
        0 : Cleaning is done
        1 : Error : Directory not found
        2 : Error : Directory not cleaned (deletion of a file not allowed)
    """
    try :
        if len(os.listdir(dir_name)) > 0 :
            for elem in os.listdir(dir_name) :
                if os.path.isdir(os.path.join(dir_name, elem)) :
                    for sub_dir_elem in os.listdir(os.path.join(dir_name, elem)) :
                        os.remove(os.path.join(dir_name, elem, sub_dir_elem))
                    os.removedirs(os.path.join(dir_name, elem))
                else : 
                    os.remove(os.path.join(dir_name, elem))
        os.removedirs(dir_name)
        return 0
    
    except :
        if not os.path.isdir(dir_name) :
            print('Error : file not found')
            return 1
        else :
            print('Error : directory ' + dir_name + ' not cleaned...')
            return 2
                     
def save_events(chain, folder='', filename_events='events.txt') :
    """
    save_events(chain, folder='', filename_events='events.txt')
    
        Save the events of the chain (stored in chain.events) in an event-file.
        
        Parameters
        ----------
        chain : chain-object
        folder : string, optional, default : ''
            Folder of the event-file
        filename_events : string, optional, default : events.txt
            Name of the event-file
    """
    events_file = open(os.path.join(folder, filename_events), 'a+')
    events_file.write(chain.events)
    events_file.close()
    chain.events = ''
    return ;
               
def save_recording(chain, filename='sim_all.dat', filename_bridges='sim_bridges.dat', folder='', chain_type='hydroosmotic', erase=True) :
    """
    save_recording(chain)
    
        Save the states of the chain
    
        Parameters
        ----------
        chain : chain-object
            Chain of the simulation
        filename : string, optional, default : 'sim_all.dat'
            File containing all data about the lumens (index, position, length)
        filename_bridges : string, optionalm default : 'sim_bridges.dat'
            File containing all data about the bridges (index, lumen1, lumen2, length)
        folder : string, optional, default : ''
            Directory where to write the files
        chain_type : string, optional, default : 'hydroosmotic'
            Type of chain
        erase : boolean, optional, default : True
            Delete entries of chain.rec, that stores all data.
    """
    try :
        os.mkdir(folder)
    except : pass
    
    # Save qties
    file_all = open(os.path.join(folder, filename), 'a+')
    file_br = open(os.path.join(folder, filename_bridges), 'a+')
    
    time_list = np.sort(list(chain.rec.keys()))

    for t in time_list :
        s = ''
        # Write the lumens
        for n in chain.rec[t].keys() :
            if n != 0 and n != -1 :
                try :
                    if chain_type == 'hydroosmotic' :
                        # Save : index, time, length, nb_ions, pos
                        s   += str(n) + '\t' + str(t) + '\t' + str(chain.rec[t][n][0]) + '\t' + str(chain.rec[t][n][1]) + '\t' + str(chain.rec[t][n][2])+ '\n'
                    elif chain_type == 'hydraulic' :
                        # Save : index, time, length, pos
                        s   += str(n) + '\t' + str(t) + '\t' + str(chain.rec[t][n][0]) + '\t' + str(chain.rec[t][n][1]) +'\n'
                except :
                    #print(chain.rec[t][n][0], chain.rec[t][n][1])
                    print(chain.rec[t].keys())
                    pass
                    #print(chain)
        
        file_all.write(s)
        
        # Write the bridges
        sb = ''
        
        for b in chain.rec_br[t].keys() :
            # Save : index, time, length, lumen1, lumen2
            sb += str(b) + '\t' + str(t) + '\t' + str(chain.rec_br[t][b][0]) + '\t'+ str(chain.rec_br[t][b][1]) + '\t'+ str(chain.rec_br[t][b][2]) +'\n'
            
        file_br.write(sb)    
        
    file_all.close()
    file_br.close()
    
    if erase :
        chain.rec = {}
        chain.rec_br = {}
    
    chain.events = ''

# ========================================================================
# ============================ EVENTS ====================================
# ========================================================================

def check_slope(index, t0, t1, rec) :
    """
    check_slope(index, t0, t1, rec)
        
        Calculate the slope of the length evolution of lumen index, given times t0 and t1
        
        Parameters
        ----------
        index : int
            Index of 
        t0 : float
            Time index
        t1 : float
            Time index
        rec : dict
            Dictionnary of recordings
    
        Returns
        -------
        slope : float
            Slope of the length L(t)
    """
    slope = rec[t1][index][0] - rec[t0][index][0] 
    return slope
    
def find_winner(chain) :
    """
    find_winner(chain)
    
        Find the index of the winner of the chain.
    """
    for k in chain.lumens_dict.keys() : 
        if k != 0 and k != -1 :
            chain.winner = k
    return chain.winner
    
def find_nmax(array, n) :
    """
    find_nmax(array, n)
    
        Returns the n-th maximum numbers of a list
        Example :
        >>> findnmax([0, 1, 2, 3, 4, 5, 6, 6, 7], 3)
        array([7, 6, 6])
        
        Parameters
        ----------
        array : array
        n : int
            n-th max
    """
    return np.sort(np.partition(array, -n)[-n:])
    
def check_ending_event(chain) :
    """
    check_ending_event(chain)
    
        Returns the last event for classification :
            Coarsening (C) : if the slope of the winner is positive
            Collapse (D)   : if the slope of the winner is negative
            Merge (M)      : if the index of the winning lumen is not in the chain, it was created by merging
    
        Parameters
        ----------
        chain : chain-object
        
        Returns
        -------
        ev : string
            Contains the event to write in events.txt
    """
    if len(chain.lumens_dict) - 2 == 0 :
        return 'Time ' + "{:4.6f}".format(chain.time) + ' : winner (' + 'D' + ') is lumen ' + str(0)
        
    else :
        w = find_winner(chain)
        tend_array = find_nmax(list(chain.rec.keys()), 3)           # finds the 3 maximum times of rec dict
        try :
            s = check_slope(w, tend_array[0], tend_array[1], chain.rec)
            
            if s > 0. :
                # Coarsening
                ev = 'Time ' + "{:4.6f}".format(chain.time) + ' : winner (C) is lumen ' + str(w)
                return ev
            elif s < 0. :
                # Disparition
                ev = 'Time ' + "{:4.6f}".format(chain.time) + ' : winner (D) is lumen ' + str(0)
                return ev
        except :
            ev = 'Time ' + "{:4.6f}".format(chain.time) + ' : winner (M) is lumen ' + str(w)
            return ev
            
def write_ending_event(end, chain, eventfilename) :
    """
    write_ending_event(end, chain, eventfilename)
    
        Write the last event in the events.txt file
        
        Parameters
        ----------
        end : int
            Type of ending events
        chain : chain-object
    
        eventfilename : float
            Name of the event filename. Usually events.txt    
        
    """
    f = open(eventfilename, 'a+')
    f.write('\n')
    print(end)
    if end == 0 :
        # Error
        f.write('Simulation stopped before the end.')
    elif end == 1 :
        print(1)
        # One lumen left
        out = check_ending_event(chain)
        f.write(out)
    elif end == 2 :
        # Zero lumen left
        f.write('Time ' + "{:4.6f}".format(chain.time) + ' : winner (' + 'D' + ') is lumen ' + str(0))
        
    elif end == 10 :
        # Max step
        f.write('Simulation stopped [End code 10 : max step reached].')
    elif end == 11 :
        # Full size
        f.write('Simulation stopped [End code 11 : full system size reached].')
    return;
    
def get_winner(filename) :
    """
    get_winner(filename)
    
        Get the winner of the chain from event-file.
    
        Parameters
        ----------
        filename : string
            Name of the event-file
        
        Returns
        -------
        w : int
            Index of the winner from events.txt
    """
    f = open(filename, 'r')
    w = int(f.readlines()[-1].split(' ')[-1])
    f.close()
    return w

# ========================================================================
# ============================ Loading ===================================
# ========================================================================

def load_file(filename, skiprows=0, max_rows=0, hydroosmotic = True) :
    """
    load_file(filename, skiprows=0, max_rows=0)
    
        Load the lumen-file and returns it as a list (to plot for instance)
        
        
        Parameters
        ----------
        filename : string
            Name of the lumen-file
        skiprows : int, optional, default : 0
            Skip the rows above it.
        max : int, optional, default :
            Skip the rows after it
        hydroosmotic : boolean, optional, default : True
            Type of the chain
    
        Returns
        -------
        L_array : dict
            Dictionnary of lengths, first indexed by time, then by lumen-index
        p_array : dict
            Dictionnary of positions, first indexed by time, then by lumen-index
        N_array : dict, if hydroosmotic
            Dictionnary of number of ions, first indexed by time, then by lumen-index
    """
    time = []
    if max_rows != 0 :
        dat = np.loadtxt(filename, skiprows=skiprows, max_rows=max_rows)
        print('Import successful !')
    else : 
        dat = np.loadtxt(filename, skiprows=skiprows)
        print('Import successful !')
    
    Nmax = int(np.max(dat[:, 0]))

    L_a = {}
    p_a = {}
    #c_a = {} ### TO CHANGE
    
    if hydroosmotic : N_a = {}
    
    for i in range(len(dat)) :
        t = dat[i, 1]
        if t not in L_a.keys() :
            L_a[t] = {}
            p_a[t] = {}
            #c_a[t] = {} ### TO CHANGE
            if hydroosmotic : N_a[t] = {}
            
        
        L_a[t][int(dat[i, 0])] = dat[i, 2]
        if hydroosmotic : 
            N_a[t][int(dat[i, 0])] = dat[i, 3]
            p_a[t][int(dat[i, 0])] = dat[i, 4]
        else :
            p_a[t][int(dat[i, 0])] = dat[i, 3]
            #c_a[t][int(dat[i, 0])] = dat[i, 4] ### TO CHANGE

    L_array = np.zeros(( len(L_a.keys()), Nmax+1 ))
    p_array = np.zeros(( len(p_a.keys()), Nmax+1 ))
    #c_array = np.zeros(( len(c_a.keys()), Nmax+1 )) ### TO CHANGE
    if hydroosmotic : N_array = np.zeros(( len(N_a.keys()), Nmax+1 ))

    step = -1
    for k in L_a.keys() :
        step += 1
        L_temp = [k]
        p_temp = [k]
        #c_temp = [k] ### TO CHANGE
        if hydroosmotic : N_temp = [k]
        for j in range(1, Nmax+1) :
            if j in L_a[k].keys() :
                L_temp += [L_a[k][j]]
                p_temp += [p_a[k][j]]
                #c_temp += [c_a[k][j]] ### TO CHANGE
                if hydroosmotic : N_temp += [N_a[k][j]]
            else :
                L_temp += [None]
                p_temp += [None]
                #c_temp += [None] ### TO CHANGE
                if hydroosmotic : N_temp += [None]
        L_array[step] = np.array(L_temp)
        if hydroosmotic : N_array[step] = np.array(N_temp)
        p_array[step] = np.array(p_temp)
        #c_array[step] = np.array(c_temp) ### TO CHANGE
    
    if hydroosmotic : 
        return L_array, N_array, p_array
    else :
        return L_array, p_array#, c_array ### TO CHANGE
    
def load_brfile(filename, skiprows=0, max_rows=0) :
    """
    load_brfile(filename, skiprows=0, max_rows=0)
    
        Load the bridge-file and returns it as a list (to plot for instance)
        
        
        Parameters
        ----------
        filename : string
            Name of the bridge-file
        skiprows : int, optional, default : 0
            Skip the rows above it.
        max : int, optional, default :
            Skip the rows after it
    
        Returns
        -------
        ell_a : list
            List of all bridges at a given time
    """
    time = []
    if max_rows != 0 :
        dat = np.loadtxt(filename, skiprows=skiprows, max_rows=max_rows)
    else : 
        dat = np.loadtxt(filename, skiprows=skiprows)
    
    Nmax = int(np.max(dat[:, 0]))
    ell_a = {}
    
    for i in range(len(dat)) :
        t = dat[i, 1]
        if t not in ell_a.keys() :
            ell_a[t] = {}
        
        #if dat[i, 3] != 0 and dat[i, 3] != -1 and dat[i, 4] != 0 dat[i, 4] != -1 :
        #    ell_a[t][int(dat[i, 0])] = dat[i, 2]
        ell_a[t][int(dat[i, 0])] = dat[i, 2]

    ell_array = np.zeros(( len(ell_a.keys()), Nmax+2 ))
    
    step = -1
    for k in ell_a.keys() :
        step += 1
        ell_temp = [k]
        
        for j in range(0, Nmax+1) :
            if j in ell_a[k].keys() :
                ell_temp += [ell_a[k][j]]
                
            else :
                ell_temp += [None]
                
        ell_array[step] = np.array(ell_temp)

    return ell_array
    
def calc_A_tot(L) :
    """
    calc_A_tot(L)
    
        Calculate the list of total area from list of lengths L, such that area is given by A = L^2 / mu
        NB : mu is fixed at mu(pi/3) = 0.6105653703843762
    
        Parameters
        ----------
        L : list
            List of lengths and times,such that 
            L = [[t1, [L1(t1), L2(t1), L3(t1), ...]], [t2, [L1(t2), ...]], ... ]
        Returns
        -------
        time : list
            List of times
        A_tot : list
            List of total area at time t
    """
    mu = 0.6105653703843762
    A_tot = []
    time = []
    for i in range(len(L)) :
        a = 0
        time += [L[i, 0]]
        for k in range(1, len(L[i])) :
            if not np.isnan(L[i][k]) :
                a += L[i, k]**2 / mu
        A_tot += [a]
    return time, A_tot
    
#
