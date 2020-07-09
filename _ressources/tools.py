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
    end_time = time.ctime()
    f = open(os.path.join(outdir, 'log.txt'), 'a')
    
    f.write('status   : ' + str(end) + '\n')
    f.write('stop     : ' + str(end_time))
    
# ========================================================================
# ============================ Functions =================================
# ========================================================================

def autocorrFFT(x) :
    """
    autocorrFFT(x)

        Compute the autocorrelation function of the vector x, using the Wiener-Khinchin theorem.

        Parameters
        ----------
        x : array
            Input signal.
        Returns
        -------
        res : array
            Autocorrelation of the signal.
    """
    N = len(x)
    F = np.fft.fft(x, n = 2*N) # 2*N because zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res = (res[:N]).real
    n = N * np.ones(N) - np.arange(0, N) # generates array of  1/(N-m), m = 0, ..., N-1
    return res / n
    
def msd_fft1d(x):
    """
    msd_fft1d(x)

        Compute the Mean Square Displacement of a trajectory x, using Fast Fourier Transform and Autocorrelation function.

        Parameters
        ----------
        x : array
            Trajectory array, in 1 dimension
        Returns
        -------
        MSD : array
            MSD of the trajectory.
    """
    N=len(x)
    D=np.square(x)
    D=np.append(D,0)
    S2=autocorrFFT(x)
    Q=2*D.sum()
    S1=np.zeros(N)
    for m in range(N):
        Q=Q-D[m-1]-D[N-m]
        S1[m]=Q/(N-m)
    return S1-2*S2
    

# ========================================================================
# ========================== Trajectories ================================
# ========================================================================

def position_lumen(index, lumens) :
    return lumens[np.argwhere(lumens[:, 0]==index), 2]
    
def positions_list(index, lumens_list) :
    pos = []
    #print lumens_list[0][1]
    for t in range(len(lumens_list)) :
        if index in lumens_list[t][1][:, 0] :
            step = lumens_list[t][0]
            pos += [[step, position_lumen(index, lumens_list[t][1])[0, 0]]]
    return np.array(pos)

def plot_trajectories(lumens_list) :
    plt.figure(figsize=(8, 8))
    for i in range(1, int(np.max(lumens_list[:, 0]))) :
        plt.scatter(positions_list(i, lumens_list)[:, 0], positions_list(i, lumens_list)[:, 1], label = i, s=1)
    plt.xlabel('Time [a.u.]')
    plt.ylabel('Positions [a.u.]')
    plt.show()
    
# ========================================================================
# ============================ Profile ===================================
# ========================================================================

def calc_height_up(x, Li, thetai, h0, xci) :
    Ri = Li/np.sin(thetai)
    yci = h0 - np.sqrt(Ri**2 - Li**2)
    return yci + np.sqrt(Ri**2 - (x-xci)**2)

def calc_height_dw(x, Li, thetai, h0, xci) :
    Ri = Li/np.sin(thetai)
    yci = - h0 + np.sqrt(Ri**2 - Li**2)
    return yci - np.sqrt(Ri**2 - (x-xci)**2)

def profile(x, chain, theta=np.pi/3., h0=0.1) :
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

def plot_profile(x, chain, theta=np.pi/3., centers = True, axis = False, savefig = False, show=True, savename = 'pic.png', picformat='png', lw = 2, contour_color='k', center_color='r', xlim=[]) :
    #fig, ax = plt.subplots(1, 1)
    
    h_u, h_d = profile(x, chain, theta=theta, h0=chain.e0)
    
    #ax[1].suptitle('t = ' + "{:5.5f}".format(chain.time))
    number=int(savename[-11:-4])
    cste = -1e-2
    if number == 0 :
        plt.plot(x-number*cste, h_d-number*cste, linewidth = lw, color = contour_color)
        plt.plot(x-number*cste, h_u-number*cste, linewidth = lw, color = contour_color)
    else :
        plt.plot(x-number*cste, h_u-number*cste, linewidth = lw, color = contour_color)
    
    ### TO REMOVE
    #def gaussian_profile(x, amp, mu, sigma, threshold) :
    #    return amp*np.exp(-(x-mu)**2/sigma**2) + threshold
    #amp, mu, sigma, threshold = 1.0, 0.4, 0.05, 1.
    #y0 = 10
    #for k in chain.lumens_dict.keys() :
    #    if k != 0 and k !=-1 :
    #        xp = chain.lumens_dict[k].pos
    #        V = -gaussian_profile(xp/chain.total_length, amp, mu, sigma, threshold)
    #        ax.plot((xp, xp), (y0, y0+V*2.), color='k')
    #ax.plot(x, 20*gaussian_profile(x, amp, mu*chain.total_length, sigma*chain.total_length, threshold))
    ###
    
    #ax.plot(x, h_d, linewidth = lw, color = contour_color)
    #ax.plot(x, h_u, linewidth = lw, color = contour_color)
    ##plt.plot(x, h_d, linewidth = lw, color = contour_color)
    ##plt.plot(x, h_u, linewidth = lw, color = contour_color)

    if len(xlim) > 0 :
        xmin, xmax = xlim[0], xlim[1]
    else :
        xmin, xmax = np.min(x), np.max(x)
    xmin, xmax = chain.total_length*0.18, chain.total_length*0.82
    
    if centers :
        for k in list(chain.lumens_dict.keys()) :
            if k == 0 or k == -1 :
                ax.scatter(chain.lumens_dict[k].pos, 0, color = 'b')
            else :
                ax.scatter(chain.lumens_dict[k].pos, 0, color = center_color)
                
    #ax.vlines(x=0., ymin=-1., ymax=1.)
    #ax.vlines(x=chain.total_length, ymin=-1., ymax=1.)
    #ax.axis('equal')
    plt.axis('equal')
    plt.axis('off')
    #ax.set_xlim(xmin, xmax)
    #plt.xlim(xmin, xmax)
    
    ##if not axis :
    ##    ax.axis('off')
    #format = 'eps'
    #savefig = 1
    
    #print(number)
    #savename=savename[:-4] + '.eps'
    #if savefig and number == 8000 :
    
    #plt.suptitle('Time = '+"{:4.4e}".format(chain.time))
    #if savefig :
    if savefig :
        plt.savefig(savename, picformat=format)

    #if show : plt.show()
    #else : plt.close()

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

# ========================================================================
# ============================ Savings ===================================
# ========================================================================
def clear_folder(dir_name) :
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

# ========================================================================
# ======================== SAVE FUNCTIONS ================================
# ========================================================================
         
def save_events(chain, folder='', filename_events='events.txt') :
    events_file = open(os.path.join(folder, filename_events), 'a+')
    events_file.write(chain.events)
    events_file.close()
    chain.events = ''
    return ;
               
def save_recording(chain, filename='sim_all.dat', filename_bridges='sim_bridges.dat', folder='', chain_type='hydroosmotic', erase=True) :
    try :
        os.mkdir(folder)
    except : pass

    # Saving events
    #save_events(chain, folder=folder, filename_events='events.txt')###
    
    # Save qties
    file_all = open(os.path.join(folder, filename), 'a+')
    file_br = open(os.path.join(folder, filename_bridges), 'a+')
    
    time_list = np.sort(list(chain.rec.keys()))

    for t in time_list :
        s = ''
        # Write the lumens
        for n in chain.rec[t].keys() :
            if n != 0 and n != -1 :
                if 1 :
                    if chain_type == 'hydroosmotic' :
                        # Save : index, time, length, nb_ions, pos
                        s   += str(n) + '\t' + str(t) + '\t' + str(chain.rec[t][n][0]) + '\t' + str(chain.rec[t][n][1]) + '\t' + str(chain.rec[t][n][2])+ '\n'
                    elif chain_type == 'hydraulic' :
                        ### TO CHANGE
                        ###s   += str(n) + '\t' + str(t) + '\t' + str(chain.rec[t][n][0]) + '\t' + str(chain.rec[t][n][1]) + '\t' + str(chain.rec[t][n][2]) +'\n'
                        ### CORRECT : 
                        # Save : index, time, length, pos
                        s   += str(n) + '\t' + str(t) + '\t' + str(chain.rec[t][n][0]) + '\t' + str(chain.rec[t][n][1]) +'\n'
                else :
                    #print(chain.rec[t][n][0], chain.rec[t][n][1])
                    print(chain.rec[t].keys())
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
    #print(rec[time2][index][0] - rec[time1][index][0])
    return rec[t1][index][0] - rec[t0][index][0]
    
def find_winner(chain) :
    for k in chain.lumens_dict.keys() : 
        if k != 0 and k != -1 :
            chain.winner = k
    return chain.winner
    
def find_nmax(array, n)  :
    """Returns the n-th maximum numbers of a list"""
    return np.sort(np.partition(array, -n)[-n:])
    
def check_ending_event(chain) :
    
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
    f = open(eventfilename, 'a+')
    if end == 0 :
        # Error
        f.write('Simulation stopped before the end.')
    elif end == 1 :
        # Zero lumen left
        f.write('Time ' + "{:4.6f}".format(chain.time) + ' : winner (' + 'D' + ') is lumen ' + str(0))
    elif end == 2 :
        # One lumen left
        out = check_ending_event(chain)
        f.write(out)
        
    elif end == 10 :
        # Max step
        f.write('Simulation stopped [End code 10 : max step reached].')
    elif end == 11 :
        # Full size
        f.write('Simulation stopped [End code 11 : full system size reached].')
    return;

# ========================================================================
# ============================ Loading ===================================
# ========================================================================

def load_file(filename, skiprows=0, max_rows=0, hydroosmotic = True) :
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
    
# ========================================================================
# ============================= Plots ====================================
# ========================================================================

def plot_evolution(L, nions, ell, show_totalarea=False, savefig=False, savename='graph.eps', figsize=(7, 7), x_logscale=False, y_logscale=False, show_meanlength = True, title='', xlim=[], nbins=0) :
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

def chemograph(L, pos, x) :
    opening = np.ones((len(L), len(x)))
    for t in range(len(L)) :
        for i in range(len(x)) :
            for k in range(1, len(pos[t, 1:])+1) :
                if np.abs(x[i]-pos[t, k]) < L[t, k] :
                    opening[t, i] = 0
    return opening

# ========================================================================
# ============================= Plots ====================================
# ========================================================================

def get_winner(filename) :
    f = open(filename, 'r')
    w = int(f.readlines()[-1].split(' ')[-1])
    f.close()
    return w
    
