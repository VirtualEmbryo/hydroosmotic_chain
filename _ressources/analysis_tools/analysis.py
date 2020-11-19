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

# =================
def lin(x, a, b) :
    return a*x+b

def fit_coeff(func, x_dat, y_dat) :
    x_dat_log, y_dat_log = x_dat, y_dat
    popt, pcov = curve_fit(func, x_dat_log, y_dat_log)
    return popt#, pcov

def average_powerlaw(a_list, k_list) :
    a_avg = np.average(a_list)
    k_avg = np.exp(np.average(k_list))
    return a_avg, k_avg

def calc_mu(theta) :
    return np.sin(theta)**2 / (2*theta - np.sin(2*theta))

def calc_chi(theta, gamma, kappa, ell0, L0) :
    mu = calc_mu(theta)
    return 0.5*mu*np.sin(theta)*gamma*kappa / (ell0*L0**3)
    #return gamma*kappa / (ell0*L0**3)

def calc_chi(theta, eps, kappa, ell0, L0) :
    mu = calc_mu(theta)
    return 0.5*mu*np.sin(theta)*eps*kappa / (ell0*L0**3)
    #return mu*np.sin(theta)*eps / (L0*ell0**3)

def fit_lin(t, N) :
    x, y = np.log(t), np.log(N)
    popt, pcov = curve_fit(lin, x, y)
    alpha, kappa = popt[0], np.exp(popt[1])
    alpha_std, kappa_std = pcov[0, 0], np.exp(pcov[1, 1])
    return kappa, alpha#, kappa_std, alpha_std

def make_path_dict(nsim, chiv, chis, main_dir, subdir, subsubdir) :
    pathdict = {}
    for n in range(nsim) :
        #subsubdir = 'chiv' + str(chiv) + '_chis' + str(chis) + s
        pathdict[n] = os.path.join(main_dir, subdir, subsubdir, 'run'+str(n).zfill(4))
    return pathdict

def import_osmotic(chiv, chis, path_list, ca = None) :
    Nt = {}
    if ca == None :
        for n in path_list[(chiv, chis)].keys() :
            Nt[n] = np.loadtxt(os.path.join(path_list[(chiv, chis)][n], 'sim_nlum.dat'))
    else :
        for n in path_list[(chiv, chis)][ca].keys() :
            Nt[n] = np.loadtxt(os.path.join(path_list[(chiv, chis)][ca][n], 'sim_nlum.dat'))
    return Nt

def plot_osmotic(chiv, chis, path_list, plot_param_list, Nt_list, rescale = False) :
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

# Averaging

def gen_log_times(tmin, tmax, npts) :
    return np.logspace(np.log10(tmin), np.log10(tmax), npts)

def find_argmins(Nt, log_times) :
    index_array = np.zeros(len(log_times), dtype=int)
    for k in range(len(log_times)) :
        index_array[k] = np.argmin(np.abs(log_times[k]-Nt[:, 0]))
    return index_array

def gen_index_array(Nt, npts) :
    size = len(Nt)
    
    tmin = np.min(Nt[1:, 0])   # skip the first time since t=0
    tmax = np.max(Nt[1:, 0])   # skip the first time since t=0
    
    log_times = gen_log_times(tmin, tmax, npts)
    index_array = find_argmins(Nt, log_times)
    
    return index_array

# Show configuration
def plot_conf(ca, chis, chiv, Nt, npts=20,  wmin=-6, wmax=3, rescaled=False, rescaled_pumping=False, ell0=10, L0=12000, tau=1, show_sim=False, ax=None) :
    global chi_dict, ca_dict
    global mu, nu, eps
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
    
    