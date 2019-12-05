import numpy as np
import os
try :
    import matplotlib.pyplot as plt
    
except : 
    pass
    


#import network as net
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

def plot_profile(x, chain, theta=np.pi/3., centers = True, axis = True, savefig = False, show=True, savename = 'pic.png', lw = 2, contour_color='k', center_color='r') :
    
    h_u, h_d = profile(x, chain, theta=theta, h0=chain.e0)
    
    plt.plot(x, h_d, linewidth = lw, color = contour_color)
    plt.plot(x, h_u, linewidth = lw, color = contour_color)
    
    if centers :
        for k in list(chain.lumens_dict.keys()) :
            if k == 0 or k == -1 :
                plt.scatter(chain.lumens_dict[k].pos, 0, color = 'b')
            else :
                plt.scatter(chain.lumens_dict[k].pos, 0, color = center_color)
    plt.axis('equal')
    
    if not axis :
        plt.axis('off')
    
    if savefig :
        plt.savefig(savename)

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
            
def save_recording(chain, filename='sim.dat', filename_events='events.log', folder='') :
    try :
        os.mkdir(folder)
    except : pass

    # Saving events
    events_file = open(os.path.join(folder, filename_events), 'a+')
    events_file.write(chain.events)
    events_file.close()
    
    # Save qties
    file_all = open(os.path.join(folder, filename), 'a+')
    time_list = np.sort(list(chain.rec.keys()))
    for t in time_list :
        s = ''
        for n in chain.rec[t].keys() :
            if n != 0 and n != -1 :
                if 1 :
                    s   += str(n) + '\t' + str(t) + '\t' + str(chain.rec[t][n][0]) + '\t' + str(chain.rec[t][n][1]) + '\t' + str(chain.rec[t][n][2])+ '\n'
                else :
                    print(chain.rec[t].keys())
                    print(chain)
        
        file_all.write(s)

    file_all.close()         
        
    chain.rec = {}
    chain.events = ''
    
# ========================================================================
# ============================ Loading ===================================
# ========================================================================

def load_file(filename) :
    time = []
    dat = np.loadtxt(filename)
    
    Nmax = int(np.max(dat[:, 0]))

    L_a = {}
    N_a = {}
    p_a = {}
    for i in range(len(dat)) :
        t = dat[i, 1]
        if t not in L_a.keys() :
            L_a[t] = {}
            N_a[t] = {}
            p_a[t] = {}
    
        L_a[t][int(dat[i, 0])] = dat[i, 2]
        N_a[t][int(dat[i, 0])] = dat[i, 3]
        p_a[t][int(dat[i, 0])] = dat[i, 4]

    L_array = np.zeros(( len(L_a.keys()), Nmax+1 ))
    N_array = np.zeros(( len(N_a.keys()), Nmax+1 ))
    p_array = np.zeros(( len(p_a.keys()), Nmax+1 ))

    step = -1
    for k in L_a.keys() :
        step += 1
        L_temp = [k]
        N_temp = [k]
        p_temp = [k]
        for j in range(1, Nmax+1) :
            if j in L_a[k].keys() :
                L_temp += [L_a[k][j]]
                N_temp += [N_a[k][j]]
                p_temp += [p_a[k][j]]
            else :
                L_temp += [None]
                N_temp += [None]
                p_temp += [None]
        L_array[step] = np.array(L_temp)
        N_array[step] = np.array(N_temp)
        p_array[step] = np.array(p_temp)
    return L_array, N_array, p_array

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
    
#
        