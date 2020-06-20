import numpy as np
import os
import sys
import warnings
import time

module_path = os.path.abspath(os.path.join('..', 'chain_lumen/'))

if module_path not in sys.path :
    sys.path.append(module_path)

try :
    import _ressources.network as net
    import _ressources.topology
    import _ressources.functions as functions
except :
    import network as net
    import topology
    import functions
    
# ========================================================================
# =========================== Warnings ===================================
# ========================================================================
# This raises RuntimeWarning as an error, that will stop the simulation.
# It typically happens in hydroosmotic fluxes, when calculating cosh(x), x > 1000
warn = False
if warn :
    warnings.filterwarnings('error')
    try : 
        warnings.warn(RuntimeWarning())
    except Warning :
        print('RuntimeWarning is now an exception.')

# ========================================================================
# ========================== Hydraulic ===================================
# ========================================================================

def func_Lj_hydraulic(index, t, L_vec, ell_vec, chain) :
    """
    Calculate the variation of the lumen index, such as :
    
    \frac{d L_index}{dt} = J_left + J_right + ca = Jjh
    
    where 
    ca is the pumping of the lumen index, 
    J_left is the flux coming from the left (lumen i), 
    J_right is the flux coming from the right (lumen k).
    
    
    
    Parameters
    ----------
    index : int
        Index of the lumen
    t : float
        Current time.
    L_vec : dict
        Dictionnary of the lengths of lumens to use to calculate the fluxes.
    N_vec : dict
        Dictionnary of the lengths to use to calculate the fluxes.
    ell_vec : dict
        Dictionnary of the lengths of bridges to use to calculate the fluxes.
    chain : chain
        The chain.
    
    Returns
    -------
    Jjh : variation
    
    """
    # CALCULATE THE FLUXES
    i, k = net.leftright_neighbors(index, chain)
    
    Jjh_left  = func_JLh(i_left = i, i_right = index, L_vec = L_vec, ell_vec = ell_vec, chain = chain)
    Jjh_right = func_JRh(i_left = index, i_right = k, L_vec = L_vec, ell_vec = ell_vec, chain = chain)
    
    # CALCULATE THE PUMPING
    try :
        func = chain.pumping_func
        args = chain.pumping_args
        L = chain.lumens_dict[index].length
        pos = chain.lumens_dict[index].pos
        x1, x2 = pos - L, pos + L
        ca = functions.integrate(func, x1, x2, args)
        chain.lumens_dict[index].ca = ca
        
    except :
        ca = chain.lumens_dict[index].ca
    
    # FINAL FLUX
    Jjh = Jjh_left + Jjh_right + ca
    
    return Jjh
    
def func_JLh(i_left, i_right, L_vec, ell_vec, chain) :
    """
    
        Flux coming from the left of a lumen, thus i_right is the lumen for which the flux is calculated, i_left designates its neighbor.
        The equation for the "incoming" flux of lumen l reads :
                JLh = |frac{1}{L_r \ell_{r,l}} ( L_l^{-1} - L_r^{-1} )
        
        NB :  that there is no time in this flux.
            
    """
    if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
        return 0.
    
    b = net.find_bridge(i_left, i_right, chain) # index of the connecting bridge
    
    ellt = ell_vec[b]
    
    L_L, L_R = L_vec[i_left], L_vec[i_right]
    
    kappa_R = chain.lumens_dict[i_right].kappa
    
    return kappa_R / (ellt*L_R)*(1./L_L - 1. / L_R)
    
def func_JRh(i_left, i_right, L_vec, ell_vec, chain) :
    """
        Flux coming from the right of a lumen, thus i_left is the lumen for which the flux is calculated, i_right designates its neighbor.
        The equation for the variation of lumen l reads :
                \frac{d L_l}{dt} = |frac{1}{L_l \ell_{l,r}} ( L_r^{-1} - L_l^{-1} )
    
        NB :  that there is no time in this flux.
        NB : ca_l depends only on the lumen l
    
    """
    if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
        return 0.
    
    b = net.find_bridge(i_left, i_right, chain) # index of the connecting bridge
    
    ellt = ell_vec[b]
    
    L_L, L_R = L_vec[i_left], L_vec[i_right]
    
    kappa_L = chain.lumens_dict[i_left].kappa
    
    return kappa_L / (ellt*L_L)*(1./L_R - 1./L_L)

# ========================================================================
# ===================== Osmotic Hydraulic ================================
# ========================================================================

def func_Lj(index, t, L_vec, N_vec, ell_vec, chain) :
    # CALCULATE THE FLUXES
    i, k = net.leftright_neighbors(index, chain)
    
    Jjv_left  = func_JLv(i_left=i, i_right=index, L_vec=L_vec, N_vec=N_vec, ell_vec=ell_vec, chain=chain)
    Jjv_right = func_JRv(i_left=index, i_right=k, L_vec=L_vec, N_vec=N_vec, ell_vec=ell_vec, chain=chain)
    
    Jjv = Jjv_left + Jjv_right
    
    # GEOMETRIC VALUES
    mu_j, nu_j, eps_j = chain.lumens_dict[index].mu, chain.lumens_dict[index].nu, chain.lumens_dict[index].eps
    
    # VARIABLES
    Lj, Nj = L_vec[index], N_vec[index]
    
    ### TO REMOVE
    # For fluxes analysis
    #if index == 2 :
        #print('Lumen 2')
    #    Jlat  = mu_j*nu_j*(mu_j * Nj / (Lj*Lj) - 1. - eps_j / Lj)
    #    Jexch = mu_j * (Jjv) / (2.*Lj)
    #    
    #    if Jexch < Jlat and 'classification' not in chain.__dict__.keys() :
            #print(Jexch, Jlat)
    #        chain.classification = 'coarsening'
    #    
    #    f = open('fluxes.dat', 'a')
    #    f.write(str(t) + '\t' + str(Jlat) + '\t' + str(Jexch) + '\n')
    #    f.close()
        
    return (mu_j*nu_j*(mu_j * Nj / (Lj*Lj) - 1. - eps_j / Lj) - mu_j * (Jjv) / (2.*Lj))/chain.tauv
    
def func_Nj(index, t, L_vec, N_vec, ell_vec, chain) :
    # CALCULATE THE FLUXES
    i, k = net.find_neighbors(index, chain.bridges_dict)

    Jjs_left  = func_JLs(i_left=i, i_right=index, L_vec=L_vec, N_vec=N_vec, ell_vec=ell_vec, chain=chain)
    Jjs_right = func_JRs(i_left=index, i_right=k, L_vec=L_vec, N_vec=N_vec, ell_vec=ell_vec, chain=chain)
    
    Jjs = Jjs_left + Jjs_right
    
    # GEOMETRIC VALUES
    mu_j, nu_j = chain.lumens_dict[index].mu, chain.lumens_dict[index].nu
    
    # VARIABLES
    Lj, Nj = L_vec[index], N_vec[index]
    
    # PUMPING
    try :
        func = chain.pumping_func
        args = chain.pumping_args
        pos = chain.lumens_dict[index].pos
        x1, x2 = pos - Lj, pos + Lj
        ca = functions.integrate(func, x1, x2, args)        # pumping of the lumen only
    except :
        ca = chain.lumens_dict[index].ca
    
        
    return (2.*nu_j*Lj*(1. + ca - mu_j*Nj/(Lj*Lj)) -  Jjs)/chain.taus
    
### OSMOTIC FLUXES
def func_JLs(i_left, i_right, L_vec, N_vec, ell_vec, chain) :
    if chain.leaks == False : 
        if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
            return 0.
    b = net.find_bridge(i_left, i_right, chain)             # index of the bridge between i and j
    ellt = ell_vec[b]                                       # length of the bridge b
    
    # Total length between i and j (for normalization of the flux)
    #L0   = L_vec[i_right] + L_vec[i_left] + ellt
    
    # SCREENING RATIOS
    chis  = chain.xis / ellt
    
    # GEOMETRICAL FACTORS
    mu_L  = chain.lumens_dict[i_left].mu
    mu_R  = chain.lumens_dict[i_right].mu
    
    # PUMPING
    try :
        func = chain.pumping_func
        args = chain.pumping_args
        x1 = chain.lumens_dict[i_left].pos + chain.lumens_dict[i_left].length
        x2 = chain.lumens_dict[i_right].pos-chain.lumens_dict[i_right].length
        ca_LR = functions.integrate(func, x1, x2, args)
    except :
        ca_LR = chain.bridges_dict[b].ca
    
    dC_L = func_deltaC_j(Nj=N_vec[i_left], Lj=L_vec[i_left], mu_j=mu_L, index=i_left)
    dC_R = func_deltaC_j(Nj=N_vec[i_right], Lj=L_vec[i_right], mu_j=mu_R, index=i_right)
    
    #return chis*ellt/L0 * ((dC_R-ca_LR)*np.cosh(1./chis) - (dC_L-ca_LR)) / np.sinh(1./chis) ### THIS IS WRONG
    return chain.xis * ((dC_R-ca_LR)*np.cosh(1./chis) - (dC_L-ca_LR)) / np.sinh(1./chis)
    
def func_JRs(i_left, i_right, L_vec, N_vec, ell_vec, chain) :
    if chain.leaks == False : 
        if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
            return 0.
    b = net.find_bridge(i_left, i_right, chain)             # index of the bridge between i and j
    ellt = ell_vec[b]                                       # length of the bridge b
    
    # Total length between i and j (for normalization of the flux)
    #L0   = L_vec[i_right] + L_vec[i_left] + ellt
    
    # SCREENING RATIO
    chis  = chain.xis / ellt
    
    # GEOMETRICAL FACTORS
    mu_L  = chain.lumens_dict[i_left].mu
    mu_R  = chain.lumens_dict[i_right].mu
    
    # PUMPING
    try :
        func = chain.pumping_func
        args = chain.pumping_args
        x1 = chain.lumens_dict[i_left].pos + chain.lumens_dict[i_left].length
        x2 = chain.lumens_dict[i_right].pos-chain.lumens_dict[i_right].length
        ca_LR = functions.integrate(func, x1, x2, args)
    except :
        ca_LR = chain.bridges_dict[b].ca
    
    dC_L = func_deltaC_j(Nj=N_vec[i_left], Lj=L_vec[i_left], mu_j=mu_L, index=i_left)   
    dC_R = func_deltaC_j(Nj=N_vec[i_right], Lj=L_vec[i_right], mu_j=mu_R, index=i_right)
    
    #return chis*ellt/L0 * ((dC_L-ca_LR)*np.cosh(1./chis) - (dC_R-ca_LR)) / np.sinh(1./chis) ### THIS IS WRONG
    return chain.xis * ((dC_L-ca_LR)*np.cosh(1./chis) - (dC_R-ca_LR)) / np.sinh(1./chis)

### HYDRAULIC FLUXES
def func_JLv(i_left, i_right, L_vec, N_vec, ell_vec, chain) :
    if chain.leaks == False : 
        if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
            return 0.
    
    b = net.find_bridge(i_left, i_right, chain)             # index of the bridge between i and j
    ellt = ell_vec[b]                                       # length of the bridge b
    
    # Total length between i and j (for normalization of the flux)
    #L0   = L_vec[i_right] + L_vec[i_left] + ellt
    
    # SCREENING RATIOS
    chis  = chain.xis / ellt
    chiv  = chain.xiv / ellt
    
    # GEOMETRICAL FACTORS
    mu_L  = chain.lumens_dict[i_left].mu
    mu_R  = chain.lumens_dict[i_right].mu
    eps_L = chain.lumens_dict[i_left].eps
    eps_R = chain.lumens_dict[i_left].eps
    
    # PUMPING
    try :
        func = chain.pumping_func
        args = chain.pumping_args
        x1 = chain.lumens_dict[i_left].pos + chain.lumens_dict[i_left].length
        x2 = chain.lumens_dict[i_right].pos-chain.lumens_dict[i_right].length
        ca_LR = functions.integrate(func, x1, x2, args)
    except :
        ca_LR = chain.bridges_dict[b].ca
    
    dC_L = func_deltaC_j(Nj=N_vec[i_left], Lj=L_vec[i_left], mu_j=mu_L, index=i_left)
    dC_R = func_deltaC_j(Nj=N_vec[i_right], Lj=L_vec[i_right], mu_j=mu_R, index=i_right)
    
    
    P_L  = func_deltaP_j(Lj=L_vec[i_left], eps_j=eps_L, index=i_left)
    P_R  = func_deltaP_j(Lj=L_vec[i_right], eps_j=eps_R, index=i_right)
    
    la = lam(-0.5, chiv, chis, dC_L, dC_R, ca_LR)*np.exp(-0.5/chiv)
    mb = mu(0.5, chiv, chis, dC_L, dC_R, ca_LR)*np.exp(-0.5/chiv)
    
    #return chiv*ellt/L0 * ((P_R-mb)*np.cosh(1./chiv)/np.sinh(1./chiv) - (P_L-la)/np.sinh(1./chiv) - mb) ### THIS IS WRONG
    return chain.xiv * ((P_R-mb)*np.cosh(1./chiv)/np.sinh(1./chiv) - (P_L-la)/np.sinh(1./chiv) - mb)
       
def func_JRv(i_left, i_right, L_vec, N_vec, ell_vec, chain) :
    if chain.leaks == False : 
        if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
            return 0.
    
    b = net.find_bridge(i_left, i_right, chain)                 # index of the bridge between i and j
    ellt = ell_vec[b]                                           # length of the bridge b
    
    # Total length between i and j (for normalization of the flux)
    #L0   = L_vec[i_right] + L_vec[i_left] + ellt
    
    # SCREENING RATIOS
    chis  = chain.xis / ellt
    chiv  = chain.xiv / ellt
    
    # GEOMETRICAL FACTORS
    mu_L  = chain.lumens_dict[i_left].mu
    mu_R  = chain.lumens_dict[i_right].mu
    eps_L = chain.lumens_dict[i_left].eps
    eps_R = chain.lumens_dict[i_left].eps
    
    # PUMPING
    try :
        func = chain.pumping_func
        args = chain.pumping_args
        x1 = chain.lumens_dict[i_left].pos + chain.lumens_dict[i_left].length
        x2 = chain.lumens_dict[i_right].pos - chain.lumens_dict[i_right].length
        ca_LR = functions.integrate(func, x1, x2, args)
    except :
        ca_LR = chain.bridges_dict[b].ca

    dC_L = func_deltaC_j(Nj=N_vec[i_left], Lj=L_vec[i_left], mu_j=mu_L, index=i_left) ###
    dC_R = func_deltaC_j(Nj=N_vec[i_right], Lj=L_vec[i_right], mu_j=mu_R, index=i_right) ###
    
    P_L  = func_deltaP_j(Lj=L_vec[i_left], eps_j=eps_L, index=i_left) ###
    P_R  = func_deltaP_j(Lj=L_vec[i_right], eps_j=eps_R, index=i_right) ###
    
    la = lam(-0.5, chiv, chis, dC_L, dC_R, ca_LR)*np.exp(-0.5/chiv)
    mb = mu(0.5, chiv, chis, dC_L, dC_R, ca_LR)*np.exp(-0.5/chiv)
    
    #return chiv*ellt/L0 * ((P_L-la)*np.cosh(1./chiv)/np.sinh(1./chiv) - (P_R-mb)/np.sinh(1./chiv) - la) ### THIS IS WRONG
    return chain.xiv * ((P_L-la)*np.cosh(1./chiv)/np.sinh(1./chiv) - (P_R-mb)/np.sinh(1./chiv) - la)

### FUNCTIONS
def lam(x, chiv, chis, dC_L, dC_R, ca_LR) :
    la = 0.5*ca_LR * (np.exp(-x/chiv) - np.exp(-0.5/chiv))
    l_L = (dC_L-ca_LR) * I_1_minus_v2(x, chiv, chis) / (2.*chiv*np.sinh(1./chis))
    l_R = (dC_R-ca_LR) * I_2_minus_v2(x, chiv, chis) / (2.*chiv*np.sinh(1./chis))
    return la + l_L - l_R

def mu(x, chiv, chis, dC_L, dC_R, ca_LR) :
    ma = 0.5*ca_LR * (np.exp(x/chiv) - np.exp(-0.5/chiv))
    m_L = (dC_L-ca_LR) * I_1_plus_v2(x, chiv, chis) / (2.*chiv*np.sinh(1./chis))
    m_R = (dC_R-ca_LR) * I_2_plus_v2(x, chiv, chis) / (2.*chiv*np.sinh(1./chis))
    return ma - m_L + m_R
    
def I_1_minus_v2(x, Xv, Xs) :
    c1, s1 = np.cosh(1./Xs), np.sinh(1./Xs)
    cmx, smx = np.cosh((x-0.5)/Xs), np.sinh((x-0.5)/Xs)
    
    if Xv == Xs :
        return 0.5*( (x-0.5)*np.exp(-0.5/Xv) + 0.5*Xv*(np.exp(-2.*x/Xv) - np.exp(-1./Xv))*np.exp(0.5/Xv) )
    
    return Xv*Xs*(np.exp(-x/Xv)*(Xv*cmx + Xs*smx) - Xv*np.exp(-0.5/Xv))/ ( (Xv-Xs)*(Xv+Xs) )

def I_2_minus_v2(x, Xv, Xs) :
    c1, s1 = np.cosh(1./Xs), np.sinh(1./Xs)
    cpx, spx = np.cosh((x+0.5)/Xs), np.sinh((x+0.5)/Xs)
    
    if Xv == Xs :
        return 0.5*( (x-0.5)*np.exp(0.5/Xv) + 0.5*Xv*(np.exp(-2.*x/Xv) - np.exp(-1./Xv))*np.exp(-0.5/Xv) )
    
    return Xv*Xs*(np.exp(-x/Xv)*(Xv*cpx + Xs*spx) - np.exp(-0.5/Xv)*(Xv*c1 + Xs*s1))/ ( (Xv-Xs)*(Xv+Xs) )

def I_1_plus_v2(x, Xv, Xs) :
    c1, s1 = np.cosh(1./Xs), np.sinh(1./Xs)
    cmx, smx = np.cosh((x-0.5)/Xs), np.sinh((x-0.5)/Xs)
    
    if Xv == Xs :
        return 0.5*( 0.5*Xv*(np.exp(2.*x/Xv) - np.exp(-1./Xv))*np.exp(-0.5/Xv) - (x+0.5)*np.exp(0.5/Xv))
    
    return Xv*Xs*(np.exp(x/Xv)*(Xv*cmx - Xs*smx) - np.exp(-0.5/Xv)*(Xv*c1 + Xs*s1))/ ( (Xv-Xs)*(Xv+Xs) )

def I_2_plus_v2(x, Xv, Xs) :
    c1, s1 = np.cosh(1./Xs), np.sinh(1./Xs)
    cpx, spx = np.cosh((x+0.5)/Xs), np.sinh((x+0.5)/Xs)
    
    if Xv == Xs :
        return 0.5*( 0.5*Xv*(np.exp(2.*x/Xv) - np.exp(-1./Xv))*np.exp(0.5/Xv) - (x+0.5)*np.exp(-0.5/Xv))
    
    return Xv*Xs*(np.exp(x/Xv)*(Xv*cpx - Xs*spx) - np.exp(-0.5/Xv)*Xv)/ ( (Xv-Xs)*(Xv+Xs) )

def func_deltaC_j(Nj, Lj, mu_j, index) :
    """
    func_deltaC_j(Nj, Lj, mu_j)
    
        Calculate the concentration of lumen j.
    
        Parameters
        ----------
        Nj : float
            Number of ions in lumen j
        Lj : float
            Length of the lumen j
        mu_j : float
            Geometrical factor
    
        Returns
        -------
        delta Cj : float
            Value of the concentration in lumen j.
    """
    if index != 0 and index != -1 :
        # delta C
        return mu_j * Nj / (Lj*Lj) - 1.
    else : 
        return 0.
    
def func_deltaP_j(Lj, eps_j, index) :
    """
    func_deltaP_j(Lj, eps_j)
    
        Calculate the pressure of lumen j.
    
        Parameters
        ----------
        Lj : float
            Length of the lumen j
        eps_j : float
            Osmotic vs hydraulic pressure term for lumen j
    
        Returns
        -------
        Pj : float
            Value of the pressure in lumen j.
    """
    if index != 0 and index != -1 :
        #return eps_j / Lj - 1.          # Correction 4/05/2020
        return eps_j / Lj               # Before the correction
    else : 
        #print('Leaks')
        return 0.