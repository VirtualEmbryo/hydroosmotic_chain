"""
flux.py library, containing analytical expressions of the functions 
used to integrate the ODE describing hydraulic or hydroosmotic lumens.


    Contains
    -------
        Hydraulic chain
    func_Lj_hydraulic   : ODE for length of hydraulic lumen
    func_JLh            : Left-hydraulic flux for hydraulic lumen
    func_JRh            : Right-hydraulic flux for hydraulic lumen

        Hydroosmotic chain
    func_Lj             : ODE for length of hydroosmotic lumen
    func_Nj             : ODE for number of ions of hydroosmotic lumen
    func_JLs            : Left-solute flux for hydroosmotic lumen
    func_JRs            : Right-solute flux for hydroosmotic lumen
    func_JLv            : Left-solvent flux for hydroosmotic lumen
    func_JRv            : Right-solvent flux for hydroosmotic lumen
    lam                 : Function lambda(x) from Supplementary information
    mu                  : Function mu(x) from Supplementary information
    I_1_minus_v2        : Function I_1^-(x) from Supplementary information
    I_2_minus_v2        : Function I_2^-(x) from Supplementary information
    I_1_plus_v2         : Function I_1^+(x) from Supplementary information
    I_2_plus_v2         : Function I_2^+(x) from Supplementary information
    func_deltaP_j       : Pressure jump of lumen j
    func_deltaC_j       : Concentration jump of lumen j

    func_Js_MF : Calculate the mean field flux of solute. Not used.
    func_Jv_MF : Calculate the mean field flux of solute. Not used.

    Requirements
    ------------
        Python libraries
    os
    sys
    warnings
    time
    numpy (np)
        
        Homemade libraries
    network
    topology
    functions
    
    
    NB : this file also inclues a shut-down of warning for RunTimeWarning.
"""


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
    func_Lj_hydraulic(index, t, L_vec, ell_vec, chain)
    
        Calculate the variation of the lumen index, such as :
    
        \frac{d L_index}{dt} = J_left + J_right + ja = Jjh
    
        where 
        ja is the pumping of the lumen index, 
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
            The chain-object.
    
        Returns
        -------
        Jjh : float 
            Variation of dL_j / dt
    
    """
    # CALCULATE THE FLUXES
    # i, k are such that :  i--j--k
    i, k = net.leftright_neighbors(index, chain)
    
    
    Jjh_left  = func_JLh(i_left = i, i_right = index, L_vec = L_vec, ell_vec = ell_vec, chain = chain)  # Flux from the left bridge
    Jjh_right = func_JRh(i_left = index, i_right = k, L_vec = L_vec, ell_vec = ell_vec, chain = chain)  # Flux from the right bridge
    
    # CALCULATE THE PUMPING
    try :
        func = chain.pumping_func                               # Search for the imposed pumping profile function
        args = chain.pumping_args                               # Search for the parameters of the profile
        L = chain.lumens_dict[index].length                     # Length of the lumen
        pos = chain.lumens_dict[index].pos                      # Center of  mass of the lumen
        x1, x2 = pos - L, pos + L                               # Positions of the borders of the lumens
        ca = functions.integrate(func, x1, x2, args)            # Calculate the total pumping for the lumen
        chain.lumens_dict[index].ca = ca
        
    except :
        ca = chain.lumens_dict[index].ca
    
    # FINAL FLUX
    Jjh = Jjh_left + Jjh_right + ca
    
    return Jjh
    
def func_JLh(i_left, i_right, L_vec, ell_vec, chain) :
    """
    func_JLh(i_left, i_right, L_vec, ell_vec, chain)
    
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
        
    return 1. / (ellt*L_R)*(1./L_L - 1. / L_R)
    
def func_JRh(i_left, i_right, L_vec, ell_vec, chain) :
    """
        Flux coming from the right of a lumen, thus i_left is the lumen for which the flux is calculated, i_right designates its neighbor.
        The equation for the variation of lumen l reads :
                \frac{d L_l}{dt} = |frac{1}{L_l \ell_{l,r}} ( L_r^{-1} - L_l^{-1} )
    
        NB :  there is no time in this flux.    
    """
    if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
        return 0.
    
    b = net.find_bridge(i_left, i_right, chain) # index of the connecting bridge
    
    ellt = ell_vec[b]
    
    L_L, L_R = L_vec[i_left], L_vec[i_right]
    
    return 1. / (ellt*L_L)*(1./L_R - 1./L_L)

# ========================================================================
# ===================== Osmotic Hydraulic ================================
# ========================================================================

def func_Lj(index, t, L_vec, N_vec, ell_vec, chain) :
    """Calculate the length variation of the lumen j, such as :
    
        $$\frac{d L_j}{dt} = \frac{1}{\tau_v}\left( \mu_j \nu_j (\mu_j  \frac{N_j}{L_j^2} - 1 - \frac{\epsilon_j}{L_j}) - \frac{\mu_j}{2 L_j} J_j^v \right)$$
    
        where
        j       : refers to the index of the considered lumen, i and k are its left and right neighbors.
        Lj, Nj  : are the length and number of ions of the lumen j
        eps_j   : physical parameter comparing hydraulic versus osmotic pressure.
        Jjv     : total flux of solvent coming from the bridges (left and right)
        tauv    : typical hydraulic relaxation time.
        mu_j, nu_j are geometrical constants    
    
    
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
            The chain-object.
    
        Returns
        -------
        dLjdt : float 
            Variation of dL_j / dt
    """
    
    # CALCULATE THE FLUXES
    i, k = net.leftright_neighbors(index, chain)
    
    Jjv_left  = func_JLv(i_left=i, i_right=index, L_vec=L_vec, N_vec=N_vec, ell_vec=ell_vec, chain=chain)
    Jjv_right = func_JRv(i_left=index, i_right=k, L_vec=L_vec, N_vec=N_vec, ell_vec=ell_vec, chain=chain)
    
    Jjv = Jjv_left + Jjv_right
    
    # GEOMETRIC VALUES
    mu_j, nu_j, eps_j = chain.lumens_dict[index].mu, chain.lumens_dict[index].nu, chain.lumens_dict[index].eps
    
    # VARIABLES
    Lj, Nj = L_vec[index], N_vec[index]
    
    return (mu_j*nu_j*(mu_j * Nj / (Lj*Lj) - 1. - eps_j / Lj) - mu_j * (Jjv) / (2.*Lj))/chain.tauv
    
def func_Nj(index, t, L_vec, N_vec, ell_vec, chain) :
    """Calculate the number of ions variation of the lumen j, such as :
    
        $$\frac{d N_j}{dt} = \frac{1}{\tau_s}\left( 2 \nu_j L_j (1 - \mu_j  \frac{N_j}{L_j^2} + j^a_i) - J_j^s \right)$$
    
        where
        j       : refers to the index of the considered lumen, i and k are its left and right neighbors.
        Lj, Nj  : are the length and number of ions of the lumen j
        j^a_i   : active pumping of the lumen j (denoted ca below)
        Jjs     : total flux of solute coming from the bridges (left and right)
        taus    : typical solute relaxation time.
        mu_j, nu_j are geometrical constants    
    
    
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
            The chain-object.
    
        Returns
        -------
        dNjdt : float 
            Variation of dN_j / dt
    """
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
        x1, x2 = pos - Lj, pos + Lj                         # positions of the left, right borders of the lumen j
        ca = functions.integrate(func, x1, x2, args)        # pumping of the lumen only
    except :
        ca = chain.lumens_dict[index].ca
    
    return (2.*nu_j*Lj*(1. + ca - mu_j*Nj/(Lj*Lj)) -  Jjs)/chain.taus
    
### OSMOTIC FLUXES
def func_JLs(i_left, i_right, L_vec, N_vec, ell_vec, chain) :
    """
    Solute flux coming from the left bridge of the lumen i_right
    
    Equation is detailed in the Supplementary Information (Eq. 20b)
    
    """
    # Fluxes from the border are considered to be zero.
    if chain.leaks == False : 
        if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
            return 0.
            
    b = net.find_bridge(i_left, i_right, chain)             # index of the bridge between i and j
    ellt = ell_vec[b]                                       # length of the bridge b
    
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
        ca_LR = functions.integrate(func, x1, x2, args)                         # Total pumping over the bridge b
    except :
        ca_LR = chain.bridges_dict[b].ca
    
    dC_L = func_deltaC_j(Nj=N_vec[i_left], Lj=L_vec[i_left], mu_j=mu_L, index=i_left)
    dC_R = func_deltaC_j(Nj=N_vec[i_right], Lj=L_vec[i_right], mu_j=mu_R, index=i_right)
    
    return chain.xis * ((dC_R-ca_LR)*np.cosh(1./chis) - (dC_L-ca_LR)) / np.sinh(1./chis)
    
def func_JRs(i_left, i_right, L_vec, N_vec, ell_vec, chain) :
    """
    Solute flux coming from the right bridge of the lumen i_left
    
    Equation is detailed in the Supplementary Information (Eq. 20a)
    
    """
    # Fluxes from the border are considered to be zero.
    if chain.leaks == False : 
        if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
            return 0.
    b = net.find_bridge(i_left, i_right, chain)             # index of the bridge between i and j
    ellt = ell_vec[b]                                       # length of the bridge b
    
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
        ca_LR = functions.integrate(func, x1, x2, args)                         # Total pumping over the bridge b
    except :
        ca_LR = chain.bridges_dict[b].ca
    
    dC_L = func_deltaC_j(Nj=N_vec[i_left], Lj=L_vec[i_left], mu_j=mu_L, index=i_left)   
    dC_R = func_deltaC_j(Nj=N_vec[i_right], Lj=L_vec[i_right], mu_j=mu_R, index=i_right)
    
    return chain.xis * ((dC_L-ca_LR)*np.cosh(1./chis) - (dC_R-ca_LR)) / np.sinh(1./chis)

### HYDRAULIC FLUXES
def func_JLv(i_left, i_right, L_vec, N_vec, ell_vec, chain) :
    """
    Solvent flux coming from the left bridge of the lumen i_right
    
    Equation is detailed in the Supplementary Information (Eq. 22b)
    
    """
    # Fluxes from the border are considered to be zero.
    if chain.leaks == False : 
        if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
            return 0.
    
    b = net.find_bridge(i_left, i_right, chain)             # index of the bridge between i and j
    ellt = ell_vec[b]                                       # length of the bridge b
    
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
        ca_LR = functions.integrate(func, x1, x2, args)                         # Total pumping over the bridge b
    except :
        ca_LR = chain.bridges_dict[b].ca
    
    dC_L = func_deltaC_j(Nj=N_vec[i_left], Lj=L_vec[i_left], mu_j=mu_L, index=i_left)
    dC_R = func_deltaC_j(Nj=N_vec[i_right], Lj=L_vec[i_right], mu_j=mu_R, index=i_right)
    
    P_L  = func_deltaP_j(Lj=L_vec[i_left], eps_j=eps_L, index=i_left)
    P_R  = func_deltaP_j(Lj=L_vec[i_right], eps_j=eps_R, index=i_right)
    
    la = lam(-0.5, chiv, chis, dC_L, dC_R, ca_LR)*np.exp(-0.5/chiv)
    mb = mu(0.5, chiv, chis, dC_L, dC_R, ca_LR)*np.exp(-0.5/chiv)
    
    return chain.xiv * ((P_R-mb)*np.cosh(1./chiv)/np.sinh(1./chiv) - (P_L-la)/np.sinh(1./chiv) - mb)
       
def func_JRv(i_left, i_right, L_vec, N_vec, ell_vec, chain) :
    """
    Solvent flux coming from the right bridge of the lumen i_left
    
    Equation is detailed in the Supplementary Information (Eq. 22b)
    
    """
    # Fluxes from the border are considered to be zero.
    if chain.leaks == False : 
        if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
            return 0.
    
    b = net.find_bridge(i_left, i_right, chain)                 # index of the bridge between i and j
    ellt = ell_vec[b]                                           # length of the bridge b
    
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
        ca_LR = functions.integrate(func, x1, x2, args)                         # Total pumping over the bridge b
    except :
        ca_LR = chain.bridges_dict[b].ca

    dC_L = func_deltaC_j(Nj=N_vec[i_left], Lj=L_vec[i_left], mu_j=mu_L, index=i_left)
    dC_R = func_deltaC_j(Nj=N_vec[i_right], Lj=L_vec[i_right], mu_j=mu_R, index=i_right)
    
    P_L  = func_deltaP_j(Lj=L_vec[i_left], eps_j=eps_L, index=i_left)
    P_R  = func_deltaP_j(Lj=L_vec[i_right], eps_j=eps_R, index=i_right)
    
    la = lam(-0.5, chiv, chis, dC_L, dC_R, ca_LR)*np.exp(-0.5/chiv)
    mb = mu(0.5, chiv, chis, dC_L, dC_R, ca_LR)*np.exp(-0.5/chiv)
    
    return chain.xiv * ((P_L-la)*np.cosh(1./chiv)/np.sinh(1./chiv) - (P_R-mb)/np.sinh(1./chiv) - la)

### FUNCTIONS
def lam(x, chiv, chis, dC_L, dC_R, ca_LR) :
    """
    Function lambda(x) as detailed in the Supplementary Information
    """
    la = 0.5*ca_LR * (np.exp(-x/chiv) - np.exp(-0.5/chiv))
    l_L = (dC_L-ca_LR) * I_1_minus_v2(x, chiv, chis) / (2.*chiv*np.sinh(1./chis))
    l_R = (dC_R-ca_LR) * I_2_minus_v2(x, chiv, chis) / (2.*chiv*np.sinh(1./chis))
    return la + l_L - l_R

def mu(x, chiv, chis, dC_L, dC_R, ca_LR) :
    """
    Function lambda(x) as detailed in the Supplementary Information
    """
    ma = 0.5*ca_LR * (np.exp(x/chiv) - np.exp(-0.5/chiv))
    m_L = (dC_L-ca_LR) * I_1_plus_v2(x, chiv, chis) / (2.*chiv*np.sinh(1./chis))
    m_R = (dC_R-ca_LR) * I_2_plus_v2(x, chiv, chis) / (2.*chiv*np.sinh(1./chis))
    return ma - m_L + m_R
    
def I_1_minus_v2(x, Xv, Xs) :
    """
    Function I_1^-(x),  as detail in the Supplementary Information.
    Solution was found using WolframAlpha
    """
    c1, s1 = np.cosh(1./Xs), np.sinh(1./Xs)
    cmx, smx = np.cosh((x-0.5)/Xs), np.sinh((x-0.5)/Xs)
    
    if Xv == Xs :
        return 0.5*( (x-0.5)*np.exp(-0.5/Xv) + 0.5*Xv*(np.exp(-2.*x/Xv) - np.exp(-1./Xv))*np.exp(0.5/Xv) )
    
    return Xv*Xs*(np.exp(-x/Xv)*(Xv*cmx + Xs*smx) - Xv*np.exp(-0.5/Xv))/ ( (Xv-Xs)*(Xv+Xs) )

def I_2_minus_v2(x, Xv, Xs) :
    """
    Function I_2^-(x),  as detail in the Supplementary Information.
    Solution was found using WolframAlpha
    """
    c1, s1 = np.cosh(1./Xs), np.sinh(1./Xs)
    cpx, spx = np.cosh((x+0.5)/Xs), np.sinh((x+0.5)/Xs)
    
    if Xv == Xs :
        return 0.5*( (x-0.5)*np.exp(0.5/Xv) + 0.5*Xv*(np.exp(-2.*x/Xv) - np.exp(-1./Xv))*np.exp(-0.5/Xv) )
    
    return Xv*Xs*(np.exp(-x/Xv)*(Xv*cpx + Xs*spx) - np.exp(-0.5/Xv)*(Xv*c1 + Xs*s1))/ ( (Xv-Xs)*(Xv+Xs) )

def I_1_plus_v2(x, Xv, Xs) :
    """
    Function I_1^+(x),  as detail in the Supplementary Information.
    Solution was found using WolframAlpha
    """
    c1, s1 = np.cosh(1./Xs), np.sinh(1./Xs)
    cmx, smx = np.cosh((x-0.5)/Xs), np.sinh((x-0.5)/Xs)
    
    if Xv == Xs :
        return 0.5*( 0.5*Xv*(np.exp(2.*x/Xv) - np.exp(-1./Xv))*np.exp(-0.5/Xv) - (x+0.5)*np.exp(0.5/Xv))
    
    return Xv*Xs*(np.exp(x/Xv)*(Xv*cmx - Xs*smx) - np.exp(-0.5/Xv)*(Xv*c1 + Xs*s1))/ ( (Xv-Xs)*(Xv+Xs) )

def I_2_plus_v2(x, Xv, Xs) :
    """
    Function I_2^+(x),  as detail in the Supplementary Information.
    Solution was found using WolframAlpha
    """
    c1, s1 = np.cosh(1./Xs), np.sinh(1./Xs)
    cpx, spx = np.cosh((x+0.5)/Xs), np.sinh((x+0.5)/Xs)
    
    if Xv == Xs :
        return 0.5*( 0.5*Xv*(np.exp(2.*x/Xv) - np.exp(-1./Xv))*np.exp(0.5/Xv) - (x+0.5)*np.exp(-0.5/Xv))
    
    return Xv*Xs*(np.exp(x/Xv)*(Xv*cpx - Xs*spx) - np.exp(-0.5/Xv)*Xv)/ ( (Xv-Xs)*(Xv+Xs) )

def func_deltaC_j(Nj, Lj, mu_j, index) :
    """
    func_deltaC_j(Nj, Lj, mu_j)
    
        Calculate the concentration jump of lumen j.
    
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
        return eps_j / Lj
    else : 
        return 0.
        
# ========================================================================
# ===================== Mean Field Fluxes ================================
# ========================================================================   
def func_Js_MF(index, L_vec, N_vec, ell_vec, chain) :
    #b = net.find_bridge(i_left, i_right, chain)             # index of the bridge between i and j
    #ellt = ell_vec[b]                                       # length of the bridge b
    
    ell_mf = 1./np.average(np.array([1./ell_vec[k] for k in ell_vec.keys() if k != 0 and k != -1])) # mean-field bridge length
    L_mf = 1./np.average(np.array([1./L_vec[k] for k in L_vec.keys() if k != 0 and k != -1]))       # mean-field lumen length
    N_mf = 1./np.average(np.array([1./N_vec[k] for k in N_vec.keys() if k != 0 and k != -1]))       # mean-field lumen nb_ions
    
    # SCREENING RATIOS
    chis  = chain.xis / ell_mf
    
    # GEOMETRICAL FACTORS
    mu_index  = chain.lumens_dict[index].mu
        
    dC_mf = 1./np.average([L_vec[k]**2 / (mu_index*N_vec[k]) for k in L_vec.keys() if k != 0 and k != -1])
    dC_index = func_deltaC_j(Nj=N_vec[index], Lj=L_vec[index], mu_j=mu_index, index=index)
    ca_LR = 0.
    
    return chain.xis * ((dC_index-ca_LR)*np.cosh(1./chis) - (dC_mf-ca_LR)) / np.sinh(1./chis)
    
def func_Jv_MF(index, L_vec, N_vec, ell_vec, chain) :
    #b = net.find_bridge(i_left, i_right, chain)             # index of the bridge between i and j
    #ellt = ell_vec[b]                                       # length of the bridge b
    ell_mf = 1./np.average(np.array([1./ell_vec[k] for k in ell_vec.keys() if k != 0 and k != -1])) # mean-field bridge length
    L_mf = 1./np.average(np.array([1./L_vec[k] for k in L_vec.keys() if k != 0 and k != -1]))       # mean-field lumen length
    N_mf = 1./np.average(np.array([1./N_vec[k] for k in N_vec.keys() if k != 0 and k != -1]))       # mean-field lumen nb_ions
    
    # SCREENING RATIOS
    chis  = chain.xis / ell_mf
    chiv  = chain.xiv / ell_mf
    
    # GEOMETRICAL FACTORS
    mu_index  = chain.lumens_dict[index].mu
    eps_index = chain.lumens_dict[index].eps
    
    dC_mf = 1./np.average([L_vec[k]**2 / (mu_index*N_vec[k]) for k in L_vec.keys() if k != 0 and k != -1])
    dC_index = func_deltaC_j(Nj=N_vec[index], Lj=L_vec[index], mu_j=mu_index, index=index)
    
    #P_mf  = eps_index / L_mf
    P_mf = eps_index / np.average(np.array([L_vec[k] for k in L_vec.keys() if k != 0 and k != -1]))
    P_index  = func_deltaP_j(Lj=L_vec[index], eps_j=eps_index, index=index)
    ca_LR = 0.
    
    la = lam(-0.5, chiv, chis, dC_mf, dC_index, ca_LR)*np.exp(-0.5/chiv)
    mb = mu(0.5, chiv, chis, dC_mf, dC_index, ca_LR)*np.exp(-0.5/chiv)
    
    return chain.xiv * ((P_index-mb)*np.cosh(1./chiv)/np.sinh(1./chiv) - (P_mf-la)/np.sinh(1./chiv) - mb)
    