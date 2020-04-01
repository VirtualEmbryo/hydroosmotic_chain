import numpy as np
import os
import sys
import warnings

module_path = os.path.abspath(os.path.join('..', 'chain_lumen/'))

if module_path not in sys.path :
    sys.path.append(module_path)

try :
    import _ressources.network as net
    import _ressources.topology
except :
    import network as net
    import topology
    
    
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
    # CALCULATE THE FLUXES
    i, k = net.leftright_neighbors(index, chain)
    
    Jjh_left  = func_JLh(i_left = i, i_right = index, L_vec = L_vec, ell_vec = ell_vec, chain = chain)
    Jjh_right = func_JRh(i_left = index, i_right = k, L_vec = L_vec, ell_vec = ell_vec, chain = chain)
    
    Jjh = Jjh_left + Jjh_right
    
    return Jjh
    
def func_JLh(i_left, i_right, L_vec, ell_vec, chain) :
    if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
        return 0.
    
    b = net.find_bridge(i_left, i_right, chain) # index of the connecting bridge
    
    ellt = ell_vec[b]
    #ellt = np.abs(chain.lumens_dict[i_left].pos - chain.lumens_dict[i_right].pos)
    #print(i_left, i_right, ellt)
    
    L_L, L_R = L_vec[i_left], L_vec[i_right]
    L0   = L_vec[i_right] + L_vec[i_left] + ellt
    
    #chis  = chain.xis / ellt
    
    mu_R  = chain.lumens_dict[i_right].mu
    
    phi_R = chain.lumens_dict[i_right].phi
    
    gamma_L = chain.lumens_dict[i_left].gamma
    gamma_R = chain.lumens_dict[i_right].gamma
    
    ca_R = chain.lumens_dict[i_right].ca
    
    return phi_R / (ellt*L_R)*(gamma_L/L_L - gamma_R / L_R) + ca_R
    
def func_JRh(i_left, i_right, L_vec, ell_vec, chain) :
    if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
        return 0.
    
    b = net.find_bridge(i_left, i_right, chain) # index of the connecting bridge
    
    ellt = ell_vec[b]
    #ellt = np.abs(chain.lumens_dict[i_left].pos - chain.lumens_dict[i_right].pos)
    #print(i_left, i_right, ellt)
    
    L_L, L_R = L_vec[i_left], L_vec[i_right]
    L0   = L_vec[i_right] + L_vec[i_left] + ellt
    
    #chis  = chain.xis / ellt
    
    mu_L  = chain.lumens_dict[i_left].mu
    
    phi_L = chain.lumens_dict[i_left].phi
    
    gamma_L = chain.lumens_dict[i_left].gamma
    gamma_R = chain.lumens_dict[i_right].gamma
    
    ca_L = chain.lumens_dict[i_left].ca 
    
    return phi_L / (ellt*L_L)*(gamma_R/L_R - gamma_L / L_L) + ca_L

# ========================================================================
# ===================== Osmotic Hydraulic ================================
# ========================================================================

def func_Lj(index, t, L_vec, N_vec, ell_vec, chain) :
    # CALCULATE THE FLUXES
    i, k = net.leftright_neighbors(index, chain)
    
    #print(index, ' has neighbors ', i, ' (left) and ', k, ' (right)')
    Jjv_left  = func_JLv(i_left=i, i_right=index, L_vec=L_vec, N_vec=N_vec, ell_vec=ell_vec, chain=chain)
    Jjv_right = func_JRv(i_left=index, i_right=k, L_vec=L_vec, N_vec=N_vec, ell_vec=ell_vec, chain=chain)
    #print('H-Flux from left is  :', Jjv_left)
    #print('H-Flux from right is :', Jjv_right)
    Jjv = Jjv_left + Jjv_right
    #Jjv = 0
    
    # GEOMETRIC VALUES
    mu_j, nu_j, eps_j = chain.lumens_dict[index].mu, chain.lumens_dict[index].nu, chain.lumens_dict[index].eps
    
    Lj, Nj = L_vec[index], N_vec[index]

    ell_ij, ell_jk = 1, 1
    save_fluxes = 0
    if save_fluxes :
        lateral_flux = mu_j*nu_j*(mu_j * Nj / (Lj*Lj) - 1. - eps_j / Lj)
        exchange_flux = mu_j * Jjv / (2.*Lj)
        print(lateral_flux, exchange_flux)
        #f = open('fluxes_L.dat', 'a')
        #f.write(str(t) + '\t' + str(index) + '\t' + str(lateral_flux) + '\t' + str(exchange_flux) + '\n')
        #f.close()
        
    return (mu_j*nu_j*(mu_j * Nj / (Lj*Lj) - 1. - eps_j / Lj) - mu_j * (Jjv) / (2.*Lj))/chain.tauv
    
def func_Nj(index, t, L_vec, N_vec, ell_vec, chain) :
    # CALCULATE THE FLUXES
    i, k = net.find_neighbors(index, chain.bridges_dict)
    #print(index, ' has neighbors ', i, ' (left) and ', k, ' (right)')
    Jjs_left  = func_JLs(i_left=i, i_right=index, L_vec=L_vec, N_vec=N_vec, ell_vec=ell_vec, chain=chain)
    Jjs_right = func_JRs(i_left=index, i_right=k, L_vec=L_vec, N_vec=N_vec, ell_vec=ell_vec, chain=chain)
    #print('O-Flux from left is  :', Jjs_left)
    #print('O-Flux from right is :', Jjs_right)
    
    Jjs = Jjs_left + Jjs_right
    #Jjs=0
    
    ell_ij, ell_jk = 1, 1
    
    mu_j, nu_j, ca = chain.lumens_dict[index].mu, chain.lumens_dict[index].nu, chain.lumens_dict[index].ca
    Lj, Nj = L_vec[index], N_vec[index]
    
    save_fluxes = 0
    if save_fluxes :
        lateral_flux = 2.*nu_j*Lj*(1. + ca - mu_j*Nj/(Lj*Lj)) 
        exchange_flux = Jjs
        print('Jjs', index, lateral_flux, exchange_flux)
        #f = open('fluxes_N.dat', 'a')
        #f.write(str(t) + '\t' + str(index) + '\t' + str(lateral_flux) + '\t' + str(exchange_flux) + '\n')
        #f.close()
        
    return (2.*nu_j*Lj*(1. + ca - mu_j*Nj/(Lj*Lj)) -  Jjs)/chain.taus
    
### OSMOTIC FLUXES
def func_JLs(i_left, i_right, L_vec, N_vec, ell_vec, chain) :
    if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
        return 0.
    
    b = net.find_bridge(i_left, i_right, chain)
    
    ellt = ell_vec[b]
    L0   = L_vec[i_right] + L_vec[i_left] + ellt
    
    chis  = chain.xis / ellt
    
    mu_L  = chain.lumens_dict[i_left].mu
    mu_R  = chain.lumens_dict[i_right].mu
    
    ca_LR = chain.bridges_dict[b].ca
    
    dC_L = func_C_j(Nj=N_vec[i_left], Lj=L_vec[i_left], mu_j=mu_L)
    dC_R = func_C_j(Nj=N_vec[i_right], Lj=L_vec[i_right], mu_j=mu_R)
    
    return chis*ellt/L0 * ((dC_R-ca_LR)*np.cosh(1./chis) - (dC_L-ca_LR)) / np.sinh(1./chis)
    
def func_JRs(i_left, i_right, L_vec, N_vec, ell_vec, chain) :
    if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
        return 0.
    
    b = net.find_bridge(i_left, i_right, chain)
    
    ellt = ell_vec[b]
    L0   = L_vec[i_right] + L_vec[i_left] + ellt
    
    chis  = chain.xis / ellt
    
    mu_L  = chain.lumens_dict[i_left].mu
    mu_R  = chain.lumens_dict[i_right].mu
    
    ca_LR = chain.bridges_dict[b].ca
    
    dC_L = func_C_j(Nj=N_vec[i_left], Lj=L_vec[i_left], mu_j=mu_L)
    dC_R = func_C_j(Nj=N_vec[i_right], Lj=L_vec[i_right], mu_j=mu_R)
    return chis*ellt/L0 * ((dC_L-ca_LR)*np.cosh(1./chis) - (dC_R-ca_LR)) / np.sinh(1./chis)

### HYDRAULIC FLUXES
def func_JLv(i_left, i_right, L_vec, N_vec, ell_vec, chain) :
    if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
        return 0.
    
    b = net.find_bridge(i_left, i_right, chain)
    
    ellt = ell_vec[b]
    L0   = L_vec[i_right] + L_vec[i_left] + ellt
    
    chis  = chain.xis / ellt
    chiv  = chain.xiv / ellt
    
    mu_L  = chain.lumens_dict[i_left].mu
    mu_R  = chain.lumens_dict[i_right].mu
    eps_L = chain.lumens_dict[i_left].eps
    eps_R = chain.lumens_dict[i_left].eps
    
    ca_LR = chain.bridges_dict[b].ca
    
    dC_L = func_C_j(Nj=N_vec[i_left], Lj=L_vec[i_left], mu_j=mu_L)
    dC_R = func_C_j(Nj=N_vec[i_right], Lj=L_vec[i_right], mu_j=mu_R)
    
    
    P_L  = func_P_j(Lj=L_vec[i_left], eps_j=eps_L, Ltot=L0)
    P_R  = func_P_j(Lj=L_vec[i_right], eps_j=eps_R, Ltot=L0)
    
    la = lam(-0.5, chiv, chis, dC_L, dC_R, ca_LR)*np.exp(-0.5/chiv)
    mb = mu(0.5, chiv, chis, dC_L, dC_R, ca_LR)*np.exp(-0.5/chiv)
    
    return chiv*ellt/L0 * ((P_R-mb)*np.cosh(1./chiv)/np.sinh(1./chiv) - (P_L-la)/np.sinh(1./chiv) - mb)
       
def func_JRv(i_left, i_right, L_vec, N_vec, ell_vec, chain) :
    if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
        return 0.
    
    b = net.find_bridge(i_left, i_right, chain)
    
    ellt = ell_vec[b]
    L0   = L_vec[i_right] + L_vec[i_left] + ellt
    
    chis  = chain.xis / ellt
    chiv  = chain.xiv / ellt
    
    mu_L  = chain.lumens_dict[i_left].mu
    mu_R  = chain.lumens_dict[i_right].mu
    eps_L = chain.lumens_dict[i_left].eps
    eps_R = chain.lumens_dict[i_left].eps
    
    ca_LR = chain.bridges_dict[b].ca
    
    dC_L = func_C_j(Nj=N_vec[i_left], Lj=L_vec[i_left], mu_j=mu_L)
    dC_R = func_C_j(Nj=N_vec[i_right], Lj=L_vec[i_right], mu_j=mu_R)
    
    
    P_L  = func_P_j(Lj=L_vec[i_left], eps_j=eps_L, Ltot=L0)
    P_R  = func_P_j(Lj=L_vec[i_right], eps_j=eps_R, Ltot=L0)
    
    la = lam(-0.5, chiv, chis, dC_L, dC_R, ca_LR)*np.exp(-0.5/chiv)
    mb = mu(0.5, chiv, chis, dC_L, dC_R, ca_LR)*np.exp(-0.5/chiv)
    return chiv*ellt/L0 * ((P_L-la)*np.cosh(1./chiv)/np.sinh(1./chiv) - (P_R-mb)/np.sinh(1./chiv) - la)

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

def func_C_j(Nj, Lj, mu_j) :
    """
    func_C_j(Nj, Lj, mu_j)
    
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
        Cj : float
            Value of the concentration in lumen j.
    """
    return mu_j * Nj / (Lj*Lj)
    
def func_P_j(Lj, eps_j, Ltot) :
    """
    func_P_j(Lj, eps_j, Ltot)
    
        Calculate the pressure of lumen j.
    
        Parameters
        ----------
        Lj : float
            Length of the lumen j
        eps_j : float
            Osmotic vs hydraulic pressure term for lumen j
        Ltot : float
            Total length of the system.
    
        Returns
        -------
        Pj : float
            Value of the pressure in lumen j.
    """
    return Ltot*eps_j / Lj