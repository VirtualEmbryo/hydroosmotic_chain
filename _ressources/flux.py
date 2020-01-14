import numpy as np
import os
import sys

module_path = os.path.abspath(os.path.join('..', 'chain_lumen/'))

if module_path not in sys.path :
    sys.path.append(module_path)

try :
    import _ressources.network as net
    import _ressources.topology
except :
    import network as net
    import topology

def calc_fluxes(chain, threshold=0.5, flux_val=1e-2, viscosity=1e-3, e0=1e-2, kappa=1) :
    fl = np.zeros(chain.nb_lumens)

    if chain.lumen_type == 'standard' :
        fl = standard_flux(chain.lumens_dict, threshold, flux_val)

    elif chain.lumen_type == 'hydraulic' :
        fl = hydraulic_flux(chain, kappa)
        
    elif chain.lumen_type == 'hydroosmotic' :
        fl = hydroosmotic_flux(chain)
               
    else :
        print('Flux is not valid')
    
    chain.fluxes = fl
    return fl
    
# ========================================================================
# ============================ Standard ==================================
# ========================================================================

def standard_flux(lumens_dict, threshold, flux_val) :
    flux = {}
    # NB : the output is LENGTH
    for l in lumens_dict.keys() :
        if l != 0 and l != -1 :
            if lumens_dict[l].length >= threshold :
                flux[l] = flux_val
            else :
                flux[l] = -flux_val
        else :
            flux[l] = 0.
            
    for b in bridges_dict.keys() :
        flux[(bridges_dict[b].lumen1, bridges_dict[b].lumen2)] = -(flux[bridges_dict[b].lumen1] + flux[bridges_dict[b].lumen2])
        
    return flux

# ========================================================================
# ========================== Hydraulic ===================================
# ========================================================================

def hydraulic_flux2(chain, kappa) :
    flux = {}
    lumens_dict, bridges_dict = chain.lumens_dict, chain.bridges_dict
    for j in lumens_dict.keys() :
        
        #phi_j = 0.5*kappa*lumens_dict[j].mu
        
        if j != 0 and j != -1 :
            i, k = net.find_neighbors(j, bridges_dict)                   # indices of the connected lumens
            br_ij, br_jk = net.connected_bridges(j, bridges_dict)        # indices of the connected bridges
            
            if i != 0 and i != -1 and k != 0 and k != -1 :
                l_ij, l_jk = bridges_dict[br_ij].length, bridges_dict[br_jk].length
                
                L_i = lumens_dict[i].length
                L_j = lumens_dict[j].length
                L_k = lumens_dict[k].length
                
                rho_i = np.sin(lumens_dict[i].theta)*chain.gamma/chain.total_length
                rho_j = np.sin(lumens_dict[j].theta)*chain.gamma/chain.total_length
                rho_k = np.sin(lumens_dict[k].theta)*chain.gamma/chain.total_length

                flux[j] = kappa/l_ij * (1./L_i-1./L_j) + kappa/l_jk * (1./L_k-1./L_j)
                
            elif (i==0 and k==-1) or (i==-1 and k==0) :
                #print 'no connection'
                flux[j] = 0.
                
            elif (i == 0 or i == -1) and (k!= 0 or k!= -1) :
                #print i, k
                #print('border left')
                L_j = lumens_dict[j].length
                L_k = lumens_dict[k].length
                
                l_jk = bridges_dict[br_jk].length
                
                rho_j = np.sin(lumens_dict[j].theta)*chain.gamma/chain.total_length
                rho_k = np.sin(lumens_dict[k].theta)*chain.gamma/chain.total_length
                
                flux[j] = 1./l_jk * (1./L_k-1./L_j) * kappa
                
            elif (k == 0 or k == -1) and (i != 0 or i != -1) :
                #print i, k
                #print('border right')
                L_i = lumens_dict[i].length
                L_j = lumens_dict[j].length
                l_ij = bridges_dict[br_ij].length
                
                rho_i = np.sin(lumens_dict[i].theta)*chain.gamma/chain.total_length
                rho_j = np.sin(lumens_dict[j].theta)*chain.gamma/chain.total_length
                
                flux[j] = 1./l_ij * (1./L_i-1./L_j) * kappa
        else :
            flux[j] = 0.
            
    for b in bridges_dict.keys() :
        flux[(bridges_dict[b].lumen1, bridges_dict[b].lumen2)] = -(np.sqrt(flux[bridges_dict[b].lumen1]) + np.sqrt(flux[bridges_dict[b].lumen2]))
    return flux

def func_Lj_hydraulic(index, t, L_vec, ell_vec, chain) :
    # CALCULATE THE FLUXES
    i, k = net.leftright_neighbors(index, chain)
    
    Jjh_left = func_JLh(i_left = i, i_right = index, L_vec = L_vec, ell_vec = ell_vec, chain = chain)
    Jjh_right = func_JRh(i_left = index, i_right = k, L_vec = L_vec, ell_vec = ell_vec, chain = chain)
    
    Jjh = Jjh_left + Jjh_right
    
    return chain.tau*Jjh
    
def func_JLh(i_left, i_right, L_vec, ell_vec, chain) :
    if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
        return 0.
    
    b = net.find_bridge(i_left, i_right, chain) # index of the connecting bridge
    
    ellt = ell_vec[b]
    L_L, L_R = L_vec[i_left], L_vec[i_right]
    L0   = L_vec[i_right] + L_vec[i_left] + ellt
    
    #chis  = chain.xis / ellt
    
    mu_R  = chain.lumens_dict[i_right].mu
    
    phi_R = 0.5*chain.kappa*mu_R / L0**2
    rho_L = chain.gamma * np.sin(chain.theta) / L0
    rho_R = chain.gamma * np.sin(chain.theta) / L0
    #ca_LR = chain.bridges_dict[b].ca    
    
    return phi_R / (ellt*L_R)*(rho_L/L_L - rho_R / L_R)
    
def func_JRh(i_left, i_right, L_vec, ell_vec, chain) :
    if i_left == 0 or i_left == -1 or i_right == 0 or i_right == -1 :
        return 0.
    
    b = net.find_bridge(i_left, i_right, chain) # index of the connecting bridge
    
    ellt = ell_vec[b]
    L_L, L_R = L_vec[i_left], L_vec[i_right]
    L0   = L_vec[i_right] + L_vec[i_left] + ellt
    
    #chis  = chain.xis / ellt
    
    mu_L  = chain.lumens_dict[i_left].mu
    
    phi_L = 0.5*chain.kappa*mu_L / L0**2
    rho_L = chain.gamma * np.sin(chain.theta) / L0
    rho_R = chain.gamma * np.sin(chain.theta) / L0
    #ca_LR = chain.bridges_dict[b].ca    
    
    return phi_L / (ellt*L_L)*(rho_R/L_R - rho_L / L_L)

# ========================================================================
# ===================== Osmotic Hydraulic ================================
# ========================================================================

def hydroosmotic_flux2(chain, left_v_flux=0., left_s_flux=0., right_v_flux=0., right_s_flux=0.) :    
    # Calculate the flux
    lumens_dict, bridges_dict = chain.lumens_dict, chain.bridges_dict
    taus, tauv = chain.taus, chain.tauv
    flux = {}
    xis, xiv = chain.xis, chain.xiv
    
    for j in lumens_dict :
        if j != 0 and j != -1 and lumens_dict[j].length >= 1e-6:
            i, k = net.leftright_neighbors(j, chain)              # indices of the connected lumens (i -- j -- k)
            br_ij, br_jk = net.connected_bridges(j, bridges_dict)       # indices of the connected bridges
            
            if (i == 0 and k == -1) :
                # The lumen j is connected to the two borders
                JLv = left_v_flux
                JRv = right_v_flux
                
                JLs = left_s_flux
                JRs = right_s_flux
                
            elif i == 0 :
                # The lumen j is connected to the left (0) border and k to the right
                JLv = left_v_flux
                JRv = func_JRv(lumen_L = lumens_dict[j], lumen_R = lumens_dict[k], bridge_LR = bridges_dict[br_jk], xis=xis, xiv=xiv, Ltot=chain.total_length)
                
                JLs = left_s_flux
                JRs = func_JRs(lumen_L = lumens_dict[j], lumen_R = lumens_dict[k], bridge_LR = bridges_dict[br_jk], xis=xis, Ltot=chain.total_length)
                
            elif k == -1 :
                # The lumen j is connected to i to the left and to the right (-1) border
                JLv = func_JLv(lumen_L = lumens_dict[i], lumen_R = lumens_dict[j], bridge_LR = bridges_dict[br_ij], xis=xis, xiv=xiv, Ltot=chain.total_length)
                JRv = right_v_flux
                
                JLs = func_JLs(lumen_L = lumens_dict[i], lumen_R = lumens_dict[j], bridge_LR = bridges_dict[br_ij], xis=xis, Ltot=chain.total_length)
                JRs = right_s_flux
                
            else :
                # The lumen j is connected to to two lumens i and k
                JLs = func_JLs(lumen_L = lumens_dict[i], lumen_R = lumens_dict[j], bridge_LR = bridges_dict[br_ij], xis=xis, Ltot=chain.total_length)
                JLv = func_JLv(lumen_L = lumens_dict[i], lumen_R = lumens_dict[j], bridge_LR = bridges_dict[br_ij], xis=xis, xiv=xiv, Ltot=chain.total_length)
                                
                JRs = func_JRs(lumen_L = lumens_dict[j], lumen_R = lumens_dict[k], bridge_LR = bridges_dict[br_jk], xis=xis, Ltot=chain.total_length)
                JRv = func_JRv(lumen_L = lumens_dict[j], lumen_R = lumens_dict[k], bridge_LR = bridges_dict[br_jk], xis=xis, xiv=xiv, Ltot=chain.total_length)
                
            
            J_jv = JLv + JRv
            J_js = JLs + JRs
            
            #print(j, JLv, JRv, JLs, JRs)
            dLdt = func_Lj(lumens_dict[j], J_jv, tauv)
            dNdt = func_Nj(lumens_dict[j], J_js, taus)
            
            flux[j] = [J_jv, J_js]
            #flux[j] = [dLdt, dNdt]
        else :
            flux[j] = [0., 0.]

    return flux

def func_Lj(index, t, L_vec, N_vec, ell_vec, chain) :
    # CALCULATE THE FLUXES
    i, k = net.leftright_neighbors(index, chain)
    
    #print(index, ' has neighbors ', i, ' (left) and ', k, ' (right)')
    Jjv_left  = func_JLv(i_left=i, i_right=index, L_vec=L_vec, N_vec=N_vec, ell_vec=ell_vec, chain=chain)
    Jjv_right = func_JRv(i_left=index, i_right=k, L_vec=L_vec, N_vec=N_vec, ell_vec=ell_vec, chain=chain)
    #print('H-Flux from left is  :', Jjv_left)
    #print('H-Flux from right is :', Jjv_right)
    Jjv = Jjv_left + Jjv_right
    
    # GEOMETRIC VALUES
    mu_j, nu_j, eps_j = chain.lumens_dict[index].mu, chain.lumens_dict[index].nu, chain.lumens_dict[index].eps
    
    Lj, Nj = L_vec[index], N_vec[index]

    ell_ij, ell_jk = 1, 1

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
    
    ell_ij, ell_jk = 1, 1
    
    mu_j, nu_j, ca = chain.lumens_dict[index].mu, chain.lumens_dict[index].nu, chain.lumens_dict[index].ca
    Lj, Nj = L_vec[index], N_vec[index]
    
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