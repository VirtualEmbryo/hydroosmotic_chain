import numpy as np

# ========================================================================
# =================== Network Generation =================================
# ========================================================================

def calc_inter_dist(L, pos) :
    """
    calc_inter_dist(L, pos)
    
        Calculate the distances between two consecutive lumens given the positions and the lengths of each lumen
    
        Parameters
        ----------
        L : array
            Array of the lengths (sizes) L_i of each lumen i
        pos : array
            Array of the positions of each lumen.
        
        Returns
        -------
        ell_list : array
            Array of the distances l_ij between lumen i and lumen j
    """
    ell_list = np.zeros(len(L)-1)
    for n in range(len(L)-1) :
        ell_list[n] = pos[n+1] - pos[n] - (L[n] + L[n+1])
    return ell_list
    
def make_network(lengths, positions, L0) :
    """
    make_network(lengths, positions, L0)
    
        Make lumens and bridges matrices for a network of connected lumens.
    
        Parameters
        ----------
        lengths : array
            Array of the lengths (sizes) L_i of each lumen i
        positions : array
            Array of the positions of each lumen.
        L0 : float
            Total length of the chain.
    
        Returns
        -------
        lumens : array
            Array of lumens, such that
            lumens = np.array([... ,
                                [i, L_i, pos_i],
                                , ...])
            where i is the lumen index. i = 0, -1 are reserved for boundaries (left and right respectively).
    
        bridges : array
            Array of bridges, such that
            bridges = np.array([... ,
                                [k, i, j, l_ij],
                                , ...])
            where k is the index of the bridge, (i, j) are indices of the connected lumens, l_ij is the length of the bridge k.
    """
    N = len(positions)

    inter_dis = calc_inter_dist(lengths, positions)
    lumens, bridges = np.zeros((N+2, 3)), np.zeros((N+1, 4))

    for i in range(1, N) :
        lumens[i] = np.array([i, lengths[i-1], positions[i-1]])
        bridges[i] = np.array([i, i, i+1, inter_dis[i-1]])
    lumens[N] = np.array([N, lengths[N-1], positions[N-1]])

    # borders
    lumens[0] = np.array([0, 0, 0])
    lumens[N+1] = np.array([-1, 0, L0])

    bridges[0] = np.array([0, 0, 1, positions[0]-lumens[1, 1]])
    bridges[N] = np.array([N, N, -1, L0-positions[N-1]-lumens[N, 1]])
    return lumens, bridges
    
def gen_random_conf(M, avg_size=0.5, std_size=0.1, avg_dist = 1., std_dist=0.1, dist_toleft=0.1, dist_toright=0.1) :
    """
    gen_random_conf(M, avg_size=0.5, std_size=0.1, avg_dist = 1., std_dist=0.1, dist_toleft=0.1, dist_toright=0.1)
    
        Generate a random configuration
        
        Parameters
        ----------
        M : int
            Number of lumens in the chain. 
            Note that the total number of lumens will be M+2 since there are also the borders, 0 and -1.
        avg_size : float, optional, default : 0.5
            Average size of a lumen. Mean of a normal distribution.
        std_size : float, optional, default : 0.1
            Standard deviation for the size of a lumen. Standard deviation of a normal distribution.
        avg_dist : float, optional, default : 1.
            Average distance between two lumens. Mean of a normal distribution.
        std_dist : float, optional, default : 0.1
            Standard deviation for the distance between two lumens. Standard deviation of a normal distribution.
    
        dist_to_left : float, optional, default : 0.1
            Distance between the left-most lumen and the left border.
        dist_to_right : float, optional, default : 0.1
            Distance between the right-most lumen and the right border.
        
        Returns
        -------
        lumens : array
            Array of lumens
        bridges : array
            Array of bridges
        L0 : float
            Total length of the chain.
    
    
    """
    lengths = np.random.normal(loc=avg_size, scale=std_size, size=M)
    
    ### TO REMOVE
    ###areas = np.random.normal(loc=avg_size, scale=std_size, size=M)
    ###theta = np.pi/3.
    ###mu = np.sin(theta)**2 / (2*theta - np.sin(2*theta))
    ###lengths = np.sqrt(mu*areas)
    ### TO REMOVE
    
    distances = np.random.normal(loc=avg_dist, scale=std_dist, size=M-1)
    positions = np.zeros(M)
    positions[0] = lengths[0] + dist_toleft
    for n in range(1, M) :
        positions[n] = positions[n-1] + lengths[n-1] + lengths[n] + distances[n-1]
        
    L0 = np.sum(distances) + 2*np.sum(lengths) + dist_toleft + dist_toright
    
    lumens, bridges = make_network(lengths, positions, L0)
    
    return lumens, bridges, L0
    
def gen_uniform_conf(M, avg_size=0.5, std_size=0.1, avg_dist = 1., std_dist=0.1, dist_toleft=0.1, dist_toright=0.1) :
    """
    gen_random_conf(M, avg_size=0.5, std_size=0.1, avg_dist = 1., std_dist=0.1, dist_toleft=0.1, dist_toright=0.1)
    
        Generate a random configuration
        
        Parameters
        ----------
        M : int
            Number of lumens in the chain. 
            Note that the total number of lumens will be M+2 since there are also the borders, 0 and -1.
        avg_size : float, optional, default : 0.5
            Average size of a lumen. Mean of a normal distribution.
        std_size : float, optional, default : 0.1
            Standard deviation for the size of a lumen. Standard deviation of a normal distribution.
        avg_dist : float, optional, default : 1.
            Average distance between two lumens. Mean of a normal distribution.
        std_dist : float, optional, default : 0.1
            Standard deviation for the distance between two lumens. Standard deviation of a normal distribution.
    
        dist_to_left : float, optional, default : 0.1
            Distance between the left-most lumen and the left border.
        dist_to_right : float, optional, default : 0.1
            Distance between the right-most lumen and the right border.
        
        Returns
        -------
        lumens : array
            Array of lumens
        bridges : array
            Array of bridges
        L0 : float
            Total length of the chain.
    
    
    """
    lengths = avg_size*np.ones(shape=M)
    distances = avg_size*np.ones(shape=M-1)
    positions = np.zeros(M)
    positions[0] = lengths[0] + dist_toleft
    for n in range(1, M) :
        positions[n] = positions[n-1] + lengths[n-1] + lengths[n] + distances[n-1]
        
    L0 = np.sum(distances) + 2*np.sum(lengths) + dist_toleft + dist_toright
    
    lumens, bridges = make_network(lengths, positions, L0)
    
    return lumens, bridges, L0
    
#=========

def calc_ell_list(chain) :
    """
    calc_ell_list(lumens, bridges)
    
        Calculate the distances between two consecutive lumens given lumens and bridges array.
    
        Parameters
        ----------
        lumens : array
            Array of the lumens
        bridges : array
            Array of the bridges
        
        Returns
        -------
        No returns, directly changes the chain.
    """
    ell_list = np.zeros(len(chain.bridges_dict))
    
    for b in chain.bridges_dict.keys() :
        i, j = chain.bridges_dict[b].lumen1, chain.bridges_dict[b].lumen2
        L_i, pos_i = chain.lumens_dict[i].length, chain.lumens_dict[i].pos
        L_j, pos_j = chain.lumens_dict[j].length, chain.lumens_dict[j].pos
        
        chain.bridges_dict[b].length = np.abs(pos_j - pos_i) - (L_i + L_j)
    
def find_neighbors(index, bridges_dict) :
    """
    
        Finds neighbors of lumen with index j
    """
    neighbors = []
    for k in bridges_dict.keys() :
        if bridges_dict[k].lumen1 == index :
            neighbors += [bridges_dict[k].lumen2]
        elif bridges_dict[k].lumen2 == index :
            neighbors += [bridges_dict[k].lumen1]
    return neighbors
    
def leftright_neighbors(j, chain) :
    """
    
        Returns indices of neighbors of lumen j, ordered by position.
    """
    i, k = find_neighbors(j, chain.bridges_dict)
    if chain.lumens_dict[i].pos <= chain.lumens_dict[k].pos :
        return [i, k]
    else :
        return [k, i]
        
def connected_bridges(index, bridges_dict) :
    """
    
    """
    bridges_index = []
    for k in bridges_dict.keys() :
        if bridges_dict[k].lumen1 == index or bridges_dict[k].lumen2 == index :
            bridges_index += [k]
    return bridges_index

def find_bridge(i, j, chain) :
    B = chain.bridges_dict
    for b in B.keys() :
        if (B[b].lumen1 == i and B[b].lumen2 == j) or (B[b].lumen1 == j and B[b].lumen2 == i) :
            br = b
    return br
    
# ========================================================================
# =================== Osmosis functions ==================================
# ========================================================================
def osmotic_equilibrium(L, nu, mu) :
    """
    
    """
    N = L**2 / mu
    return N
    
def gen_osmotic_equilibrium(L_list, M, nu_list, mu_list) :
    """
    gen_osmotic_equilibrium(L_list, M, nu_list, mu_list)
    
        Calculate an osmotic equilibrium
    """
    N_list = np.zeros(M)
    for m in range(M) :
        #N_list[m] = L_list[m]**2 / (2.*nu_list[m]*L_list[m]*mu_list[m])
        N_list[m] = osmotic_equilibrium(L_list[m], nu_list[m], mu_list[m])
        #print(m, L_list[m])
    #print(N_list)
    return N_list

def gen_ion_array(lumens, theta=np.pi/3., equilibrium = True, avg_nion=0.5, std_nion=1e-2) :
    """
    gen_ion_array(lumens, theta=np.pi/3., equilibrium = True)
    
        Generate an array for the ions given the lumens array
        
        Parameters
        ----------
        lumens : array
            Array of lumens
        theta : optional, float, default : np.pi/3.
            Angle of the lumens. Used to compute the factor nu and nu.
        equilibrium : optional, boolean, default : True
            If true, the configuration will be at osmotic equilibrium
            Else, it will be random.
        avg_nion : optional, float, default : 0.5
            Average number of ions in a lumen for random configuration
        std_nion : optional, float, default : 0.01
            Standard deviation of the number of ions in a lumen for random configuration
    """
    
    mat = np.zeros((len(lumens), 2))
    
    M = len(lumens)-2
    L_list = lumens[1:-1, 1]
    
    nu_list, mu_list = calc_nuj_list(np.ones(M)*theta), calc_muj_list(np.ones(M)*theta)        ### TO BE FIXED  !!!
    
    if equilibrium :
        Nions = gen_osmotic_equilibrium(L_list, M, nu_list, mu_list)
        mat[:, 0] = lumens[:, 0]
        mat[1:-1, 1] = Nions
    else :
        mat[:, 0] = lumens[:, 0]
        mat[1:-1, 1] = np.random.normal(avg_nion, std_nion, size=M)
    return mat
    
def calc_nuj_list(theta_list) :
    return theta_list / np.sin(2*theta_list)

def calc_muj_list(theta_list) :
    return np.sin(theta_list)**2 / (2*theta_list - np.sin(2*theta_list))

#