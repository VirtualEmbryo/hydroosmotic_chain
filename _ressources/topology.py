import numpy as np
import os
import sys
module_path = os.path.abspath(os.path.join('..', 'chain_lumen/'))

if module_path not in sys.path :
    sys.path.append(module_path)
try :
    import _ressources.network as net
    import _ressources.lumenclass as lc
except :
    import network as net
    import lumenclass as lc

def topology(chain) :
    # Delete lumens
    delete(chain)
    
    count = 0    
    while borders(chain) or collisions(chain) :
        count += 1
        
        if count > 10 :
            print('More than 10 topological events... The time step is maybe too big.')
    
    chain.nb_lumens = len(chain.lumens_dict)-2
    #print(chain.time, chain.nb_lumens)

# ========================================================================
# =========================== Delete =====================================
# ========================================================================
def check_emptylumens(chain) :
    lumens_dict, l_dis = chain.lumens_dict, chain.l_dis
    indices_lumens = []
    
    for k in lumens_dict.keys() :
        if k != 0 and k != -1 and lumens_dict[k].length < l_dis :
            indices_lumens += [k]
    return indices_lumens

def delete_lumen(index, chain, tolerance = 1e-6) :
    #print('Delete Lumen', index)
    lumens_dict, bridges_dict = chain.lumens_dict, chain.bridges_dict
    neighbors = net.find_neighbors(index, bridges_dict)
    connect_br = net.connected_bridges(index, bridges_dict)
        
    if len(neighbors) < 2 :
        print('Border is not allow to disappear !')
    else :
        new_bridge_index = max(bridges_dict.keys()) + 1
        new_bridge_length = bridges_dict[connect_br[0]].length + bridges_dict[connect_br[1]].length + 2*lumens_dict[index].length
        
        if chain.lumen_type == 'hydraulic' :
            new_bridge = lc.Bridge(index = new_bridge_index, lumen1=neighbors[0], lumen2=neighbors[1], length=new_bridge_length)
            chain.bridges_dict[new_bridge_index] = new_bridge
        elif chain.lumen_type == 'hydroosmotic' :
            new_ca = 0.5*(bridges_dict[connect_br[0]].ca + bridges_dict[connect_br[1]].ca)      # pumping of the bridge is the average of the two neighboring bridges
            new_bridge = lc.Osmotic_Bridge(index = new_bridge_index, lumen1=neighbors[0], lumen2=neighbors[1], length=new_bridge_length, ca=new_ca)
            chain.bridges_dict[new_bridge_index] = new_bridge
    
        # Delete the lumen and the bridges
        del lumens_dict[index]
        del bridges_dict[connect_br[0]], bridges_dict[connect_br[1]]
        
        # Check
        old_L = chain.total_length
        new_L = chain.__calc_total_length__()
        
        if np.abs(old_L - new_L) > tolerance :
            print('Error for the total length !')
            print(old_L, new_L)
    
def delete(chain) :
    """
    Returns 1 if delete lumens
    Returns 0 else
    """
    lumens_dict, bridges_dict = chain.lumens_dict, chain.bridges_dict
    l_dis = chain.l_dis
    
    dis_lumens = check_emptylumens(chain)
    
    if len(dis_lumens) > 0 :
        disappearing = True
        for lum_index in dis_lumens :
            chain.events += 'Time ' + str(chain.time) + ' : lumen ' + str(lum_index) + ' disappears\n'
            delete_lumen(lum_index, chain)
        return 1
    return 0

# ========================================================================
# ======================== Collisions ====================================
# ========================================================================
def check_merginglumens(chain) :
    bridges_dict, l_merge = chain.bridges_dict, chain.l_merge
    to_merge = []
    for b in bridges_dict.keys() :
        if bridges_dict[b].lumen1 != 0 and bridges_dict[b].lumen2 != 0 and bridges_dict[b].lumen1 != -1 and bridges_dict[b].lumen2 != -1 and bridges_dict[b].length < l_merge :
            to_merge +=  [[b, bridges_dict[b].length]]
    to_merge = np.array(to_merge)
    if len(to_merge) > 0 :
        to_merge = to_merge[to_merge[:, 1].argsort()]  # sort the merging array from lowest to highest merging bridge length
    return to_merge

def merge_lumens(i, j, b_index, chain) :
    lumens_dict, bridges_dict = chain.lumens_dict, chain.bridges_dict
    if i == 0 or j == 0 or i == -1 or j == -1 :
        print('No merge with borders !')
        return
    # Connected bridges :
    br_i = net.connected_bridges(i, bridges_dict)
    br_j = net.connected_bridges(j, bridges_dict)
    for elem in br_i :
        if elem in br_j : # Remove the common bridge because of the merge !
            br_i.remove(elem)
            br_j.remove(elem)
    
    if bridges_dict[br_i[0]].lumen1 == i :
        a = bridges_dict[br_i[0]].lumen2
    else :
        a = bridges_dict[br_i[0]].lumen1
    pos_a = lumens_dict[a].pos
    L_a = lumens_dict[a].length
    
    if bridges_dict[br_j[0]].lumen1 == j :
        b = bridges_dict[br_j[0]].lumen2 
    else :
        b = bridges_dict[br_j[0]].lumen1
    pos_b = lumens_dict[b].pos
    L_b = lumens_dict[b].length
        
        
    # New lumen from merging
    k = chain.nmax+1
    
    # Merging procedure
    if chain.lumen_type == 'hydraulic' :
        pos_k, L_k = chain.__merge__(k, i, j)
        
    elif chain.lumen_type == 'hydroosmotic' :        
        pos_k, L_k = chain.__merge__(k, i, j)
        
    # plus others...
    
    # New bridges
    bmax = max(bridges_dict.keys())

    l_ak = np.abs(pos_a - pos_k) - (L_k+L_a)
    l_bk = np.abs(pos_b - pos_k) - (L_k+L_b)
    
    if chain.lumen_type == 'hydraulic' :
        bridges_dict[bmax+1] = lc.Bridge(index=bmax+1, lumen1=min(a, k), lumen2=max(a, k), length=l_ak)
        bridges_dict[bmax+2] = lc.Bridge(index=bmax+2, lumen1=min(b, k), lumen2=max(b, k), length=l_bk)

    elif chain.lumen_type == 'hydroosmotic' :
        bridges_dict[bmax+1] = lc.Osmotic_Bridge(index=bmax+1, lumen1=min(a, k), lumen2=max(a, k), length=l_ak, ca=bridges_dict[br_i[0]].ca)
        bridges_dict[bmax+2] = lc.Osmotic_Bridge(index=bmax+2, lumen1=min(b, k), lumen2=max(b, k), length=l_bk, ca=bridges_dict[br_j[0]].ca)
    
    del bridges_dict[br_i[0]], bridges_dict[br_j[0]], bridges_dict[b_index]
    del lumens_dict[i], lumens_dict[j]
    
    return k, 'Time ' + str(chain.time) + ' : ' + str(i) + ' and '+ str(j) + ' merge to give ' + str(k) + ' \n'
    
def collisions(chain, n_max_collisions=10) :
    """
    Returns 1 if merge lumens
    Returns 0 else
    """
    bridges_dict = chain.bridges_dict
    merging_lumens = check_merginglumens(chain)
    
    n_collisions = 0
    if len(merging_lumens) > 0 :
        while len(merging_lumens) > 0 :
            n_collisions += 1
            b_index = int(merging_lumens[0, 0])
            index1, index2 = bridges_dict[b_index].lumen1, bridges_dict[b_index].lumen2
            new_lumen_index, e = merge_lumens(index1, index2, b_index, chain)
            chain.events += 'Time ' + str(chain.time) + ' : ' + str(index1) + ' and '+ str(index2) + ' merge to give ' + str(new_lumen_index) + ' \n'
            chain.nmax += 1
            merging_lumens = check_merginglumens(chain)
        
            if n_collisions > n_max_collisions :
                print('More than '+str(n_max_collisions)+' collisions detected. Maybe the time step is too large !')
        return 1
    else :
        return 0
    
# ========================================================================
# ============================ Borders ===================================
# ========================================================================
def detect_lumens_borders(index, chain) :
    lumens_dict = chain.lumens_dict
    l_merge = chain.l_merge
    length_to_move = 0.
    
    l_ia = np.abs(lumens_dict[index].pos-lumens_dict[0].pos) - (lumens_dict[index].length)
    l_ib = np.abs(lumens_dict[index].pos-lumens_dict[-1].pos) - (lumens_dict[index].length)
    
    if l_ia < 0 or lumens_dict[index].pos < lumens_dict[0].pos:
        #print 'Must move right !'
        length_to_move = lumens_dict[index].length - (lumens_dict[index].pos-lumens_dict[0].pos) + l_merge#*0.1
        return  length_to_move
    
    elif l_ib < 0 or lumens_dict[index].pos > lumens_dict[-1].pos :
        #print 'Must move left !'
        length_to_move = -lumens_dict[index].length + (lumens_dict[-1].pos-lumens_dict[index].pos) - l_merge#*0.1
        return  length_to_move
    
    else : return length_to_move
    
def check_borders(chain) :
    lumens_dict = chain.lumens_dict
    to_move = []
    for i in lumens_dict.keys() :
        if i != 0 or i != -1 :
            displacement = detect_lumens_borders(i, chain)
            if displacement != 0. :
                to_move += [[i, displacement]]
    return to_move

def move_lumens_borders(index, displacement, chain) :
    lumens_dict, bridges_dict = chain.lumens_dict, chain.bridges_dict
    
    # new position
    lumens_dict[index].pos += displacement
    net.calc_ell_list(chain)

def borders(chain) :
    
    lumens_dict, bridges_dict = chain.lumens_dict, chain.bridges_dict
    
    moving_lumensborders = check_borders(chain)
    if len(moving_lumensborders) > 0 :
        for elem in moving_lumensborders :
            move_lumens_borders(elem[0], elem[1], chain)
        return 1
    else :
        return 0

# ========================================================================
# ============================ No ions ===================================
# ========================================================================
def empty_ions(chain) :
    empty_lumens = []
    for k in chain.lumens_dict.keys() :
        if k != 0 and k != -1 :
            if chain.lumens_dict[k].nb_ions < 0 :
                empty_lumens += [k]
                chain.lumens_dict[k].nb_ions = 0.
    return empty_lumens
    

#