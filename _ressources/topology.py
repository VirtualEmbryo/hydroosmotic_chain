"""
topology.py is a library of functions used to manage the topological events 
that happen on a chain of micro-lumens.

    Contains
    --------
    topology                :  General function that handles topological events.

        Collapse
    check_empty_lumens      : Check for empty lumens in the  chain.
    delete_lumen            : Delete a given lumen from the chain.
    delete                  : Delete all empty lumens from the chain.

        Coalescence
    check_merginglumens     : Check for lumens that coalese/merge.
    merge_lumens            : Merge two given lumens.
    collisions              : Merge all coalescing lumens.
        
        Borders
    detect_lumens_borders   : Calculate the distances of a lumen to the borders.
    check_borders           : Check for lumens leaving the chain at the borders.
    move_lumens_borders     : Move a given lumen out of the borders.
    borders                 : Move all lumens from the borders.
    
        No ions
    empty_ions              : Avoids a negative number of ions within a lumen.
    
    Requirements
    ------------
        Python libraries
    os
    sys
    numpy (np)
    
        Homemade libraries
    network (net)
    lumenclass (lc)
    
Mathieu Le Verge--Serandour, 2020
"""


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
    """
    topology(chain)
        
        Handles the topological events of the chain. Proceeds as :
        i) Delete empty lumens
        ii) Merge coalescing lumens
        iii) Move the lumens overlapping the borders
        iv) Decrease the total number of lumens in the chain according to topological events.
    
        If an event is detected, the chain is immediatlely changed.
    
        A counter (count=10) is set such that a too large number of coalescence or border events is flagged and printed. 
        This does not stop the simulation.
    
        Parameters
        ----------
        chain : chain-object
    
    """
    # Delete lumens
    delete(chain)
    
    count = 0    
    if chain.merge :
        while borders(chain) or collisions(chain) :
            count += 1
        
            if count == 10 :
                print('More than 10 topological events (borders)... The time step is maybe too big.')
                
    else :
        #print('Merging not allowed')
        while borders(chain) :
            count += 1
        
            if count == 10 :
                print('More than 10 topological events (merge)... The time step is maybe too big.')
                
    
    chain.nb_lumens = len(chain.lumens_dict)-2

# ========================================================================
# ============================= Collapse =================================
# ========================================================================
def check_emptylumens(chain) :
    """
    check_emptylumens(chain)
    
        Check if lumens in the chain must disappear (collapse). 
        The rule for a collapsing lumen is determined with its length :
            L_i < l_dis
        where  l_dis is the disappearance length , and set inn the configuration.
        Usually, l_dis = 0.1
    
        Parameters
        ----------
        chain : chain-object
    
        Returns
        -------
        indices_lumens : list
            List of indices of the lumens that collapse.
    """
    lumens_dict, l_dis = chain.lumens_dict, chain.l_dis
    indices_lumens = []
    
    for k in lumens_dict.keys() :
        if k != 0 and k != -1 and lumens_dict[k].length < l_dis :
            indices_lumens += [k]
    return indices_lumens

def delete_lumen(index, chain, tolerance = 1e-6) :
    """
    delete_lumen(index, chain, tolerance = 1e-6)
    
        Delete a lumen with given index (j) from the chain. Borders are not allowed to disappear.
        Proceeds as :
        i) Neigbors (i, k) and connecting bridges (i-j) & (j-k) of the lumen are found.
        ii) Creates a new bridge between the neighbors (i-k) of the collapsing lumen.
        iii) Delete the collapsing lumen j and the two old connecting bridges (i-j) & (j-k)
        iv) Check whether the total length of the chain is conserved. If not, print a flag.
        
    
        Parameters
        ----------
        index : int
            Index of the lumen to delete
        chain : chain-object
            The chain from which to remove the lumen
        tolerance : float, optional, default : 1e-6
            Tolerance used to check the conservation of the total length ot the chain.
    """
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
    delete(chain)
    
        Delete all collapsing lumens from the chain. 
        Flag the events with time stamp for the events.txt file.
    
        Parameters
        ----------
        chain : chain-object
        
        Returns
        -------
        1 if deletion(s) happened
        0 otherwise
    
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
# ======================== Coalescence ===================================
# ========================================================================
def check_merginglumens(chain) :
    """
    check_merginglumens(chain)
    
        Check wether two lumens coalesce. Borders are not allowed to merge with lumens.
        Detects the lengths of bridges shorter than a length l_merge, set in the configuration file.
        Usually, l_merge = 0.001
    
        Parameters
        ----------
        chain : chain-object
        
        Returns
        -------
        to_merge : list
            List of bridges shorter than l_merge.
    """
    bridges_dict, l_merge = chain.bridges_dict, chain.l_merge
    to_merge = []
    
    for b in bridges_dict.keys() :
        # Prevents borders from merge
        if bridges_dict[b].lumen1 != 0 and bridges_dict[b].lumen2 != 0 and bridges_dict[b].lumen1 != -1 and bridges_dict[b].lumen2 != -1 and bridges_dict[b].length < l_merge :
            to_merge +=  [[b, bridges_dict[b].length]]
    to_merge = np.array(to_merge)
    
    if len(to_merge) > 0 :
        to_merge = to_merge[to_merge[:, 1].argsort()]  # sort the merging array from lowest to highest merging bridge length
    return to_merge

def merge_lumens(i, j, b_index, chain) :
    """
    merge_lumens(i, j, b_index, chain)
    
        Merge the lumens i and j together into a new lumen k, such that the total mass is conserved.
        Borders cannot merge with lumens.
    
        Proceeds as
        i) Find bridges connected to i (br_ia) and j (br_jb), and neighbors of i (a) and j (b)
        ii) Creates the new lumen from the merge, with index k, and calculates its position, length and nb of ions
        iii) Deletes old bridge and old lumenss, connect the new lumen k with neighors of i and j
    
        Rules for merge (mass conservation, new position, nb of ions conservation etc.) are defined in lumenclass.
    
        Parameters
        ----------
        i, j : int
            Indices of lumens to merge.
        b_index : int
            Index of the bridge connecting i and j
        chain : chain-object
    
        Returns
        -------
        k : int
            Index of the the lum
        event : str
            Line to print in the events.txt file.
        
    """
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
    
    # Neighbors of i and  j
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
    collisions(chain, n_max_collisions=10)
    
        Merge all coalescing lumens.
    
        Parameters
        ----------
        chain : chain-object
        n_max_collisions : int, optional, default : 10
            Max number of collisions in a time-step
    
        Returns
        -------
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
            
            # Check for coalescence that remain
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
    """
    detect_lumens_borders(index, chain)
    
        Calculate the length by which a lumen overlaps the border, to move later.
    
        Parameters
        ----------
        index : int
            Index of the lumen to check the overlap
        chain : chain-object
        
        Returns
        -------
        length_to_move : float
            Length by which a lumen overlaps the border.
    """
    lumens_dict = chain.lumens_dict
    l_merge = chain.l_merge
    length_to_move = 0.
    
    l_ia = np.abs(lumens_dict[index].pos-lumens_dict[0].pos) - (lumens_dict[index].length)
    l_ib = np.abs(lumens_dict[index].pos-lumens_dict[-1].pos) - (lumens_dict[index].length)
    
    if l_ia < 0 or lumens_dict[index].pos < lumens_dict[0].pos:
        #print('Must move right !')
        length_to_move = lumens_dict[index].length - (lumens_dict[index].pos-lumens_dict[0].pos) + l_merge
        return  length_to_move
    
    elif l_ib < 0 or lumens_dict[index].pos > lumens_dict[-1].pos :
        #print('Must move left !')
        length_to_move = -lumens_dict[index].length + (lumens_dict[-1].pos-lumens_dict[index].pos) - l_merge
        return  length_to_move
    
    else : return length_to_move
    
def check_borders(chain) :
    """
    check_borders(chain)
    
        Find all lumen that overlap the borders.
        
        Parameters
        ----------
        chain : chain-object
    
        Returns
        -------
        to_move : list
            List of lumens that overlap the border, such that it contains their index and 
            displacement by which to move them.
            
    """
    lumens_dict = chain.lumens_dict
    to_move = []
    for i in lumens_dict.keys() :
        if i != 0 or i != -1 :
            displacement = detect_lumens_borders(i, chain)
            if displacement != 0. :
                to_move += [[i, displacement]]
    return to_move

def move_lumens_borders(index, displacement, chain) :
    """
    move_lumens_borders(index, displacement, chain)
        
        Move a lumen by displacement.
        
        Parameters
        ----------
        index : int
            Index of the lumen to move
        displacement : float
            Value by which the lumen must move.
        chain : chain_object
    
    """
    lumens_dict, bridges_dict = chain.lumens_dict, chain.bridges_dict
    
    # new position
    lumens_dict[index].pos += displacement
    net.calc_ell_list(chain)

def borders(chain) :
    """
    borders(chain)
        Check the lumens that need to be displaced and move them.
    
        Parameters
        ----------
        chain : chain-object
    
        Returns
        -------
        1 if lumens were displaced
        0 otherwise
    """
    
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
    """
    empty_ions(chain)
    
        Forbid negative number of ions due to numerical errors.
    
        Parameters
        ----------
        chain : chain-object
    
        Returns
        -------
        empty_lumens : list
            List of the lumens where negative number of ions was  detected
    """
    empty_lumens = []
    for k in chain.lumens_dict.keys() :
        if k != 0 and k != -1 :
            if chain.lumens_dict[k].nb_ions < 0 :
                empty_lumens += [k]
                chain.lumens_dict[k].nb_ions = 0.
    return empty_lumens
    

#