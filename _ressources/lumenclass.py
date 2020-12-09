"""
lumenclass.py is a library containing the chain-object and associated sub-classes,
such as Lumen, Bridge, Osmotic_lumen, Osmotic_bridge and Osmotic_chain.

    Contains
    --------
        Hydraulic chain
    Chain           : Class defining the hydraulic chain.
    Lumen           : Class defining a (hydraulic) lumen.
    Bridge          : Class defining a (hydraulic) bridge.
    
        Osmotic Chain
    Osmotic_Chain   : Inherited class defining the hydro-osmotic chain.
    Osmotic_Lumen   : Inherited class defining a hydro-osmotic lumen.
    Osmotic_Bridge  : Inherited class defining a hydro-osmotic bridge.

    Requirements
    ------------
        Python Libraries
    os
    sys
    numpy (np)

        Homemade Libraries
    network (net)
    functions

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
    import _ressources.functions as functions
except :
    import network as net
    import functions
# ============================================================
# =================  Hydraulic Chain  ========================
# ============================================================

class Chain :
    def __init__(self, nb_lumens, theta=np.pi/3., l_dis=1e-2, l_merge=1e-2, pbc=False) :
        """
        Chain(nb_lumens, theta=np.pi/3., l_dis=1e-2, l_merge=1e-2, pbc=False)
        
            Object representing the chain. Contains attributes of the chain 
            
            Attributes
            ----------
                General parameters
            time : current time of the chain. Default : 0
            nb_lumens : number of lumens
            nmax : maximum number of lumens. Used for indexation
            lumen_type : 'hydraulic'
            total_length : total length of the chain, summing lumen and bridge lengths
            events : string that contains the events of the chain
                
                Contents
            lumens_dict : dictionnary, that contains the Lumen-objects. Each is indexed by an int. 0 and -1 are reserved only for borders
            bridges_dict : dictionnary, that contains the Bridges-objects. Each is indexed by an int
            
                Geometrical parameters
            theta : contact angle
            l_dis : length of collapse/disappearance
            l_merge : length of coalescence/merge
        
                Boundary conditions
            pbc : boolean, False for no periodic boundary conditions
            
                Recording
            rec : recording of the lumens
            rec_br : recording of the bridges
        
        """
        # General parameters
        self.time = 0
        self.nb_lumens = nb_lumens
        self.nmax = nb_lumens
        self.lumen_type = 'hydraulic'
        self.total_length = 0
        self.events = ''
        
        # Contents
        self.lumens_dict = {}
        self.bridges_dict = {}
        
        # Geometrical parameters
        self.theta = theta
        self.l_dis = l_dis
        self.l_merge = l_merge
        
        # Boundary conditions
        self.pbc = pbc
        
        # Recording
        self.rec = {} 
        self.rec_br = {}
        
    def __gen_network_lumen_object__(self, avg_size=0.5, std_size=0.1, avg_dist = 1., std_dist=0.1, dist_toleft=0.1, dist_toright=0.1, ca_lumen_list=[], ca_bridge_list=[]) :
        """
        self.__gen_network_lumen_object__(avg_size=0.5, std_size=0.1, avg_dist = 1., std_dist=0.1, 
                                          dist_toleft=0.1, dist_toright=0.1, ca_lumen_list=[], ca_bridge_list=[])
        
            Generate a network of hydraulic lumens.
        
            Parameters
            ----------
            avg_size : float, optional, default : 0.5
                Average size of the lumens
            std_size : float, optional, default : 0.1
                Standard deviation of the size of the lumens
            avg_dist : float, optional, default : 1.
                Average distance between the lumens
            std_dist : float, optional, default : 0.1
                Standard deviation of the distances between the lumens
            dist_toleft : float, optional, default : 0.1
                Distance of the left-most lumen from the left border
            dist_toright : float, optional, default : 0.1
                Distance of the right-most lumen from the right border
            ca_lumen_list : list, optional, default : []
                List of active pumping constant for the lumens
            ca_bridge_list :  list, optional, default : []
                List of active pumping constant for the bridges
        """
        # Generate random configuration of the lumens and the bridges
        lumens, bridges, self.total_length = net.gen_random_conf(self.nb_lumens, avg_size=avg_size, std_size=std_size, avg_dist=avg_dist, std_dist=std_dist, dist_toleft=dist_toleft, dist_toright=dist_toright)
        
        for b in range(len(bridges)) :
            self.bridges_dict[int(bridges[b, 0])] = Bridge(index=int(bridges[b, 0]), lumen1=bridges[b, 1], lumen2=bridges[b, 2], length=bridges[b, 3], ca=0)
        
        for m in range(self.nb_lumens+2) :
            if int(lumens[m, 0]) == 0 or int(lumens[m, 0]) == -1 or len(ca_lumen_list)==0 :
                ca = 0.
            else :
                ca = ca_lumen_list[m]
            self.lumens_dict[int(lumens[m, 0])] = Lumen(index = int(lumens[m, 0]), init_length=lumens[m,1], init_pos=lumens[m,2], theta=self.theta, ca = 0)
                    
        self.nmax = max(self.lumens_dict.keys())
        
        self.__record_state__()
    
    def __import_config__(self, lumens, bridges, eps=1e-3) :
        """
        Deprecated
        
        self.__import_config__(lumens, bridges, eps=1e-3)
            
            Import a configuration. 
        """
        self.nb_lumens = len(lumens)-2
        for m in range(self.nb_lumens+2) :
            self.lumens_dict[int(lumens[m, 0])] = Lumen(index = int(lumens[m, 0]), init_length = lumens[m, 2], init_pos = lumens[m, 1], theta = self.theta)
        
        for b in range(len(bridges)) :
            self.bridges_dict[int(bridges[b, 0])] = Bridge(index=int(bridges[b, 0]), lumen1=bridges[b, 1], lumen2=bridges[b, 2], length=bridges[b, 3])
            
        net.calc_ell_list(self)
        self.total_length = 2*np.sum(lumens[:, 2]) + np.sum(bridges[:, 3])
        
        if abs(self.lumens_dict[-1].pos - self.total_length) > eps :
            print('The right border position ('+str(self.lumens_dict[-1].pos)+') does not correspond to total length ('+str(self.total_length)+')')
    
    def __calc_total_length__(self) :
        """
        self.__calc_total_length__()
        
            Calculate the total length of the chain
        """
        total_length = 0.
        for b in self.bridges_dict.keys() :
            total_length += self.bridges_dict[b].length
        for k in self.lumens_dict.keys() :
            total_length += 2*self.lumens_dict[k].length
        return total_length
           
    def __give_lengths__(self) :
        """
        self.__give_lengths__()
            
            Gives the lengths of all lumens in the chain at the current step.
            Output is non-indexed.
        """
        lengths = np.zeros(self.nb_lumens)
        n = 0
        for k in self.lumens_dict.keys() :
            if k != 0 and k != -1 :
                lengths[n] = self.lumens_dict[k].length
                n+=1
        return lengths
    
    def __give_positions__(self) :
        """
        self.__give_positions__()
            
            Gives the positions of all lumens in the chain at the current step.
            Output is non-indexed.
        """
        positions = np.zeros((self.nb_lumens+2, 2))
        n = 0
        for k in self.lumens_dict.keys() :
            positions[n] = [self.lumens_dict[k].index, self.lumens_dict[k].pos]
            n += 1
        positions = positions[positions[:, 1].argsort()]
        return positions
        
    def __update_chain__(self) :
        """
        Deprecated 
        
        self.__update_chain__()
        
            Update the chain given the fluxes
        """
        for k in self.lumens_dict.keys() :
            self.lumens_dict[k].length += self.fluxes[k]
            
            if k != 0 and k != -1 :
                self.lumens_dict[k].__update__()
                
        for b in self.bridges_dict.keys() :
            self.bridges_dict[b].length += self.fluxes[(self.bridges_dict[b].lumen1, self.bridges_dict[b].lumen2)]
                
    def __record_state__(self) :
        """
        self.__record_state__()
        
            Record the current state of the chain in self.rec and self.rec_br
            Useful to save every time step and calculate end events.
        
        """
        self.rec[self.time] = {}
        self.rec_br[self.time] = {}
        
        for k in self.lumens_dict.keys() :
            self.rec[self.time][k] = [self.lumens_dict[k].length, self.lumens_dict[k].pos, self.lumens_dict[k].ca]
            #self.rec[self.time][k] = [self.lumens_dict[k].length, self.lumens_dict[k].pos]
        for b in self.bridges_dict.keys() :
            self.rec_br[self.time][b] = [self.bridges_dict[b].length, self.bridges_dict[b].lumen1, self.bridges_dict[b].lumen2]
           
    def __calc_ell_avg__(self) :
        """
        self.__calc_ell_avg__()
        
            Calculate the average of all the lengths of the bridges, 
            if there is at least one non-border lumen left (3 bridges or more).
            Otherwise, return None.
        
        """
        L = []
        if len(self.bridges_dict.keys()) > 2 :
            for b in self.bridges_dict.keys() :
                if self.bridges_dict[b].lumen1 != 0 and self.bridges_dict[b].lumen1 != -1 and self.bridges_dict[b].lumen2 != 0 and self.bridges_dict[b].lumen2 != -1 :
                    L += [self.bridges_dict[b].length]
        
            ellt_avg = np.average(L)
        
            return ellt_avg
        else : return None
        
    def __calc_L_avg__(self) :
        """
        self.__calc_L_avg__()
        
            Calculate the average of all the lengths of the lumens, 
            if there is at least one non-border lumen left (3 lumens or more).
            Otherwise, return None.
        
        """
        L = []
        if len(self.lumens_dict.keys()) > 2 :
            for l in self.lumens_dict.keys() :
                if l != 0 and l != -1 : 
                    L += [self.lumens_dict[l].length]
        
            Lt_avg = np.average(L)
        
            return Lt_avg
        else : return None
        
    def __calc_L_mf__(self) :
        """
        self.__calc_L_mfg__()
        
            Calculate the mean-field length of all the lengths of the lumens, 
            if there is at least one non-border lumen left (3 lumens or more).
            Otherwise, return None.
        
            The mean-field average is defined as the harmonic average.
        
        """
        L = []
        if len(self.lumens_dict.keys()) > 2 :
            for l in self.lumens_dict.keys() :
                if l != 0 and l != -1 : 
                    L += [1./self.lumens_dict[l].length]
        
            L_mf = 1./np.average(L)
        
            return L_mf
        else : return None
             
    def __L_vec__(self) :
        """
        self.__L_vec__()
            
            Returns the lengths of all lumens (including borders) as a column vector.
            Used for integration.
        
        """
        L_vec = {}
        for j in self.lumens_dict.keys() :
            L_vec[j] = float(self.lumens_dict[j].length)
        return L_vec
        
    def __ell_vec__(self) :
        """
        self.__ell_vec__()
            
            Returns the lengths of all bridge (including borders) as a column vector.
            Used for integration.
        
        """
        ell_vec = {}
        for b in self.bridges_dict.keys() :
            ell_vec[b] = float(self.bridges_dict[b].length)
        return ell_vec
        
    def __bridge_pos__(self, b) :
        """
        self.__bridge_pos__(b)
        
            Calculate the position of the center of the bridge b
        """
        lumen1 = self.bridges_dict[b].lumen1
        lumen2 = self.bridges_dict[b].lumen2
        pos1 = self.lumens_dict[lumen1].pos
        pos2 = self.lumens_dict[lumen2].pos
        return 0.5*(pos1+pos2)
        
    def __bridges_pos_list__(self) :
        """
        self.__bridges_pos_list__()
            
            Give the bridges positions
        """
        pos = {}
        for b in self.bridges_dict.keys() :
            pos[b] = self.__bridge_pos__(b)
        return pos
    
    def __merge__(self, k, i, j) :
        """
        self.merge(k, i, j)
            
            Merge lumens i and j into lumen k.
        
            The total mass is conserved.
        
            Returns
            -------
            pos_k : position of the new lumen k
            L_k : length of the new lumen k
        """
        # Lumen
        pos_i, pos_j = self.lumens_dict[i].pos, self.lumens_dict[j].pos
        L_i, L_j = self.lumens_dict[i].length, self.lumens_dict[j].length
        mu_i, mu_j = self.lumens_dict[i].mu, self.lumens_dict[j].mu
        area_i, area_j = L_i**2 / mu_i, L_j**2 / mu_j
        
        ca_i, ca_j = self.lumens_dict[i].ca, self.lumens_dict[j].ca
        
        # New lumen size
        mu_k = 0.5*(mu_i + mu_j)
        area_k = area_i+area_j
        L_k = np.sqrt(area_k * mu_k)
        
        # New lumen position
        # Given by the center of mass of the merging lumens
        pos_k = (pos_i*area_i + pos_j*area_j) / area_k
        
        #if len(self.pumping_args) > 0 :
        try :
            x1, x2 = pos_k - L_, pos_k + L_k
            ca_k = functions.integrate(self.pumping_func, x1, x2, self.pumping_args)
        #else :
        except :
            ca_k = 0.5*(ca_i+ca_j)
        
        self.lumens_dict[k] = Lumen(index=k, init_pos=pos_k, init_length=L_k, theta=self.theta, ca = ca_k)
        
        return pos_k, L_k
        
    def __str__(self) :
        """
        print(chain)
        
            Print the chain : its main attributes and the state of the lumens and bridges.
        """
        print('======= CHAIN =======')
        print('Type         : '+str(self.lumen_type))
        print('Total length : '+str(self.total_length))
        print('Current Time : '+str(self.time))
        print('======= LUMENS =======')
        print('Nb lumens : '+str(self.nb_lumens))
        for k in list(self.lumens_dict.keys()) :
            print(self.lumens_dict[k])
        print('======= BRIDGES ======')
        for b in self.bridges_dict.keys() :
            print(self.bridges_dict[b])
        return ''
    
    def __copy__(self) :
        """
        self.copy()
            
            Copy the chain in its current state.
        """
        attributes = self.__dict__.keys()
        
        cp_chain = Chain(self.nb_lumens, self.theta, self.l_dis, self.l_merge, self.pbc)
        
        cp_chain.lumens_dict = {k : self.lumens_dict[k].__copy__() for k in self.lumens_dict.keys()}
        cp_chain.bridges_dict = {b : self.bridges_dict[b].__copy__() for b in self.bridges_dict.keys()}
        
        # other attributes
        cp_chain.total_length = self.total_length
        cp_chain.time = self.time
        cp_chain.rec = self.rec
        cp_chain.rec_br = self.rec_br
        cp_chain.nmax = self.nmax
        cp_chain.events = self.events
        cp_chain.merge = self.merge
        
        try :
            cp_chain.pumping_args = self.pumping_args
            cp_chain.pumping_func = self.pumping_func
        except : pass
        
        return cp_chain
    
    def __save__(self, filename) :
        """
        self.save(filename)
        
            Save the state of the chain at filename.
        """
        s = ''
        s += '======= CHAIN =======' + '\n'
        s += 'Type         : '+str(self.lumen_type) + '\n'
        s += 'Total length : '+str(self.total_length) + '\n'
        s += 'Current Time : '+str(self.time) + '\n'
        
        s += '======= LUMENS =======' + '\n'
        s += 'Nb lumens : '+str(self.nb_lumens) + '\n'
        for k in list(self.lumens_dict.keys()) :
            s += self.lumens_dict[k].__save__() + '\n'
        s += '======= BRIDGES ======' + '\n'
        for b in self.bridges_dict.keys() :
            s += self.bridges_dict[b].__save__() + '\n'
        f = open(filename, 'w')
        f.write(s)
        f.close()
        return s
        
class Lumen :
    def __init__(self, index, init_pos, init_length, theta, ca=0.) :
        """
        Lumen(index, init_pos, init_length, theta, ca=0.)
        
            Object representing a hydraulic lumen.
        
            Attributes
            ----------
                Generation
            index : index of the lumen.
            init_pos : initial position of the lumen
            init_length : initial length of the lumen. Note that this is its half length
            theta : contact angle of the lumen
            ca : active pumping constant of the lumen
                
                Calculated quantities
            mu : constant depending on the contact angle theta. See Supplementary Information.
            nu : constant depending on the contact angle theta. See Supplementary Information.
            area : area of the lumen. Defined as A_i = L_i**2 / mu_i
            
                Dynamic quantities
            length : current length of the lumen
            pos : current position of the lumen
        """
        # Generation
        self.index = index
        self.init_pos = init_pos
        self.theta = theta
        self.ca = ca
        
        # Dynamic quantities
        self.length = init_length
        self.pos = self.init_pos
        
        # Calculated quantities
        self.__calc_geom_mu__()
        self.__calc_geom_nu__()
        self.__calc_area__()
        
    def __update__(self) :
        """
        Deprecated
        
        self.__update__()
                
            Update the chain by computing its area.
        """
        self.__calc_area__()
        
    def __calc_radius__(self) :
        """
        self.__calc_radius__()
            
            Return the radius of the lumen
        """
        self.radius = self.length / np.sin(self.theta)
        return self.radius
        
    def __calc_contact_angle__(self) :
        """
        Deprecated 
        
        self.__contact_angle__()
            
            Return the contact angle
        """
        return np.arccos(0.5 * self.contact_tension / self.tension)
    
    def __calc_gammafactor__(self) :
        """
        Deprecated
        
        self.__calc_gamma_factor()
        
            Return a gamma factor, similar to the one we defined in Dumortier et al., Science, 2019
        
        """
        return self.tension*np.sqrt(2*self.theta - np.sin(2*self.theta))
    
    def __calc_geom_mu__(self) :
        """
        self.__calc_geom_mu__()
        
            Returns the constant mu for the lumen
        """
        self.mu = np.sin(self.theta)**2 / (2*self.theta - np.sin(2*self.theta))
        return self.mu
    
    def __calc_geom_nu__(self) :
        """
        self.__calc_geom_nu__()
        
            Returns the constant nu for the lumen
        """
        self.nu = self.theta / np.sin(self.theta)
        return self.nu
        
    def __calc_area__(self) :
        """
        self.__calc_geom_mu__()
        
            Returns the area of the lumen
        """
        self.area = self.length**2 / self.mu
        return self.area
        
    def __str__(self) :
        """
        print(self)
        
            Print the state of the lumen
        """
        return "Lumen {0} is at position {1:.5f} with length {2:.5f} and pumping {3:.5f}".format(self.index, self.pos, self.length, self.ca)
        
    def __save__(self) :
        """
        self.save()
        
            Save the state of the lumen. Used to sve the chain state.
        """
        return "Lumen {0} is at position {1:.5f} with length {2:.5f} and pumping {3:.5f}".format(self.index, self.pos, self.length, self.ca)
    
    def __copy__(self) :
        """
        self.copy()
        
            Returns a copy of the lumen state.
        """
        return Lumen(self.index, self.init_pos, self.length, self.theta, self.ca)
        
class Bridge :
    def __init__(self, index, lumen1, lumen2, length, ca=0.) :
        """
        Bridge(index, lumen1, lumen2, length, ca=0.)
        
            Object representing a hydraulic bridge.
        
            Attributes
            ----------
            index : index of the bridge.
            lumen1, lumen2 : indices of the lumens connected via the bridges
            length : length of the bridge, delimited by the left and right borders of the lumens.
            ca : active pumping constant of the bridge. M
        
        """
        self.index = index
        self.lumen1 = int(lumen1)
        self.lumen2 = int(lumen2)
        self.length = length
        
        self.ca = ca

    def __str__(self) :
        """
        print(self)
        
            Print the state of the bridge.
        """
        return "Bridge {0} : ({1}, {2}) has length {3:.5f} with pumping {4:.5f}".format(self.index, self.lumen1, self.lumen2, self.length, self.ca)
    
    def __save__(self) :
        """
        self.__save__()
        
            Returns the current state of the bridge as a string.
            Used to save the chain.
        """
        return "Bridge {0} : ({1}, {2}) has length {3:.5} with pumping {4:.5f}".format(self.index, self.lumen1, self.lumen2, self.length, self.ca)
    
    def __copy__(self) :
        """
        self.copy()
            
            Returns a copy of the Bridge.
        """
        return Bridge(self.index, self.lumen1, self.lumen2, self.length, self.ca)
        
# ============================================================
# ==================  Osmotic Chain  =========================
# ============================================================

class Osmotic_Chain(Chain):
    
    def __init__(self, nb_lumens, theta=np.pi/3., l_dis=1e-2, l_merge=1e-2, pbc=False) :
        """
        Osmotic_Chain(nb_lumens, theta=np.pi/3., l_dis=1e-2, l_merge=1e-2, pbc=False)
        
            Object representing the chain. Contains attributes of the chain.
            Inherited from the class Chain.
            
            Attributes
            ----------
                General parameters
            time : current time of the chain. Default : 0
            nb_lumens : number of lumens
            nmax : maximum number of lumens. Used for indexation
            lumen_type : 'hydroosmotic'
            total_length : total length of the chain, summing lumen and bridge lengths
            events : string that contains the events of the chain
                
                Contents
            lumens_dict : dictionnary, that contains the Lumen-objects. Each is indexed by an int. 0 and -1 are reserved only for borders
            bridges_dict : dictionnary, that contains the Bridges-objects. Each is indexed by an int
            
                Geometrical parameters
            theta : contact angle
            l_dis : length of collapse/disappearance
            l_merge : length of coalescence/merge
        
                Boundary conditions
            pbc : boolean, False for no periodic boundary conditions
            
                Recording
            rec : recording of the lumens
            rec_br : recording of the bridges
            
                Screening lengths
            xis : solute screening length. Defined in the script chain.loadconfig from configuration file.
            xiv : solvent screening length. Defined in the script chain.loadconfig from configuration file.
            
        """
        
        Chain.__init__(self, nb_lumens, theta=theta, l_dis=l_dis, l_merge=l_merge, pbc=pbc)
        self.lumen_type = 'hydroosmotic'
        
    def __gen_network_lumen_object__(self, avg_size=0.5, std_size=0.1, avg_dist = 1., std_dist=0.1, dist_toleft=0.1, dist_toright=0.1, eps = 1e-3, ca_lumen_list=[], ca_bridge_list=[], pattern='normal', equilibrium=True, nions_avg=2, nions_std=1.) :
        """
        self.__gen_network_lumen_object__(avg_size=0.5, std_size=0.1, avg_dist = 1., std_dist=0.1, 
                                          dist_toleft=0.1, dist_toright=0.1, ca_lumen_list=[], ca_bridge_list=[])
        
            Generate a network of hydroosmotic lumens.
        
            Parameters
            ----------
            avg_size : float, optional, default : 0.5
                Average size of the lumens
            std_size : float, optional, default : 0.1
                Standard deviation of the size of the lumens
            avg_dist : float, optional, default : 1.
                Average distance between the lumens
            std_dist : float, optional, default : 0.1
                Standard deviation of the distances between the lumens
            dist_toleft : float, optional, default : 0.1
                Distance of the left-most lumen from the left border
            dist_toright : float, optional, default : 0.1
                Distance of the right-most lumen from the right border
            ca_lumen_list : list, optional, default : []
                List of active pumping constant for the lumens
            ca_bridge_list :  list, optional, default : []
                List of active pumping constant for the bridges
            pattern : string, optional, default : 'normal'
                May either be 'normal' or 'uniform'. 
                normal : the lumens are generated using a normal distribution with average area = avg_size and std deviation = atd_size
                uniform : the lumens all have the same size = avg_area and bridges have the same length = avg_dist
            equilibrium : boolean, optional, default : True
                If True, the lumens will be at osmotic equilibrium, the number of ions is calculated according to N = L**2 / mu                
                If False, the lumens will have a number of ions generated according to a normal distribution of ions.
            nions_avg : float, optional, default : 2
                Average number of ions in the lumens. Used only if equilibrium = False.
            nions_std : float, optional, default : 1.
                Standard deviation of the number of ions in the lumens. Used only if equilibrium = False. Standard deviation of a gaussian distribution.
        """
        # Generate the lumens and bridges array according to the pattern
        if pattern == 'normal' :
            lumens, bridges, self.total_length = net.gen_random_conf(self.nb_lumens, avg_size=avg_size, std_size=std_size, avg_dist=avg_dist, std_dist=std_dist, dist_toleft=dist_toleft, dist_toright=dist_toright)
        elif pattern == 'uniform' :           
            lumens, bridges, self.total_length = net.gen_uniform_conf(self.nb_lumens, avg_size=avg_size, std_size=std_size, avg_dist=avg_dist, std_dist=std_dist, dist_toleft=dist_toleft, dist_toright=dist_toright)
        
        for b in range(len(bridges)) :
            self.bridges_dict[int(bridges[b, 0])] = Osmotic_Bridge(index=int(bridges[b, 0]), lumen1=bridges[b, 1], lumen2=bridges[b, 2], length=bridges[b, 3], ca=0)
                
        for m in range(self.nb_lumens+2) :
            nu, mu = net.calc_nuj_list(self.theta), net.calc_muj_list(self.theta)
            if equilibrium == 'True' :
                nion = net.osmotic_equilibrium(L=lumens[m, 1], nu=nu, mu=mu)
            else :
                if int(lumens[m, 0]) != 0 and int(lumens[m, 0]) != -1 :
                    nion = -1
                    while nion <= 0 :
                        nion = np.random.normal(nions_avg, nions_std)
                else : 
                    nion = 0.
                    
            #self.lumens_dict[int(lumens[m, 0])] = Osmotic_Lumen(index = int(lumens[m, 0]), init_length=lumens[m,1], init_pos=lumens[m,2], init_nb_ions=nion, theta=self.theta, eps=eps, ca = ca_lumen_list[m])
            self.lumens_dict[int(lumens[m, 0])] = Osmotic_Lumen(index = int(lumens[m, 0]), init_length=lumens[m,1], init_pos=lumens[m,2], init_nb_ions=nion, theta=self.theta, eps=eps, ca = 0)
                    
        self.nmax = max(self.lumens_dict.keys())
        
        self.__record_state__()
        
    def __record_state__(self) :
        """
        self.__record_state__()
        
            Record the current state of the chain in self.rec and self.rec_br
            Useful to save every time step and calculate end events.
        
        """
        self.rec[self.time] = {}
        self.rec_br[self.time] = {}
        
        for k in self.lumens_dict.keys() :
            self.rec[self.time][k] = [self.lumens_dict[k].length, self.lumens_dict[k].nb_ions, self.lumens_dict[k].pos]
        for b in self.bridges_dict.keys() :
            self.rec_br[self.time][b] = [self.bridges_dict[b].length, self.bridges_dict[b].lumen1, self.bridges_dict[b].lumen2]
        
    def __import_config__(self, lumens, bridges, eps=1e-3) :
        """
        __import_config__(self, lumens, bridges, eps=1e-3)
        
            Import a configuration from pre-existing lumens and bridges arrays.
        
            Parameters
            ----------
            lumens :
                Array of the lumens.
            bridges : array
                Array of the bridges.
            eps : float, optional, default : 1e-3
                Laplace/Osmotic pressures ratio.
        """
        self.nb_lumens = len(lumens)-2
        for m in range(self.nb_lumens+2) :
            nu, mu = net.calc_nuj_list(self.theta), net.calc_muj_list(self.theta)
            self.lumens_dict[int(lumens[m, 0])] = Osmotic_Lumen(index = int(lumens[m, 0]), init_length = lumens[m, 2], init_pos = lumens[m, 1], init_nb_ions = lumens[m, 3], theta = self.theta, eps = eps, ca = lumens[m, 4])
        
        for b in range(len(bridges)) :
            self.bridges_dict[int(bridges[b, 0])] = Osmotic_Bridge(index=int(bridges[b, 0]), lumen1=bridges[b, 1], lumen2=bridges[b, 2], length=bridges[b, 3], ca=bridges[b, 4])
            
        net.calc_ell_list(self)
        self.total_length = 2*np.sum(lumens[:, 2]) + np.sum(bridges[:, 3])
        
        if abs(self.lumens_dict[-1].pos - self.total_length) > 1e-3 :
            print('The right border position ('+str(self.lumens_dict[-1].pos)+') does not correspond to total length ('+str(self.total_length)+')')
    
    def __L_vec__(self) :
        """
        self.__L_vec__()
            
            Returns the lengths of all lumens (including borders) as a column vector.
            Used for integration.
        
        """
        L_vec = {}
        for j in self.lumens_dict.keys() :
            L_vec[j] = float(self.lumens_dict[j].length)
        return L_vec
        
    def __N_vec__(self) :
        """
        self.__N_vec__()
            
            Returns the number of ions of all lumens (including borders) as a column vector.
            Used for integration.
        
        """
        N_vec = {}
        for j in self.lumens_dict.keys() :
            N_vec[j] = float(self.lumens_dict[j].nb_ions)
        return N_vec
        
    def __ell_vec__(self) :
        """
        self.__ell_vec__()
            
            Returns the lengths of all bridge (including borders) as a column vector.
            Used for integration.
        
        """
        ell_vec = {}
        for b in self.bridges_dict.keys() :
            ell_vec[b] = float(self.bridges_dict[b].length)
        return ell_vec
    
    def __calc_ell_avg__(self) :
        """
        self.__calc_ell_avg__()
        
            Calculate the average of all the lengths of the bridges, 
            if there is at least one non-border lumen left (3 bridges or more).
            Otherwise, return None.
        
        """
        L = []
        if len(self.bridges_dict.keys()) > 2 :
            for b in self.bridges_dict.keys() :
                if self.bridges_dict[b].lumen1 != 0 and self.bridges_dict[b].lumen1 != -1 and self.bridges_dict[b].lumen2 != 0 and self.bridges_dict[b].lumen2 != -1 :
                    L += [self.bridges_dict[b].length]
        
            ellt_avg = np.average(L)
        
            return ellt_avg
        else : return None
        
    def __calc_L_avg__(self) :
        """
        self.__calc_L_avg__()
        
            Calculate the average of all the lengths of the lumens, 
            if there is at least one non-border lumen left (3 lumens or more).
            Otherwise, return None.
        
        """
        L = []
        if len(self.lumens_dict.keys()) > 2 :
            for l in self.lumens_dict.keys() :
                if l != 0 and l != -1 : 
                    L += [self.lumens_dict[l].length]
        
            ellt_avg = np.average(L)
        
            return ellt_avg
        else : return None
                
    def __merge__(self, k, i, j) :
        """
        self.merge(k, i, j)
            
            Merge lumens i and j into lumen k.
        
            The total mass/area and number of ions are conserved.
        
            Returns
            -------
            pos_k : position of the new lumen k
            L_k : length of the new lumen k
        """
        # Lumen
        pos_i, pos_j = self.lumens_dict[i].pos, self.lumens_dict[j].pos
        L_i, L_j = self.lumens_dict[i].length, self.lumens_dict[j].length
        mu_i, mu_j = self.lumens_dict[i].mu, self.lumens_dict[j].mu
        mu_k = 0.5*(mu_i + mu_j)
        area_i, area_j = L_i**2 / mu_i, L_j**2 / mu_j
        
        area_k = area_i+area_j
        L_k = np.sqrt(area_k * mu_k)
        pos_k = (pos_i*area_i + pos_j*area_j) / area_k
        
        nb_ions_k = self.lumens_dict[i].nb_ions + self.lumens_dict[j].nb_ions
        eps_k = 0.5*(self.lumens_dict[i].eps + self.lumens_dict[j].eps)
        
        #if len(self.pumping_args) > 0 :
        try :
            x1, x2 = pos_k - L_, pos_k + L_k
            ca_k = functions.integrate(self.pumping_func, x1, x2, self.pumping_args)
        #else :
        except :
            ca_k = 0.5*(self.lumens_dict[i].ca+self.lumens_dict[j].ca)
        
        self.lumens_dict[k] = Osmotic_Lumen(index=k, init_pos=pos_k, init_length=L_k, theta=self.theta, init_nb_ions=nb_ions_k, eps=eps_k, ca=ca_k)
        
        # Bridge
        
        
        return pos_k, L_k

    def __copy__(self) :
        """
        self.__copy__()
            
            (Deep)-Copy the chain. Used to save the chain during time-integration.
        
            Returns
            -------
            cp_chain : chain-object
                Copy of self
        """
        attributes = self.__dict__.keys()
        
        cp_chain = Osmotic_Chain(self.nb_lumens, self.theta, self.l_dis, self.l_merge, self.pbc)
        
        cp_chain.lumens_dict = {k : self.lumens_dict[k].__copy__() for k in self.lumens_dict.keys()}
        cp_chain.bridges_dict = {b : self.bridges_dict[b].__copy__() for b in self.bridges_dict.keys()}
        
        # other attributes
        cp_chain.total_length = self.total_length
        cp_chain.time = self.time
        cp_chain.rec = self.rec
        cp_chain.rec_br = self.rec_br
        cp_chain.nmax = self.nmax
        cp_chain.events = self.events
        cp_chain.merge = self.merge
        
        try :
            cp_chain.pumping_args = self.pumping_args
            cp_chain.pumping_func = self.pumping_func
        except : pass
        
        cp_chain.taus = self.taus
        cp_chain.tauv = self.tauv
        cp_chain.xis = self.xis
        cp_chain.xiv = self.xiv
        cp_chain.leaks = self.leaks
        
        return cp_chain
    
    def __str__(self) :
        """
        print(self)
            
            Print the state of the chain.
        """
        print('======= CHAIN =======')
        print('Type         : '+str(self.lumen_type))
        print('Total length : '+str(self.total_length))
        print('Current Time : '+str(self.time))
        print('Screening lengths : ')
        print('        xi_s = '+str(self.xis))
        print('        xi_v = '+str(self.xiv))
        print('Permeation times :')
        print('       tau_s = '+str(self.taus))
        print('       tau_v = '+str(self.tauv))
        print('Pumping : ' + str(self.pumping))
        print('======= LUMENS =======')
        print('Nb lumens : '+str(self.nb_lumens))
        for k in list(self.lumens_dict.keys()) :
            print(self.lumens_dict[k])
        print('======= BRIDGES ======')
        for b in self.bridges_dict.keys() :
            print(self.bridges_dict[b])
        return ''
        
    def __save__(self, filename) :
        """
        print(self)
            
            Save the chain in a file
            
            Parameters
            ----------
            filename : string
                Name of the file.
        
            Returns
            -------
            s : string  
                String that contains the saved informations of the chain.
        """
        s = ''
        s += '======= CHAIN =======' + '\n'
        s += 'Type         : '+str(self.lumen_type) + '\n'
        s += 'Total length : '+str(self.total_length) + '\n'
        s += 'Current Time : '+str(self.time) + '\n'
        s += 'Screening lengths : ' + '\n'
        s += '        xi_s = '+str(self.xis) + '\n'
        s += '        xi_v = '+str(self.xiv) + '\n'
        s += 'Permeation times :' + '\n'
        s += '       tau_s = '+str(self.taus) + '\n'
        s += '       tau_v = '+str(self.tauv) + '\n'
        
        s += '======= LUMENS =======' + '\n'
        s += 'Nb lumens : '+str(self.nb_lumens) + '\n'
        for k in list(self.lumens_dict.keys()) :
            s += self.lumens_dict[k].__save__() + '\n'
        s += '======= BRIDGES ======' + '\n'
        for b in self.bridges_dict.keys() :
            s += self.bridges_dict[b].__save__() + '\n'
        f = open(filename, 'w')
        f.write(s)
        f.close()
        return s
        
class Osmotic_Lumen(Lumen) :
    def __init__(self, index, init_pos, init_length, init_nb_ions, theta, eps, ca) :
        """
        Osmotic_Lumen(index, init_pos, init_length, init_nb_ions, theta, eps, ca)
        
            Object representing a hydro-osmotic lumen. Class inherited from Lumen.
        
            Attributes
            ----------
                Generation
            index           : index of the lumen.
            init_pos        : initial position of the lumen
            init_length     : initial length of the lumen. Note that this is its half length
            init_nb_ions    : initial number of ions of the lumen
            theta           : contact angle of the lumen
            eps             : Laplace/Osmotic pressure ratio. Depends on the local tension, hence is lumen-dependent.
            ca              : active pumping constant of the lumen
                
                Calculated quantities
            mu              : constant depending on the contact angle theta. See Supplementary Information.
            nu              : constant depending on the contact angle theta. See Supplementary Information.
            area            : area of the lumen. Defined as A_i = L_i**2 / mu_i
            
                Dynamic quantities
            length          : current length of the lumen
            pos             : current position of the lumen
            nb_ions         : current number of ions of the lumen
        """
        #Lumen.__init__(self, index, init_pos, init_length, theta, eps)
        self.index = index
        self.init_length = init_length
        self.init_pos = init_pos
        self.init_nb_ions = init_nb_ions
        self.theta = theta
        
        self.length = init_length
        self.pos = init_pos
        self.nb_ions = init_nb_ions
        
        self.__calc_geom_mu__()
        self.__calc_geom_nu__()
        self.__calc_area__()
        
        self.eps = eps
        self.ca = ca
        
    def __calc_dconcentration__(self) :
        """
        self.__calc_dconcentration__()
            
            Calculate the concentration jump delta_C for lumen i such that
            $$\delta C = \mu \frac{N}{L^2} - 1$$            
        """
        self.dconcentration = self.mu *self.nb_ions / self.length**2 - 1.
        return self.dconcentration
    
    def __str__(self) :
        """
        print(self)
        
            Print the state of the lumen
        """
        return "Lumen {0} is at position {1:.5f} with length {2:.5f} and {3:.5f} ions with pumping {4:.5f}".format(self.index, self.pos, self.length, self.nb_ions, self.ca)
        
    def __copy__(self) :
        """
        self.copy()
        
            Returns a copy of the lumen state.
        """
        return Osmotic_Lumen(self.index, self.pos, self.length, self.nb_ions, self.theta, self.eps, self.ca)

    def __save__(self) :
        """
        self.save()
        
            Save the state of the lumen. Used to sve the chain state.
        """
        return "Lumen {0} is at position {1:.3f} with length {2:.3f} and {3:.3f} ions with pumping {4:.3f}".format(self.index, self.pos, self.length, self.nb_ions, self.ca)
        
class Osmotic_Bridge(Bridge) :
    def __init__(self, index, lumen1, lumen2, length, ca=0.) :
        """
        Osmotic_Bridge(index, lumen1, lumen2, length, ca=0.)
        
            Object representing a hydro-osmotic bridge. Class inherited from Bridge
        
            Attributes
            ----------
            index : index of the bridge.
            lumen1, lumen2 : indices of the lumens connected via the bridges
            length : length of the bridge, delimited by the left and right borders of the lumens.
            ca : active pumping constant of the bridge.
        
        """
        Bridge.__init__(self, index, lumen1, lumen2, length)
        self.ca = ca
        
    def __str__(self) :
        """
        print(self)
        
            Print the state of the bridge.
        """
        return "Bridge {0} : ({1}, {2}) has length {3:.5f} with pumping {4:.5f}".format(self.index, self.lumen1, self.lumen2, self.length, self.ca)
        
    def __save__(self) :
        """
        self.__save__()
        
            Returns the current state of the bridge as a string.
            Used to save the chain.
        """
        return "Bridge {0} : ({1}, {2}) has length {3:.5} with pumping {4:.3f}".format(self.index, self.lumen1, self.lumen2, self.length, self.ca)
        
    def __copy__(self) :
        """
        self.copy()
            
            Returns a copy of the Bridge.
        """
        return Osmotic_Bridge(self.index, self.lumen1, self.lumen2, self.length, self.ca)

#