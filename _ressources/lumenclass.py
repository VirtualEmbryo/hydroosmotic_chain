import numpy as np
import os
import sys
module_path = os.path.abspath(os.path.join('..', 'chain_lumen/'))

if module_path not in sys.path :
    sys.path.append(module_path)

import _ressources.network as net

# ============================================================
# =================  Hydraulic Chain  ========================
# ============================================================

class Chain :
    def __init__(self, nb_lumens, e0 = 0.1, theta=np.pi/3., l_dis=1e-2, l_merge=1e-2, pbc=True) :
        self.time = 0
        
        self.nb_lumens = nb_lumens
        self.nmax = nb_lumens
        self.lumen_type = 'hydraulic'
        self.lumens_dict = {}
        self.bridges_dict = {}
        self.total_length = 0
        self.events = ''
        
        # General parameters
        self.theta = theta
        self.e0 = e0
        self.l_dis = l_dis
        self.l_merge = l_merge
        
        # Boundary conditions
        self.pbc = pbc
        
        # Recording
        self.rec = {} 
        self.rec_br = {}
        
    def __gen_network_lumen_object__(self, avg_size=0.5, std_size=0.1, avg_dist = 1., std_dist=0.1, dist_toleft=0.1, dist_toright=0.1, eps = 1e-3) :
        lumens, bridges, self.total_length = net.gen_random_conf(self.nb_lumens, avg_size=avg_size, std_size=std_size, avg_dist=avg_dist, std_dist=std_dist, dist_toleft=dist_toleft, dist_toright=dist_toright)
        
        for b in range(len(bridges)) :
            self.bridges_dict[int(bridges[b, 0])] = Bridge(index=int(bridges[b, 0]), lumen1=bridges[b, 1], lumen2=bridges[b, 2], length=bridges[b, 3])
        
        for m in range(self.nb_lumens+2) :
            self.lumens_dict[int(lumens[m, 0])] = Lumen(index = int(lumens[m, 0]), init_length=lumens[m,1], init_pos=lumens[m,2], theta=self.theta)
                    
        self.nmax = max(self.lumens_dict.keys())
        
        self.__record_state__()
        
    def __calc_total_length__(self) :
        total_length = 0.
        for b in self.bridges_dict.keys() :
            total_length += self.bridges_dict[b].length
        for k in self.lumens_dict.keys() :
            total_length += 2*self.lumens_dict[k].length
        #self.total_length = total_length
        return total_length
        
    def __str__(self) :
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
    
    def __give_positions__(self) :
        positions = np.zeros((self.nb_lumens+2, 2))
        n = 0
        for k in self.lumens_dict.keys() :
            positions[n] = [self.lumens_dict[k].index, self.lumens_dict[k].pos]
            n += 1
        positions = positions[positions[:, 1].argsort()]
        return positions
        
    def __update_chain__(self) :
        for k in self.lumens_dict.keys() :
            self.lumens_dict[k].length += self.fluxes[k]
            
            if k != 0 and k != -1 :
                self.lumens_dict[k].__update__()
                
        for b in self.bridges_dict.keys() :
            self.bridges_dict[b].length += self.fluxes[(self.bridges_dict[b].lumen1, self.bridges_dict[b].lumen2)]
                
    def __record_state__(self) :
        self.rec[self.time] = {}
        self.rec_br[self.time] = {}
        
        for k in self.lumens_dict.keys() :
            self.rec[self.time][k] = [self.lumens_dict[k].length, self.lumens_dict[k].nb_ions, self.lumens_dict[k].pos]
        for b in self.bridges_dict.keys() :
            self.rec_br[self.time][b] = [self.bridges_dict[b].length]
    
class Lumen :
    def __init__(self, index, init_pos, init_length, theta) :
        self.index = index
        self.init_pos = init_pos
        self.init_length = init_length
        self.theta = theta
        
        self.length = self.init_length
        self.pos = self.init_pos
        
        self.length_list = [init_length]
        self.pos_list = [init_pos]
        
        self.__calc_geom_mu__()
        self.__calc_geom_nu__()
        self.__calc_area__()
        
    def __update__(self) :
        self.__calc_area__()
        
    def __calc_radius__(self) :
        self.radius = self.length / np.sin(self.theta)
        return self.radius
        
    def __calc_contact_angle__(self) :
        return np.arccos(0.5 * self.contact_tension / self.tension)
    
    def __calc_gammafactor__(self) :
        return self.tension*np.sqrt(2*self.theta - np.sin(2*self.theta))
    
    def __calc_geom_mu__(self) :
        self.mu = np.sin(self.theta)**2 / (2*self.theta - np.sin(2*self.theta))
        return self.mu
    
    def __calc_geom_nu__(self) :
        self.nu = self.theta / np.sin(self.theta)
        return self.nu
        
    def __calc_area__(self) :
        self.area = self.length**2 / self.mu
        return self.area
        
    def __str__(self) :
        return "Lumen {0} is at position {1:.5f} with length {2:.5f}".format(self.index, self.pos, self.length)
        
    def __save__(self, time) :
        self.length_list += [self.length]
        self.pos_list += [self.pos]
        #return 1
        
class Bridge :
    def __init__(self, index, lumen1, lumen2, length) :
        self.index = index
        self.lumen1 = int(lumen1)
        self.lumen2 = int(lumen2)
        self.length = length

    def __str__(self) :
        return "Bridge {0} : ({1}, {2}) has length {3:.5}".format(self.index, self.lumen1, self.lumen2, self.length)

# ============================================================
# ==================  Osmotic Chain  =========================
# ============================================================

class Osmotic_Chain(Chain):
    def __init__(self, nb_lumens, taus, tauv, xis=1., xiv=1., e0 = 0.1, theta=np.pi/3., l_dis=1e-2, l_merge=1e-2, pbc=True) :
        Chain.__init__(self, nb_lumens, e0 = e0, theta=theta, l_dis=l_dis, l_merge=l_merge, pbc=pbc)
        self.lumen_type = 'hydroosmotic'
        self.xis  = xis
        self.xiv  = xiv
        self.taus = taus
        self.tauv = tauv
        
    def __gen_network_lumen_object__(self, avg_size=0.5, std_size=0.1, avg_dist = 1., std_dist=0.1, dist_toleft=0.1, dist_toright=0.1, eps = 1e-3, ca_lumen_list=[], ca_bridge_list=[], equilibrium=True, nions_avg=2, nions_std=1.) :
        
        lumens, bridges, self.total_length = net.gen_random_conf(self.nb_lumens, avg_size=avg_size, std_size=std_size, avg_dist=avg_dist, std_dist=std_dist, dist_toleft=dist_toleft, dist_toright=dist_toright)
        
        for b in range(len(bridges)) :
            self.bridges_dict[int(bridges[b, 0])] = Osmotic_Bridge(index=int(bridges[b, 0]), lumen1=bridges[b, 1], lumen2=bridges[b, 2], length=bridges[b, 3], ca=ca_bridge_list[b])
                
        for m in range(self.nb_lumens+2) :
            nu, mu = net.calc_nuj_list(self.theta), net.calc_muj_list(self.theta)
            if equilibrium :
                nion = net.osmotic_equilibrium(L=lumens[m, 1], nu=nu, mu=mu)
            else : 
                nion = -1
                while nion <= 0 :
                    nion = np.random.normal(nions_avg, nions_std)
            self.lumens_dict[int(lumens[m, 0])] = Osmotic_Lumen(index = int(lumens[m, 0]), init_length=lumens[m,1], init_pos=lumens[m,2], init_nb_ions=nion, theta=self.theta, eps=eps, ca = ca_lumen_list[m])
                    
        self.nmax = max(self.lumens_dict.keys())
        
        self.__record_state__()
        
    def __record_state__(self) :
        
        self.rec[self.time] = {}
        self.rec_br[self.time] = {}
        
        for k in self.lumens_dict.keys() :
            self.rec[self.time][k] = [self.lumens_dict[k].length, self.lumens_dict[k].nb_ions, self.lumens_dict[k].pos]
        for b in self.bridges_dict.keys() :
            self.rec_br[self.time][b] = [self.bridges_dict[b].length]
             
    def __make_flux_vec__(self) :
        vec = []
        pos_vec = []
        k = 0
        for elem in self.fluxes.keys() :
            if type(elem) == int :
                vec += [self.fluxes[elem][0], self.fluxes[elem][1]]
                pos_vec += [[k, elem, 'length'], [k+1, elem, 'nb_ions']]
                k += 2
            elif type(elem) == tuple :
                vec += [self.fluxes[elem]]
                pos_vec += [[k, elem, 'bridge'], ]
                k += 1
        self.flux_vec = vec
        self.flux_vec_correspondances = pos_vec
        
    def __make_Y0__(self) :
        vec = []
        pos_vec = []
        i = 0
        
        for k in self.lumens_dict.keys() :
            vec += [self.lumens_dict[k].length, self.lumens_dict[k].nb_ions]
            pos_vec += [[i, k, 'length'], [i, k, 'nb_ions']]
            i += 2

        #for b in self.bridges_dict.keys() :
        #    vec += [self.bridges_dict[b].length]
        #    pos_vec += [[i, b, 'bridge']]
        #    i+= 1
            
        self.y0 = vec
        self.y0_correspondances = pos_vec
        
    def __update_chain__(self) :
        for k in self.lumens_dict.keys() :
            self.lumens_dict[k].length += self.fluxes[k][0]
            self.lumens_dict[k].nb_ions += self.fluxes[k][1]
            
            if k != 0 and k != -1 :
                self.lumens_dict[k].__update__()
                
        for b in self.bridges_dict.keys() :
            self.bridges_dict[b].length += self.fluxes[(self.bridges_dict[b].lumen1, self.bridges_dict[b].lumen2)]
    
    def __copy__(self) :
        attributes = self.__dict__.keys()
        
        cp_chain = Osmotic_Chain(self.nb_lumens, self.xis, self.xiv, self.taus, self.tauv, self.e0, self.theta, self.l_dis, self.l_merge, self.pbc)
        
        cp_chain.lumens_dict = {k : self.lumens_dict[k].__copy__() for k in self.lumens_dict.keys()}
        cp_chain.bridges_dict = {b : self.bridges_dict[b].__copy__() for b in self.bridges_dict.keys()}
        
        # other attributes
        cp_chain.total_length = self.total_length
        cp_chain.time = self.time
        cp_chain.rec = self.rec
        cp_chain.rec_br = self.rec_br
        cp_chain.nmax = self.nmax
        cp_chain.events = self.events
        
        
        return cp_chain
        
    def __import_config__(self, lumens, bridges, eps=1e-3) :
            
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
            
    def __str__(self) :
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
    
    def __L_vec__(self) :
        L_vec = {}
        for j in self.lumens_dict.keys() :
            L_vec[j] = self.lumens_dict[j].length
        return L_vec
        
    def __N_vec__(self) :
        N_vec = {}
        for j in self.lumens_dict.keys() :
            N_vec[j] = self.lumens_dict[j].nb_ions
        return N_vec
        
    def __ell_vec__(self) :
        ell_vec = {}
        for b in self.bridges_dict.keys() :
            ell_vec[b] = self.bridges_dict[b].length
        return ell_vec
    
class Osmotic_Lumen(Lumen) :
    def __init__(self, index, init_pos, init_length, init_nb_ions, theta, eps, ca) :
        Lumen.__init__(self, index, init_pos, init_length, theta)
        self.init_nb_ions = init_nb_ions
        self.nb_ions = init_nb_ions
        self.eps = eps
        self.ca = ca
        
    def __calc_dconcentration__(self) :
        self.dconcentration = self.mu *self.nb_ions / self.length**2 - 1.
        return self.dconcentration
    
    def __str__(self) :
        return "Lumen {0} is at position {1:.3f} with length {2:.3f} and {3:.3f} ions with pumping {4:.3f}".format(self.index, self.pos, self.length, self.nb_ions, self.ca)
        
    def __copy__(self) :
        return Osmotic_Lumen(self.index, self.init_pos, self.init_length, self.init_nb_ions, self.theta, self.eps, self.ca)

    def __save__(self) :
        return "Lumen {0} is at position {1:.3f} with length {2:.3f} and {3:.3f} ions with pumping {4:.3f}".format(self.index, self.pos, self.length, self.nb_ions, self.ca)
        
class Osmotic_Bridge(Bridge) :
    def __init__(self, index, lumen1, lumen2, length, ca=0.) :
        Bridge.__init__(self, index, lumen1, lumen2, length)
        self.ca = ca
        
    def __str__(self) :
        return "Bridge {0} : ({1}, {2}) has length {3:.5} with pumping {4:.3f}".format(self.index, self.lumen1, self.lumen2, self.length, self.ca)
        
    def __save__(self) :
        return "Bridge {0} : ({1}, {2}) has length {3:.5} with pumping {4:.3f}".format(self.index, self.lumen1, self.lumen2, self.length, self.ca)
        
    def __copy__(self) :
        return Osmotic_Bridge(self.index, self.lumen1, self.lumen2, self.length, self.ca)

#