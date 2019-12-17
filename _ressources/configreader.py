import sys
import os

module_path = os.path.abspath(os.path.join('..', 'chain_lumen/'))

if module_path not in sys.path :
    sys.path.append(module_path)

class Config() :
    def __init__(self) :
        self.config = {}
        self.categories = []
    
    def read_file(self, filename) :
        f = open(filename, 'r')
        s = f.readlines()
        f.close()
        return s
    

    def read_config(self, string) :
        config = {}
        cat_list = []
        for line in string :
            if line != '\n' :
                if line[0].startswith('[') and line[-2].startswith(']') :
                    cat = line[1:-2]
                    config[cat] = {}
                    cat_list += [str(cat)]
                else :
                    res = line.split('=')
                    arg = res[0].replace(" ", "")
                    val = res[1].replace("\n", "").replace(" ", "")
                    config[cat][arg] = val
        return config, cat_list

    def read(self, filename) :
        string_file = self.read_file(filename)
        self.config, self.categories = self.read_config(string_file)
        return self.config
    
    def write(self, filename) :
        f = open(filename, 'w')
        for key in self.config.keys() :
            f.write('['+str(key)+']\n')
            for sub_key in self.config[key].keys() :
                f.write(str(sub_key) + ' = ' + str(self.config[key][sub_key]) + '\n')
            f.write('\n')
        f.close()
    
    def add(self, category, name, value) :
        if category in self.categories :
            self.config[category][name] = str(value)
        else :
            self.categories += [str(name)]
            self.config[category] = {}
            self.config[category][name] = str(value)
            
    def __str__(self) :
        for key in self.config.keys() :
            print('['+str(key)+']')
            for sub_key in self.config[key].keys() :
                print(str(sub_key) + ' = ' + str(self.config[key][sub_key]) )
            print('')
        return ''