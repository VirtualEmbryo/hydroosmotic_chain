from numpy import exp, where, ones
from scipy.integrate import quad


def linear(x, slope=1., offset=0.) :
    return slope*x+offset
    
def gauss(x, amp=1., mu=0., sigma=1.) :
    return amp*exp(-(x-mu)**2/sigma**2)
    
def rectangular(x, fmax=1., fmin=0., start=0., stop=1.) :
    """
        rectangular(x, fmax=1., fmin=0., start=0., stop=1.)
    """
    return where(abs(x-0.5*(start+stop))<=abs(stop-start)/2., fmax, fmin)
    
def sigmoide(x, fmin=0., fmax=1., slope=50., inflexion_point=0.5) :
    return fmin + (fmax-fmin) / (1+exp(-slope*(x-inflexion_point)))
    
def get_parameters(args, func) :
    if func == '' or func == None :
        return None
    
    elif func == 'constant' :
        try : value = args['value']
        except : value = 1.
            
        return value
        
    elif func == 'linear' :
        try : slope = args['slope']
        except : slope = 1.
        
        try : offset = args['offset']
        except : offset = 0.
            
        return slope, offset
    
    elif func == 'gaussian' :
        try : amp = args['amp']
        except : amp = 1.
        
        try : mu = args['mu']
        except : mu = 1.
        
        try : sigma = args['sigma']
        except : sigma = 1.
        
        return amp, mu, sigma
        
        
    elif func == 'rectangular' :
        try : fmin = args['fmin']
        except : fmin = 0.
        
        try : fmax = args['fmax']
        except : fmax = 1.
        
        try : start = args['start']
        except : start = 0.
        
        try : stop = args['stop']
        except : stop = 1.
        
        return fmin, fmax, start, stop
    
    elif func == 'sigmoide' :
        try : fmin = args['fmin']
        except : fmin = 0.
        
        try : fmax = args['fmax']
        except : fmax = 1.
        
        try : slope = args['slope']
        except : slope = 10.
        
        try : inflexion_point = args['inflexion_point']
        except : inflexion_point = 0.5
        
        return fmin, fmax, slope, inflexion_point
       
def calc_func(x, func, args, Ltot=0) :
    if func == '' or func == None or len(x) <= 0 :
        return None
        
    elif func == 'constant' :
        value = args
        return value*ones(len(x))
    
    elif func == 'linear' :
        slope, offset = args
        return linear(x, slope, offset)
    
    elif func == 'gaussian' :
        amp, mu, sigma = args
        return gauss(x, amp, mu, sigma)
        
    elif func == 'rectangular' :
        fmin, fmax, start, stop = args
        return rectangular(x, fmax, fmin, start, stop)
        
    elif func == 'sigmoide' :
        fmin, fmax, slope, inflexion_point = args
        return sigmoide(x, fmin, fmax, slope, inflexion_point)
 
def integrate(func, x1, x2, args={}) :
    if func == '' or func == None :
        return None
        
    elif func == 'constant' :
        value = get_parameters(args, func)
        res = value*abs(x2-x1)
        
    elif func == 'linear' :
        slope, offset = get_parameters(args, func)
        res, err = quad(func=lambda x, slope=slope, offset=offset : slope*x+offset, a=x1, b=x2)
        
    elif func == 'gaussian' :
        amp, mu, sigma = get_parameters(args, func)
        res, err = quad(func=lambda x, amp=amp, mu=mu, sigma=sigma : amp*exp(-(x-mu)**2/sigma**2), a=x1, b=x2)
        
    elif func == 'rectangular' :
        fmin, fmax, start, stop = get_parameters(args, func)
        res, err = quad(func=lambda x, fmin=fmin, fmax=fmax, start=start, stop=stop : where(abs(x-0.5*(start+stop))<=abs(stop-start)/2., fmax, fmin), a=x1, b=x2)
    
    elif func == 'sigmoide' :
        fmin, fmax, slope, inflexion_point = get_parameters(args, func)
        res, err = quad(func=lambda x, fmin=fmin, fmax=fmax, slope=slope, inflexion_point=inflexion_point : fmin + (fmax-fmin) / (1+exp(-slope*(x-inflexion_point))), a=x1, b=x2)
        
    return res / abs(x2-x1)
    
    
    
