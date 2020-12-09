"""
functions.py contains mathematical functions used essentially for
active pumping profile.

    Contains
    --------
    linear          : Linear profile with slope and offset.
    gauss           : Gaussian profile with average, standard deviation, amplitude and offset.
    rectangular     : Rectangular profile with min, max values.
    sigmoide        : Sigmoide profile with min, max, slope and inflexion point.
    get_parameters  : Extract parameters from configuration.
    calc_func       : Calculate the profile along a chain.
    integrate       : Integrate the profile on a interval [x1, x2].

    Requirements
    ------------
    numpy
    scipy.integrate.quad

Mathieu Le Verge--Serandour, 2020
"""

from numpy import exp, where, ones
from scipy.integrate import quad


def linear(x, slope=1., offset=0.) :
    """
    linear(x, slope=1., offset=0.)
    
        Linear profile of active pumping, such that
            p(x) = slope * x + offset
        
        Parameters
        ----------
        x : float
            Position of the profile along the chain
        slope : float, optional, default : 1.
            Slope of the linear profile
        offset : float, optional, default : 0.
            Offset of the linear profile.
    """
    return slope*x+offset
    
def gauss(x, amp=1., mu=0., sigma=1., fmin=0.) :
    """
    gauss(x, amp=1., mu=0., sigma=1., fmin=0.)
        
        Linear profile of active pumping, such that
            p(x) = amp * exp(-(x-mu)**2/sigma**2) + fmin
        
        Parameters
        ----------
        x : float
            Position of the profile along the chain
        amp : float, optional, default : 1.
            Amplitude of the gaussian profile
        mu : float, optional, default : 0.
            Average of the gaussian profile
        sigma : float, optional, default : 1.
            Standard deviation of the gaussian profile
        fmin : float, optional, default : 0.
            Offset of the gaussian profile
    """
    return amp*exp(-(x-mu)**2/sigma**2) + fmin
    
def rectangular(x, fmax=1., fmin=0., start=0., stop=1.) :
    """
    rectangular(x, fmax=1., fmin=0., start=0., stop=1.)
        
        Rectangular profile of active pumping, such that
            p(x) = fmax, if start <= x <= stop
            p(x) = fmin, otherwise
        
        Parameters
        ----------
        x : float
            Position of the profile along the chain
        fmax : float, optional, default : 1.
            Maximum value of the rectangular function
        fmin : float, optional, default : 0.
            Minimal value of the rectangular function
        start : float, optional, default : 0.
            Position of the left side of the rectangle
        stop : float, optional, default : 1.
            Position of the right side of the rectangle
    """
    return where(abs(x-0.5*(start+stop))<=abs(stop-start)/2., fmax, fmin)
    
def sigmoide(x, fmin=0., fmax=1., slope=50., inflexion_point=0.5) :
    """
    sigmoide(x, fmin=0., fmax=1., slope=50., inflexion_point=0.5)
        
        Sigmoidal profile of active pumping, such that
            p(x) = fmin + (fmax-fmin) / (1+exp(-slope*(x-inflexion_point)))
        
        Parameters
        ----------
        x : float
            Position of the profile along the chain
        fmax : float, optional, default : 1.
            Maximum value of the sigmoid function
        fmin : float, optional, default : 0.
            Maximum value of the sigmoid function
        slope : float, optional, default : 50.
            Slope of the sigmoide function
        inflexion_point : float, optional, default : .5
            Position of the inflexion point along the chain
    """
    return fmin + (fmax-fmin) / (1+exp(-slope*(x-inflexion_point)))
    
def get_parameters(args, func) :
    """
    get_parameters(args, func)
    
        Parameters
        ----------
        args : dict
            Dictionnary that contains the arguments/values of the puming  profile.
            Arguments depends on the profile (default values are indicated below) :
                constant    : {'value' : 1.}
                linear      : {'slope' : 1., 'offset' : 0.}
                gaussian    : {'amp' : 1., 'mu' : 0., 'sigma' : 1., 'fmin' : 0.}
                rectangular : {'fmin' : 0., 'fmax': 1., 'start' : 0., 'stop' : 1.}
                sigmoide    : {'fmin' : 0., 'fmax': 1., 'slope': 10., 'inflexion_point': 0.5}
                
        func : str
            Name of the function (constant, linear, gaussian, rectangular, sigmoide)
    
        Returns
        -------
        Returns the values of the parameters of the profile.
    """
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
        
        try : fmin = args['fmin']
        except : fmin = 0.
        
        return amp, mu, sigma, fmin
                
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
    """
    calc_func(x, func, args, Ltot=0)
    
        Calculate the value of the pumping profile p(x) along the chain
    
        Parameters
        ----------
        x : float or array
            Position(s) of the chain (from 0 to 1.)
        func : string
            Name of the profile function
        args : dict
            Arguments for the pumping profile
        Ltot : float
            Total length of the chain.
    
        Returns
        -------
        p(x) : float or array
            Value of the pumping profile.        
    """
    if func == '' or func == None or len(x) <= 0 :
        return None
        
    elif func == 'constant' :
        value = args
        return value*ones(len(x))
    
    elif func == 'linear' :
        slope, offset = args
        return linear(x, slope, offset)
    
    elif func == 'gaussian' :
        amp, mu, sigma, fmin = args
        return gauss(x, amp, mu, sigma, fmin)
        
    elif func == 'rectangular' :
        fmin, fmax, start, stop = args
        return rectangular(x, fmax, fmin, start, stop)
        
    elif func == 'sigmoide' :
        fmin, fmax, slope, inflexion_point = args
        return sigmoide(x, fmin, fmax, slope, inflexion_point)
 
def integrate(func, x1, x2, args={}) :
    """
    integrate(func, x1, x2, args={})
    
        Integrate the profile p(x) between [x1, x2].
        
        Parameters
        ----------
        func : function
            Function that calculates the pumping profile.
        x1, x2 : floats
            Borders (left, right) of the integral
        args : dict
            Arguments of the function
    
        Returns
        -------
        res : float or array
            Normalized value (integral) of the pumping profile.
    """
    if func == '' or func == None :
        return None
        
    elif func == 'constant' :
        value = get_parameters(args, func)
        res = value*abs(x2-x1)
        
    elif func == 'linear' :
        slope, offset = get_parameters(args, func)
        res, err = quad(func=lambda x, slope=slope, offset=offset : slope*x+offset, a=x1, b=x2)
        
    elif func == 'gaussian' :
        amp, mu, sigma, fmin = get_parameters(args, func)
        res, err = quad(func=lambda x, amp=amp, mu=mu, sigma=sigma, fmin=fmin : amp*exp(-(x-mu)**2/sigma**2)+fmin, a=x1, b=x2)
        
    elif func == 'rectangular' :
        fmin, fmax, start, stop = get_parameters(args, func)
        res, err = quad(func=lambda x, fmin=fmin, fmax=fmax, start=start, stop=stop : where(abs(x-0.5*(start+stop))<=abs(stop-start)/2., fmax, fmin), a=x1, b=x2)
    
    elif func == 'sigmoide' :
        fmin, fmax, slope, inflexion_point = get_parameters(args, func)
        res, err = quad(func=lambda x, fmin=fmin, fmax=fmax, slope=slope, inflexion_point=inflexion_point : fmin + (fmax-fmin) / (1+exp(-slope*(x-inflexion_point))), a=x1, b=x2)
        
    return res / abs(x2-x1)
    
    
    
