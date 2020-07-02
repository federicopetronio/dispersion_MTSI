

import numpy as np



def generalized_Maxw_f(v, a,b,c,d):
    """general distribution function used for the fits"""
    return a * np.exp( - (v**2)**(c/2)/d)


def maxw_f(v, vt):
    return np.exp(- (v/vt)**2)/np.sqrt(np.pi)/vt
