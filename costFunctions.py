import numpy as np


def kinematic_cost_function(u_FEM, u_DIC, weights=None):
    """
    Compute the cost kinematic cost function from displacement fields 
    
    Parameters
    ----------
    u_DIC : numpy array
        Array containing the displacements, given by DIC measurements
    u_FEM : numpy array
        Array containing the displacements, given by FEM results
    w : numpy array
        Weights on displacement error function 

    Returns
    -------
    float
        Cost function
    """
    delta_u = u_DIC-u_FEM[:, :2]
    if weights is None:
        weights=np.ones(len(delta_u))
    else:
        weights=1/weights
    K = np.dot(np.sum(u_DIC**2, axis=1), weights)
    return np.dot(weights, np.sum(delta_u**2, axis=1))/K    

def static_cost_function(eps_exp, sigma_exp, eps_FEM, sigma_FEM):
    """
    Compute the static cost function, by comparing two tensile curves.

    Parameters
    ----------
    eps_exp : numpy array
        Elongation measured experimentally
    sigma_exp : TYPE
        True stress measured experimentally
    eps_FEM : TYPE
        Simulated elongation
    sigma_FEM : TYPE
        Simulated tensile stress

    Returns
    -------
    float
        Cost function

    """
    sigma_exp_interp = np.interp(eps_FEM, eps_exp, sigma_exp)
    K = np.sum(sigma_exp_interp**2)
    delta_sigma = sigma_FEM-sigma_exp_interp
    return np.sum(delta_sigma**2)/K

def weighted_cost_function(f1, f2, w1=0.5):
    """
    Compute the weighted mixture of two cost functions

    Parameters
    ----------
    f1 : float
        1st cost function
    f2 : float
        2nd cost function
    w1 : float, optional
        Weight associated to f1. The weight associated to f2 is 1-w1. 
        The default is None.

    Returns
    -------
    float
        Weighted mixture of cost functions

    """
    return w1*f1 + (1-w1)*f2