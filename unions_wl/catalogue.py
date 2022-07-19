# -*- coding: utf-8 -*-                                                         
                                                                                
"""CATALOGUE MODULE.                                                                 
                                                                                
This module provides functions to handle samples of objects and
catalogues.

"""                                        

import numpy as np


def y_equi(cdf, n):
    """ Y EQUI
    
    Split sample into n equi-populated bins, return bin boundaries.
    
    Parameters
    ----------
    cdf : statsmodels.distributions.empirical_distribution.ECDF
        distribution data
    n : int
        number of splits
        
    Returns
    -------
    list
        bin boundaries of samples
        
    """
    x_list = []
    for denom in range(1, n):
        idx =  np.where(cdf.y >= denom/n)[0][0]
        x_list.append(cdf.x[idx])
    
    return x_list


def bin_edges2centers(bin_edges):                                               
    """BIN EDGES TO CENTERS                                                     
                                                                                
    Transform bin edge values to central values                                 
                                                                                
    Parameters                                                                  
    ----------                                                                  
    bin_edges : list                                                            
        bin edge values                                                         
                                                                                
    Returns                                                                     
    -------                                                                     
    list                                                                        
        bin central values                                                      
                                                                                
    """                                                                         
    bin_means = 0.5 * (                                                         
        bin_edges[1:] + bin_edges[:-1]                                          
    )                                                                           
                                                                                
    return bin_means 
