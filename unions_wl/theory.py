# -*- coding: utf-8 -*-                                                         
                                                                                
"""THEORY MODULE.                                                                 
                                                                                
This module provides theoretical predictions of weak-lensing observables.
                                                                                
""" 

import numpy as np

import pyccl as ccl                                                             
import pyccl.nl_pt as pt                                                        
import pyccl.ccllib as lib


def pk_gm_theo(cosmo, bias_1, log10k_min=-4, log10k_max=2, nk_per_decade=20):
    """PK GM THEO.

    3D galaxy-matter power spectrum

    Parameters
    ----------
    cosmo : pyccl.core.Cosmology
        Cosmological parameter
    bias_1 : float
        linear bias
    log10k_min : float, optional
        minimum 3D Fourier scale (log-10), default=-4
    log10k_max : float, optional
        maximum 3D Fourier scale (log-10), default=2
    nk_per_decade : int
        number of k-modes per log-10  interval in k, default=20

    Returns
    -------
     array_like
        3D power spectrum on a grid in (k, z)

    """                                                                         
    # Tracers                                                                   
    # Galaxies with constant linear bias                                        
    ptt_g = pt.PTNumberCountsTracer(b1=bias_1)                                  
                                                                                
    # Dark matter                                                               
    ptt_m = pt.PTMatterTracer()                                                 

    # Power spectrum pre-computation                                            
    ptc = pt.PTCalculator(
        with_NC=True,
        with_IA=False,
        log10k_min=log10k_min,
        log10k_max=log10k_max,
        nk_per_decade=nk_per_decade,
    )                                                                           
                                                                                
    # 3D galaxy - dark-matter cross power spectrum                              
    pk_gm = pt.get_pt_pk2d(cosmo, ptt_g, tracer2=ptt_m, ptc=ptc)                
                                                                                
    return pk_gm


def gamma_t_theo(                                                               
        theta_deg,                                                              
        cosmo,                                                                  
        dndz_lens,                                                              
        dndz_source,                                                            
        bias_1,                                                                 
        ell=None,                                                               
        p_of_k=None,                                                            
        integr_method='FFTlog',                                                 
):                                                                              
    """GAMMA T THEO.                                                            
                                                                                
    Theoretical prediction of the tangential shear of a source                  
    population around lenses using the ccl library.                             
                                                                                
    Parameters                                                                  
    ----------                                                                  
    theta_deg : array                                                           
        Angular scales in degrees                                               
    cosmo : pyccl.core.Cosmology                                                
        Cosmological parameters                                                 
    dndz_lens : tuple of arrays                                                 
        Lens redshift distribution (z, n(z))                                    
    dndz_source : tuple of arrays                                               
        Source redshift distribution (z, n(z))                                  
    bias_1 : float                                                              
        linear bias                                                             
    ell : array, optional                                                       
        2D Fourier mode, default is                                             
        np.geomspace(2, 10000, 1000)                                            
    p_of_k : array_like, optional                                               
        3D power spectrum on a grid in (k, z). If not given,                    
        the function ``pk_gm_theo`` is called                                   
    integr_method : str, optional                                               
        Method of integration over the Bessel function times                    
        the angular power spectrum, default is 'FFT_log'                        
                                                                                
    Returns                                                                     
    -------                                                                     
    array :                                                                     
        Tangential shear at scales ``theta``                                    
                                                                                
    """                                                                         
    z_lens = dndz_lens[0]                                                       
                                                                                
    # 2D tracers 

    # Galaxies (lenses)                                                         
    bias_g = np.ones_like(z_lens) * bias_1                                      
    tracer_g = ccl.NumberCountsTracer(                                          
            cosmo,                                                              
            False,                                                              
            dndz=dndz_lens,                                                     
            bias=(z_lens, bias_g),                                              
    )                                                                           
                                                                                
    # Weak lensing (sources)                                                    
    n_nz = len(dndz_source[0])                                                  
    tracer_l = ccl.WeakLensingTracer(                                           
        cosmo,                                                                  
        dndz=dndz_source,                                                       
        n_samples=n_nz,                                                         
    )                                                                           
                                                                                
    # Angular cross-power spectrum                                              
    if ell is None:                                                             
        ell_min = 2                                                             
        ell_max = 10000                                                         
        n_ell = 1000                                                            
        ell = np.geomspace(ell_min, ell_max, num=n_ell)                         
                                                                                
    if not p_of_k:                                                              
        pk_gm = pk_gm_theo(cosmo, bias_1)                                       
    else:                                                                       
        pk_gm = p_of_k                                                          
    cls_gG = ccl.angular_cl(cosmo, tracer_g, tracer_l, ell, p_of_k_a=pk_gm)     
                                                                                
    # Tangential shear                                                          
    gt = ccl.correlation(                                                       
        cosmo,                                                                  
        ell,                                                                    
        cls_gG,                                                                 
        theta_deg,                                                              
        type='NG',                                                              
        method=integr_method,                                                   
    )                                                                           
                                                                                
    return gt
