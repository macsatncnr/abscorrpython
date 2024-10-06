#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 29 14:13:54 2024

@author: Jose A. Rodriguez-Rivera
jose.rodriguez@nist.gov
joarodri@gmail.com
University of Maryland
NIST Center for Neutron Research
"""

import numpy as np
import abscorr as abcr
import time


def single_angle_trans():
    angle=35     # twotheta value in degress
    r=1.0           # cylinder radius in cm
    mu_ei=2.0       # mu ei in cm^-1
    mu_ef=2.0       # mu ef in cm^-1
    print(abcr.sphere_smpl_transm(angle,r,mu_ei,mu_ef))

def trans_sample_sequential(angles,**kwargs):
    
    """
    Calculate the transmission for a given twotheta list of angles
    
    Parameters
    ----------

    angles:  
                angles :    One column list [twotheta]
                            One columns list of floats
                            [theta, twotheta] values to compute transmission.
                            
                            String
                            if  angles is a "string" then angles is a file name 
                            with the angles theta_twotheta.
                            The file  must have two columns list of floats  [theta, twotheta]
                                                                             
                         
    **kwargs :  file_out :  file to save the calculated transmission with 3 columns:
                            [theta, twotheta, transmission]

    Returns
    -------
                A list with with two columns [twotheta, transmission]            
    """ 
    r=1.0               #   sphere radius in cm
    mu_ei=2.0           #   sphere mu_ei in cm^-1
    mu_ef=2.0           #   sphere mu_ef in cm^-1
    startime=time.time()
    if type(angles)==str:
        file2read=angles
        ff=open(file2read, 'r')
        angles=np.loadtxt(ff)
        ff.close()
        
    if 'file_out' in kwargs:
       file2write=kwargs["file_out"]
       ff=open(file2write, 'w')
       ff.close()
    
    transmi=[abcr.sphere_smpl_transm(angle,r,mu_ei,mu_ef,**kwargs) for angle in angles]
        
    print('done.', time.time()-startime)
    return transmi


def trans_test():
    """
    
    Returns
    -------
    None.
    
    Calculate the transmission sequencially and save it in the file test_cuboid.txt'   
    For multiprocessign open file cuboid_mpt.py
    

    """
    angles_twotheta=np.arange(0,190,10)
    file_out='test_sphere.txt'
    trans_sample_sequential(angles_twotheta,file_out=file_out)
    
     