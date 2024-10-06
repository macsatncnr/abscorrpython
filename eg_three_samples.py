#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 10:22:10 2024

@author: yayo
"""

import numpy as np
import abscorr as abcr
import time

def three_samples():
    """
    example of how to define three samples
    
    returns a list of 3 dictionaries with the samples information
    cube sample: 2x2x2cm cube located at (2,2,2)
    shpere sample: radius=2cm  located at (-2,-2,-2)
    cylinder sample: raius=1.5cm lenght=1.5cm located at (-2,2,0)

    """
    
    cubsmpl = abcr.smpl_trans(abcr.boxsample(2.0,2.0,2.0),[2.0,2.0,2.0])
    sphsmpl = abcr.smpl_trans(abcr.sphere_sample(2.0),[-2.0,-2.0,-2.0])
    cylsmpl = abcr.smpl_trans(abcr.cylinder_sample(1.5,1.5),[-2.0,2.0,0.0])
    samples=[cylsmpl, sphsmpl, cubsmpl]
    return samples

def three_samples_transm(theta,twotheta,**kwargs):
    """
    Example program to calculate the transmission of multiple samples
    
    Parameters
    ----------
    theta : float
        theta value
    twotheta : fluat
        twotheta value.

    Returns
    -------
    absorption : fluat
        Calculated transmission given a theta and twotheta value for a sample define as a box.

    """

    s0=[0.0,1.0,0.0]    #   Incident beam
    s=[0.0,1.0,0.0]     #   Scattered beam

    mu_ei=2.0
    mu_ef=2.0
   
    # We must rotate s and s0 -theta degrees and then s twotheta degrees
    # The vectors s and s0 rotates in oposite direction of the sample.
    
    ss0=abcr.srotxy(s0,-1*theta)          # We rotate incident beam s0 along xy axis -theta deg
    ss=abcr.srotxy(abcr.srotxy(s,-1*theta),twotheta)  # We We rotate scattered beam s along xy axis -theta deg and another twotheta deg along xy axis
    
    transm=abcr.integ(ss0, ss, three_samples(),mu_ei,mu_ef)
    print(theta,twotheta,transm)
    if 'file_out' in kwargs:
       file2write=kwargs["file_out"]
       with open(file2write, 'a') as ff:
           np.savetxt(ff,[[theta,twotheta,transm]],fmt='%1.3f',delimiter='\t')
       
    if 'return_list' in kwargs:
        return_list=kwargs["return_list"] 
        return_list.append([theta,twotheta,transm])       
    return transm


def trans_sample_sequential(angles,**kwargs):
    """
    Parameters
    ----------

    angles:  
                angles :    Two columns list [theta, twotheta]
                            List with two columns list of floats
                            [theta, twotheta] values to compute transmission.
                            
                            String
                            if  angles is a "string" then angles is a file name 
                            with the angles theta_twotheta.
                            The file  must have two columns list of floats  [theta, twotheta]
                                                                             
                         
    **kwargs :  file_out :  file to save the calculated transmission with 3 columns:
                            [theta, twotheta, transmission]

    Returns
    -------
                A list with with three columns [theta, twotheta, transmission]
    """ 
    startime=time.time()
    if type(angles)==str:
        file2read=angles
        with open(file2read, 'r') as ff:
            angles=np.loadtxt(ff)
        
    if 'file_out' in kwargs:
       file2write=kwargs["file_out"]
       ff=open(file2write, 'w')
       ff.close()
    
    transmi=[three_samples_transm(angle[0],angle[1],**kwargs) for angle in angles]
      
    print('done.', time.time()-startime)
    return transmi
    
              

def trans_test():
    """
    
    Returns
    -------
    None.
    
    Calculate the transmission sequencially and save it to test_cylinder90.txt'   
    For multiprocessign see file cylinder_90deg_mpt.py
    

    """
    range_theta=[0,30,10]
    range_twotheta=[0,30,10]
    angles=abcr.generate_Sample_theta_2theta(range_theta, range_twotheta)
    print('angles',angles)
    file_out='test_three_smpl.txt'
    trans_sample_sequential(angles,file_out=file_out)
    
    
def plot_transmission():
    abcr.generate_contour_plot("test_three_smpl.txt") 

