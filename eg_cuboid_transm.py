#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 11:04:10 2024

@author: Jose A. Rodriguez-Rivera
jose.rodriguez@nist.gov
joarodri@gmail.com
University of Maryland
NIST Center for Neutron Research
"""
# we calculate the absorption of the box given the theta-twotheta position.
# The incident beam s0 is along the y axis s0=[0,1,0] s=[0,1,0] for theta=2theta=0


import numpy as np
import abscorr as abcr
import time



def boxsmpltransm(theta,twotheta,**kwargs):
    """
    Example program to calculate the transmission of a cuboid
    
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
    mu_ei=4.0           #   Linear absorption factor for the incident beam in cm^-1
    mu_ef=4.0           #   Linear absorption factor for the scattered beam in cm^-1
    lx=2.0              #   cuboid width lenght and height
    ly=0.5
    lz=2.0
    samplebox=abcr.boxsample(lx,ly,lz)
    ss0=abcr.srotxy(s0,-1*theta)
    ss=abcr.srotxy(abcr.srotxy(s,-1*theta),twotheta)
    transm=abcr.integ(ss0, ss, samplebox,mu_ei,mu_ef)
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
    Calculate the transmssion sequentially giving a list of [theta,twotheta] angles
    
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
    
    transmi=[boxsmpltransm(angle[0],angle[1],**kwargs) for angle in angles]
        
    print('done.', time.time()-startime)
    return transmi
    
              

def trans_test():
    """
    Calculate the transmission sequencially and save it in the file test_cuboid.txt'   
    For multiprocessign open file cuboid_mpt.py
    """
    range_theta=[0,10,5]        # Range of theta angles [min, max, steps]
    range_twotheta=[0,10,5]     # Range of twotheta angles [min, max, steps]
    angles=abcr.generate_Sample_theta_2theta(range_theta, range_twotheta)
    file_out='test_cuboid.txt'  # file to save thge results
    trans_sample_sequential(angles,file_out=file_out)
    
    
def plot_transmission():
    abcr.generate_contour_plot("test_cuboid_mtp2.txt")    

    

