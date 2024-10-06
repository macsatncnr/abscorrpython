
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

def cylinder90_smpl_transm(theta,twotheta,**kwargs):
    """
    Transmission of a long cylinder laying in the y-z plane.

    Parameters
    ----------
    theta : float
        theta value
    twotheta : fluat
        twotheta value.

    Returns
    -------
    absorption : float
        Calculated transmission given a theta and twotheta value 
        for a sample define as a cylinder The cylinder height is along k directon
        and the incident beam is also along k direction the scattering plane is y-z.
        
        
        
        The Circumference of the cylinder is defined in the x-y plane
        the height of the cylinder is defined along the z plan.
        
       the value o sample theta is 0 when the incident beam s0 
       is pointing along the z direction.
        
    """

    s0=[0.0,0.0,1.0]    #   Incident beam
    s=[0.0,0.0,1.0]     #   Scattered beam
    r=0.5
    l=4.0
    mu_ei=4.0
    mu_ef=4.0
    sphere_smpl=abcr.cylinder_sample(r,l)
    ss0=abcr.srotyz(s0,-1*theta)
    ss=abcr.srotyz(abcr.srotyz(s,-1*theta),twotheta)
    transm=abcr.integ(ss0, ss, sphere_smpl,mu_ei,mu_ef)
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
    
    transmi=[cylinder90_smpl_transm(angle[0],angle[1],**kwargs) for angle in angles]
      
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
    file_out='test_cylinder90.txt'
    trans_sample_sequential(angles,file_out=file_out)
    
    
def plot_transmission():
    abcr.generate_contour_plot("test_cylinder90.txt") 
    
    