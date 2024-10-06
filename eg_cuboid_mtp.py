#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 15:37:50 2024

@author: Jose A. Rodriguez-Rivera
jose.rodriguez@nist.gov
joarodri@gmail.com
University of Maryland
NIST Center for Neutron Research

Calculate the transmssion of a cuboid for a range of theta,twotheta angles.
Returns
-------
None.

Calculate the transmission sequencially and save it in the file test_cuboid.txt'   
For multiprocessign open file cuboid_mpt.py
in a multitasking loop
 
The multitasking loop must always run in the main program. 
"""


import multiprocessing as multip
import cuboid_transm as cb
import abscorr as abcr
import time

if __name__ == '__main__':
    range_theta=[0,90,10]       # Range of theta angles [min, max,steps]
    
    range_twotheta=[0,90,10]    # Range of twotheta angles [min, max,steps]
    angles=abcr.generate_Sample_theta_2theta(range_theta, range_twotheta)
    
    manager = multip.Manager()
    return_list = manager.list()
    file2write='test_cuboid_mtp2.txt'   #file to save
    ff=open(file2write, 'w')
    ff.close()
    kwargs={"file_out":file2write}
    
    
    jobs = []
    startime=time.time()
    for tt in angles:
         p = multip.Process(target=cb.boxsmpltransm, args=(tt[0], tt[1]),kwargs=kwargs)
         jobs.append(p)
         p.start()
    for proc in jobs:
         proc.join()
    print('done.', time.time()-startime)
