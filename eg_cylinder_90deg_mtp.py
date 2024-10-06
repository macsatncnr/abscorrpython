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
# we calculate the absorption of the cylinder rotated 90deg  given the theta-twotheta position.
# The incident beam s0 is along the y axis s0=[0,0,1] s=[0,0,1] for theta=2theta=0


import multiprocessing as multip
import cylinder_90deg as cy90
import abscorr as abcr
import time

if __name__ == '__main__':
    range_theta=[0,360,5]
    range_twotheta=[0,180,5]
    angles=abcr.generate_Sample_theta_2theta(range_theta, range_twotheta)
    
    manager = multip.Manager()
    return_list = manager.list()
    file2write='test_cylinder90_mtp.txt'
    ff=open(file2write, 'w')
    ff.close()
    kwargs={"file_out":file2write}
    
    jobs = []
    startime=time.time()
    for tt in angles:
         p = multip.Process(target=cy90.cylinder90_smpl_transm, args=(tt[0], tt[1]),kwargs=kwargs)
         jobs.append(p)
         p.start()
    for proc in jobs:
         proc.join()
    print('done.', time.time()-startime)

