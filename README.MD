
Python routines to calculate the neutron and xray transmission for arbitrry shape samples:

By 	Jose. A Rodriguez-Rivera
     	University of Maryland
	NIST Center for Neutron Research
	joarodri@nist.gov
	jrodrig9@umd.edu

If you use these routines please cite the web site and include it in the acknowledgments.

The sample boundaries  must be defined as quadratic surfaces of the form 

F= Ax^2 + By^2 + Cz^2+ Dxy + Eyz + Gzx +Hx +Iy +Jz.

Eg.  A cyiinder with a radius r and length l can be defined with three surfaces:

x^2+y^2=r^2  		A=1 B=1 F=r^2
 z=l/2   			J=1  F=1/2
z=-l/2  			J=2  F=-l/2


programs:
abscorr.py              Absorption correction program.
eg_cuboid_trasm         Example program to calculate the transmission of a cuboid
eg_cuboid_mtp           Example  multiprocess program verison of cuboid_trasm 
eg_sphere_transm.py     Example sphere transmission
eg_cylinder_transm.py   Example cylinder transmission
eg_half_cylinder        Example half transmission
eg_three_samples        Example multiple samples transmission
eg_cylinder_90deg       Example lying cylinder along scattering plane


#########################################################################
Program files:  abscorr.py:
methods in abscorr.py:

cylinder_smpl_transm(twotheta,r,l,mu_ei,mu_ef,**kwargs):
    Calculate the transmission for a cylinder sample giving the radius and
    the linear absorption coefficients.
    theta : float
        theta value
    twotheta : float
        twotheta value.
    r : float
        sphere radius in cm
    l : float
        sphere length in cm    
    mu_ei : float
        inear absorption coefficient for ei in cm^-1
    mu_ef : float
        inear absorption coefficient for ef in cm^-1
    Returns
    absorption : float
        Calculated transmission given a theta and twotheta value for a sample define as a box.

sphere_smpl_transm(twotheta,r,mu_ei,mu_ef,**kwargs):
    Calculate the transmission for a sphere sample giving the radius and
    the linear absorption coefficients.
    theta : float
        theta value
    twotheta : float
        twotheta value.
    r : float
        sphere radius in cm
    mu_ei : float
        inear absorption coefficient for ei in cm^-1
    mu_ef : float
        inear absorption coefficient for ef in cm^-1
    Returns
    absorption : float
        Calculated transmission given a theta and twotheta value for a sample define as a box.
    
def cylinder_sample(r,l)
            Cylinder sample
    r   :   cylinder radius     (x-y plane)
    l   :   cylinder length     (z plane) 

def sphere_sample(r)
            Sphere sample
    r   :   sphere radius
    
def boxsample(lx,ly,lz)
            Cuboid sample.
    lx,ly,lz    :  x, y and z cuboid lenghts
    
def tdis(fu,s,dd):
    This procedure returns [t1,t2] distances using the s vector direction
    at the dd(x,y,z) volume element positions for the Fu defined boundary
    function F= Ax^2 + By^2 + Cz^2+ Dxy + Eyz + Gzx +Hx +Iy +Jz.
    fu  :   Boundary functuon  float(10)
    s   :   scattered beam beam unitary vector float(3)
    dd  :   Volume element position float(3)
    return: [t1,t2] distances  

def todis(fu,s,dd):
    This procedure returns [t1,t2] distances using the so vector direction
    at the dd(x,y,z) volume element positions for the Fu defined boundary
    function F= Ax^2 + By^2 + Cz^2+ Dxy + Eyz + Gzx +Hx +Iy +Jz.
    fu  :   Boundary functuon  float(10)
    s   :   Incident beam unitary vector float(3)
    dd  :   Volume element position float(3)
    return: [t1,t2] distances  


get_t_to_distances(F,dd,s0,s)
    Return t and t0 distances from sample volume element dd to the nearest sample 
    boundary function in the direction of incident beam s0 and scattered beam S    
    returns [tn0,tp0,tn,tp] where:
    tn,tp  Negative and positive distances in the direction of s
    tn0,tp0 Negative and positive distances in the direction of s0.
     F : Sample dictionary
         Set N boundary functions of a sample
     dd : float(3)
          volume element coordinates
     s0 : float(3)
         Incident beam unitary vector
     s : float(3)
         Scattered beam unitary vector
         Returns
         list  [tn0,tp0,tn,tp]
         tn,tp  Negative and positive distances in the direction of s
         tn0,tp0 Negative and positive distances in the direction of s0.
     
    
get_t_t0_distances_VEI_out(smpl,ddl,s0l,sl)
    Return t and t0 distances traveled in other M samples in the direction of
    the incident beam s0 land scattered beam sl from the volume element ddl 
    of the N sample (N!=M)
    returns [tn0,tp0,tn,tp] where:
    tn,tp  Negative and positive distances in the direction of s
    tn0,tp0 Negative and positive distances in the direction of s0.
    smpl : Sample dictionary
        Set N boundary functions of a sample + more sample info
    ddl : float(3)
         volume element coordinates
    s0l : float(3)
        Incident beam unitary vector
    sl : float(3)
        Scattered beam unitary vector
    Returns
        list  [tn0,tp0,tn,tp]
        tn,tp  Negative and positive distances in the direction of s
        tn0,tp0 Negative and positive distances in the direction of s0.


is_inside(sample,s00,ss,dd1)
    This function check if the volume element is inside the sample.
    The routine check if there are boundaries between the dd volume position
    and the origin.
    sample :  Sample dictionary
        Set N boundary functions of a sample + more sample info
    s00 : float(3)
        Incident beam unitary vector
    ss : float(3)
        Scattered beam unitary vector
    dd1 : float(3)
        Volume element coordinates
    Returns
    list  [insid,t0,t]
        insid : 0 if the volume element is outside the sample.
                1 if is inside the sample
        t0    : t0 distance to the nearest boundary function and dd1 in the
                direction of the incident beam s00
        t    : t distance to the nearest boundary function and dd1 in the
                direction of the scattered beam ss

equal_samples(smpl1,smpl2)
    Check if smpl1==smpl2
    smpl1 : Sample dictionary
        Set N boundary functions of a sample + more sample info
    smpl2 : Sample dictionary
        Set N boundary functions of a sample + more sample info
    Returns
    cmpr : TRUE if equal, FALSE if different.

integ(s0,s,samples,mu_Ei,mu_Ef)
    This function calculates the total transmission of a collection of samples.
    given the direction of the incident and scattered beam
    s0 : float(3)
        Unitary vector incident beam
    s : float(3)
        Unitary vector scattered beam
    samples : matrix of sample dictionaries.
    mu_Ei : Ei linear absorption coefficient in cm
    mu_Ef : Ef linear absorption coefficient in cm
    Returns
        Transmission float

srotxy(s,alpha)
    Rotate the s vector around xy plane alpha degrees
   s :  float(3)
       beam direction (unitary vector)
   alpha : float
       rotation in degrees.
    Returns:
       float(3)
           rotated unitary vector

srotyz(s,alpha)
    Rotate the s vector around yz plane alpha degrees
   s :  float(3)
       beam direction (unitary vector)
   alpha : float
       rotation in degrees.
    Returns:
       float(3)
           rotated unitary vector
srotzx(s,alpha)
    Rotate the s vector around zx plane alpha degrees
   s :  float(3)
       beam direction (unitary vector)
   alpha : float
       rotation in degrees.
    Returns:
       float(3)
           rotated unitary vector

smpl_trans(sample,new_center):
    apply transformations to translate the sample and volume element
    New_center is at x0,y0,z0 position
    sample F= Ax^2 + By^2 + Cz^2+ Dxy + Eyz + Gzx +Hx +Iy +Jz = fu(x,y,z)
    A=smpl(1,1), B=smpl(1,2), C=smpl(1,3), D=smpl(1,4), E=smpl(1,5)
    G=smpl(1,6), H=smpl(1,7), I=smpl(1,8), J=smpl(1,9), F=smpl(1,10)
    H_new=H+2Ax0-Dy0-Gz0
    I_new=I+2By0-Dx0-Ez0
    J_new=J+2*Cz0-Ey0-Gx0
    F_new=F-(Ax0^2+By0^2+Cz0^2+Dx0y0+Ez0y0+Gz0x0-Hx0-Iy0-Jz0)

generate_Sample_theta_2theta(range_theta,range_twotheta,**kwargs):
    Generate a list of theta and twotheta angles.
    range_theta : float
        Sample theta range list  [Min, Max, Steps]
    range_twotheta : float
        Sample twotheta range  ist  [Min, Max, Steps]
    **kwargs : File name
        String file name to save the file with generated theta-twotheta angles.
    Return a list of theta-twotheta angles.
    if file='file' the theta-twotheta values will be saved in a file

generate_contour_plot(file2read):

    file2read : ascii file with the transmission data in 3 colums: theta, 2theta, transmission
    Returns
    Plots a countour map of the data.



#########################################################################
Program files:  eg_cuboid_transm.py:
methods in cuboid_transm.py:

boxsmpltransm(theta,twotheta,**kwargs):
    Example program to calculate the transmission of a cuboid
    theta : float
        theta value
    twotheta : float
        twotheta value.
    Returns
    absorption : float
        Calculated transmission given a theta and twotheta value for a sample define as a box.

trans_sample_sequential(angles,**kwargs):
    Calculate the transmssion sequentially giving a list of [theta,twotheta] angles
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
                A list with with three columns [theta, twotheta, transmission]

trans_test():
    Calculate the transmission sequencially and save it in the file test_cuboid.txt'   
    For multiprocessign open file cuboid_mpt.py
    
    
#########################################################################
Program files:  eg_sphere_transm.py and eg_cylinder_transm.py

trans_sample_sequential(angles,**kwargs):
    Calculate the transmission for a given twotheta list of angles
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
                A list with with two columns [twotheta, transmission]   
                
                
#########################################################################
Program files: eg_half_cylinder:    

half_cylinder_sample(r,l):
    Defining a half cylinder sample radius r height l
    4 boundary surfaces, 
    IMPORTANT:
        Boundaries are defined as
        F= Ax^2 + By^2 + Cz^2+ Dxy + Eyz + Gzx +Hx +Iy +Jz = fu(x,y,z)
        The sample must be located at (0,0,0). There is not need to be centered
        The sample must be defined as a set of convex boundaries.
        Cylinder radius r, height l        
        x^2+y^2=r^2  A=1 B=1 F=r^2      boundary eq 1
        z=l/2   J=1  F=1/2              boundary eq 2
        z=-l/2  J=2  F=-l/2             boundary eq 3
        y=0     I=1  F=0                boundary eq 4
        inside value is a defined volume element inside the volume of interest.
        all the lenght units are in cm

half_cylinder_sample_transm(theta,twotheta,r,l,mu_ei,mu_ef,**kwargs):
    Calculate the transmission for a half cylinder given theta and twotheta
    theta=0 when the flat side of the cylider is perpendicular to the incindent
    beam.  s0=[0,1,0],     flat side plane  y=0
    theta : float
        theta value
    twotheta : float
        twotheta value.
    r : float
        sphere radius in cm
    l : float
        sphere length in cm    
    mu_ei : float
        inear absorption coefficient for ei in cm^-1
    mu_ef : float
        inear absorption coefficient for ef in cm^-1
    absorption : fluat
        Calculated transmission given a theta and twotheta value for a sample define as a box.

trans_sample_sequential(angles,**kwargs):
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
                A list with with three columns [theta, twotheta, transmission]  
                
#########################################################################    
Program files: eg_three_samples:

three_samples():
    example of how to define three samples   
    returns a list of 3 dictionaries with the samples information
    cube sample: 2x2x2cm cube located at (2,2,2)
    shpere sample: radius=2cm  located at (-2,-2,-2)
    cylinder sample: raius=1.5cm lenght=1.5cm located at (-2,2,0)

three_samples_transm(theta,twotheta,**kwargs):
    Example program to calculate the transmission of multiple samples
    theta : float
        theta value
    twotheta : fluat
        twotheta value.
    Returns
    absorption : fluat
        Calculated transmission given a theta and twotheta value for a sample define as a box.

trans_sample_sequential(angles,**kwargs):
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
                A list with with three columns [theta, twotheta, transmission]
                
                
#########################################################################    
Program files: eg_cylinder_90deg

cylinder90_smpl_transm(theta,twotheta,**kwargs):
    Transmission of a long cylinder laying in the y-z plane.
    theta : float
        theta value
    twotheta : fluat
        twotheta value.
    Returns
    absorption : float
        Calculated transmission given a theta and twotheta value 
        for a sample define as a cylinder The cylinder height is along k directon
        and the incident beam is also along k direction the scattering plane is y-z.
        The Circumference of the cylinder is defined in the x-y plane
        the height of the cylinder is defined along the z plan.      
       the value o sample theta is 0 when the incident beam s0 
       is pointing along the z direction.
       
trans_sample_sequential(angles,**kwargs):
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
                A list with with three columns [theta, twotheta, transmission]

