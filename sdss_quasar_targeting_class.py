#!/usr/bin/env python

'''
DESCRIPTION:

A class containing functions that will determine if an object(s) would be targeted by the Sloan Digital Sky Survey (SDSS) QSO targeting algorithim (Richards et al 2002) given ugriz photometry and PSF fraction. Compatible with observed or simulated data.


CONTAINS:

1) not_in_locus(): determines if a source lies > 4 \sigma from the stellar locus in 3D color-color-color space (ugri or  griz)

2) is_unresolved(): determines if a source is resloved in sdss imaging based on its calculated PSF fraction

3) color_mag_cuts_sdss(): determines if a source lies in the allowed i-band magnitude range and if it passes the specified color cuts

4) is_targeted(): determines if a source passes all of the criteria necessary to have been targeted my the SDSS quasar targeting algorithm 


USAGE:

sdssQSO = sdss_quasar_selection(psf_mags, psf_frac, locus_name=None)


INPUT:

1) psf_mags: (5) x ngal array containing the ugriz(ugriz) apparent PSF magnitudes for each galaxy in sample.

2) psf_frac: 1D, length ngal array containing the i-band PSF fraction for each galaxy in sample.

3) locus_name: String. The particular color cube being used. Options are 'ugri' (low-z) or 'griz' (high-z). Default value is 'ugri'.


OUTPUT:

All functions in here output 1D, length ngal boolean arrays. True means would be targeted.

'''



class sdss_quasar_selection:
    
    
    
    
    def __init__(self, psf_mags, psf_frac, locus_name=None):
        
        self.psf_mags = psf_mags
        self.psf_frac = psf_frac
        
        if locus_name==None:
            self.locus_name = 'ugri'
        else:
            self.locus_name = locus_name
            
    
    
    
    
    
    
    def not_in_locus(self):

        ### import all necessary packages
        from astropy.cosmology import FlatLambdaCDM
        import numpy as np
        
        ### unpack psf_mags
        u_psf_mag = self.psf_mags[0]
        g_psf_mag = self.psf_mags[1]
        r_psf_mag = self.psf_mags[2]
        i_psf_mag = self.psf_mags[3]
        z_psf_mag = self.psf_mags[4]
        #bol_psf_mag = self.psf_mags[5]

        ### set colors
        ug_color = u_psf_mag - g_psf_mag
        gr_color = g_psf_mag - r_psf_mag
        ri_color = r_psf_mag - i_psf_mag
        iz_color = i_psf_mag - z_psf_mag
    
        if self.locus_name == 'ugri':

            ### import richards et al table 3 values
            ug_locus = np.array([0.855, 1.002, 1.136, 1.262, 1.382, 1.499, 1.609, 1.717, 1.808, 1.904, 2.007, 2.117, 2.234, 2.361, 2.478, 2.518, 2.510])
            gr_locus = np.array([0.259, 0.344, 0.410, 0.466, 0.517, 0.565, 0.611, 0.655, 0.700, 0.748, 0.802, 0.866, 0.945, 1.047, 1.191, 1.327, 1.355])    
            ri_locus = np.array([0.094, 0.126, 0.150, 0.170, 0.189, 0.207, 0.223, 0.238, 0.255, 0.273, 0.293, 0.317, 0.351, 0.403, 0.502, 0.707, 1.068])

            mjax = np.array([0.282, 0.247, 0.221, 0.219, 0.216, 0.217, 0.224, 0.227, 0.233, 0.248, 0.266, 0.278, 0.309, 0.382, 0.463, 0.484, 0.569])
            mnax = np.array([0.135, 0.129, 0.124, 0.126, 0.125, 0.129, 0.131, 0.127, 0.132, 0.129, 0.134, 0.136, 0.136, 0.145, 0.156, 0.180, 0.212])

            theta = np.array([-1.165, -1.147, -1.075, -1.026, -0.977, -0.983, -0.986, -0.989, -1.040, -1.002, -1.017, -1.023, -1.033, -1.051, -1.108, -1.244, -1.669])
            theta_deg = theta*180/np.pi

            k_ug = np.array([0.851, 0.869, 0.893, 0.907, 0.915, 0.915, 0.912, 0.904, 0.888, 0.878, 0.860, 0.827, 0.773, 0.646, 0.355, 0.053, -0.22])        
            k_gr = np.array([0.492, 0.467, 0.422, 0.396, 0.379, 0.379, 0.387, 0.402, 0.430, 0.449, 0.478, 0.521, 0.573, 0.650, 0.634, 0.278, 0.076])
            k_ri = np.array([0.182, 0.172, 0.154, 0.145, 0.140, 0.135, 0.132, 0.146, 0.165, 0.166, 0.178, 0.213, 0.271, 0.400, 0.688, 0.959, 0.997])      

            x_locus = ug_locus
            y_locus = gr_locus
            z_locus = ri_locus

            k_x = k_ug
            k_y = k_gr
            k_z = k_ri
            
            data = np.array([ug_color, gr_color, ri_color])


        elif self.locus_name == 'griz': 
            ###import richards et al table 4
            gr_locus = np.array([0.204, 0.304, 0.382, 0.454, 0.525, 0.594, 0.659, 0.723, 0.787, 0.853, 0.922, 0.991, 1.063, 1.132, 1.202, 1.262, 1.313, 1.343, 1.355, 1.352, 1.347, 1.350, 1.361])
            ri_locus = np.array([0.071, 0.110, 0.137, 0.166, 0.194, 0.219, 0.242, 0.265, 0.288, 0.313, 0.341, 0.371, 0.409, 0.454, 0.507, 0.569, 0.651, 0.754, 0.874, 0.996, 1.116, 1.240, 1.385])
            iz_locus = np.array([0.003, 0.027, 0.044, 0.066, 0.087, 0.105, 0.123, 0.140, 0.155, 0.171, 0.188, 0.206, 0.227, 0.251, 0.280, 0.314, 0.356, 0.408, 0.465, 0.525, 0.583, 0.646, 0.729])      

            mjax = np.array([0.207, 0.165, 0.154, 0.159, 0.164, 0.162, 0.151, 0.150, 0.153, 0.157, 0.160, 0.163, 0.171, 0.175, 0.178, 0.185, 0.193, 0.213, 0.246, 0.250, 0.265, 0.246, 0.300])
            mnax = np.array([0.146, 0.126, 0.128, 0.134, 0.133, 0.133, 0.133, 0.127, 0.128, 0.124, 0.125, 0.123, 0.123, 0.125, 0.127, 0.135, 0.129, 0.131, 0.137, 0.135, 0.133, 0.121, 0.139])

            theta = np.array([0.067, -2.907, -2.990, -0.029, -0.194, -0.315, -0.610, -0.858, -0.917, -0.921, -0.898, -0.949, -1.033, -1.127, -1.323, -1.423, -1.554, -1.628, -1.667, -1.647, -1.652, -1.530])     
            theta_deg = theta*180/np.pi

            k_gr = np.array([0.911, 0.916, 0.910, 0.895, 0.905, 0.913, 0.911, 0.915, 0.916, 0.906, 0.897, 0.876, 0.832, 0.778, 0.704, 0.566, 0.362, 0.168, 0.035, -0.031, -0.008, 0.047, 0.067])
            k_ri = np.array([0.351, 0.339, 0.340, 0.356, 0.342, 0.325, 0.330, 0.332, 0.335, 0.360, 0.380, 0.420, 0.485, 0.551, 0.626, 0.729, 0.832, 0.885, 0.900, 0.899, 0.895, 0.879, 0.868])
            k_iz = np.array([0.218, 0.213, 0.237, 0.268, 0.253, 0.246, 0.246, 0.231, 0.220, 0.222, 0.227, 0.237, 0.267, 0.301, 0.342, 0.386, 0.420, 0.434, 0.435, 0.438, 0.446, 0.475, 0.493])

            x_locus = gr_locus
            y_locus = ri_locus
            z_locus = iz_locus

            k_x = k_gr
            k_y = k_ri
            k_z = k_iz
            
            data = np.array([gr_color, ri_color, iz_color])

        else: 
            print('pick a valid locus')

        ### set up unit vectors
        loci = np.transpose(np.array([x_locus, y_locus, z_locus]))

        k_vecs = np.transpose(np.array([k_x, k_y, k_z])) ### roughly unit vectors
        j_vecs = np.transpose(np.array([k_y, -1*k_x, np.zeros(len(k_y))])/(np.sqrt(k_y**2 + k_x**2))) ### from Newberg & Yanny
        i_vecs = np.transpose(np.array([-1*k_x*k_z, -1*k_y*k_z, (k_x**2 + k_y**2)])/(np.sqrt(k_y**2 + k_x**2)))

        theta_matrix = np.ones(np.shape(i_vecs))
        mjax_matrix = np.ones(np.shape(i_vecs))

        for i in range(len(theta)):
            theta_matrix[i] = theta_matrix[i]*theta[i] ### so we can broadcast with j and i vecs

        l_vecs = (i_vecs*np.cos(theta_matrix) + j_vecs*np.sin(theta_matrix))
        m_vecs = -1*i_vecs*np.sin(theta_matrix) + j_vecs*np.cos(theta_matrix)

        p_vecs = np.zeros(np.shape(k_vecs))
        p_mdpt = np.zeros(np.shape(k_vecs))

        for j in range(len(p_vecs)):
            if ((j!=0) & (j!=(len(p_vecs))-1)):
                p_vecs[j] = loci[j+1] - loci[j]
                p_mdpt[j] = loci[j] + p_vecs[j]/2

        ### this is the bulk of the code

        locusFlag = np.zeros(len(np.transpose(data)), dtype=bool) ### to initialize a while loop

        for i in range(len(k_vecs)):

            ### Determine if data points are in infinite cylinder

            ### rotate k along the z axis into the x-y plane
            rot_z_hyp = np.sqrt(k_vecs[i][0]**2 + k_vecs[i][1]**2)
            rot_z = np.array([[k_vecs[i][1]/rot_z_hyp, -1*k_vecs[i][0]/rot_z_hyp, 0], [k_vecs[i][0]/rot_z_hyp, k_vecs[i][1]/rot_z_hyp, 0], [0, 0, 1]])

            k_rotated_z = np.matmul(rot_z, k_vecs[i])
            l_rotated_z = np.matmul(rot_z, l_vecs[i])
            m_rotated_z = np.matmul(rot_z, m_vecs[i])

            ### rotate k along the x axis into the z-axis
            rot_x_hyp = np.sqrt(k_vecs[i][2]**2 + k_vecs[i][1]**2)
            rot_x = np.array([[1, 0, 0], [0, k_rotated_z[2]/rot_x_hyp, -1*k_rotated_z[1]/rot_x_hyp], [0, k_rotated_z[1]/rot_x_hyp, k_rotated_z[2]/rot_x_hyp]])

            k_rotated = np.matmul(rot_x, k_rotated_z)
            l_rotated = np.matmul(rot_x, l_rotated_z)
            m_rotated = np.matmul(rot_x, m_rotated_z)

            ### rotate m into the x axis
            rot_final_hyp = np.sqrt(m_rotated[0]**2 + m_rotated[1]**2)
            rot_final = np.array([[m_rotated[0]/rot_final_hyp, m_rotated[1]/rot_final_hyp, 0], [-1*m_rotated[1]/rot_final_hyp, m_rotated[0]/rot_final_hyp, 0], [0,0,1]])

            ### rotate all of the relevant vectors relative to the ellipse using these transformations

            lvec_final = np.matmul(rot_final, l_rotated)/np.linalg.norm(np.matmul(rot_final, l_rotated))
            mvec_final = np.matmul(rot_final, m_rotated)/np.linalg.norm(np.matmul(rot_final, m_rotated))
            kvec_final = np.matmul(rot_final, k_rotated)/np.linalg.norm(np.matmul(rot_final, k_rotated))

            ### want to translate data based on central locus point
            translation_mat = -1*loci[i]
            data_trans = np.array([data[0]+translation_mat[0], data[1]+translation_mat[1], data[2]+translation_mat[2]])

            ### rotate the translated data
            data_rotated = np.matmul(rot_final, np.matmul(rot_x , np.matmul(rot_z, data_trans)))

            ### determine if the data falls within the radius of ellipse 
            rad_cc = (data_rotated[0]**2/mnax[i]**2) + (data_rotated[1]**2/mjax[i]**2)
            ellipseFlag = rad_cc < 1.


            ### Cut off each cylinder at the plane specified in Richards et al.

            if ((i!=0) & (i!=(len(k_vecs))-1)): ### this is for the body of the structure- we need to treat the ends seperately
                ### for our specific project, none of our sources would fall in the end region so we will update later
                
                ### bottom part of cylinder
                ### rotation around y axis to get p into the y-z plane
                b_hyp = np.sqrt(p_vecs[i][2]**2 + p_vecs[i][0]**2)
                bottom_rot1 = np.array([[p_vecs[i][2]/b_hyp, 0, -1*p_vecs[i][0]/b_hyp], [0,1,0], [p_vecs[i][0]/b_hyp, 0, p_vecs[i][2]/b_hyp]])

                rotated_pvec_i = np.matmul(bottom_rot1, p_vecs[i])

                ### rotation around x axis to get p into the z axis
                b_hyp2 = np.sqrt(p_vecs[i][2]**2 + p_vecs[i][1]**2)
                bottom_rot2 = np.array([[1,0,0], [0, rotated_pvec_i [2]/b_hyp2, -1*rotated_pvec_i [1]/b_hyp2], [0,rotated_pvec_i [1]/b_hyp2,rotated_pvec_i [2]/b_hyp2]])

                ### translation from p midpoint to origin
                bottom_trans = -1*p_mdpt[i]


                ### top part of cylinder

                ### rotation around y axis to get p into the y-z plane
                t_hyp = np.sqrt(p_vecs[i+1][2]**2 + p_vecs[i+1][0]**2)
                top_rot1 = np.array([[p_vecs[i+1][2]/t_hyp, 0, -1*p_vecs[i+1][0]/t_hyp], [0,1,0], [p_vecs[i+1][0]/t_hyp, 0, p_vecs[i+1][2]/t_hyp]])

                rotated_pvec_i1 = np.matmul(top_rot1, p_vecs[i+1])

                ### rotation around x axis ro get p into the z axis
                t_hyp2 = np.sqrt(p_vecs[i+1][2]**2 + p_vecs[i+1][1]**2)
                top_rot2 = np.array([[1,0,0], [0, rotated_pvec_i1 [2]/t_hyp2, -1*rotated_pvec_i1 [1]/t_hyp2], [0,rotated_pvec_i1 [1]/t_hyp2,rotated_pvec_i1 [2]/t_hyp2]])

                ### translation from p midpoint to origin
                top_trans = -1*p_mdpt[i+1]

                ### transform the data to determine planes for cylinders
                data_bot_trans = np.array([data[0]+bottom_trans[0], data[1]+bottom_trans[1], data[2]+bottom_trans[2]])
                data_top_trans = np.array([data[0]+top_trans[0], data[1]+top_trans[1], data[2]+top_trans[2]])

                data_bot = np.matmul(bottom_rot2, np.matmul(bottom_rot1, (data_bot_trans)))
                data_top = np.matmul(top_rot2, np.matmul(top_rot1, (data_top_trans)))

                
                planeFlag = ((data_bot[2]>0)*(data_top[2]<0)) ### plane flag and locus flag can be taken out of the if statement after I define the ends

                locusFlag = np.logical_or(locusFlag, (planeFlag*ellipseFlag))
                not_in_locusFlag = np.logical_not(locusFlag)

        return not_in_locusFlag
    
    
    
    
    
    
    
    def is_unresolved(self): 

        ### import all necessary packages
        from astropy.cosmology import FlatLambdaCDM
        import numpy as np
        from matplotlib import pyplot as plt
        
        delta_mag = np.abs(-2.5*np.log10(1/self.psf_frac))
        unresolved_Flag = delta_mag < 0.145

        return unresolved_Flag
      
    
    
    
    
    
    def color_mag_cuts_sdss(self):

        ### import all necessary packages
        from astropy.cosmology import FlatLambdaCDM
        import numpy as np
        from matplotlib import pyplot as plt
        
        ### unpack psf mags
        u_psf_mag = self.psf_mags[0]
        g_psf_mag = self.psf_mags[1]
        r_psf_mag = self.psf_mags[2]
        i_psf_mag = self.psf_mags[3]
        z_psf_mag = self.psf_mags[4]

        ### set colors
        ug_color = u_psf_mag - g_psf_mag
        gr_color = g_psf_mag - r_psf_mag
        ri_color = r_psf_mag - i_psf_mag
        iz_color = i_psf_mag - z_psf_mag

        ### define extended sources
        unresolved_Flag = self.is_unresolved()

        extendedFlag = np.logical_not(unresolved_Flag)

        ### low redshift sources
        if self.locus_name == 'ugri':

            ### okay I'm going to have to rotate things again because there are l and k cuts.


            ### cut and paste part of the in_locus() function

            ### import richards et al table 3 values
            ug_locus = np.array([0.855, 1.002, 1.136, 1.262, 1.382, 1.499, 1.609, 1.717, 1.808, 1.904, 2.007, 2.117, 2.234, 2.361, 2.478, 2.518, 2.510])
            gr_locus = np.array([0.259, 0.344, 0.410, 0.466, 0.517, 0.565, 0.611, 0.655, 0.700, 0.748, 0.802, 0.866, 0.945, 1.047, 1.191, 1.327, 1.355])    
            ri_locus = np.array([0.094, 0.126, 0.150, 0.170, 0.189, 0.207, 0.223, 0.238, 0.255, 0.273, 0.293, 0.317, 0.351, 0.403, 0.502, 0.707, 1.068])

            mjax = np.array([0.282, 0.247, 0.221, 0.219, 0.216, 0.217, 0.224, 0.227, 0.233, 0.248, 0.266, 0.278, 0.309, 0.382, 0.463, 0.484, 0.569])
            mnax = np.array([0.135, 0.129, 0.124, 0.126, 0.125, 0.129, 0.131, 0.127, 0.132, 0.129, 0.134, 0.136, 0.136, 0.145, 0.156, 0.180, 0.212])

            theta = np.array([-1.165, -1.147, -1.075, -1.026, -0.977, -0.983, -0.986, -0.989, -1.040, -1.002, -1.017, -1.023, -1.033, -1.051, -1.108, -1.244, -1.669])
            theta_deg = theta*180/np.pi

            k_ug = np.array([0.851, 0.869, 0.893, 0.907, 0.915, 0.915, 0.912, 0.904, 0.888, 0.878, 0.860, 0.827, 0.773, 0.646, 0.355, 0.053, -0.22])        
            k_gr = np.array([0.492, 0.467, 0.422, 0.396, 0.379, 0.379, 0.387, 0.402, 0.430, 0.449, 0.478, 0.521, 0.573, 0.650, 0.634, 0.278, 0.076])
            k_ri = np.array([0.182, 0.172, 0.154, 0.145, 0.140, 0.135, 0.132, 0.146, 0.165, 0.166, 0.178, 0.213, 0.271, 0.400, 0.688, 0.959, 0.997])      

            x_locus = ug_locus
            y_locus = gr_locus
            z_locus = ri_locus

            k_x = k_ug
            k_y = k_gr
            k_z = k_ri

            ### make sure the data is the right colors 
            data = np.array([ug_color, gr_color, ri_color])

        elif self.locus_name == 'griz':

            ### set high-z magnitude limits 
            magFlag = (i_psf_mag > 15) * (i_psf_mag < 20.2)


            ###import richards et al table 4
            gr_locus = np.array([0.204, 0.304, 0.382, 0.454, 0.525, 0.594, 0.659, 0.723, 0.787, 0.853, 0.922, 0.991, 1.063, 1.132, 1.202, 1.262, 1.313, 1.343, 1.355, 1.352, 1.347, 1.350, 1.361])
            ri_locus = np.array([0.071, 0.110, 0.137, 0.166, 0.194, 0.219, 0.242, 0.265, 0.288, 0.313, 0.341, 0.371, 0.409, 0.454, 0.507, 0.569, 0.651, 0.754, 0.874, 0.996, 1.116, 1.240, 1.385])
            iz_locus = np.array([0.003, 0.027, 0.044, 0.066, 0.087, 0.105, 0.123, 0.140, 0.155, 0.171, 0.188, 0.206, 0.227, 0.251, 0.280, 0.314, 0.356, 0.408, 0.465, 0.525, 0.583, 0.646, 0.729])      

            mjax = np.array([0.207, 0.165, 0.154, 0.159, 0.164, 0.162, 0.151, 0.150, 0.153, 0.157, 0.160, 0.163, 0.171, 0.175, 0.178, 0.185, 0.193, 0.213, 0.246, 0.250, 0.265, 0.246, 0.300])
            mnax = np.array([0.146, 0.126, 0.128, 0.134, 0.133, 0.133, 0.133, 0.127, 0.128, 0.124, 0.125, 0.123, 0.123, 0.125, 0.127, 0.135, 0.129, 0.131, 0.137, 0.135, 0.133, 0.121, 0.139])

            theta = np.array([0.067, -2.907, -2.990, -0.029, -0.194, -0.315, -0.610, -0.858, -0.917, -0.921, -0.898, -0.949, -1.033, -1.127, -1.323, -1.423, -1.554, -1.628, -1.667, -1.647, -1.652, -1.530])     
            theta_deg = theta*180/np.pi

            k_gr = np.array([0.911, 0.916, 0.910, 0.895, 0.905, 0.913, 0.911, 0.915, 0.916, 0.906, 0.897, 0.876, 0.832, 0.778, 0.704, 0.566, 0.362, 0.168, 0.035, -0.031, -0.008, 0.047, 0.067])
            k_ri = np.array([0.351, 0.339, 0.340, 0.356, 0.342, 0.325, 0.330, 0.332, 0.335, 0.360, 0.380, 0.420, 0.485, 0.551, 0.626, 0.729, 0.832, 0.885, 0.900, 0.899, 0.895, 0.879, 0.868])
            k_iz = np.array([0.218, 0.213, 0.237, 0.268, 0.253, 0.246, 0.246, 0.231, 0.220, 0.222, 0.227, 0.237, 0.267, 0.301, 0.342, 0.386, 0.420, 0.434, 0.435, 0.438, 0.446, 0.475, 0.493])

            x_locus = gr_locus
            y_locus = ri_locus
            z_locus = iz_locus

            k_x = k_gr
            k_y = k_ri
            k_z = k_iz

            data = np.array([gr_color, ri_color, iz_color])

        ### do all the transformations to the data to get l and k values for cuts
        ### set up unit vectors
        loci = np.transpose(np.array([x_locus, y_locus, z_locus]))

        k_vecs = np.transpose(np.array([k_x, k_y, k_z])) ### roughly unit vectors
        j_vecs = np.transpose(np.array([k_y, -1*k_x, np.zeros(len(k_y))])/(np.sqrt(k_y**2 + k_x**2))) ### from Newberg & Yanny
        i_vecs = np.transpose(np.array([-1*k_x*k_z, -1*k_y*k_z, (k_x**2 + k_y**2)])/(np.sqrt(k_y**2 + k_x**2)))

        theta_matrix = np.ones(np.shape(i_vecs))
        mjax_matrix = np.ones(np.shape(i_vecs))

        for i in range(len(theta)):
            theta_matrix[i] = theta_matrix[i]*theta[i] ### so we can broadcast with j and i vecs

        l_vecs = (i_vecs*np.cos(theta_matrix) + j_vecs*np.sin(theta_matrix))
        m_vecs = -1*i_vecs*np.sin(theta_matrix) + j_vecs*np.cos(theta_matrix)

        ### set up l and k values for data since I will cut on this
        ### these arrays are the length of v_vec: will determine which to pull based on nearest locus point to data
        l_dat_all = np.zeros((len(k_vecs), len(data[0])))
        k_dat_all = np.zeros((len(k_vecs), len(data[0])))
        m_dat_all = np.zeros((len(k_vecs), len(data[0])))
        dist_euc = np.zeros((len(k_vecs), len(data[0])))

        ### rotations
        for i in range(len(k_vecs)):

            ### Determine if data points are in infinite cylinder

            ### rotate k along the z axis into the x-y plane
            rot_z_hyp = np.sqrt(k_vecs[i][0]**2 + k_vecs[i][1]**2)
            rot_z = np.array([[k_vecs[i][1]/rot_z_hyp, -1*k_vecs[i][0]/rot_z_hyp, 0], [k_vecs[i][0]/rot_z_hyp, k_vecs[i][1]/rot_z_hyp, 0], [0, 0, 1]])

            k_rotated_z = np.matmul(rot_z, k_vecs[i])
            l_rotated_z = np.matmul(rot_z, l_vecs[i])
            m_rotated_z = np.matmul(rot_z, m_vecs[i])

            ### rotate k along the x axis into the z-axis
            rot_x_hyp = np.sqrt(k_vecs[i][2]**2 + k_vecs[i][1]**2)
            rot_x = np.array([[1, 0, 0], [0, k_rotated_z[2]/rot_x_hyp, -1*k_rotated_z[1]/rot_x_hyp], [0, k_rotated_z[1]/rot_x_hyp, k_rotated_z[2]/rot_x_hyp]])

            k_rotated = np.matmul(rot_x, k_rotated_z)
            l_rotated = np.matmul(rot_x, l_rotated_z)
            m_rotated = np.matmul(rot_x, m_rotated_z)

            ### rotate m into the x axis
            rot_final_hyp = np.sqrt(m_rotated[0]**2 + m_rotated[1]**2)
            rot_final = np.array([[m_rotated[0]/rot_final_hyp, m_rotated[1]/rot_final_hyp, 0], [-1*m_rotated[1]/rot_final_hyp, m_rotated[0]/rot_final_hyp, 0], [0,0,1]])

            ### rotate all of the relevant vectors relative to the ellipse using these transformations

            lvec_final = np.matmul(rot_final, l_rotated)/np.linalg.norm(np.matmul(rot_final, l_rotated))
            mvec_final = np.matmul(rot_final, m_rotated)/np.linalg.norm(np.matmul(rot_final, m_rotated))
            kvec_final = np.matmul(rot_final, k_rotated)/np.linalg.norm(np.matmul(rot_final, k_rotated))

            ### want to translate data based on central locus point
            translation_mat = -1*loci[i]
            data_trans = np.array([data[0]+translation_mat[0], data[1]+translation_mat[1], data[2]+translation_mat[2]])

            ### rotate the translated data
            data_rotated = np.matmul(rot_final, np.matmul(rot_x , np.matmul(rot_z, data_trans)))

            l_dat_all[i] = data_rotated[1]
            k_dat_all[i] = data_rotated[2]
            m_dat_all[i] = data_rotated[0]



        ### determine which locus points are attributed with each galaxy

            dist_euc[i] = np.sqrt((data[0] -  x_locus[i])**2 + (data[1] -  y_locus[i])**2  + (data[2] -  z_locus[i])**2)

        dist_euc = np.transpose(dist_euc)

        l_dat_all = np.transpose(l_dat_all)
        k_dat_all = np.transpose(k_dat_all)
        m_dat_all = np.transpose(m_dat_all)


        l_dat = np.zeros(len(self.psf_mags[0]))
        k_dat = np.zeros(len(self.psf_mags[0]))
        m_dat = np.zeros(len(self.psf_mags[0]))



        for i in range(len(self.psf_mags[0])):


            if sum(np.isnan(dist_euc[i])) != len(dist_euc[i]):
            
                l_dat[i] =l_dat_all[i][np.where(dist_euc[i] == min(dist_euc[i]))[0]]
                k_dat[i] = k_dat_all[i][np.where(dist_euc[i] == min(dist_euc[i]))[0]]
                m_dat[i] = m_dat_all[i][np.where(dist_euc[i] == min(dist_euc[i]))[0]]


        ### low redshift sources
        if self.locus_name == 'ugri':

            ### magnitude limits
            magFlag = (i_psf_mag > 15) * (i_psf_mag < 19.1)  ###19.1 note this is for SDSS DR4 

            ### color cuts (this defines keep for extended objects)
            u_minus_gFlag = ug_color < 0.9 ### originally 0.9 in richards+2002 ###0.7 in DR4

            ### exclusion region (Not currently being used)

            ### white dwarfs
            
            wdFlag = (ug_color > -0.8)*(ug_color < 0.7)*(gr_color > -0.8)*(gr_color < -0.1)*(ri_color > -0.6)*(ri_color < -0.1)*(iz_color > -1)*(iz_color < -0.1)
            not_wdFlag = np.logical_not(wdFlag)


            ### A stars
            
            astarFlag = (ug_color > 0.7)*(ug_color < 1.4)*(gr_color > -0.5)*(gr_color < 0)*(ri_color > -0.5)*(ri_color < 0.2)*(iz_color > -0.4)*(iz_color < -0.2)
            not_astarFlag = np.logical_not(astarFlag)

            no_astar_wd_Flag = not_astarFlag*not_wdFlag

            ### red blue pairs

            redblueFlag = (gr_color > -0.3)*(gr_color < 1.25)*(ri_color > 0.6)*(ri_color < 2)*(iz_color > 0.4)*(iz_color < 1.2)
            not_redblueFlag = np.logical_not(redblueFlag)

            ### inclusion region
            includeFlag_midz = (ug_color > 0.6)*(ug_color < 1.5)*(gr_color > 0.0)*(gr_color < 0.2)*(ri_color > -0.1)*(ri_color < 0.4)*(iz_color > -0.1)*(iz_color < 0.4)*unresolved_Flag

            include_z_lt_02_Flag = (ug_color < 0.6)*magFlag*not_wdFlag



            ### for extended sources l < 0 and k < 0
            lk_flag = (l_dat <= 0)*(m_dat <= 0)

            mag_color_Flag_gen = magFlag*(np.logical_or(unresolved_Flag, (u_minus_gFlag*lk_flag)))*not_astarFlag*not_wdFlag*not_redblueFlag

            ### or statement between inclusion regions
            args = np.array((mag_color_Flag_gen, includeFlag_midz,include_z_lt_02_Flag*not_wdFlag))
            
            mag_color_Flag = np.logical_or.reduce(args) 


        ### high redshift sources
        elif self.locus_name == 'griz':

            ### magnitude limits
            magFlag = (i_psf_mag > 15) * (i_psf_mag < 20.2) ### note this is for SDSS DR4

            ### color cuts

            g_minus_rFlag = gr_color < 1.0
            u_minus_gFlag = ug_color >= 0.8
            

            ### general mag/color flag
            mag_color_Flag_gen = np.logical_not(g_minus_rFlag*u_minus_gFlag*np.logical_or((i_psf_mag >= 19.1), (ug_color < 2.5)))*magFlag


            ### exclusion regions

            ### white dwarfs
            
            wdFlag = (ug_color > -0.8)*(ug_color < 0.7)*(gr_color > -0.8)*(gr_color < -0.1)*(ri_color > -0.6)*(ri_color < -0.1)*(iz_color > -1)*(iz_color < -0.1)
            not_wdFlag = np.logical_not(wdFlag)


            ### A stars
            
            astarFlag = (ug_color > 0.7)*(ug_color < 1.4)*(gr_color > -0.5)*(gr_color < 0)*(ri_color > -0.5)*(ri_color < 0.2)*(iz_color > -0.4)*(iz_color < -0.2)
            not_astarFlag = np.logical_not(astarFlag)

            no_astar_wd_Flag = not_astarFlag*not_wdFlag

            ### red blue pairs

            redblueFlag = (gr_color > -0.3)*(gr_color < 1.25)*(ri_color > 0.6)*(ri_color < 2)*(iz_color > 0.4)*(iz_color < 1.2)
            not_redblueFlag = np.logical_not(redblueFlag)


            ### inclusion regions
            incl_1b_Flag = np.logical_or((ug_color > 1.5), (u_psf_mag > 20.6))
            incl_1d_Flag = np.logical_or((gr_color > 2.1), (ri_color < (0.44*gr_color - 0.358)))

            incl_1_Flag = incl_1b_Flag*(gr_color > 0.7)*(incl_1d_Flag)*(iz_color < 0.25)*(iz_color > -1)

            
            incl_2f_Flag = (iz_color < (0.52*ri_color - 0.412))

            incl_2_Flag = (u_psf_mag > 21.5)*(g_psf_mag > 21)*(ri_color > 0.6)*(iz_color > -1)*incl_2f_Flag


            incl_3f_Flag = (gr_color < (0.44*ug_color - 0.56))

            incl_3_Flag = (u_psf_mag > 20.6)*(ug_color > 1.5)*(gr_color < 1.2)*(ri_color < 0.3)*(iz_color > -1)*incl_3f_Flag

            args = np.array((mag_color_Flag_gen, incl_1_Flag, incl_2_Flag, incl_3_Flag))
            mag_color_Flag = np.logical_or.reduce(args)*not_astarFlag*not_wdFlag*not_redblueFlag
            

        return mag_color_Flag
          
    
    
    
    
    
    def is_targeted(self):
        
        ### import all necessary packages
        from astropy.cosmology import FlatLambdaCDM
        import numpy as np
        from matplotlib import pyplot as plt
        
        not_in_locus_Flag = self.not_in_locus()
        color_mag_Flag = self.color_mag_cuts_sdss()
        
        targeted_Flag = not_in_locus_Flag*color_mag_Flag
        
        return targeted_Flag
        
        






