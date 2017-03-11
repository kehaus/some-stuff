# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import os
import numpy as np
import matplotlib.pyplot as plt


#==============================================================================
#   constants
#==============================================================================

speed_of_light = 2.997 * 10**8          # [m/s]

#==============================================================================
#   exercise sheet 2
#==============================================================================

def straight_line(x, m, q):
    """ """
    return x*m + q


def sellmeier_eq(lam, Bs, Cs):
    """returns square of refration index calculated by sellmeier equation"""
    
    n2 = 1    
    for i in range(len(Bs)):
        n2 += Bs[i]*lam**2 / (lam**2 - Cs[i])
        
    return n2 #np.array(n2)
    

class Glas(object):
    """basic Glas-class """
    
    def __init__(self, lam_min, lam_max, Bs, Cs):
        """ """
        self.lam_min = lam_min
        self.lam_max = lam_max
        self.Bs = Bs
        self.Cs = Cs        
        return


    def calc_n2(self, lam):
        """returns square of refraction index n at wavelenth lam"""
        return sellmeier_eq(lam, self.Bs, self.Cs)


    def calc_n(self, lam):
        """returns refraction index n at the wavelength lam"""        
        return np.sqrt(self.calc_n2(lam))
        
        
    def get_n_derivative(self, lam, d_lam=0.00001):
        """returns derivative of n at position lam
        
        derivaive is calculated with centeral difference method with a 
        stepsize equal to d_lam
        
        """
        return (self.calc_n(lam + d_lam) - self.calc_n(lam - d_lam))/ (2*d_lam)
        
        
    def get_n_II_derivative(self, lam, d_lam=0.00001):
        """returns second derivative of n at position lam """
        
        return (self.calc_n(lam + 2*d_lam) - 2*self.calc_n(lam) + self.calc_n(lam - 2*d_lam)) / (4*d_lam**2)
        
    
    def get_straight_line_at_n(self, lam, lam_0):
        """ """
        
        m = self.get_n_derivative(lam_0)
        q = self.calc_n(lam_0) - lam_0*m
        return straight_line(lam, m, q)
        
        
    def get_phase_velocity(self, lam):
        """returns phase velocity [m/s] at given wavelength. eq.: (2.46)"""
        return speed_of_light / self.calc_n(lam)
        

    def get_group_velocity(self, lam):
        """returns group velocity [m/s] at given wavelength. eq.: (2.57)"""
        return speed_of_light / (self.calc_n(lam) - lam*self.get_n_derivative(lam))

        
class BK_7(Glas):
    """BK-7 class derived from glas class  """
    
    
    def __init__(self):
        """ """ 
        
        lam_min = 0.25     # [um]
        lam_max = 3.00     # [um]
        
        Bs = [1.03961212,
              0.231792344,
              1.01046945]
              
        Cs = [6.00069867*10**-3,
              2.00179144*10**-2,
              1.03560653*10**2]
        
        super().__init__(lam_min, lam_max, Bs, Cs)
        return
        



def solve_ex_sheet2():
    """ """
    
    ## ex1
    # do calculations
    bk_7 = BK_7() 
    lam = np.arange(0.25, 3.05, 0.01)
    
    # do plot
    fig1, ax1 = do_plot_es2_ex1(bk_7, lam)
    

    # print results    
    print('refraction index is for wavelengths larger than 1.5 um anormal')
    print('this is visible because the derivation is negativ in this region')
    print('')
    
    print('n: ', bk_7.calc_n(0.4))
    print('deriv n: ', bk_7.get_n_derivative(0.4))
    print('')

    print('phase_velocity at 800nm = {:.4f}c'.format(bk_7.get_phase_velocity(0.8)/speed_of_light))
    print('group_velocity at 800nm = {:.4f}c'.format(bk_7.get_group_velocity(0.8)/speed_of_light))
    print('')    
    
    print('phase_velocity at 400nm = {:.4f}c'.format(bk_7.get_phase_velocity(0.4)/speed_of_light))
    print('group_velocity at 400nm = {:.4f}c'.format(bk_7.get_group_velocity(0.4)/speed_of_light))
    print('')    

    print('distance at which wavepackets will be by t= 1ps')
    print('s = {:.4f}um'.format(10**6 * get_distance(10**-12, 
                                                     bk_7.get_group_velocity(0.4),
                                                     bk_7.get_group_velocity(0.8))))    
    print('\n \n \n')
    
      
    ## ex2
    lam2, delta, beta = import_data_es2_ex2()  
     
    fig2, ax2 = do_plot_es2_ex2(lam2, delta, beta)
    
    
    
    return


def get_distance(dt, v1, v2):
    """returns distance at which wavepacktes with velocities v1, v2 are 
    displaced by time dt
    
    """
    
    return dt * (v2-v1)    
#    return dt * v1*v2/(v1-v2)


def do_plot_es2_ex1(bk_7, lam, saving_path=None):
    """ """

    fig = plt.figure(figsize=(15,6))
    ax = fig.add_subplot(111)

    ax.plot(lam, bk_7.calc_n(lam), '-r', label='refraction index for BK7')
    ax.plot(lam, bk_7.get_straight_line_at_n(lam, 0.8), '--g', label='tangent at n(0.8)')
    ax.plot(lam, bk_7.get_straight_line_at_n(lam, 0.4), ':g', label='tangent at n(0.4)')
    ax.set_title('refraction index for BK-7 glas calculated with sellmeier' + 
                 'equation')
    ax.set_xlabel('wavelength [um]')
    ax.set_ylabel('refraction index n')
    ax.set_ylim([1.45, 1.60])
    ax.legend(loc=0)
    ax.grid()
    
    if saving_path == None:
        fig.savefig(os.getcwd() +'//'+ 'qe_es2_ex1_plot.png')
    
    return fig, ax


def import_data_es2_ex2():
    """ """
    
    fn = 'bk_7_refractive index'
    
    data = np.array(np.loadtxt(fn, skiprows=2))
    return data[:,0], data[:,1], data[:,2]


def do_plot_es2_ex2(lam, delta, beta, saving_path=None):
    """ """
    
    fig = plt.figure(figsize=(15,6))
    ax = fig.add_subplot(111)
    
    ax.loglog(lam, delta, '-r', label='delta')
    ax.loglog(lam, beta, '-g', label='beta')
    ax.set_title('index of refraction of Si64B24Na20Al18O201')
    ax.set_ylabel('delta beta,')
    ax.set_xlabel('wavelength [nm]')
    ax.set_xlim([0.1,10])
    ax.legend(loc=0)
    ax.grid()
    
    if saving_path == None:
        fig.savefig(os.getcwd() +'//'+ 'qe_es2_ex2_plot.png')
    
    return fig, ax


if __name__ == "__main__":
    
    solve_ex_sheet2()    
    
    
    
    pass