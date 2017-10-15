# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import os
import numpy as np
import math
import matplotlib.pyplot as plt


#==============================================================================
#   constants
#==============================================================================

speed_of_light = 2.997 * 10**8          # [m/s]


#==============================================================================
#   functions
#==============================================================================

def straight_line(x, m, q):
    """ """
    return x*m + q


#==============================================================================
#   Glas class
#==============================================================================

def E_gauss_z0(t, t_0, omega):
    """gaussian shaped electric field pulse at z=0 """
    return np.exp(-(t / t_0)**2) * np.cos(t*omega)


def spectral_intensity(f, omega, t_0):
    """returns values of spectral intensity"""
    return np.exp(-2*np.pi**2 *(f - omega)**2 * t_0**2)
    

#==============================================================================
#   Glas class
#==============================================================================


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


class SiO2(Glas):
    """SiO2 class derived from glas class"""

    def __init__(self):
        """ """
        
        lam_min = 0.25      # [um]
        lam_max = 0.80      # [um]
        
        Bs = [0.696166300,
              0.407942600,
              0.897479400]
        
        Cs = [4.67914826*10**-3,
              1.35120631*10**-2,
              9.79340025*10]
              
        super().__init__(lam_min, lam_max, Bs, Cs)
        return


#==============================================================================
#   medium-interface class
#==============================================================================

def calc_rs(theta1, theta2, n1, n2):
    """returns fresnel coefficient r_s of reflection of s-polarized beam"""
    return (n1*np.cos(theta1) - n2*np.cos(theta2))/(n1*np.cos(theta1) + n2*np.cos(theta2)) 

def calc_ts(theta1, theta2, n1, n2):
    """returns fresnel coefficient t_s of transmission of s-polarized beam"""
    return 1 + calc_rs(theta1, theta2, n1, n2)

def calc_Rs(theta1, theta2, n1, n2):
    """returns beams reflection coefficient for s-polarized part"""
    return calc_rs(theta1, theta2, n1, n2)**2
    
def calc_Ts(theta1, theta2, n1, n2):
    """returns beams transmission coefficent for s-polarized part"""
    return 1 - calc_Rs(theta1, theta2, n1, n2)


def calc_rp(theta1, theta2, n1, n2):
    """returns fresnel coefficient r_p of reflection of p-polarized beam"""
    return (n1*np.cos(theta2) - n2*np.cos(theta1))/(n1*np.cos(theta2) + n1*np.cos(theta1))
    
def calc_tp(theta1, theta2, n1, n2):
    """returns fresnel coefficient t_p of reflection of p-polarized beam"""
    return (1 + calc_rp(theta1, theta2, n1, n2)) * (np.cos(theta1)/np.cos(theta2))

def calc_Rp(theta1, theta2, n1, n2):
    """returns beams reflection coefficient for p-polarized part"""
    return calc_rp(theta1, theta2, n1, n2)**2
    
def calc_Tp(theta1, theta2, n1, n2):
    """returns beams transmission coefficent for p-polarized part"""
    return 1 - calc_Rp(theta1, theta2, n1, n2)



class MediumInterface(object):
    """ """
    
    def __init__(self, n2, n1=1):
        
        self.n1 = n1
        self.n2 = n2
        return
    
    def theta2(self, theta1):
        """returns theta2 calculated with snells law and takes total reflection into account"""
        
        if not hasattr(theta1, '__iter__'): # check if theta1 is iterable
            theta1 = [theta1]        
        
        sin_theta2 = self.n1/ self.n2 * np.sin(theta1)
        theta2 = np.zeros(len(list(theta1)))
        
        # take total reflection into account
        for idx,val in enumerate(list(sin_theta2)):
            if math.isnan(np.arcsin(val)):
                theta2[idx] = np.pi/2
            else:
                theta2[idx] = np.arcsin(val)
        
        # check if theta2 is list/array or single number
        if len(theta2) > 1:
            return theta2
        else:
            return theta2[0]
                

    def theta2_(self, theta1):
        """theta2 with neglecting total reflection"""
        return np.arcsin(self.n1/ self.n2 * np.sin(theta1))
    
    def get_Rs(self, theta1):
        """returns reflection coeff of s-polarized beam at incident angle theta1"""
        return calc_Rs(theta1, self.theta2(theta1), self.n1, self.n2)
        
    def get_Ts(self, theta1):
        """returns transmission coeff of s-polarized beam at incident angle theta1"""
        return calc_Ts(theta1, self.theta2(theta1), self.n1, self.n2)
    
    def get_Rp(self, theta1):
        """returns reflection coeff of p-polarized beam at incident angle theta1"""
        return calc_Rp(theta1, self.theta2(theta1), self.n1, self.n2)
    
    def get_Tp(self, theta1):
        """returns transmission coeff of p-polarized beam at incident angle theta1"""
        return calc_Tp(theta1, self.theta2(theta1), self.n1, self.n2)
        
    

#==============================================================================
#   exercise sheet 2
#==============================================================================

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


#==============================================================================
#   exercise sheet 3
#==============================================================================

def solve_ex_sheet3():
    """runs all functions related to exercise sheet3 """
    
#    do_plot_es3_ex1()
    ax1, fig1 = do_plot_es3_ex2()
    
    

def do_plot_es3_ex1():
    """ """
    
    fig = plt.figure(figsize=(24, 12))
    
    
    t = np.linspace(-5, 5, 1001)* 10**-14
    t_0 = 10 * 10**-15 / (2*np.sqrt(np.log(2)))          # [s]
    omega = 5*10**15
    
    ax1 = fig.add_subplot(121)
    ax1.plot(t, E_gauss_z0(t, t_0, omega), '-r')    
    ax1.grid()
    
    
    f = omega + np.linspace(-5, 5, 1001) * 10**14
    I_w = spectral_intensity(f, omega, t_0)        
    
    ax2 = fig.add_subplot(122)
    ax2.plot(f, spectral_intensity(f, omega, t_0), '-g')
    ax2.grid()
    
    
    # ex1.f.
    
    
    sio2 = SiO2()
    
    return



def do_plot_es3_ex2():
    """ """
    
    # initialize variables
    mi = MediumInterface(1.5)    
    theta1 = np.linspace(0,0.5, 1001) * np.pi

    # plot
    fig1 = plt.figure(figsize=(14,6))
    
    ax1 = fig1.add_subplot(111)
    ax1.plot(theta1, mi.get_Rs(theta1), ':g', label='Rs')
    ax1.plot(theta1, mi.get_Ts(theta1), '--g', label='Ts')
    ax1.plot(theta1, mi.get_Rp(theta1), ':r', label='Rp')
    ax1.plot(theta1, mi.get_Tp(theta1), '--r', label='Tp')
    
    
    ax1.set_xlabel(r'$ \theta\ [rad] $')
    ax1.legend(loc=0)
    ax1.grid()
    

    return ax1, fig1
    
    
def multilayer(theta1):
    """ """
    
    n = 1.5
    nc = np.sqrt(n)
    
    # assume 3 different media with n1, n2, n3
    # points at which beam hits interfaces
    mi1 = MediumInterface(n1=1, n2=nc)
    mi2 = MediumInterface(n1=nc, n2=n)
    mi3 = MediumInterface(n1=nc, n2=1)
    
#    R_s1 = mi1.calc_rs(theta1)**2 + (mi3.calc_ts(theta1)*
#                                     mi2.calc_rs(theta1)*
#                                     mi1.calc_ts(tehta1))**2
    
    return R_s1

#==============================================================================
#   main function
#==============================================================================


if __name__ == "__main__":
    
#    solve_ex_sheet2()    
    solve_ex_sheet3()   
    
    
    pass