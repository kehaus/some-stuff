# -*- coding: utf-8 -*-
""" 
Tight binding calculations for *d Fermi surface of LSCO

follows the caluclations presented in *Three-Dimensional Fermi Surface of Overdoped 
La-base Cuprates



TODO:
    * kx,ky,kz input variables are neglected if one of those is None in __init__(). 
    change this implementation to a more accurate bahavoir
    * implement test which check if prm kwargs match to kwargs of tb_func in __init__()
    * include functionality to plot multiple dos-datasets in the same figure/axes


"""
__author__ = "kha"
__version__ = "0.0.1"


import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


# =====================================
# Tight binding Functions
# =====================================
def e_3d(kx, ky, kz, t=1, t_=1, t__=1, tz=1, c=1, c_=1, a=1, e0=1):
    """ 2d tight binding model """
#    print('t  :', t)
#    print('t_ :', t_)
#    print('t__:', t__)
#    print('tz :', tz)
#    print('c  :', c)
#    print('a  :', a)
#    print('e0 :', e0)
    res = e_2d(kx, ky, t=t, t_=t_, t__=t__, a=a, e0=e0)
    res += e_z(kx, ky, kz, tz=tz, c_=c_, a=a)
    return res

def e_2d(kx, ky, t=1, t_=1, t__=1, a=1, e0=0):
    """ 2d tight binding model"""
    res = e0 + 2*t * (np.cos(kx*a) + np.cos(ky*a))
    res += 4*t_ * np.cos(kx*a) * np.cos(ky*a)
    res += 2*t__ * (np.cos(2*kx*a) + np.cos(2*ky*a)) 
    return res

def e_z(kx, ky, kz, tz=1, c_=1, a=1):
    """ kz.dispersion tight binding model """
    res = 2*tz * np.cos(kz*c_) * (np.cos(kx*a) - np.cos(ky*a))**2
    res *= np.cos(kx*a/2) * np.cos(ky*a/2)
    return res

def e_simple(kx ,ky, kz, t=1, a=1, e0=0):
    """3d, nearest neighbour hopping"""
    res = e0 + 2*t*(np.cos(kx*a) + np.cos(ky*a) + np.cos(kz*a))
    return res

def e_simple_2d(kx, ky, t=1, a=1, e0=0):
    """2d, nearest neighbout hopping"""
    return e_simple(kx, ky, kz=np.pi/(2*a), a=a, e0=e0)

# =====================================
# parameter
# =====================================
t_lsco = 0.25	# not kown correct value !!
lsco_prm = {
    'a':	3.76,		# [Angstrom]
    'c':	13.22,		# [Angstrom]
    'c_':	13.22/2,	# [Angstrom]
    't':	t_lsco,
    't_':	-0.12*t_lsco,		
    't__':	0.06*t_lsco,
    'e0':	-0.93*t_lsco,
    'tz':	0.07*t_lsco,      # correct parameter
}

t_leusco = 1 	# not kown correct value !!
Leusco_prm = {
    'a':	3.76,		# [Angstrom]
    'c':	13.14,		# [Angstrom]
    'c_':	13.14/2,	# [Angstrom]
    't':	t_leusco,
    't_':	-0.14*t_leusco,		
    't__':	0.07*t_leusco,
    'e0':	-0.95*t_leusco,
    'tz':	0.07*t_leusco,
}

simple_prm = {
    'a': 3.76,
    't': 1,
    'e0': 0
}

# =====================================
# plot function
# =====================================
def plot_contour(kz_cut_lst, num, levels):
    """ """
    for idx, kz_cut in enumerate(kz_cut_lst):
        fig1 = plt.figure(num=num+idx); fig1.clf()
        plt.contourf(kX[:,:,kz_cut],kY[:,:,kz_cut], tb[:,:,kz_cut], levels=levels)
        plt.show()
    return


# =====================================
# Plot settings
# =====================================
AX_DCT = {
    'xlabel':   r'$ k_x\ [\pi/a] $',
    'ylabel':   r'$ k_y\ [\pi/b] $',
    'zlabel':   r'$ E(k_x, k_y)\ [eV] $',
}

AX_DCT_DOS = {
    'xlabel':   r'$ E\ [eV] $',
    'ylabel':   r'$ DOS $',
}

LINE_DCT1 = {   # matplotlib.lines.Line2D
#    'color':        'r',
    'linewidth':    1.5,
}



# =====================================
# TB - Exception
# =====================================
class TightBindingError(Exception):
    """ """
    pass


# =====================================
# TB - class
# =====================================
class BandStructure():
    """ """
    def __init__(self):
        """ """

        self.type = None


class TightBindingBS(BandStructure):
    """ """

    def __init__(self, prm, tb_func=None, verbose=False):
        """

        TODO: implement test which check if prm kwargs match to kwargs of tb_func

         """
        self.type = 'tb'
        self.verbose = verbose

        self.kspace = prm.pop('kspace')
        self.dim = len(self.kspace)
        self.kspace = np.meshgrid(*self.kspace)
        self.prm = prm

        if tb_func != None:
            self.set_tb_model(tb_func)

    def _get_default_kspace(self):
        """

        **not used**

         """
        kx = np.linspace(-1.5, 1.5, 101)*np.pi
        ky = np.linspace(-1.5, 1.5, 101)*np.pi
        kz = np.linspace(-2.5, 2.5, 101)*np.pi
        kX, kY, kZ = np.meshgrid(kx, ky, kz)
        return kX, kY, kZ

    def tb_model(self):
        """tight-binding model. Must be defined"""
        pass

    def set_tb_model(self, tb_func):
        """ sets tb_model and calls calc_tb function"""
        self.tb_model = tb_func
        self.calc_tb(verbose=self.verbose)

    def calc_tb(self, verbose=False):
        """calculates tight binding dset with self.tb_model"""
        if verbose: 
            print('calculate tight binding dset ...')
            print('kspace dimension: ', self.kspace[0].shape)

        self.tb_dset = self.tb_model(*self.kspace, **self.prm)
        if verbose: print('tight binding dset calculated.')


    def plot_tb_3Drelief(self, kz_idx, num=1, figsize=(14,12), E_F=2.2):
        """ """
        TightBindingBS.plot_3Drelief(*self.kspace, self.tb_dset, 
            kz_idx, num=num, figsize=figsize, E_F=E_F)

    @staticmethod
    def calc_dos_old(tb_dset, E_lst, E_step, verbose=False):
        """calculates normalized DOS by summing over 3D iso-energy surface"""
        dos = []
        for idx, E in enumerate(E_lst):
            bool_arr1 = np.ones(tb_dset.shape)*E <= tb_dset
            bool_arr2 = tb_dset < np.ones(tb_dset.shape)*(E+E_step)
            dos.append(np.logical_and(bool_arr1, bool_arr2).sum())
            del(bool_arr1, bool_arr2)
            if verbose: print('idx: {:d}'.format(idx))
        return dos/np.sum(dos)

    @staticmethod
    def calc_dos(tb_dset, E_lst, E_step, verbose=False):
        """calculates normalized DOS by using ``np.histogram`` """
        dset_ = tb_dset.reshape(1, tb_dset.size)
        E_bin = np.concatenate((np.array([E_lst[0]]),
                                np.array(E_lst)+E_step/2.), axis=0)
        dos, b = np.histogram(dset_, bins=E_bin)
        return dos/np.sum(dos)


    @staticmethod
    def calc_dos_kz_average(tb_dset, E_lst, E_step, verbose=False):
        """calc DOS by averaging over kz and summing over 2D iso-energy surface"""
        if len(tb_dset.shape) != 3:
            raise TightBindingError("tb_dset dimension need to be 3!")
        
        dos_lst = []       # fill with DOS for every kz value
        for idx in range(tb_dset.shape[2]):
            if verbose: print('    ', end='')
            dos = TightBindingBS.calc_dos(tb_dset[:,:,idx], E_lst, E_step, verbose=verbose)
            dos_lst.append(dos)
        if verbose: print('idx: {:d}'.format(idx))
        return np.mean(dos_lst, axis=0)

    @staticmethod
    def plot_3Drelief(kX, kY, kZ, tb_dset, kz_idx, num=1, E_F=2.2, figsize=(14,12), ax_dct=None):
        """ """

        if ax_dct == None:
            ax_dct = AX_DCT

        fig3 = plt.figure(num=num, figsize=figsize); fig3.clf()
        ax = fig3.gca(projection='3d', **ax_dct)
        cont = ax.contour(kX[:,:,kz_idx],kY[:,:,kz_idx], tb_dset[:,:,kz_idx], 
        levels=np.linspace(-4.5,4.5,151), cmap=cm.coolwarm)

        cont = ax.contour(kX[:,:,kz_idx],kY[:,:,kz_idx], tb_dset[:,:,kz_idx], 
        levels=[E_F], color='k')
        ax.set_title(r'$ k_z = {:.1f}\ \pi/c $'.format(kZ[0,0,kz_idx]/np.pi))

        ax.view_init(elev=50.,azim=90)
        
    @staticmethod
    def plot_dos(E_lst, dos, num=1, figsize=(8,6), ax_dct=None, line_dct=None):
        """ """
        if ax_dct == None:
            ax_dct = AX_DCT_DOS

        if line_dct == None:
            line_dct = LINE_DCT1

        fig = plt.figure(num=num, figsize=figsize)
        ax = fig.add_axes([0.1, 0.1, 0.85, 0.85], **ax_dct)
        ax.plot(E_lst, dos, **line_dct)
        ax.grid()



# =====================================
# main
# =====================================
if __name__ == "__main__":
#    if 'tb1' not in locals().keys():
#        tb1 = TightBindingBS(lsco_prm, tb_func=e_3d)
#    if not hasattr(tb1, 'dos1'):
#        tb1.dos1 = TightBindingBS.calc_dos_3d(tb1.tb_dset, E_lst, E_step, verbose=True)
#    if not hasattr(tb1, 'dos2'):
#        tb1.dos2 = TightBindingBS.calc_dos_2(tb1.tb_dset, E_lst, E_step, verbose=True)

    if 'tb2' not in locals().keys():
        tb2_3d  = TightBindingBS(simple_prm, tb_func=e_simple)
    if not hasattr(tb2_3d, 'dos1'):
        tb2_3d.dos1 = TightBindingBS.calc_dos_3d(tb2_3d.tb_dset, E_lst, E_step, verbose=True)
    if not hasattr(tb2_3d, 'dos2'):
        tb2_3d.dos2 = TightBindingBS.calc_dos_2(tb2_3d.tb_dset, E_lst, E_step, verbose=True)


    plt.plot(E_lst, tb2_3d.dos1)

#   pass
