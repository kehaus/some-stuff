# -*- coding: utf-8 -*-
""" 
Contains tight binding examples


"""
__author__ = "kha"
__version__ = "0.0.1"



import sys
import os
import time


if 'fermi_3d' in sys.modules:
	print('fermi_3d.py reloaded')
	import importlib
	importlib.reload(fermi_3d)
else:
	import fermi_3d


from fermi_3d import *


# =====================================
# Settings
# =====================================
kx = np.linspace(-1.5, 1.0, 2001)*np.pi
ky = np.linspace(-1.5, 1.0, 2001)*np.pi
kz = np.linspace(-1.5, 1.0, 81)*np.pi

kx2 = np.linspace(-1.5, 1.0, 801)*np.pi
ky2 = np.linspace(-1.5, 1.0, 801)*np.pi

E_min = -1.5; E_max = 1.0; E_step = 0.001
E_lst = np.arange(E_min, E_max, E_step)


#simple_prm_3d = simple_prm.copy(); simple_prm_3d['kspace'] = [kx, ky, kz]
#simple_prm_2d = simple_prm.copy(); simple_prm_2d['kspace'] = [kx2, ky2]

lsco_prm_3d = lsco_prm.copy(); lsco_prm_3d['kspace'] = [kx/lsco_prm['a'], ky/lsco_prm['a'], kz/lsco_prm['c_']]
lsco_prm_2d = lsco_prm.copy(); lsco_prm_2d['kspace'] = [kx2/lsco_prm['a'], ky2/lsco_prm['a']]
[lsco_prm_2d.pop(key) for key in ['c', 'c_', 'tz']]


#lsco_prm_3d['tz'] = 0.14*0.25
lsco_prm_4 = lsco_prm_3d.copy(); lsco_prm_4['tz'] = 0

# =====================================
# TB - simple
# =====================================
#tb_3d = TightBindingBS(simple_prm_3d, tb_func=e_simple)
#tb_2d = TightBindingBS(simple_prm_2d, tb_func=e_simple_2d)
#
## dos
#tb_3d.dos = TightBindingBS.calc_dos(tb_3d.tb_dset, E_lst, E_step, verbose=True)
#tb_3d.dos = tb_3d.dos / np.sum(tb_3d.dos)
#TightBindingBS.plot_dos(E_lst, tb_3d.dos, num=2)
#
#tb_2d.dos = TightBindingBS.calc_dos(tb_2d.tb_dset, E_lst, E_step, verbose=True)
#tb_2d.dos = tb_2d.dos / np.sum(tb_2d.dos)
#TightBindingBS.plot_dos(E_lst, tb_2d.dos, num=2)


# =====================================
# TB - paper
# =====================================
time0 = time.time()

tb2_3d = TightBindingBS(lsco_prm_3d, tb_func=e_3d, verbose=True)
print('time: {:.2f}'.format(time.time()-time0))
tb2_3d.dos = TightBindingBS.calc_dos(tb2_3d.tb_dset, E_lst, E_step, verbose=False)
tb2_3d.dos2 = TightBindingBS.calc_dos_kz_average(tb2_3d.tb_dset, E_lst, E_step, verbose=False)
TightBindingBS.plot_dos(E_lst, tb2_3d.dos2, num=4)

print('time: {:.2f}'.format(time.time()-time0))

tb2_4 = TightBindingBS(lsco_prm_4, tb_func=e_3d)
print('time: {:.2f}'.format(time.time()-time0))
tb2_4.dos = TightBindingBS.calc_dos(tb2_4.tb_dset, E_lst, E_step, verbose=False)
TightBindingBS.plot_dos(E_lst, tb2_4.dos, num=4)

print('time: {:.2f}'.format(time.time()-time0))

#tb2_2d = TightBindingBS(lsco_prm_2d, tb_func=e_2d)
#tb2_2d.dos = TightBindingBS.calc_dos(tb2_2d.tb_dset, E_lst, E_step, verbose=True)
#TightBindingBS.plot_dos(E_lst, tb2_2d.dos, num=3)


from scipy.integrate import trapz
from scipy.interpolate import interp1d
# transfrom from energy E_lst to doping n
tb2_3d_xx = [trapz(tb2_3d.dos[:idx]) for idx in range(len(tb2_3d.dos))]
#tb2_3d_f = interp1d(E_lst, tb2_3d_xx, kind='cubic')
#tb2_3d_n = tb2_3d_f(E_lst)

tb2_4_xx = [trapz(tb2_4.dos[:idx]) for idx in range(len(tb2_4.dos))]
#tb2_4_f = interp1d(E_lst, tb2_4_xx, kind='cubic')
#tb2_4_n = tb2_4_f(E_lst)








# Sr2Ru4O7 TB model

def e_Sr2Ru4O7(ky, kx, t_1, t_2, t_3, t_4, t_5, t_6, mu, lam):
	"""

	**check formulas again**
	 """
	e_yz = -2*t_2*np.cos(kx) - 2*t_1*np.cos(ky)
	e_xz = -2*t_1*np.cos(kx) - 2*t_2*np.cos(ky)
	e_xy = -2*t_3*(np.cos(kx) + np.cos(ky)) - 4*t_4*np.cos(kx)*np.cos(ky) 
	e_xy += -2*t_5*(np.cos(kx) + np.cos(ky))
	e_off = -4*t_6*np.sin(kx)*np.sin(ky)

	# build matrix
	A = np.matrix([
		[e_yz-mu,      e_off + 1j*lam, -lam   ],
		[e_off-1j*lam, e_xz-mu,        1j*lam ],
		[-lam,         -1j*lam,        e_xy-mu]
	])


	return



