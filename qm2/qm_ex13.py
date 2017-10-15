# -*- coding: utf-8 -*-
"""
Created on Mon May 29 13:37:03 2017

Calculations to exercise sheet 13 from quantum mechanics 2


@author: kh
"""

import scipy.constants as const
import numpy as np


# =============================================================================
# constans
# =============================================================================

kB = const.k            # [J/K]
hbar = const.hbar       # [Js]
m_p = const.m_p         # [kg]
m_n = const.m_n
m_deuteron = const.physical_constants['deuteron mass'][0]       # [kg]


# =============================================================================
# ex1
# =============================================================================

## parameter
n = 0.9983
T = 20          # [K]

## calc ex1.b.
a = hbar / np.sqrt(np.log(9*n) * 2*m_p*kB*T)


## calc ex1.c. 
n = 1./9 * np.exp(hbar**2 / (a**2*2*(m_deuteron)*kB*T))


## print result ##
print("*** ex1.b. ***")
print('Distance between protons a = {:.2e}m'.format(a))
print('covalent radius of hydrogen molecules (from wikipedia): {:.2e}m'.format(25*10**-12))
print("\n")

print("*** ex1.c. ***")
print('mixing ratio for Deuterium n = {:.2f}'.format(n))
print("")
